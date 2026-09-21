# Rebuild presentation outputs without entering analysis/cache invalidation paths.
# Core result RDS remains unchanged. C1 companion plots need selected-feature
# data and matching; C4 independently validates/updates its imaging-only cache.
le8_figure_result <- function(path, fields = character()) {
  required_file(path, "result for output regeneration")
  z <- tryCatch(readRDS(path), error = function(e)
    stop("Cannot read result for output regeneration: ", path, ": ", conditionMessage(e), call. = FALSE))
  if (!is.list(z) || length(setdiff(fields, names(z))))
    stop("Incomplete result for output regeneration: ", path,
      "; required fields: ", paste(fields, collapse = ", "), call. = FALSE)
  z
}

le8_figure_csv <- function(rawdir, file) {
  as_tibble(data.table::fread(required_file(file.path(rawdir, file),
    "table for output regeneration"), data.table = FALSE, check.names = FALSE))
}

# Older C1 results did not embed this table. Recover it from the same plotting
# cohort and cached analysis metadata when the optional CSV has been removed.
le8_figure_c1_cohort <- function(obj, dat, layer, rawdir) {
  if (!is.null(obj$cohort)) return(obj$cohort)
  path <- file.path(rawdir, "c1.cohort.csv")
  if (file.exists(path)) return(le8_figure_csv(rawdir, "c1.cohort.csv"))
  message("C1: rebuilding missing cohort summary from trajectory inputs and cached metadata")
  meta <- obj$meta
  outcome <- meta$trait
  events <- sum(dat[[paste0(outcome, ".Yt2e")]] == 1, na.rm = TRUE)
  if ((!is.null(meta$N) && nrow(dat) != meta$N) ||
      (!is.null(meta$events) && events != meta$events))
    stop("C1 cohort inputs no longer match the cached sample/event counts; rerun C1 with --replace TRUE", call. = FALSE)
  dat <- add_attained_age_time(dat, outcome)
  audit <- obj$input_feature_annotation_audit
  removed <- meta$le8_covariates_removed %||% character()
  cohort <- tibble(layer = layer, N_omics = nrow(dat), incident_events = events,
    prevalent_cases = sum(make_prevalent_status(dat, outcome) == 1, na.rm = TRUE),
    features = length(audit$feature %||% unique(obj$association_adj2$term)),
    annotation_matched = if (is.null(audit$annotation_matched)) NA_integer_ else sum(audit$annotation_matched, na.rm = TRUE),
    annotation_unmatched = if (is.null(audit$annotation_matched)) NA_integer_ else sum(!audit$annotation_matched, na.rm = TRUE),
    attained_age_N = sum(is.finite(dat$.attained_entry) & is.finite(dat$.attained_exit)),
    bi2e_available = sum(is.finite(dat[[paste0(outcome, ".bi2e")]])),
    PGS_matched = meta$pgs_matched %||% NA_integer_,
    PGS_file_signature = meta$pgs_signature %||% NA_character_,
    LE8_components = length(intersect(vars.le8, names(dat))),
    covs_use = meta$covs_use %||% NA_character_,
    LE4_components = length(meta$le4_covariates),
    le8_covariates_removed = paste(removed, collapse = ";"),
    overlap_covariates_removed = paste(meta$overlap_covariates_removed %||% removed, collapse = ";"),
    treatment_covariates = paste(meta$treatment_covariates, collapse = ";"))
  write_raw_csv(cohort, "c1.cohort.csv", rawdir)
  cohort
}

# Explicit workbook mappings retain the original sheet names and avoid exporting
# fitted model objects or individual-level data that were never in the workbook.
le8_figure_tables <- function(obj, mapping) {
  lapply(mapping, function(path) {
    z <- obj
    for (key in strsplit(path, "$", fixed = TRUE)[[1]]) z <- z[[key]]
    z %||% tibble()
  })
}

le8_figure_workbook <- function(tables, file, review = list()) {
  # Write once: reopening a 200+ MB workbook just to append review sheets can
  # cost more time and memory than all of the figures combined.
  if(length(review)) {
    names(review)<-substr(paste0("review_",names(review)),1,31)
    tables<-c(tables[!names(tables)%in%names(review)],review)
  }
  message("Writing workbook: ",file)
  write_xlsx2(tables,file)
}

le8_restore_outputs <- function(layer, module) {
  outdir <- if (layer == "protein") out.prot else out.met
  rawdir <- le8_job_dir(outdir, module)
  oldwd <- getwd(); on.exit(setwd(oldwd), add = TRUE); setwd2(outdir)
  cache <- file.path(rawdir, paste0(sub("_.*$", "", module), ".res.rds"))
  fields <- switch(module,
    c1_correlate=c("association", "clusters", "cluster_selection", "enrichment"),
    c2_cause=c("MR", "MR_reverse", "MR_best", "observational"),
    c3_coloc=c("summary", "regional", "variants"),
    c4_connect=c("scan", "membership", "modules", "mediation"),
    c5_consolidate=c("scores", "prediction", "summary", "score_vs_topN"))
  message("Rebuilding ", module, "/", layer, " PNG/XLSX; reading cached result: ", cache)
  obj <- le8_figure_result(cache, c("meta", fields))
  le8_check_options(obj)
  if (!is.null(obj$meta$trait) && !identical(obj$meta$trait, Y))
    stop("Result trait does not match --trait: ", cache, call. = FALSE)
  if (!is.null(obj$meta$layer) && !identical(obj$meta$layer, layer))
    stop("Result layer does not match --biom: ", cache, call. = FALSE)
  # Load only definitions, without revision initialization/on-exit analyses.
  revisions <- switch(module,
    c2_cause = c("c2_revision_dandelion.R", "c2_revision_evidence.R"), character())
  for (f in revisions) sys.source(file.path(Sys.getenv("LE8_FDIR"), f), envir = .GlobalEnv)
  renderer <- get(paste0("le8_render_", sub("_.*$", "", module)), mode = "function")
  renderer(obj, layer, outdir, rawdir)
  finalize_outputs(module, outdir)
  invisible(obj)
}

le8_render_c1 <- function(obj, layer, outdir, rawdir) {
  needed <- c("association_adj2", "prevalent_adj2", "birthline_adj2", "pgs_incident", "pgs_prevalent", "pgs_attained_age")
  if (length(setdiff(needed, names(obj)))) stop("C1 result lacks fields needed for figures: ", paste(setdiff(needed, names(obj)), collapse=", "))
  aliases <- c(assoc="association", assoc_prevalent="prevalent", assoc_basic="association_basic",
    assoc_adj2="association_adj2", prevalent_basic="prevalent_basic", prevalent_adj2="prevalent_adj2",
    birthline_basic="birthline_basic", birthline_adj2="birthline_adj2",
    assoc_adj2_full_le8="association_adj2_full_le8_sensitivity", prevalent_adj2_full_le8="prevalent_adj2_full_le8_sensitivity",
    birthline_adj2_full_le8="birthline_adj2_full_le8_sensitivity", reverse_adj2="reverse_prevalent",
    duration_adj2="prevalent_duration", landmark_adj2="landmark_incident", risk_window_adj2="diagnosis_window_riskset",
    attenuation_sameN="attenuation", enrich="enrichment", enrich_prev="enrichment_prevalent",
    pgs_incident="pgs_incident", pgs_prevalent="pgs_prevalent", pgs_attained_age="pgs_attained_age",
    pgs_incident_same="pgs_incident_same_omic", pgs_prevalent_same="pgs_prevalent_same_omic",
    pgs_attained_age_same="pgs_attained_age_same_omic", pgs_concordance="pgs_actual_concordance",
    input_feature_audit="input_feature_annotation_audit")
  list2env(le8_figure_tables(obj, aliases), envir=environment())
  features_all <- input_feature_audit$feature %||% unique(assoc_adj2$term)
  covs_adj2 <- obj$meta$covariates
  if (is.null(covs_adj2)) stop("C1 cache does not record plotting covariates; rerun C1 explicitly first")
  fig_features <- unique(c(head(obj$top_features, YY_TOP), obj$clusters$feature,
    assoc$term[order(assoc$p.value)][seq_len(min(nrow(assoc), max(TOP_N, CLUSTER_TOP)))],
    assoc_adj2$term[order(assoc_adj2$p.value)][seq_len(min(nrow(assoc_adj2), GRADIENT_TOP))],
    prevalent_adj2$term[order(prevalent_adj2$p.value)][seq_len(min(nrow(prevalent_adj2), max(GRADIENT_TOP, CLUSTER_TOP)))]))
  message("C1: reading selected-feature trajectory inputs; association/PGS scans and enrichment are reused")
  biom <- if (layer=="protein") read_prot() else read_met()
  fig_features<-unique(c(fig_features,assoc_adj2$term[is.finite(assoc_adj2$p.value)&assoc_adj2$p.value<.05/nrow(assoc_adj2)]))
  missing <- setdiff(fig_features, names(biom))
  if(length(missing)) stop("C1 plotting inputs lack cached features: ",paste(missing,collapse=", "))
  biom <- biom[,unique(c("eid",fig_features)),drop=FALSE]; invisible(gc())
  need <- unique(c("eid","ethnic.c",vars.basic,vars.le8,covs_adj2,"birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y)))
  dat <- read_all(need) |> filter_analysis_cohort() |> inner_join(biom,by="eid") |> make_outcome(Y)
  rm(biom); invisible(gc())
  if(length(setdiff(covs_adj2,names(dat)))) stop("C1 plotting inputs lack recorded covariates")
  cohort <- le8_figure_c1_cohort(obj, dat, layer, rawdir)
  covs_basic <- intersect(vars.basic,names(dat))
  tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e");bvar<-paste0(Y,".b2e")
  pgs<-list(status=obj$pgs_status)
  vldl_deep <- if(layer=="metabolite") build_vldl_tg_deep_dive(
    list(incident=pgs_incident,prevalent=pgs_prevalent,attained_age=pgs_attained_age),
    list(incident=pgs_incident_same,prevalent=pgs_prevalent_same,attained_age=pgs_attained_age_same),
    list(incident=assoc_adj2,prevalent=prevalent_adj2,attained_age=birthline_adj2),
    list(incident=assoc_basic,prevalent=prevalent_basic,attained_age=birthline_basic),
    list(incident=assoc_adj2_full_le8,prevalent=prevalent_adj2_full_le8,attained_age=birthline_adj2_full_le8),
    risk_window_adj2,obj$L_VLDL_TG_pct_pgs_conditional,obj$L_VLDL_TG_pct_measured_conditional) else
      list(data=tibble(),trajectory=tibble(),pgs_conditional=tibble(),measured_conditional=tibble())
  top<-assoc|>filter(is.finite(p.value))|>slice_min(p.value,n=TOP_N,with_ties=FALSE)|>pull(term)
  top6<-head(top,YY_TOP);top_inc<-assoc_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=GRADIENT_TOP,with_ties=FALSE)|>pull(term)
  top_prev<-prevalent_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=GRADIENT_TOP,with_ties=FALSE)|>pull(term)
  cluster_features<-obj$clusters|>filter(analysis=="incident_behavioral_LE4")|>pull(feature)

  # Fig1 + Fig2: the third row is an attained-age sensitivity analysis with
  # delayed entry at baseline. The left column is a fixed-at-conception omic
  # PGS in the full genetic cohort; the right column is the measured adult omic
  # level with behavioral-LE4 adjustment. Basic observed-level analyses remain in the
  # workbook and raw files but are no longer displayed here.
  plot_c1_fig12(layer,features_all,pgs_incident,pgs_prevalent,pgs_attained_age,
    assoc_adj2,prevalent_adj2,birthline_adj2,outdir)

  # Fig3: same participants in basic and behavioral-LE4 residualized YY panels.
  yy_need<-unique(c(bvar,top6,covs_adj2));yy_dat<-dat[complete.cases(dat[,intersect(yy_need,names(dat)),drop=FALSE]) & is.finite(dat[[bvar]]),intersect(yy_need,names(dat)),drop=FALSE]
  omic_name<-ifelse(layer=="protein","protein","metabolite")
  p3a<-make_yy_panel(yy_dat,top6,bvar,covs_basic,paste0("a. Yin–Yang ",omic_name," patterns — basic adjustment"),
                     ifelse(layer=="protein","Mean protein z-score","Mean metabolite z-score"),smooth_lines=FALSE)
  adj2_label<-if(le8_custom_adjustment())paste0("adjusted for ",paste(le8_custom_covars,collapse=", ")) else "basic + behavioral LE4 adjustment"
  p3b<-make_yy_panel(yy_dat,top6,bvar,covs_adj2,paste0("b. Yin–Yang trajectories — ",adj2_label),
                     ifelse(layer=="protein","Mean protein z-score","Mean metabolite z-score"),smooth_lines=TRUE)
  save_plot(p3a/p3b+plot_layout(heights=c(1,1)),"c1.Fig3.yy_top.png",16,13.5,outdir=outdir)

  # Fig4 shows adjusted diagnosis-timed cross-sectional omics gradients.
  p4i<-make_gradient_panel(dat,top_inc,bvar,"incident",GRADIENT_STEP,MIN_BIN_N,
    paste0("a. Incident (Yin): top ",length(top_inc),if(le8_custom_adjustment())" adjusted biomarkers" else " behavioral-LE4 biomarkers"),covars=covs_adj2)
  p4p<-make_gradient_panel(dat,top_prev,bvar,"prevalent",GRADIENT_STEP,MIN_BIN_N,
    paste0("b. Baseline-prevalent (Yang): top ",length(top_prev),if(le8_custom_adjustment())" adjusted biomarkers" else " behavioral-LE4 biomarkers"),covars=covs_adj2)
  save_plot((p4i/p4p+plot_layout(guides="collect"))&theme(legend.position="right"),
            "c1.Fig4.gradient.png",17,9.6,outdir=outdir)
  gradient_rank<-bind_rows(
    assoc_adj2|>filter(term%in%top_inc)|>transmute(analysis="incident_behavioral_LE4",feature=term,p_value=p.value)|>arrange(p_value)|>mutate(rank=row_number()),
    prevalent_adj2|>filter(term%in%top_prev)|>transmute(analysis="prevalent_behavioral_LE4",feature=term,p_value=p.value)|>arrange(p_value)|>mutate(rank=row_number()))

  # Fig5: separately selected and clustered incident and prevalent adj2 sets.
  cluster_prev<-obj$clusters|>filter(analysis=="prevalent_behavioral_LE4")|>pull(feature)
  cli<-make_cluster_figure(dat,cluster_features,bvar,CLUSTER_STEP,YY_MAX_YEAR,"incident","a",
    saved_membership=filter(obj$clusters,analysis=="incident_behavioral_LE4"),
    saved_metrics=filter(obj$cluster_selection,analysis=="incident_behavioral_LE4"))
  clp<-make_cluster_figure(dat,cluster_prev,bvar,CLUSTER_STEP,YY_MAX_YEAR,"prevalent","b",
    saved_membership=filter(obj$clusters,analysis=="prevalent_behavioral_LE4"),
    saved_metrics=filter(obj$cluster_selection,analysis=="prevalent_behavioral_LE4"))
  # Keep the five top-level patchwork rows explicit.  Nesting the three-row
  # cluster_main object here is flattened by patchwork, so a three-entry outer
  # heights vector leaves room for only three of the resulting five panels.
  cluster_figure<-cli$cluster_plot/plot_spacer()/clp$cluster_plot/plot_spacer()/
    (cli$diagnostic_plot|clp$diagnostic_plot)+
    plot_layout(heights=c(1,.10,1,.05,.38))
  save_plot(cluster_figure,
            "c1.Fig5.gradient_cluster.png",18,10.5,outdir=outdir)
  cl_members<-bind_rows(cli$cluster|>mutate(analysis="incident_behavioral_LE4"),clp$cluster|>mutate(analysis="prevalent_behavioral_LE4"))
  cl_metrics<-bind_rows(cli$metrics|>mutate(analysis="incident_behavioral_LE4"),clp$metrics|>mutate(analysis="prevalent_behavioral_LE4"))

  # Fig6: adjusted Q5-vs-Q1 HRs on the same adj2 complete-case sample; no subtitle.
  qplots<-plot_quantile_top(dat,top6,tvar,evar,covs_adj2,paste0("Adjusted for ",paste(covs_adj2,collapse=", ")))
  save_plot(wrap_plots(qplots,ncol=3),"c1.Fig6.quantile_top.png",16,10,outdir=outdir)

  # Fig8: functional coherence of incident and baseline-prevalent adj2 scans.
  n_sig<-sum(is.finite(assoc_adj2$p.value)&assoc_adj2$p.value*nrow(assoc_adj2)<.05)
  n_sig_prev<-sum(is.finite(prevalent_adj2$p.value)&prevalent_adj2$p.value*nrow(prevalent_adj2)<.05)
  pe_i<-plot_functional_enrichment(enrich,n_sig,assoc_adj2,layer,if(le8_custom_adjustment())"a. Incident (Yin), selected covariates" else "a. Incident (Yin), behavioral LE4")
  pe_p<-plot_functional_enrichment(enrich_prev,n_sig_prev,prevalent_adj2,layer,if(le8_custom_adjustment())"b. Baseline-prevalent (Yang), selected covariates" else "b. Baseline-prevalent (Yang), behavioral LE4")
  if(layer != "protein") save_plot((pe_i/pe_p+plot_layout(guides="collect"))&theme(legend.position="right"),
            "c1.Fig8.enrich_sig.png",18,10.5,outdir=outdir)

  # Fig9 separates distal antecedent prediction from disease-state evidence.
  # GDF15 and PCSK9 are anchors by default, but the table is generated for all
  # assayed proteins/metabolites.
  directionality<-build_directionality_table(assoc_adj2,prevalent_adj2,duration_adj2,landmark_adj2,birthline_adj2,reverse_adj2)
  save_plot(plot_directionality_triage(directionality,C1_DIRECTION_ANCHORS),
            "c1.Fig9.directionality_triage.png",16,9,outdir=outdir)
  save_plot(plot_landmark_and_attained_age(landmark_adj2,assoc_adj2,birthline_adj2,C1_DIRECTION_ANCHORS),
            "c1.Fig10.landmark_birthline_sensitivity.png",16,9,outdir=outdir)
  save_plot(plot_risk_window_scan(risk_window_adj2,C1_DIRECTION_ANCHORS),
            "c1.Fig11.diagnosis_window_riskset.png",18,11,outdir=outdir)
  save_plot(plot_directionality_supplement(directionality,C1_DIRECTION_ANCHORS),
            "c1.Fig12.directionality_detail.png",14,7.5,outdir=outdir)
  save_plot(plot_reverse_time_exploratory(reverse_adj2,prevalent_adj2,duration_adj2,C1_DIRECTION_ANCHORS),
            "c1.Fig13.reverse_time_exploratory.png",15.5,8,outdir=outdir)
  save_plot(plot_pgs_actual_concordance(pgs_concordance,C1_DIRECTION_ANCHORS),
            "c1.Fig14.pgs_actual_concordance.png",18,10.5,outdir=outdir)
  if(layer=="metabolite")save_plot(vldl_deep$figure,
            "c1.Fig15.L_VLDL_TG_pct_deep_dive.png",19,26,outdir=outdir)

  le8_mock_c1(dat,assoc_adj2,enrich,layer,covs_adj2,tvar,evar,outdir,enrich_prev)
  gradient_rank<-obj$gradient_top10_provenance;cl_members<-obj$clusters;cl_metrics<-obj$cluster_selection
  le8_figure_workbook(list(cohort=cohort,input_feature_audit=input_feature_audit,association=assoc,prevalent=assoc_prevalent,incident_basic=assoc_basic,incident_adj2=assoc_adj2,
                   birthline_basic=birthline_basic,birthline_adj2=birthline_adj2,
                   prevalent_basic=prevalent_basic,prevalent_adj2=prevalent_adj2,
                   association_adj2_full_le8_sensitivity=assoc_adj2_full_le8,
                   prevalent_adj2_full_le8_sensitivity=prevalent_adj2_full_le8,
                   birthline_adj2_full_le8_sensitivity=birthline_adj2_full_le8,
                   attenuation_sameN=attenuation_sameN,
                   reverse_prevalent=reverse_adj2,prevalent_duration=duration_adj2,incident_landmark=landmark_adj2,
                   diagnosis_window_riskset=risk_window_adj2,directionality=directionality,
                   pgs_status=pgs$status,pgs_incident=pgs_incident,pgs_prevalent=pgs_prevalent,
                   pgs_attained_age=pgs_attained_age,pgs_incident_same_omic=pgs_incident_same,
                   pgs_prevalent_same_omic=pgs_prevalent_same,pgs_attained_same_omic=pgs_attained_age_same,
                   pgs_actual_concordance=pgs_concordance,
                   L_VLDL_TG_pct_deep_dive=vldl_deep$data,
                   L_VLDL_TG_pct_riskset_trajectory=vldl_deep$trajectory,
                   L_VLDL_TG_pct_pgs_conditional=vldl_deep$pgs_conditional,
                   L_VLDL_TG_pct_measured_conditional=vldl_deep$measured_conditional,
                   gradient_top10_provenance=gradient_rank,
                   cluster_membership=cl_members,cluster_selection=cl_metrics,
                   enrichment_incident=enrich,enrichment_prevalent=enrich_prev),"c1.out.xlsx",obj$review %||% list())
  sys.source(file.path(Sys.getenv("LE8_FDIR"),"c1_revision_temporal.R"),envir=.GlobalEnv)
  review<-obj$review %||% list()
  if(length(review)) le8_plot_c1_temporal(review$paired_associations %||% tibble(),
    review$time_heterogeneity %||% tibble(),review$paired_status %||% tibble(status="Unavailable in saved result"),
    features_all,layer)
}

le8_render_c2 <- function(obj,layer,outdir,rawdir) {
  aliases<-c(mr="MR",reverse_mr="MR_reverse",reverse_audit="MR_reverse_audit",assoc="observational",
    method_scope="method_scope",evidence_grades="evidence_grades",top_candidates="top_candidates",
    genetic_score_manifest="genetic_score_manifest",best="MR_best",availability="instrument_availability",
    top_r2="R2_QTL",arch="architecture",directionality_integration="directionality_integration",
    jobs_audit="MRLink2_job_audit",jobs="MRLink2_jobs")
  list2env(le8_figure_tables(obj,aliases),envir=environment())
  individual_decomposition<-obj$individual_decomposition %||% list()
  dandelion<-obj$DANDELION %||% list();mrlink2<-obj$MRLink2 %||% read_mrlink2_results(rawdir)
  fig1<-plot_c2_fig1(mr,assoc,layer);save_c2_plot(fig1$plot,"c2.Fig1.prots.top.png",17,13.5,outdir=outdir)
  save_c2_plot(plot_qtl_variance(mr,layer),"c2.Fig2.pQTL_R2.png",15.5,10.5,outdir=outdir)
  fig3<-plot_c2_fig2(mr,assoc,layer);save_plot(fig3$plot,"c2.Fig3.effect_concordance.png",16,12.5,outdir=outdir)
  write_raw_csv(fig3$wide,"c2.cis_trans_comparison.csv",le8_job_dir(outdir,"c2_cause"))
  le8_emit_restored_c2(mr,assoc,layer,outdir)
  save_plot(plot_c2_fig4(mr,layer),"c2.Fig4.sensitivity_architecture.png",15.5,11.5,outdir=outdir)
  plot_mrlink2_results(mrlink2,outdir)
  if(layer=="protein") {
    le8_dandelion_plot_bundle(dandelion,outdir)
    plot_dandelion_mr_integration(dandelion,mr,assoc,outdir)
  } else {
    suffix<-c("dandelion","dandelion_evidence","dandelion_mr_integration","dandelion_network")
    for(i in seq_along(suffix)) save_plot(blank_plot(paste0("C2 Figure ",i+5),
      "DANDELION is a gene/protein regulatory-network analysis and is not defined for metabolites"),
      paste0("c2.Fig",i+5,".",suffix[i],".png"),10,6,outdir=outdir)
  }
  save_plot(plot_c2_directionality(directionality_integration),"c2.Fig10.directionality_causal.png",16,12,outdir=outdir)
  save_plot(plot_bidirectional_mr(mr,reverse_mr,layer),"c2.Fig11.bidirectional_mr.png",16,8.5,outdir=outdir)
  save_plot(plot_individual_genetic_decomposition(individual_decomposition$summary %||% tibble()),"c2.Fig12.genetic_decomposition.png",17,13,outdir=outdir)
  save_plot(plot_component_leadtime(individual_decomposition$trajectory %||% tibble()),"c2.Fig13.genetic_component_leadtime.png",17,11,outdir=outdir)
  save_plot(plot_c2_evidence_grades(evidence_grades,layer),"c2.Fig14.evidence_grades.png",17,9,outdir=outdir)
  le8_figure_workbook(list(MR_all=mr,MR_reverse=reverse_mr,MR_reverse_audit=reverse_audit,
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
                   MRLink2_status=mrlink2$status%||%tibble()),"c2.out.xlsx",obj$review %||% list())
}

le8_render_c3 <- function(obj,layer,outdir,rawdir) {
  # The shell's completed-result route bypasses run_c3_layer; CIGMA must also
  # run here so adding a manifest works without deleting the expensive coloc cache.
  le8_c3_cigma(layer,outdir)
  plot_coloc_results(obj$summary,obj$regional,obj$variants,layer,outdir)
  gpu<-obj$GPU_coloc %||% read_gpu_coloc_results(rawdir)
  plot_gpu_coloc_validation(gpu,obj$summary,outdir)
  aud<-obj$credible_set_audit %||% credible_set_audit(obj$summary,obj$variants)
  tri<-obj$pgs_triangulation %||% read_c3_pgs_integration(layer,outdir,obj$summary)
  plot_c3_pgs_integration(tri,outdir)
  lists<-obj$causal_lists %||% list()
  le8_figure_workbook(list(coloc_summary=obj$summary,credible_set_audit=aud$overall,credible_set_by_locus=aud$by_locus,
    variant_posteriors=obj$variants,regional=obj$regional,GPU_results=gpu$results,GPU_status=gpu$status,
    GPU_manifest=obj$manifest,causal_sets=if(length(lists))stack(lists)else tibble(),pgs_observed_coloc=tri),"c3.out.xlsx",obj$review %||% list())
}

le8_render_c4 <- function(obj,layer,outdir,rawdir) {
  if(nrow(obj$modules$membership %||% tibble()) && nrow(obj$modules$metrics %||% tibble())) {
    obj$modules$metrics$selected_k <- n_distinct(obj$modules$membership$module)
    write_raw_csv(obj$modules$metrics,"c4.supervised_module_selection.csv",rawdir)
  }
  sets<-list(primary=obj$primary,membership=obj$membership,YS_edges=obj$YS_edges)
  if(identical(obj$meta$status,"unavailable")) {
    suffix<-c("proxy_heatmap","group_pillar_flow","connection_bridge","connection_evidence","mediation_forest",
      "mediation_diagnostics","supervised_atlas","network_globe","selection_mediation","state_network_remodeling","state_network_edges")
    if(layer=="metabolite")suffix[3:5]<-c("module_network","relationship_globe","mediation_wheel")
    reason<-paste(obj$status$detail %||% "Unavailable in saved result",collapse="; ")
    for(i in seq_along(suffix))save_plot(blank_plot(paste0("C4 Figure ",i),reason),paste0("c4.Fig",i,".",suffix[i],".png"),10,6,outdir=outdir)
  } else {
    plot_c4(obj$scan,sets,obj$mediation,obj$genetic_edges,layer,outdir)
    plot_supervised_atlas(obj$modules,sets,obj$mediation,layer,outdir)
    plot_module_globe(obj$modules,outdir)
    plot_c4_state_network(obj$state_network,layer,outdir)
  }
  tables<-le8_figure_tables(obj,c(primary_assignment="primary",proxy_membership="membership",YS_edges="YS_edges",
    supervised_modules="modules$membership",module_selection="modules$metrics",all_LE8_associations="scan",
    PRS_associations="genetic_scan",genetic_omic_bridges="genetic_edges",mediation="mediation",
    state_network_status="state_network$status",state_counts="state_network$state_counts",
    state_network_edges="state_network$edges",state_network_hubs="state_network$hubs"))
  if(!is.null(obj$status))tables$status<-obj$status
  # Reuse C4 core fits. The separate imaging cache checks its own input/option signature.
  imaging<-run_c4_imaging(layer,outdir=outdir)
  if(!is.null(imaging)) {
    associations<-imaging$associations %||% tibble()
    if(nrow(associations))plot_c4_imaging(associations,unique(associations$feature),outdir)
    tables<-c(tables,le8_figure_tables(list(imaging=imaging),c(imaging_status="imaging$status",imaging_associations="imaging$associations",imaging_fields="imaging$fields")))
  }
  le8_figure_workbook(tables,"c4.out.xlsx",obj$review %||% list())
}

le8_render_c5 <- function(obj,layer,outdir,rawdir) {
  le8_mock_c5(obj,outdir)
  all_glmnet_label<-if(layer=="protein")"Pradeep-style / glmnet"else"All-metabolite / glmnet"
  lightgbm_label<-if(layer=="protein")"Yu-style / LightGBM"else"MWAS-ranked / LightGBM"
  # The final model cache can predate prevalent-row augmentation. The exported
  # person_scores table is the input actually used for the completed figures.
  score_rows<-le8_figure_csv(rawdir,"c5.person_scores.csv")
  pred<-obj$prediction;psum<-obj$summary
  pairs_native<-obj$preclinical_pairs %||% list();boot_by_method<-obj$yy_auc_bootstrap %||% list()
  orders<-list(c(all_glmnet_label,lightgbm_label,"User specified"),
    c("Distal antecedent","Genetic-region evidence","Hybrid triangulated","Evidence-selected compact"),
    c("C4 NS","C4 YS","C4 YSplus"))
  titles<-c("Prediction paradigms","Evidence-screened prediction","C4 connection-guided prediction")
  files<-c("c5.Fig1.simple_pred.png","c5.Fig2.screened_pred.png","c5.Fig3.connection_pred.png")
  for(i in seq_along(orders))save_plot(c5_make_5row_grid(orders[[i]],pred,score_rows,pairs_native,boot_by_method,titles[i]),files[i],24,13.5,outdir=outdir)
  save_plot(plot_leadtime_prediction(obj$leadtime,obj$leadtime_windows),"c5.Fig4.leadtime_prediction.png",20,11.25,outdir=outdir)
  save_plot(plot_topn_mechanism(obj$score_vs_topN),"c5.Fig5.score_vs_topN_mechanism.png",19,10.7,outdir=outdir)
  save_plot(plot_evidence_matrix(obj$evidence)|plot_causal_reactive_map(obj$evidence),"c5.Fig6.evidence_matrix.png",18,10.5,outdir=outdir)
  save_plot(plot_performance_benchmark(psum),"c5.Fig7.performance_benchmark.png",13,9,outdir=outdir)
  save_plot(plot_incremental_performance(psum),"c5.Fig8.incremental_performance.png",13,9,outdir=outdir)
  save_plot(plot_complexity_performance(psum),"c5.Fig9.parsimony_performance.png",10.5,8,outdir=outdir)
  save_plot(plot_score_correlation(obj$score_correlations),"c5.Fig10.score_concordance.png",11,9.5,outdir=outdir)
  save_plot(plot_subgroup_auc(obj$subgroup_AUC),"c5.Fig11.subgroup_discrimination.png",12,10,outdir=outdir)
  save_plot(plot_attained_age_score_sensitivity(obj$attained_age_sensitivity),"c5.Fig12.attained_age_sensitivity.png",14,8.5,outdir=outdir)
  review<-obj$review %||% list()
  if(nrow(review$budget_metrics %||% tibble())) {
    sys.source(file.path(Sys.getenv("LE8_FDIR"),"c5_revision_validation.R"),envir=.GlobalEnv)
    le8_plot_c5_review(review$budget_metrics,review$calibration,review$decision_curves,
      unique(review$budget_metrics$horizon),unique(review$budget_metrics$budget),outdir)
  }
  tables<-le8_figure_tables(obj,c(upstream_availability="upstream_availability",prediction_summary="summary",
    input_set_sizes="set_sizes",candidate_sets="candidate_sets",mechanism_weight_prior="mechanism_weight_prior",
    distal_landmark_screen="distal_screen",score_coverage_audit="score_coverage_audit",sequential_forward="sequential$log",
    leadtime_discrimination="leadtime",leadtime_headline="leadtime_headline",leadtime_windows="leadtime_windows",
    score_vs_topN="score_vs_topN$summary",topN_individuals="score_vs_topN$individuals",topN_individual_range="score_vs_topN$individual_range",
    attained_age_sensitivity="attained_age_sensitivity",evidence_consolidation="evidence",score_correlations="score_correlations",
    subgroup_AUC="subgroup_AUC",nested_cv_summary="nested_cv$summary"))
  tables$training_screen<-le8_figure_csv(rawdir,"c5.training_association_screen.csv")
  tables$score_method_summary<-score_rows|>group_by(method,kind,split)|>summarise(rows=n(),finite_scores=sum(is.finite(score_z)),mean_score=mean(score_z,na.rm=TRUE),sd_score=sd(score_z,na.rm=TRUE),.groups="drop")
  tables$preclinical_pairs<-bind_rows(imap(pairs_native,~.x|>mutate(method=.y)))
  tables$yy_auc_bootstrap<-bind_rows(imap(boot_by_method,~.x|>mutate(method=.y)))
  le8_figure_workbook(tables,"c5.prediction_panels.xlsx",obj$review %||% list())
}
