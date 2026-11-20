# Assemble the publication report from matching cached phases.
clk<-load_cache("clock");asc<-load_cache("assoc");mark<-load_cache("marker");enr<-load_cache("enrichment")
clock_summary<-if(is.null(clk))list(notes="Clock has not been run for this phenotype/configuration.") else clk$data[setdiff(names(clk$data),c("predictions","folds"))]
if(!is.null(mark))clock_summary$markers<-mark$data$marker
assoc_summary<-if(is.null(asc))list() else asc$data[setdiff(names(asc$data),"settings")]
if(!is.null(asc)) {
  assoc_summary$global<-asc$data$assoc
  assoc_summary$pooled<-asc$data$pooled
}
payload<-list(title=if(cfg$demo)"CKM - synthetic data" else "CKM - continuous staging and stage-specific factors",
  version=cfg$version,generated=format(Sys.time(),"%Y-%m-%d %H:%M UTC",tz="UTC"),
  phase=if(is.null(clk)&&is.null(asc))"phenotype_only" else "analysis",
  input_notes=c(if(cfg$demo)"SYNTHETIC TEST DATA. All numbers and gene sets are simulated, not UK Biobank findings." else "UK Biobank observational research outputs.",
    paste("Phenotype:",ph_meta$definition,"/",ph_meta$stage3_method),
    paste("Follow-up end:",ph_meta$date_follow_end,"from 0f/phenotype.R::date_follow_end; verified earlier coverage/loss/death may shorten individual follow-up."),
    "Stage comparison reference model includes age spline, sex and categorical stage.",
    "Clinical CVD includes CHD, heart failure, stroke, peripheral arterial disease and atrial fibrillation.",
    "Clock accuracy must be judged from held-out metrics and calibration; no accuracy gain is guaranteed.",
    ph_meta$score2_scope),
  phenotype=ckm_phenotype_summary(d),clock=clock_summary,associations=assoc_summary,
  enrichment=if(is.null(enr))data.frame(status="not_run") else enr$data,
  validation=list(phenotype_schema=ph_meta$schema,staging_definition=ph_meta$definition,
    date_follow_end=ph_meta$date_follow_end,followup_end_source="0f/phenotype.R::date_follow_end",shared_helpers=cfg$shared,
    clinical_assay_source=ph_meta$clinical_assay_source,statin_source=ph_meta$statin_source,
    R=R.version.string,packages=as.list(stats::setNames(vapply(c("data.table","survival","jsonlite"),function(p)as.character(utils::packageVersion(p)),character(1)),c("data.table","survival","jsonlite"))),
    source_reference_commit="971a61b44b84500ab46e2f5592478505a0eae50e",
    missing_incident_domains=ph_meta$missing_incident_domains,
    temporal_scope=ph_meta$temporal_scope,stage_uncertainty=ph_meta$stage_uncertainty,
    clock_status=if(is.null(clk))"not_run_or_changed_inputs" else "out_of_fold_results_available"))
report_json<-file.path(cfg$out_dir,"tmp","report.json")
jsonlite::write_json(payload,report_json,auto_unbox=TRUE,dataframe="rows",null="null",na="null",digits=10)
python<-Sys.getenv("CKM_PYTHON","python3")
exit<-system2(python,c(shQuote(file.path(fdir,"report.py")),"--input",shQuote(report_json),"--output",shQuote(file.path(cfg$out_dir,"ckm.html"))))
if(exit!=0L)stop("HTML report generation failed")
