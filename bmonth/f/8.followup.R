source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
engine <- file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R")
source(engine)
source(file.path(Sys.getenv("BMONTH_FDIR"),"0.followup_utils.R"))

run_followup <- function() {
  dat <- load_cohort()
  components <- csv_env("BMONTH_FOLLOWUP_COMPONENTS","targeted,cycle")
  if (length(setdiff(components,c("targeted","cycle"))))
    stop("BMONTH_FOLLOWUP_COMPONENTS must contain targeted and/or cycle")
  manifest <- list()
  record <- function(module,status,detail) {
    manifest[[length(manifest)+1L]] <<- data.frame(module=module,status=status,detail=detail,
      timestamp_utc=format(Sys.time(),tz="UTC",usetz=TRUE),stringsAsFactors=FALSE)
    bm_write_csv_atomic(bind_rows(manifest),file.path(rawdir,"followup_manifest.csv"))
  }
  save_table <- function(x,name) {
    if (!is.data.frame(x)) x <- bind_rows(x)
    bm_write_csv_atomic(x,file.path(rawdir,paste0(name,".csv")))
    bm_atomic_rds(x,file.path(rawdir,paste0(name,".rds")))
  }
  on.exit({ if (length(manifest)) bm_write_csv_atomic(bind_rows(manifest),file.path(rawdir,"followup_manifest.csv")) },add=TRUE)

  # Descriptive recruitment/availability audit. A missing derived trait can
  # mean an intentional exclusion (e.g. not a two-child family), not nonresponse.
  availability <- list()
  avail_traits <- intersect(c("cnt_older",EXPLORATORY_TRAITS),names(dat))
  avail_variables <- intersect(c("sex","center","assessment_year","ethnic.c","ethnicity"),names(dat))
  if ("age" %in% names(dat)) {
    dat$.age_band_followup <- cut(dat$age,breaks=c(-Inf,45,50,55,60,65,70,Inf),right=FALSE)
    avail_variables <- c(avail_variables,".age_band_followup")
  }
  for (trait in avail_traits) for (variable in avail_variables) {
    level <- as.character(dat[[variable]]); level[is.na(level)] <- "<missing>"
    for (lev in sort(unique(level))) {
      index <- level==lev; n <- sum(index); available <- sum(index & !is.na(dat[[trait]]))
      availability[[length(availability)+1L]] <- data.frame(exposure_or_source=trait,variable=variable,
        level=lev,N=n,available=available,not_available=n-available,proportion_available=available/n,
        note=if(trait=="cnt_older")"Recorded older-sibling count"else"Derived trait availability includes intentional eligibility restrictions",
        stringsAsFactors=FALSE)
    }
  }
  save_table(bind_rows(availability),"followup_birthorder_availability")

  if ("targeted" %in% components) {
    candidate_file <- Sys.getenv("BMONTH_FOLLOWUP_CANDIDATES",
      file.path(dirname(Sys.getenv("BMONTH_FDIR")),"config","followup_candidates.csv"))
    candidates <- read.csv(candidate_file,colClasses="character",check.names=FALSE)
    if (!all(c("phenotype","description") %in% names(candidates)) ||
        anyNA(candidates$phenotype) || anyDuplicated(candidates$phenotype))
      stop("Candidate file requires unique nonmissing phenotype codes and descriptions")
    candidates$selection_status <- "exploratory: locked after 2026-10-20 audit; not independent validation"
    bm_write_csv_atomic(candidates,file.path(rawdir,"followup_candidate_registry.csv"))

    input <- phewas_input_path(validate=FALSE)
    provenance <- phewas_provenance(input)
    panel <- phewas_packed_panel(provenance)
    take_y <- function(phenotype) {
      if(!phenotype %in% panel$phenotypes) return(NULL)
      bm_panel_column(panel,phenotype)
    }
    traits <- intersect(EXPLORATORY_ALL,EXPLORATORY_TRAITS)
    if (!length(traits)) record("targeted","skipped","Birth-order traits disabled by configuration")
    results <- strata <- stratum_counts <- adoption <- age_sensitivity <- extras <- selection <- comparisons <- list()
    extra_covars <- csv_env("BMONTH_FOLLOWUP_EXTRA_COVARS",Sys.getenv("BMONTH_EXTRA_COVARS",""))
    missing_extra <- setdiff(extra_covars,names(dat))
    if (length(missing_extra)) record("extra_covariates","skipped",
      paste("Requested variables missing; model not weakened:",paste(missing_extra,collapse=", ")))
    if (!length(extra_covars)) record("extra_covariates","not_requested",
      "Supply explicit BMONTH_FOLLOWUP_EXTRA_COVARS; no ancestry/smoking field is guessed")
    record("targeted","started",paste(length(candidates$phenotype),"locked candidates;",length(traits),"birth-order definitions"))
    for (phenotype in candidates$phenotype) {
      yy <- take_y(phenotype)
      if (is.null(yy)) { record(paste0("phecode_",phenotype),"skipped","Absent from the selected phenotype panel"); next }
      z <- inner_join(dat,yy,by="eid")
      for (exposure in traits) {
        family <- if (exposure=="border12") character() else "sibship_size_cat"
        d2 <- unique(c(covars_stage$D2_adultSES,exposure,"eid",".y",family))
        full <- fu_complete(z,d2)
        common <- fu_complete(full,covars_stage$D5_migration)
        for (adjust_family in c(FALSE,TRUE)) {
          fs <- if (adjust_family) family else character()
          r <- fu_scalar(full,exposure,phenotype,c(covars_stage$D2_adultSES,fs),
            if (adjust_family) "full_D2_family_size" else "full_D2", "full_D2",
            family_adjusted=adjust_family,
            extra_note=if (exposure=="border12") "Family size is fixed at two; no family-size coefficient is fitted" else "")
          results[[length(results)+1L]] <- r
          # All five stages use the exact same outcome/exposure/covariate-complete participants.
          for (stage in names(covars_stage)) {
            results[[length(results)+1L]] <- fu_scalar(common,exposure,phenotype,
              c(covars_stage[[stage]],fs),paste0("common_",stage,if(adjust_family)"_family_size"else""),
              "common_D5",family_adjusted=adjust_family)
          }
        }
        stratum_label <- fu_sibship_stratum(full$sibship_size)
        for (s in c("2","3","4","5","6+")) {
          ds <- full[which(stratum_label==s),,drop=FALSE]
          stratum_counts[[length(stratum_counts)+1L]] <- data.frame(phenotype=phenotype,exposure=exposure,
            sibship_stratum=s,N=nrow(ds),cases=sum(ds$.y==1),controls=sum(ds$.y==0),
            fit_planned=s!="6+" && (exposure!="border12" || s=="2"),
            note=if(s=="6+")"Descriptive only: pooled 6+ does not hold family size constant"else"Exact family size",stringsAsFactors=FALSE)
          if (s=="6+" || (exposure=="border12" && s!="2")) next
          strata[[length(strata)+1L]] <- fu_scalar(ds,exposure,phenotype,covars_stage$D2_adultSES,
            "exact_sibship_size",paste0("size_",s),family_adjusted=TRUE,stratum=s,
            extra_note="Exact-size stratum; overlapping birth-order definitions are correlated tests")
        }
        if ("adopted_num" %in% names(full)) {
          known <- full[which(!is.na(full$adopted_num) & full$adopted_num==0),,drop=FALSE]
          adoption[[length(adoption)+1L]] <- fu_scalar(known,exposure,phenotype,
            c(covars_stage$D2_adultSES,family),"known_not_adopted_D2_family_size","known_not_adopted",
            family_adjusted=TRUE,extra_note="Compared with full_D2_family_size: restriction changes sample, not the adjustment set")
        }
        age <- fu_age_basis(full)
        if (age$status=="ok") {
          age_data <- fu_complete(age$data,c("age",age$terms))
          age_data$.age_linear <- age_data[[age$source]]
          # Three matched rows distinguish changed age source from changed shape.
          a <- fu_scalar(age_data,exposure,phenotype,c(covars_stage$D2_adultSES,family),
            "age_sensitivity_matched_current_linear","age_complete",TRUE)
          a2 <- fu_scalar(age_data,exposure,phenotype,
            c(setdiff(covars_stage$D2_adultSES,"age"),".age_linear",family),
            "age_sensitivity_matched_source_linear","age_complete",TRUE)
          b <- fu_scalar(age_data,exposure,phenotype,
            c(setdiff(covars_stage$D2_adultSES,"age"),age$terms,family),
            "age_sensitivity_ns4","age_complete",TRUE)
          a$age_source <- "age"; a2$age_source <- b$age_source <- age$source
          a$age_source_note <- a2$age_source_note <- b$age_source_note <- if(age$source=="age_fractional")
            "Fractional age uses month-15 birth-date approximation; not exact DOB" else "Integer baseline-age spline sensitivity"
          age_sensitivity[[length(age_sensitivity)+1L]] <- bind_rows(a,a2,b)
        }
        if (length(extra_covars) && !length(missing_extra)) {
          extra_data <- fu_complete(full,extra_covars)
          a <- fu_scalar(extra_data,exposure,phenotype,c(covars_stage$D2_adultSES,family),
            "extra_covariates_matched_base","extra_complete",TRUE)
          b <- fu_scalar(extra_data,exposure,phenotype,c(covars_stage$D2_adultSES,family,extra_covars),
            "extra_covariates_adjusted","extra_complete",TRUE,
            extra_note="Explicit user-specified variables; attenuation is descriptive and does not establish mediation")
          extras[[length(extras)+1L]] <- bind_rows(a,b)
        }
      }
      message("[followup] completed targeted phenotype ",phenotype)
      # Publish completed candidates: interruption does not hide all finished work.
      save_table(bind_rows(results),"followup_targeted_partial")
    }
    result <- bind_rows(results)
    if (nrow(result)) {
      # One family per planned model specification, across all candidates AND
      # all four rank definitions; this remains exploratory after selection.
      result <- fu_adjust_family(result,c("analysis","family_adjusted"))
      groups <- split(result[result$sample_set=="common_D5",],
        interaction(result$phenotype[result$sample_set=="common_D5"],result$exposure[result$sample_set=="common_D5"],drop=TRUE))
      if (any(vapply(groups,function(g)length(unique(g$sample_key))!=1L,logical(1))))
        stop("Internal error: D1-D5 sample changed within a targeted contrast")
      # Absolute log-odds changes are retained even when a relative attenuation
      # ratio would be unstable or crosses zero. These are not causal fractions.
      pairs <- split(result,interaction(result$phenotype,result$exposure,drop=TRUE))
      compared <- lapply(pairs,fu_sample_adjustment)
      selection <- bind_rows(lapply(compared,`[[`,"summary"))
      comparisons <- bind_rows(lapply(compared,`[[`,"comparisons"))
    }
    strata <- bind_rows(strata); if(nrow(strata)) strata <- fu_adjust_family(strata,c("analysis","sibship_stratum"))
    save_table(result,"followup_targeted")
    save_table(selection,"followup_sample_vs_adjustment")
    save_table(comparisons,"followup_sample_vs_adjustment_comparisons")
    save_table(strata,"followup_exact_family_size")
    save_table(bind_rows(stratum_counts),"followup_family_size_counts")
    save_table(bind_rows(adoption),"followup_known_adoption")
    save_table(bind_rows(age_sensitivity),"followup_age_sensitivity")
    save_table(bind_rows(extras),"followup_extra_covariates")
    record("targeted","completed",paste(nrow(result),"model rows;",nrow(strata),"exact-size rows; failures retained"))
    record("known_adoption",if(length(adoption))"completed"else"skipped",
      if(length(adoption))"Explicit adopted_num==0 restriction; failed/empty models retained"else"No usable adopted_num field or no targeted contrasts")
    record("targeted_age_sensitivity",if(length(age_sensitivity))"completed"else"skipped",
      if(length(age_sensitivity))"Matched current-age linear, selected-age linear, and 4-df natural spline models"else"No age source with enough variation or no targeted contrasts")
    if (length(extra_covars) && !length(missing_extra)) record("extra_covariates","completed",paste(extra_covars,collapse=", "))
    record("sibship_size_6plus","descriptive_only","No pooled 6+ coefficient is presented as an exact-family-size contrast")
  } else record("targeted","not_requested","BMONTH_FOLLOWUP_COMPONENTS does not include targeted")

  if ("cycle" %in% components && length(OUTCOMES)) {
    summary <- curves <- list()
    for (outcome in OUTCOMES) {
      fit <- fu_joint_cox(dat,outcome,covars,"D4_joint_cycle")
      summary[[length(summary)+1L]] <- fit$summary; curves[[length(curves)+1L]] <- fit$curve
      age <- fu_age_basis(dat)
      if (age$status=="ok") {
        matched <- fu_complete(age$data,c("age",age$terms,covars))
        matched$.age_linear <- matched[[age$source]]
        for (model in c("D4_joint_cycle_age_matched_current_linear",
                        "D4_joint_cycle_age_matched_source_linear","D4_joint_cycle_age_ns4")) {
          cv <- if (grepl("ns4$",model)) c(setdiff(covars,"age"),age$terms) else
            if (grepl("source_linear$",model)) c(setdiff(covars,"age"),".age_linear") else covars
          s <- fu_joint_cox(matched,outcome,cv,model)
          s$summary$age_source <- if(grepl("current_linear$",model))"age"else age$source
          summary[[length(summary)+1L]] <- s$summary; curves[[length(curves)+1L]] <- s$curve
        }
      } else record(paste0("cycle_age_",outcome),"skipped",age$status)
      # Modest time sensitivities, each compared on the same participants.
      # Never include age + exact birth year + recruitment year simultaneously.
      time_terms <- list(
        birth_cohort_5y=if("birth_year_bm"%in%names(dat) &&
          !isTRUE(attr(dat,"bmonth_metadata")$birth_year_approx))
            factor(5*floor(dat$birth_year_bm/5)) else NULL,
        recruitment_year=if("assessment_year"%in%names(dat))
            factor(dat$assessment_year) else NULL)
      for(label in names(time_terms)) {
        term <- time_terms[[label]]
        if(is.null(term) || length(unique(term[!is.na(term)]))<2L) {
          record(paste0("cycle_",label,"_",outcome),"skipped","Recorded field absent or insufficient variation")
          next
        }
        dtime <- dat; dtime$.time_sensitivity <- term
        dtime <- fu_complete(dtime,c(covars,".time_sensitivity"))
        fits <- lapply(c(FALSE,TRUE),function(adjusted)
          fu_joint_cox(dtime,outcome,c(covars,if(adjusted)".time_sensitivity"),
            paste0("D4_joint_cycle_",label,if(adjusted)"_adjusted"else"_matched_base")))
        if(!identical(fits[[1]]$summary$sample_key,fits[[2]]$summary$sample_key))
          stop("Internal error: calendar sensitivity samples differ")
        for(s in fits) {
          s$summary$time_sensitivity <- label
          summary[[length(summary)+1L]] <- s$summary
          curves[[length(curves)+1L]] <- s$curve
        }
      }
    }
    s <- bind_rows(summary)
    if(nrow(s)) s <- fu_adjust_family(s,"analysis")
    save_table(s,"followup_cycle_cox")
    save_table(bind_rows(curves),"followup_cycle_monthly_curves")
    record("cycle","completed","Joint cos+sin coefficient covariance, 2-df Wald test, pointwise January contrasts, PH diagnostics; no splicing of scalar models")
  } else if(!length(OUTCOMES)) record("cycle","skipped_no_outcomes","Survival outcomes disabled; old cycle outputs are not current")
  else record("cycle","not_requested","BMONTH_FOLLOWUP_COMPONENTS does not include cycle")

  for (item in c("targeted_incident_phecodes","smoking_pathway","immune_cell_composition",
                 "school_cutoff","within_sibship_replication","weather_anomaly"))
    record(item,"requires_additional_design_or_inputs",
      "See 0.README.md; not executed or represented as a confirmed finding by this module")
  record("run","completed","Completed means output writing succeeded; inspect per-model status and valid/eligible counts")
}
run_followup()
