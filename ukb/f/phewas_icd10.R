# UK ICD-10 -> phecode panel with an explicitly observable population.
# This file has no side effects when sourced. The main entry point is
# ukb_rebuild_phewas(); ukb_build_phewas() also accepts synthetic mapping tables.

ukb_phewas_hash <- function(x) {
  f <- tempfile(); on.exit(unlink(f)); saveRDS(x, f, compress=FALSE)
  unname(tools::md5sum(f))
}
ukb_phewas_atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(".phewas-", tmpdir=dirname(path)); on.exit(unlink(tmp))
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
  invisible(path)
}
ukb_icd10_normalize <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z[is.na(x) | !nzchar(z)] <- NA_character_
  z <- gsub(".", "", z, fixed=TRUE)
  valid <- !is.na(z) & grepl("^[A-Z][0-9]{2}[A-Z0-9]{0,4}$", z)
  out <- rep(NA_character_, length(z))
  out[valid] <- ifelse(nchar(z[valid])>3L,
    paste0(substr(z[valid],1,3),".",substring(z[valid],4)),z[valid])
  out
}
ukb_phewas_date <- function(x) {
  if (inherits(x,"Date")) return(as.Date(x))
  if (inherits(x,"POSIXt")) return(as.Date(x))
  if (is.numeric(x))
    stop("Unclassed numeric diagnosis dates require explicit upstream conversion to Date")
  suppressWarnings(as.Date(as.character(x), format="%Y-%m-%d"))
}
ukb_phewas_missing <- function(x) {
  # The existing RAP text export and code0 parser preserve literal "NA".
  # It is a missing-value token, not a malformed ICD10 diagnosis.
  z <- trimws(as.character(x))
  is.na(x) | !nzchar(z) | toupper(z)=="NA"
}
ukb_phewas_check_ids <- function(x, id, label) {
  if (!is.data.frame(x) || !id %in% names(x)) stop(label, " requires ", id)
  x[[id]] <- as.character(x[[id]])
  if (anyNA(x[[id]]) || any(!nzchar(x[[id]])) || anyDuplicated(x[[id]]))
    stop(label, " IDs must be nonmissing and unique")
  x
}
ukb_phewas_roster <- function(roster, coverage_definition=NULL) {
  if (is.null(coverage_definition)) coverage_definition <- attr(roster,"coverage_definition")
  roster <- as.data.frame(roster)
  if ("eid" %in% names(roster) && !"id" %in% names(roster)) names(roster)[names(roster)=="eid"] <- "id"
  roster <- ukb_phewas_check_ids(roster,"id","PheWAS roster")
  req <- c("id","phewas_observable","sex")
  if (length(setdiff(req,names(roster))))
    stop("PheWAS roster requires id, phewas_observable (logical), sex (F/M/NA)")
  if (!is.logical(roster$phewas_observable))
    stop("phewas_observable must be logical; convert documented coverage flags explicitly")
  if (!is.character(coverage_definition) || length(coverage_definition)!=1L ||
      is.na(coverage_definition) || !nzchar(trimws(coverage_definition)))
    stop("Supply coverage_definition documenting the participant population, source snapshot and verified extraction completeness; incident follow-up dates are a separate contract")
  roster$sex <- toupper(trimws(as.character(roster$sex)))
  roster$sex[is.na(roster$sex) | roster$sex==""] <- NA_character_
  if (any(!is.na(roster$sex) & !roster$sex %in% c("F","M")))
    stop("Roster sex must be F/M/NA; explicitly map UKB field 31 as 0=F, 1=M")
  # Optional coverage windows are retained for future incident analyses. They do
  # not retroactively turn a lifetime first-diagnosis date into a visit history.
  for (nm in intersect(c("observation_start","observation_end"),names(roster)))
    roster[[nm]] <- ukb_phewas_date(roster[[nm]])
  if (all(c("observation_start","observation_end") %in% names(roster))) {
    bad <- !is.na(roster$observation_start) & !is.na(roster$observation_end) &
      roster$observation_start > roster$observation_end
    if (any(bad)) stop("Roster observation_start exceeds observation_end")
  }
  attr(roster,"coverage_definition") <- coverage_definition
  roster
}
ukb_phewas_maps <- function() {
  if (!requireNamespace("PheWAS",quietly=TRUE)) stop("Install the PheWAS package before rebuilding")
  list(vocabulary=as.data.frame(PheWAS::phecode_map_icd10),
       rollup=as.data.frame(PheWAS::phecode_rollup_map),
       exclusions=as.data.frame(PheWAS::phecode_exclude),
       sex=as.data.frame(PheWAS::sex_restriction),
       package_version=as.character(utils::packageVersion("PheWAS")),
       mapping_name="phecode_map_icd10",mapping_release="1.2b1")
}
ukb_phewas_validate_maps <- function(maps) {
  required <- list(vocabulary=c("vocabulary_id","code","phecode"),
    rollup=c("code","phecode_unrolled"),exclusions=c("code","exclusion_criteria"),
    sex=c("phecode","male_only","female_only"))
  for (nm in names(required)) {
    if (!is.data.frame(maps[[nm]]) || length(setdiff(required[[nm]],names(maps[[nm]]))))
      stop("Invalid PheWAS mapping table: ",nm)
    maps[[nm]] <- as.data.frame(maps[[nm]][,required[[nm]],drop=FALSE])
  }
  maps$source_table_md5 <- vapply(maps[c("vocabulary","rollup","exclusions","sex")],ukb_phewas_hash,character(1))
  v <- maps$vocabulary
  if (anyNA(v$vocabulary_id) || any(as.character(v$vocabulary_id)!="ICD10"))
    stop("Only the UK/WHO ICD10 vocabulary is allowed; ICD10CM is not interchangeable")
  v$code <- ukb_icd10_normalize(v$code); v$phecode <- as.character(v$phecode)
  if (anyNA(v$code)) stop("Invalid codes in ICD10 mapping table")
  # The official 1.2b1 table includes rows with an empty phecode: these are
  # explicitly unmapped ICD10 codes, not an empty-name phenotype.
  absent_phecode <- is.na(v$phecode) | !nzchar(trimws(v$phecode))
  maps$unmapped_dictionary_rows <- v[absent_phecode,,drop=FALSE]
  maps$vocabulary <- unique(v[!absent_phecode,,drop=FALSE])
  v <- maps$vocabulary
  for (nm in c("rollup","exclusions"))
    maps[[nm]][] <- lapply(maps[[nm]],as.character)
  missing_self <- setdiff(unique(v$phecode),
    maps$rollup$code[maps$rollup$code==maps$rollup$phecode_unrolled])
  if (length(missing_self)) stop("Rollup must include self mappings for every direct phecode")
  maps$sex$phecode <- as.character(maps$sex$phecode)
  if (!is.logical(maps$sex$male_only) || !is.logical(maps$sex$female_only) ||
      anyNA(maps$sex) || any(maps$sex$male_only & maps$sex$female_only) || anyDuplicated(maps$sex$phecode))
    stop("Invalid authoritative sex restriction table")
  maps
}
ukb_icd10_to_phecode_lookup <- function(codes, maps) {
  # Small dictionary expansion, equivalent to mapCodesToPhecodes with explicit
  # vocabulary.map=phecode_map_icd10 and rollup.map=phecode_rollup_map. Joining
  # the expanded dictionary to records avoids copying the full UKB long table
  # through the package's dplyr pipeline multiple times.
  v <- maps$vocabulary[,c("code","phecode"),drop=FALSE]
  z <- merge(data.frame(code=unique(codes[!is.na(codes)]),stringsAsFactors=FALSE),v,by="code")
  names(z)[names(z)=="code"] <- "source_code"
  z <- merge(z,maps$rollup,by.x="phecode",by.y="code")
  unique(data.frame(source_code=z$source_code,phecode=z$phecode_unrolled,stringsAsFactors=FALSE))
}
ukb_icd10_long <- function(code0,date0) {
  if (!requireNamespace("data.table",quietly=TRUE)) stop("Install data.table")
  code0 <- ukb_phewas_check_ids(as.data.frame(code0),"eid","ICD10 codes")
  date0 <- ukb_phewas_check_ids(as.data.frame(date0),"eid","ICD10 dates")
  if (!identical(code0$eid,date0$eid)) stop("ICD10 code and date IDs/order differ")
  cn <- setdiff(names(code0),"eid")
  if (!length(cn) || any(!grepl("^i[0-9]+$",cn))) stop("Code columns must be i1,i2,...")
  dn <- paste0("d",substring(cn,2))
  if (!all(dn %in% names(date0))) stop("Every ICD10 i-position requires its d-position")
  rows <- lapply(seq_along(cn),function(j) {
    raw <- as.character(code0[[cn[j]]]); keep <- !ukb_phewas_missing(raw)
    if (!any(keep)) return(NULL)
    raw_date <- date0[[dn[j]]][keep]
    parsed_date <- ukb_phewas_date(raw_date)
    data.table::data.table(id=code0$eid[keep],raw_code=raw[keep],
      source_code=ukb_icd10_normalize(raw[keep]),first_date=parsed_date,
      date_parse_failed=!ukb_phewas_missing(raw_date)&is.na(parsed_date))
  })
  if (!any(lengths(rows))) return(data.table::data.table(id=character(),raw_code=character(),
    source_code=character(),first_date=as.Date(character()),date_parse_failed=logical()))
  data.table::rbindlist(rows,use.names=TRUE)
}

ukb_build_phewas <- function(code0,date0,roster,maps=NULL,coverage_definition=NULL,
                            validate_asthma=TRUE,uncertain_phenotypes=NULL) {
  if (!requireNamespace("data.table",quietly=TRUE)) stop("Install data.table")
  official_maps <- is.null(maps)
  maps <- ukb_phewas_validate_maps(if (official_maps) ukb_phewas_maps() else maps)
  roster <- ukb_phewas_roster(roster,coverage_definition)
  long <- ukb_icd10_long(code0,date0)
  # An explicit roster may contain eligible people with zero inpatient diagnoses.
  # Its IDs must be in the complete ICD10 export, even when every code is blank.
  # This prevents an incomplete export from silently creating disease-free rows.
  absent <- setdiff(roster$id[!is.na(roster$phewas_observable)&roster$phewas_observable],as.character(code0$eid))
  if (length(absent)) stop(length(absent)," observable roster IDs absent from complete ICD10 export")
  long <- long[id %in% roster$id]
  lookup <- ukb_icd10_to_phecode_lookup(unique(long$source_code),maps)
  sentinel_codes <- c("J45.0","J45.1","J45.8","J45.9","J46")
  sentinel_map <- ukb_icd10_to_phecode_lookup(sentinel_codes,maps)
  sentinel <- data.frame(icd10=sentinel_codes,mapped_to_asthma=vapply(sentinel_codes,function(x)
    any(grepl("^495($|\\.)",sentinel_map$phecode[sentinel_map$source_code==x])),logical(1)))
  if (validate_asthma && !all(sentinel$mapped_to_asthma))
    stop("ICD10 asthma sentinel failed: J45.0/J45.1/J45.8/J45.9/J46 must map to asthma (495 family)")

  mapped_codes <- unique(lookup$source_code)
  universe_all <- sort(unique(maps$rollup$phecode_unrolled[maps$rollup$code %in% maps$vocabulary$phecode]))
  unsupported <- setdiff(universe_all,maps$sex$phecode)
  # Do not guess restrictions for newer/custom unmapped metadata codes (e.g.
  # 1000-series placeholders). Exclude those phenotypes and report them.
  universe <- setdiff(universe_all,unsupported)
  if (!length(universe)) stop("No phecodes have authoritative risk-set metadata")
  unsupported_phecodes <- data.frame(phecode=unsupported,
    reason=rep("authoritative_sex_metadata_unavailable",length(unsupported)))
  code_qc <- long[,list(n_records=.N,n_participants=data.table::uniqueN(id),
                       n_missing_date=sum(is.na(first_date)),n_date_parse_failed=sum(date_parse_failed)),
                    by=list(raw_code,source_code)]
  code_qc[,mapping_status:=ifelse(is.na(source_code),"invalid_code",
    ifelse(source_code %in% mapped_codes,"mapped","unmapped"))]
  data.table::setorder(code_qc,mapping_status,-n_records,raw_code)
  raw_counts <- long[,list(n_raw_records=.N,n_invalid_codes=sum(is.na(source_code))),by=id]
  mapped <- merge(long[,list(id,source_code,first_date,date_parse_failed)],
    data.table::as.data.table(lookup),by="source_code",allow.cartesian=TRUE,sort=FALSE)
  mapped <- unique(mapped[,list(id,phecode,first_date,date_parse_failed)])
  mapped_counts <- mapped[,list(n_mapped_phecodes=data.table::uniqueN(phecode),
    n_analyzable_phecodes=data.table::uniqueN(phecode[phecode %in% universe])),by=id]
  dates <- mapped[,list(first_date=if(all(is.na(first_date))) as.Date(NA) else min(first_date,na.rm=TRUE),
    has_missing_date=any(is.na(first_date)),has_date_parse_failure=any(date_parse_failed),
    n_distinct_first_dates=data.table::uniqueN(first_date[!is.na(first_date)])),by=list(id,phecode)]
  # This count is a count of distinct *first-diagnosis dates*, never visits.
  roster_qc <- merge(data.table::as.data.table(roster),raw_counts,by="id",all.x=TRUE,sort=FALSE)
  roster_qc <- merge(roster_qc,mapped_counts,by="id",all.x=TRUE,sort=FALSE)
  for (nm in c("n_raw_records","n_invalid_codes","n_mapped_phecodes","n_analyzable_phecodes"))
    data.table::set(roster_qc,which(is.na(roster_qc[[nm]])),nm,0L)
  # Observability describes data acquisition, not mapping success. A malformed
  # or unmapped diagnosis must not erase unrelated, definite diagnoses.
  roster_qc[,panel_observable:=!is.na(phewas_observable)&phewas_observable]
  roster_qc[,panel_exclusion_reason:=ifelse(is.na(phewas_observable)|!phewas_observable,
    "not_declared_observable","")]
  roster_qc[,mapping_audit:=ifelse(n_invalid_codes>0L,"invalid_codes_present",
    ifelse(n_raw_records>0L&n_analyzable_phecodes==0L,"no_analyzable_mapped_diagnoses",
      ifelse(n_raw_records==0L,"no_recorded_diagnoses","mapped")))]
  roster_qc <- roster_qc[match(roster$id,id)]
  if(is.null(uncertain_phenotypes))
    uncertain_phenotypes <- data.frame(id=character(),phecode=character(),reason=character())
  if(!is.data.frame(uncertain_phenotypes) ||
     !all(c("id","phecode","reason") %in% names(uncertain_phenotypes)))
    stop("Phenotype uncertainty requires id,phecode,reason columns")
  uncertain_phenotypes[] <- lapply(uncertain_phenotypes,as.character)
  if(anyNA(uncertain_phenotypes[,c("id","phecode","reason")]) ||
     any(!nzchar(trimws(uncertain_phenotypes$reason))) ||
     any(!uncertain_phenotypes$id %in% roster$id) ||
     any(!uncertain_phenotypes$phecode %in% universe))
    stop("Phenotype uncertainty must name roster participants, supported phecodes and nonempty reasons")
  uncertain_by_phecode <- split(uncertain_phenotypes$id,uncertain_phenotypes$phecode)

  # Include the map's prespecified universe, not only codes observed as cases.
  # Eligibility based on case/control counts is decided later by the model code.
  by_phecode <- split(mapped$id,mapped$phecode)
  excluded_sources <- split(maps$exclusions$exclusion_criteria,maps$exclusions$code)
  idx <- stats::setNames(seq_len(nrow(roster_qc)),roster_qc$id)
  panel <- data.frame(id=roster_qc$id,stringsAsFactors=FALSE)
  risk_counts <- vector("list",length(universe))
  for (j in seq_along(universe)) {
    phe <- universe[j]
    y <- rep.int(0L,nrow(roster_qc))
    y[!roster_qc$panel_observable] <- NA_integer_
    exclusion_ids <- unique(unlist(by_phecode[excluded_sources[[phe]]],use.names=FALSE))
    ix_exclusion <- unname(idx[exclusion_ids]); ix_exclusion <- ix_exclusion[!is.na(ix_exclusion)]
    y[ix_exclusion] <- NA_integer_
    ix_uncertain <- unname(idx[uncertain_by_phecode[[phe]]])
    y[ix_uncertain] <- NA_integer_
    case_ids <- unique(by_phecode[[phe]])
    ix_case <- unname(idx[case_ids]); ix_case <- ix_case[!is.na(ix_case)]
    ix_case <- ix_case[roster_qc$panel_observable[ix_case]]
    y[ix_case] <- 1L  # A definite case takes precedence over a control exclusion.
    rule <- maps$sex[match(phe,maps$sex$phecode),,drop=FALSE]
    sex_ineligible <- rep(FALSE,length(y))
    if (rule$male_only) sex_ineligible <- is.na(roster_qc$sex)|roster_qc$sex!="M"
    if (rule$female_only) sex_ineligible <- is.na(roster_qc$sex)|roster_qc$sex!="F"
    sex_conflicting_cases <- sum(y==1L & sex_ineligible,na.rm=TRUE)
    y[sex_ineligible] <- NA_integer_
    panel[[phe]] <- y
    risk_counts[[j]] <- data.frame(phecode=phe,n_cases=sum(y==1L,na.rm=TRUE),
      n_controls=sum(y==0L,na.rm=TRUE),n_excluded=sum(is.na(y)),
      n_reviewed_uncertain=length(unique(ix_uncertain)),
      n_sex_ineligible=sum(sex_ineligible),n_sex_conflicting_cases=sex_conflicting_cases,
      male_only=rule$male_only,female_only=rule$female_only)
  }
  qc_summary <- data.frame(metric=c("roster_n","declared_observable_n","panel_observable_n",
    "raw_diagnosis_records","mapped_diagnosis_records","unmapped_diagnosis_records",
    "invalid_diagnosis_records","unknown_diagnosis_date_records","people_without_analyzable_mapping",
    "mapped_record_fraction","unmapped_record_fraction","invalid_record_fraction","unsupported_phecodes",
    "literal_NA_code_cells"),
    value=c(nrow(roster_qc),sum(roster_qc$phewas_observable,na.rm=TRUE),sum(roster_qc$panel_observable),
      nrow(long),sum(long$source_code %in% mapped_codes),sum(!is.na(long$source_code)&!long$source_code %in% mapped_codes),
      sum(is.na(long$source_code)),sum(is.na(long$first_date)),
      sum(roster_qc$n_raw_records>0L & roster_qc$n_analyzable_phecodes==0L),
      if(nrow(long)) mean(long$source_code %in% mapped_codes) else NA_real_,
      if(nrow(long)) mean(!is.na(long$source_code)&!long$source_code %in% mapped_codes) else NA_real_,
      if(nrow(long)) mean(is.na(long$source_code)) else NA_real_,length(unsupported),
      sum(vapply(code0[setdiff(names(code0),"eid")],function(x)
        sum(toupper(trimws(as.character(x)))=="NA",na.rm=TRUE),numeric(1)))))
  list(panel=panel,dates=as.data.frame(dates),roster=as.data.frame(roster_qc),
    qc_summary=qc_summary,code_qc=as.data.frame(code_qc),risk_sets=do.call(rbind,risk_counts),
    sentinel=sentinel,unsupported_phecodes=unsupported_phecodes,uncertain_phenotypes=uncertain_phenotypes,
    maps=maps,official_maps=official_maps,
    coverage_definition=attr(roster,"coverage_definition"))
}

ukb_rebuild_phewas <- function(code_file,date_file,roster_file,output_dir,
                              coverage_definition=NULL,overwrite=FALSE,uncertainty_file=NULL) {
  if (!all(file.exists(c(code_file,date_file,roster_file))))
    stop("Code/date/explicit observable roster files must all exist")
  output_panel <- file.path(output_dir,"PheWAS.rds")
  output_manifest <- file.path(output_dir,"PheWAS.manifest.rds")
  if (file.exists(output_panel) && !isTRUE(overwrite))
    stop("PheWAS.rds already exists. Use a new PHEWAS_REBUILD_OUTDIR or explicitly set PHEWAS_REPLACE=TRUE after archiving old outputs")
  roster <- if (grepl("\\.rds$",roster_file,ignore.case=TRUE)) readRDS(roster_file) else {
    columns <- names(utils::read.csv(roster_file,nrows=0L))
    classes <- stats::setNames(rep(NA_character_,length(columns)),columns)
    classes[intersect(c("id","eid"),columns)] <- "character"
    utils::read.csv(roster_file,stringsAsFactors=FALSE,colClasses=classes)
  }
  result <- ukb_build_phewas(readRDS(code_file),readRDS(date_file),roster,
    coverage_definition=coverage_definition,uncertain_phenotypes=if(is.null(uncertainty_file))NULL else
      utils::read.csv(uncertainty_file,colClasses="character"))
  # Publish the validated manifest last. Readers verify both panel and sidecar
  # hashes, including a crash before replacing a panel whose values did not change.
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  ukb_phewas_atomic_rds(result$dates,file.path(output_dir,"PheWAS.dates.rds"))
  ukb_phewas_atomic_rds(result$roster,file.path(output_dir,"PheWAS.roster.rds"))
  ukb_phewas_atomic_rds(result$uncertain_phenotypes,file.path(output_dir,"PheWAS.uncertainty.rds"))
  for (nm in c("qc_summary","code_qc","risk_sets","sentinel","unsupported_phecodes"))
    utils::write.csv(result[[nm]],file.path(output_dir,paste0("PheWAS.",nm,".csv")),row.names=FALSE,na="")
  ukb_phewas_atomic_rds(result$panel,output_panel)
  map_hashes <- vapply(result$maps[c("vocabulary","rollup","exclusions","sex")],ukb_phewas_hash,character(1))
  manifest <- list(format_version=2L,vocabulary_id="ICD10",mapping_name="phecode_map_icd10",
    mapping_release=result$maps$mapping_release,package_version=result$maps$package_version,
    panel_md5=unname(tools::md5sum(output_panel)),roster_validated=TRUE,
    control_exclusions_applied=TRUE,sex_restrictions_applied=TRUE,
    created_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),
    input_files=c(code=normalizePath(code_file),date=normalizePath(date_file),roster=normalizePath(roster_file)),
    input_md5=tools::md5sum(c(code_file,date_file,roster_file)),mapping_md5=map_hashes,
    mapping_source_md5=result$maps$source_table_md5,
    mapping_empty_phecode_rows=nrow(result$maps$unmapped_dictionary_rows),
    unsupported_metadata_phecodes=result$unsupported_phecodes,
    sidecar_md5=tools::md5sum(file.path(output_dir,c("PheWAS.dates.rds","PheWAS.roster.rds","PheWAS.uncertainty.rds"))),
    coverage_definition=result$coverage_definition,
    roster_validation="Structural checks plus explicit supplied linkage/export eligibility; no inference from diagnosis presence",
    phenotype_definition="ever recorded ICD10-derived phecode, >=1 mapped source code; exclusions and sex risk sets applied",
    unknown_dates="Retained as ever-diagnosed cases; has_missing_date flags uncertainty for any incident/prevalence derivation",
    missing_tokens="R NA, blank cells and literal NA (case-insensitive) are missing, not invalid diagnosis codes",
    unmapped_policy="Mapping quality audited independently; missing/unmapped/invalid codes never change participant observability or erase definite cases",
    uncertainty_policy="Reviewed id/phecode/reason exclusions apply only to named phenotypes; definite cases override uncertainty and related-control exclusions, subject to coverage/sex risk sets",
    uncertainty_input=uncertainty_file,
    uncertainty_input_md5=if(is.null(uncertainty_file))NA_character_ else unname(tools::md5sum(uncertainty_file)),
    roster_evidence=attr(roster,"extraction_evidence"),
    first_date_count_note="n_distinct_first_dates is NOT independent encounters and MUST NOT be used for a two-visit phenotype",
    n_roster=nrow(result$panel),n_phenotypes=ncol(result$panel)-1L,
    builder_version="2026-10-20-v2-local1",
    builder_code_md5=ukb_phewas_hash(lapply(list(ukb_icd10_normalize,ukb_phewas_missing,ukb_phewas_roster,
      ukb_phewas_validate_maps,ukb_icd10_to_phecode_lookup,ukb_icd10_long,ukb_build_phewas,
      ukb_rebuild_phewas),function(f)paste(deparse(body(f)),collapse="\n"))),
    R=R.version.string,data_table_version=as.character(utils::packageVersion("data.table")))
  ukb_phewas_atomic_rds(manifest,output_manifest)
  invisible(manifest)
}
