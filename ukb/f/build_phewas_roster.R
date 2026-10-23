#!/usr/bin/env Rscript
# Derive the hospital ever-recorded population from verified complete exports.
# This audits extraction coverage, not incident follow-up or clinical absence.
args <- commandArgs(TRUE)
opt <- list("dir-in"=Sys.getenv("BMONTH_INDIR","/mnt/d/data/ukb"),
            "dir-out"="",replace=FALSE)
if(any(args %in% c("-h","--help"))) {
  cat("Usage: build_phewas_roster.R --dir-in DIR --dir-out DIR [--replace]\n",
      "Verifies raw participant -> ICD10/date exports -> code0/date0 RDS, including blank rows.\n",
      "Produces a derived roster, extraction evidence and coverage definition; no models are fitted.\n")
  quit(save="no")
}
i <- 1L
while(i<=length(args)) {
  if(args[i]=="--replace") {opt$replace<-TRUE;i<-i+1L;next}
  if(!args[i] %in% c("--dir-in","--dir-out") || i==length(args))stop("Unknown/incomplete argument: ",args[i])
  opt[[substring(args[i],3)]] <- args[i+1L];i<-i+2L
}
if(!nzchar(opt[["dir-out"]]))opt[["dir-out"]]<-file.path(
  Sys.getenv("BMONTH_OUTDIR","/mnt/d/analysis/bmonth"),"phewas")
script <- sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])[1]
fdir <- dirname(normalizePath(script,mustWork=TRUE))
source(file.path(fdir,"phewas_icd10.R"))
if(!requireNamespace("data.table",quietly=TRUE))stop("Install data.table")
indir <- normalizePath(opt[["dir-in"]],mustWork=TRUE)
outdir <- opt[["dir-out"]]
paths <- c(participant=file.path(indir,"rap/raw/pheno.tab.gz"),
  code_export=file.path(indir,"rap/icd10.code.tab.gz"),
  date_export=file.path(indir,"rap/icd10.date.tab.gz"),
  code_rds=file.path(indir,"Rdata/icd10.code0.rds"),
  date_rds=file.path(indir,"Rdata/icd10.date0.rds"),
  cohort=file.path(indir,"Rdata/all.rds"),
  withdrawals=file.path(indir,"rap/raw/withdrawal.csv"),
  extraction_script=file.path(dirname(fdir),"phe.sh"),rds_script=file.path(fdir,"phe.R"))
if(any(!file.exists(paths)))stop("Missing actual extraction inputs: ",paste(paths[!file.exists(paths)],collapse=", "))
target <- file.path(outdir,"phewas_observable_roster.rds")
if(file.exists(target)&&!opt$replace)stop("Roster exists; use --replace after reviewing extraction changes")
clean <- function(x) {x<-as.character(x);x[ukb_phewas_missing(x)]<-NA_character_;x}
read_export <- function(path) data.table::fread(path,colClasses="character",na.strings=c("","NA"),nThread=2L)
header <- function(path) {con<-gzfile(path,"rt");on.exit(close(con));strsplit(readLines(con,n=1L),"\t",fixed=TRUE)[[1]]}
code_names <- header(paths[["code_export"]]);date_names<-header(paths[["date_export"]])
raw_names <- header(paths[["participant"]])
if(code_names[1]!="eid" || date_names[1]!="eid" ||
   length(code_names)!=2L || !grepl("^p41270($|_)",code_names[2]) ||
   any(!grepl("^p41280($|_)",date_names[-1])))
  stop("Expected eid + the complete p41270 array and corresponding p41280 date fields")
sex_name <- intersect(c("p31_i0","p31"),raw_names)
if(length(sex_name)!=1L)stop("Cannot identify one UKB field 31 in the participant export")
selected <- unique(c("eid",sex_name,code_names[-1],date_names[-1]))
if(any(!selected %in% raw_names))stop("Required ICD10/date/sex columns absent from the raw participant export")

# cut selects columns only; the complete gzip stream is consumed and pipefail
# rejects corruption/truncation. No phenotype or participant filtering occurs.
indices <- sort(match(selected,raw_names))
pipeline <- paste("gzip -dc --",shQuote(paths[["participant"]]),"| cut -f",paste(indices,collapse=","))
command <- paste("bash -o pipefail -c",shQuote(pipeline))
message("[roster] Verify complete participant snapshot against extracted ICD10/date fields")
raw <- data.table::fread(cmd=command,colClasses="character",na.strings=c("","NA"),nThread=2L)
invisible(ukb_phewas_check_ids(as.data.frame(raw[,"eid",with=FALSE]),"eid","raw participant export"))
code <- read_export(paths[["code_export"]])
date <- read_export(paths[["date_export"]])
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
gaps <- list()
for(label in c("code","date")) {
  exported <- get(label)
  invisible(ukb_phewas_check_ids(as.data.frame(exported[,"eid",with=FALSE]),"eid",paste(label,"export")))
  absent<-setdiff(raw$eid,exported$eid);extra<-setdiff(exported$eid,raw$eid)
  if(length(absent)||length(extra)) {
    gaps[[label]]<-list(missing_participants=absent,unexpected_participants=extra)
    next
  }
  exported<-exported[match(raw$eid,exported$eid)]
  for(field in setdiff(names(exported),"eid")) {
    a<-clean(raw[[field]]);b<-clean(exported[[field]])
    mismatch<-which(xor(is.na(a),is.na(b)) | (!is.na(a)&!is.na(b)&a!=b))
    if(length(mismatch))gaps[[paste(label,field,sep="/")]]<-raw$eid[mismatch]
  }
}
if(length(gaps)) {
  ukb_phewas_atomic_rds(gaps,file.path(outdir,"extraction_gaps.rds"))
  stop("Raw-to-export mismatch; affected participants/fields saved in extraction_gaps.rds")
}
code<-code[match(raw$eid,code$eid)]
date<-date[match(raw$eid,date$eid)]
message("[roster] Raw snapshot matches both complete exports; verify derived RDS positions")
code0<-ukb_phewas_check_ids(readRDS(paths[["code_rds"]]),"eid","code0")
date0<-ukb_phewas_check_ids(readRDS(paths[["date_rds"]]),"eid","date0")
if(!setequal(code0$eid,raw$eid)||!setequal(date0$eid,raw$eid)||
   !identical(code0$eid,date0$eid))stop("RDS participant sets/order differ from complete raw export")
code0<-code0[match(raw$eid,code0$eid),,drop=FALSE]
date0<-date0[match(raw$eid,date0$eid),,drop=FALSE]
code_columns<-setdiff(names(code0),"eid")
if(!identical(code_columns,paste0("i",seq_along(code_columns))) ||
   !identical(names(date0),c("eid",paste0("d",seq_along(code_columns)))) ||
   length(date_names)-1L!=length(code_columns))
  stop("ICD10/date positional schema differs between text and RDS")
tokens<-data.table::tstrsplit(code[[code_names[2]]],"|",fixed=TRUE)
if(length(tokens)>length(code_columns))stop("RDS has truncated ICD10 array positions")
for(j in seq_along(code_columns)) {
  expected<-if(j<=length(tokens))clean(tokens[[j]])else rep(NA_character_,nrow(raw))
  if(!identical(clean(code0[[code_columns[j]]]),expected))stop("code0 mismatch at position ",j)
  if(!identical(clean(date0[[j+1L]]),clean(date[[date_names[j+1L]]])))stop("date0 mismatch at position ",j)
}
raw_ids<-raw$eid
sex<-clean(raw[[sex_name]])
if(any(!is.na(sex)&!sex%in%c("0","1")))stop("Unexpected field-31 sex coding")
sex<-ifelse(is.na(sex),NA_character_,ifelse(sex=="0","F","M"))
empty_n<-sum(ukb_phewas_missing(code[[code_names[2]]]))
rm(raw,code,date,code0,date0,tokens);invisible(gc())
message("[roster] All ICD10/date values agree; reconcile analysis IDs and withdrawals")
cohort<-ukb_phewas_check_ids(readRDS(paths[["cohort"]]),"eid","analysis cohort")
ids<-cohort$eid;rm(cohort);invisible(gc())
withdrawn<-trimws(readLines(paths[["withdrawals"]],warn=FALSE))
withdrawn<-unique(withdrawn[nzchar(withdrawn)])
if(any(!grepl("^[0-9]+$",withdrawn)))stop("Withdrawal file must be the verified one-column ID list without a header")
present<-ids%in%raw_ids;removed<-ids%in%withdrawn
roster<-data.frame(id=ids,phewas_observable=ifelse(removed,FALSE,ifelse(present,TRUE,NA)),
  sex=sex[match(ids,raw_ids)],coverage_source="UKB hospital ICD10 fields 41270/41280; verified local participant snapshot",
  in_participant_export=present,withdrawn=removed,
  coverage_reason=ifelse(removed,"withdrawal_list",ifelse(present,"complete_participant_export","absent_from_source_snapshot")),
  stringsAsFactors=FALSE)
gaps<-list(analysis_ids_absent_from_source=ids[!present],
           source_ids_absent_from_analysis=setdiff(raw_ids,ids),
           withdrawn_analysis_ids=ids[removed])
ukb_phewas_atomic_rds(gaps,file.path(outdir,"extraction_gaps.rds"))
message("[roster] Record source versions and completeness evidence")
file_info<-file.info(paths)
evidence<-list(contract_version="hospital-ever-recorded-v1",
  files=data.frame(role=names(paths),path=unname(normalizePath(paths)),
    bytes=file_info$size,mtime_utc=format(file_info$mtime,tz="UTC",usetz=TRUE),
    md5=unname(tools::md5sum(paths)),stringsAsFactors=FALSE),
  checks=c(raw_to_code_all_values=TRUE,raw_to_date_all_values=TRUE,
    code_export_to_rds_all_values=TRUE,date_export_to_rds_all_values=TRUE,
    preserves_blank_rows=TRUE,no_diagnosis_selection=TRUE),
  counts=c(analysis_n=length(ids),source_n=length(raw_ids),source_no_diagnoses=empty_n,
    analysis_absent_from_source=sum(!present),withdrawal_list_n=length(withdrawn),
    withdrawn_analysis_n=sum(removed),observable_n=sum(roster$phewas_observable,na.rm=TRUE)),
  extraction_rule="phe.sh prep_raw substitutes blank cells with NA; extract_multi performs cut on full participant rows for p41270/p41280; all source and derived values independently compared",
  snapshot_note="Local source hash/mtime identify this snapshot; filesystem dates are not clinical or regional censoring dates",
  incident_followup_validated=FALSE)
coverage<-paste(
  "Hospital ever-recorded PheWAS from the local UKB participant snapshot, fields 41270/41280.",
  "The full raw participant export was compared with both ICD10/date text extracts and every derived RDS array position; no diagnosis-based row restriction.",
  sprintf("Analysis population: %d; raw snapshot: %d; raw empty-diagnosis rows retained: %d.",
          length(ids),length(raw_ids),empty_n),
  sprintf("Explicit withdrawal list excludes %d analysis participants; %d analysis IDs absent from the raw snapshot remain unknown rather than controls.",
          sum(removed),sum(!present & !removed)),
  "TRUE means inclusion in this verified hospital-record snapshot under the complete participant-extract assumption, not proof of lifelong healthcare coverage.",
  "Zero means no mapped record of that phenotype in this source snapshot, subject to phenotype-specific exclusions and sex risk sets; mapping quality is audited separately.",
  "No individual or regional follow-up endpoints are asserted. Incident/Cox follow-up requires a separate validated source/version-specific contract.",
  "Source versions, hashes and extraction checks are stored in the roster extraction_evidence attribute and phewas_roster_evidence.rds.",
  sep="\n")
attr(roster,"coverage_definition")<-coverage
attr(roster,"extraction_evidence")<-evidence
ukb_phewas_atomic_rds(evidence,file.path(outdir,"phewas_roster_evidence.rds"))
writeLines(coverage,file.path(outdir,"phewas_coverage_definition.txt"))
ukb_phewas_atomic_rds(roster,target)
cat("Derived roster:",target,"\n")
print(evidence$counts)
