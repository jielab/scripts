# Core preflight deliberately excludes omics and HadUK-Grid inputs.
required <- c("dplyr","tidyr","tibble","survival","PheWAS","brglm2",
              "data.table","ggplot2","ggrepel","digest","patchwork","writexl")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing core R packages: ",paste(missing,collapse=", "))
source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
files <- file.path(indir,"Rdata",c("all.rds","move.rds"))
if(any(!file.exists(files))) stop("Missing core input files: ",paste(files[!file.exists(files)],collapse=", "))
panel <- phewas_input_path(validate=FALSE)
provenance <- phewas_provenance(panel)
cat("PheWAS panel:",panel,"; MD5:",provenance$panel_md5,"\n")
cat("PheWAS source:",provenance$validation_status,";",provenance$source_note,"\n")
if(length(provenance$issues))cat("Optional provenance not used:",paste(provenance$issues,collapse="; "),"\n")
cat("PheWAS structural check:",provenance$n_participants,"participants;",
    length(provenance$phenotype_counts),"binary phenotype columns\n")

dat <- read_all() %>% merge_supplemental_covariates() %>% merge_move(required=TRUE) %>% make_birth_month()
bys <- attr(dat,"birth_year_source"); if(is.null(bys)) bys <- "unknown"
dat <- add_retained_covariates(prep_covars(dat)) %>% make_birth_order()
requested <- unique(c(csv_env("BMONTH_EXTRA_COVARS"),csv_env("BMONTH_FOLLOWUP_EXTRA_COVARS"),
                      csv_env("BMONTH_OMICS_COVARIATES")))
absent <- setdiff(requested,names(dat))
if(length(absent)) stop("Requested supplemental/retained covariates are missing: ",paste(absent,collapse=", "))
cat("Input cohort:",nrow(dat),"\n")
cat("tdi.birth nonmissing:",sum(!is.na(dat$tdi.birth)),"\n")
cat("move nonmissing:",sum(!is.na(dat$move_km)),"\n")
cat("birth-year source:",bys,"\n")
if(bm_truth(Sys.getenv("BMONTH_CHECK_SURVIVAL","FALSE"))) {
  for(y in OUTCOMES) {
    dat <- ensure_outcome(dat,y)
    cat("Cox input:",y,"; valid time/status rows:",
        sum(complete.cases(dat[,paste0(y,c(".t2e",".Yt2e"))])),"\n")
  }
} else cat("PheWAS preflight: survival inputs are checked by the requested Cox stages.\n")
for(x in c(BASE_EXPOSURES,EXPLORATORY_TRAITS))cat(x,sum(!is.na(dat[[x]])),"nonmissing\n")
print(birth_order_qc(dat),n=Inf)
versions <- data.frame(package=required,
  version=vapply(required,function(p)as.character(packageVersion(p)),character(1)),stringsAsFactors=FALSE)
bm_write_csv_atomic(versions,file.path(rawdir,"package_versions.csv"))
writeLines(c(paste("Code",BMONTH_CODE_VERSION),capture.output(sessionInfo())),file.path(rawdir,"session_info.txt"))
cat("Core preflight passed. Omics and temperature inputs are checked only by their requested modules.\n")
