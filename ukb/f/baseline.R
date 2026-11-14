# Baseline diabetes; future diagnoses are retained only in dm.ever.
ukb_baseline_diabetes <- function(data) {
    n <- nrow(data)
    number <- function(name) if (name %in% names(data)) suppressWarnings(as.numeric(as.character(data[[name]]))) else rep(NA_real_, n)
    date <- function(x) {
        if (inherits(x, "Date")) return(x)
        if (is.numeric(x)) return(as.Date(x, origin="1970-01-01"))
        suppressWarnings(as.Date(as.character(x), format="%Y-%m-%d"))
    }
    baseline <- date(data$date_attend)
    sentinels <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))
    diagnosed <- ever <- uncertain <- rep(FALSE, n)
    for (name in intersect(c("fod_srd_t2dm", "fod_icd10_t2dm", "fod_ref_t2dm"), names(data))) {
        d <- date(data[[name]])
        known <- !is.na(d) & !(d %in% sentinels)
        ever <- ever | !is.na(d)
        diagnosed <- diagnosed | (known & !is.na(baseline) & d <= baseline)
        uncertain <- uncertain | (!is.na(d) & !known)
    }
    positive <- rep(FALSE,n)
    for (name in unique(c(intersect(c("icd10Ct_t2dm", "srdCt_t2dm", "case_status_t2dm"), names(data)),grep("^cnt_.*t2dm$",names(data),value=TRUE)))) positive <- positive | (!is.na(number(name)) & number(name)>0)
    dated <- rep(FALSE,n)
    for (name in intersect(c("fod_srd_t2dm", "fod_icd10_t2dm", "fod_ref_t2dm"),names(data))) {
        d <- date(data[[name]]); dated <- dated | (!is.na(d) & !(d %in% sentinels))
    }
    uncertain <- uncertain | (positive & !dated) | is.na(baseline)
    medcols <- grep("^(drug[.]big3|(p)?(6153|6177)_i0(_a[0-9]+)?$)",names(data),value=TRUE)
    medcols <- medcols[!grepl("[_.]i[1-9]",medcols)]
    preserved <- grep("^[.]baseline_med_",names(data),value=TRUE)
    if (length(preserved)) medcols <- preserved
    med <- rep(NA_integer_,n)
    if (length(medcols)) {
        text <- do.call(paste,c(lapply(data[medcols],as.character),sep=" "))
        med <- vapply(text,function(z) {
            tokens <- suppressWarnings(as.integer(regmatches(z,gregexpr("-?[0-9]+",z))[[1L]]))
            if (3L %in% tokens) return(1L)
            known <- any(tokens %in% c(-7L,1:5)) && !any(tokens %in% c(-1L,-3L))
            if (known) 0L else NA_integer_
        },integer(1))
    }
    h <- number("bb_HBA1C")*.0915+2.15
    h[!is.finite(h) | number("bb_HBA1C")<=0] <- NA_real_
    lab <- !is.na(h) & h>=6.5
    drug <- !is.na(med) & med==1
    status <- ifelse(diagnosed | drug | lab,1L,ifelse(uncertain | is.na(med),NA_integer_,0L))
    data.frame(dm.baseline=status, dm.ever=as.integer(ever|positive|drug|lab),
        dm.baseline_source=ifelse(diagnosed,"dated_at_or_before_baseline",ifelse(drug,"baseline_medication",ifelse(lab,"baseline_hba1c",ifelse(is.na(status),"uncertain","no_baseline_evidence")))),
        drug.dm.baseline=med, dm.baseline_version="2026-11-11-v1")
}

# Fixed mapping/validity only: no cohort-dependent filtering or imputation.
ukb_save_raw_assay <- function(data, path, nonnegative=TRUE) {
    data <- as.data.frame(data); data$eid <- as.character(data$eid)
    if (anyDuplicated(data$eid) || anyNA(data$eid)) stop("Raw assay requires unique nonmissing eid")
    for (name in setdiff(names(data),"eid")) {
        x <- suppressWarnings(as.numeric(as.character(data[[name]])))
        x[!is.finite(x) | (nonnegative & x<0)] <- NA_real_
        data[[name]] <- x
    }
    lineage <- list(schema="ukb-named-raw-assay-v1", preprocessing="fixed_mapping_numeric_validity_only",
        imputation="none", cohort_filtering="none", signed_scale=!nonnegative,
        rows=nrow(data), columns=names(data), source_version="2026-11-11")
    attr(data,"preprocessing_lineage") <- lineage
    tmp <- paste0(path,".tmp"); saveRDS(data,tmp)
    if (!file.rename(tmp,path)) stop("Cannot publish raw assay: ",path)
    lineage$md5 <- unname(tools::md5sum(path))
    jsonlite::write_json(lineage,paste0(path,".lineage.json"),auto_unbox=TRUE,pretty=TRUE)
    invisible(data)
}
