#!/usr/bin/env Rscript
# Synthetic publication contracts: score filtering must not destroy raw inputs.
this <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1L])
root <- normalizePath(file.path(dirname(this),"../.."))
source(file.path(root,"f/common.R"))
check <- function() {
    tmp <- tempfile("diet-disco-"); dir.create(tmp)
    on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
    a <- data.frame(eid=c("001","002"),diet.dash.sum=c(20,30),
                    diet.upf.sum=c(NA_real_,NA_real_),energy.kJ=c(8000,9000),
                    diet.dash.component=c(2,3))
    b <- data.frame(eid=c("002","003"),diet.dash.sum=c(10,15),
                    diet.dose.sum=c(4,5),dose=c(4,5),dose_component=c(1,2))
    attr(a,"recipe") <- "maha"; attr(b,"recipe") <- "dose"
    paths <- c(maha=file.path(tmp,"maha.rds"),dose=file.path(tmp,"dose.rds"))
    saveRDS(a,paths[[1]]); saveRDS(b,paths[[2]])
    log <- capture.output(x <- diet_publish_disco("ukb",tmp,paths,"dash,upf,chdi"),type="message")
    stopifnot(identical(names(x),c("id","eid","diet.dash.sum","dose.diet.dash.sum")),
              any(grepl("SKIP ukb: chdi",log,fixed=TRUE)),
              any(grepl("SKIP ukb: diet.upf.sum",log,fixed=TRUE)))
    raw <- readRDS(file.path(tmp,"disco/ukb.raw.rds"))
    stopifnot(all(c("energy.kJ","diet.dash.component","dose_component","dose","diet.dose.sum") %in% names(raw)),
              !any(c("energy.kJ","dose_component","diet.dose.sum") %in% names(x)))
    for(profile in names(paths)) {
        original <- readRDS(paths[[profile]])
        restored <- diet_profile(x,profile,"ukb")
        stopifnot(identical(restored$eid,original$eid),
                  identical(attr(restored,"recipe"),attr(original,"recipe")))
        for(nm in names(original)) stopifnot(identical(restored[[nm]],original[[nm]]))
    }
    catalog <- read.delim(file.path(tmp,"disco/ukb.scores.tsv"))
    stopifnot(setequal(catalog$variable,c("diet.dash.sum","dose.diet.dash.sum")))
    # Selection can change while reusing the same completed source parts.
    x <- diet_publish_disco("ukb",tmp,paths,"dose")
    stopifnot(identical(names(x),c("id","eid","diet.dose.sum")))
    # No hidden fallback to stale source caches if the paired raw file changes.
    raw$energy.kJ[match("001",raw$id)] <- 1; saveRDS(raw,file.path(tmp,"disco/ukb.raw.rds"))
    stopifnot(inherits(try(diet_profile(x,"dose","ukb"),silent=TRUE),"try-error"))
    # Existing schema-1 full tables remain readable.
    stopifnot(identical(diet_profile(raw,"maha","ukb")$energy.kJ,c(1,9000)))
}
check()
cat("PASS selected totals, unavailable-score logs, full raw preservation, profile reconstruction and stale-raw rejection\n")
