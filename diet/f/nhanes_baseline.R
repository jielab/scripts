# Recompute DOSE Supplementary Table 13 from the fitted NHANES primary sample.
args <- commandArgs(trailingOnly=TRUE)
root <- file.path(args[1], "main/dose", "nhanes")
bundle <- readRDS(file.path(root, "NHANES_analysis.rds"))
d <- bundle$primary
if (is.null(d) || !nrow(d)) stop("Missing two-day 2003-2018 primary NHANES sample")
options(survey.lonely.psu="adjust")
s <- survey::svydesign(ids=~psu_id, strata=~stratum_id, weights=~weight, nest=TRUE, data=d)
rows <- list(data.frame(variable="N",level="",n=nrow(d),mean=NA_real_,sd=NA_real_,percent=NA_real_))
for (v in c("dose","dash","energy_kJ","age","income","bmi")) {
    f <- reformulate(v)
    rows[[length(rows)+1L]] <- data.frame(variable=v, level="", n=sum(!is.na(d[[v]])),
        mean=as.numeric(coef(survey::svymean(f,s,na.rm=TRUE))),sd=sqrt(as.numeric(coef(survey::svyvar(f,s,na.rm=TRUE)))),percent=NA_real_)
}
for (v in c("sex","ethnicity_ygb","education_ygb","smoking","alcohol_status","pa_moderate","diabetes","hypertension","stroke_prevalent")) {
    if (!v %in% names(d)) stop("Missing baseline variable: ", v)
    for (lev in unique(as.character(na.omit(d[[v]])))) {
        s$variables$indicator <- as.numeric(as.character(d[[v]])==lev)
        rows[[length(rows)+1L]] <- data.frame(variable=v,level=lev,n=sum(s$variables$indicator,na.rm=TRUE),
            mean=NA_real_,sd=NA_real_,percent=100*as.numeric(coef(survey::svymean(~indicator,s,na.rm=TRUE))))
    }
}
write.csv(do.call(rbind,rows), file.path(root,"NHANES_baseline.csv"), row.names=FALSE,na="")
