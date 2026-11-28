source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
#!/usr/bin/env Rscript
# FLT MRI: observed outcomes, covs.all plus ICV for volumes/WMH, separate test families.
# Defaults are absolute and do not depend on the working directory.
# Rewrites outputs in analysis/01lm/mri; archive old results manually if needed.
options(stringsAsFactors=FALSE)
PROJECT <- Sys.getenv("FLT_OUT")
YGB <- Sys.getenv("FLT_YGB")
RAW <- Sys.getenv("FLT_MRI_RAW", unset=file.path(YGB,"rap/img.tab.gz"))
OUT <- file.path(PROJECT,"01lm/mri")
# Default reads original MRI table. FLT_MRI_RAW may select a verified equivalent TSV.
DICT <- Sys.getenv("FLT_MRI_DICTIONARY")
MODEL <- "covs.all"
MIN_N <- 100L
# Match s1_3_otherphe_lm.R: heteroskedasticity-robust HC1 SE.
if (!requireNamespace("data.table", quietly=TRUE)) stop("Install package: data.table")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write_tab <- function(x,n) write.csv(x,file.path(OUT,n),row.names=FALSE,na="")
ids <- function(x,label) {
    x <- as.data.frame(x)
    if (!"eid" %in% names(x) || anyDuplicated(names(x))) stop("Invalid columns: ",label)
    x$eid <- as.character(x$eid)
    if (anyNA(x$eid)||any(!nzchar(x$eid))||anyDuplicated(x$eid)) stop("Invalid/duplicate eid: ",label)
    x
}
numeric_clean <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    z <- suppressWarnings(as.numeric(x)); z[!is.finite(z)] <- NA_real_; z
}
standardize <- function(x) {
    s <- sd(x); if(!is.finite(s)||s<=0) stop("Zero or invalid SD")
    as.numeric(scale(x, center=TRUE, scale=TRUE))
}
complete <- function(d) {
    ok <- complete.cases(d)
    for(n in names(d)) if(is.numeric(d[[n]])) ok <- ok & is.finite(d[[n]])
    ok
}
D <- read.csv(DICT,check.names=FALSE)
required <- c("field_id","family","category","label","side","atlas_key","transform","adjust_icv")
if(!all(required %in% names(D))) stop("Missing dictionary columns")
D$field_id <- as.character(D$field_id)
if(anyNA(D$adjust_icv) || any(!tolower(as.character(D$adjust_icv)) %in% c("true","false"))) stop("Invalid adjust_icv flag")
D$dictionary_adjust_icv <- D$adjust_icv
D$adjust_icv <- D$family %in% c("Volume139","WMH")
expected <- c(FA27=27L,MD27=27L,ThicknessDK=66L,Volume139=139L,WMH=1L)
stopifnot(!anyNA(D$field_id),!anyDuplicated(D$field_id),identical(table(factor(D$family,levels=names(expected))) |> as.integer(),unname(expected)),nrow(D)==260L)
# Latest agreed analysis: all outcomes directly standardized, including WMH.
# Preserve the supplied dictionary setting rather than silently editing its source.
D$dictionary_transform <- D$transform
D$transform <- "none"
write_tab(D,"00_dictionary_effective.csv")
cfg <- readRDS(file.path(Sys.getenv("FLT_COHORT"),"analysis_config.rds"))
dat <- ids(readRDS(file.path(Sys.getenv("FLT_COHORT"),cfg$primary_data_file)),"cohort")
liking_selected <- c(`Butter on bread`="p20628", `Total salad vegetable`="F-Salad_vegetables", Chicken="p20638")
exposures <- unname(liking_selected)
covs.all <- c("age", "sex", "bmi", "edu", "inc", "tdi", "smoke_status", "alcohol_status", "pa_met_tot")
covs <- covs.all
if (!all(c("eid",exposures,covs.all) %in% names(dat))) stop("Missing FLT/cohort covariates")
if (!all(vapply(dat[exposures],is.numeric,logical(1)))) stop("Non-numeric FLT")
# Same primary cohort as s1_3; no disease-specific eligibility restriction.
write_tab(data.frame(variable=covs.all,retained=covs.all %in% covs),"00_inherited_covariate_audit.csv")

# Read only explicitly selected instance-2 fields; never load all modalities for KNN.
if(!file.exists(RAW)) stop("Raw MRI table not found: ",RAW)
con <- if(grepl("\\.gz$",RAW)) gzfile(RAW,"rt") else file(RAW,"rt")
header <- strsplit(readLines(con,n=1L),"\t",fixed=TRUE)[[1]]; close(con)
D$column <- paste0("p",D$field_id,"_i2")
D$available <- D$column %in% header
write_tab(D,"00_MRI_field_mapping.csv")
needed <- unique(c(D$column,"p26521_i2"))
img <- ids(data.table::fread(RAW,select=c("eid",intersect(needed,header))),"raw MRI")
for(n in setdiff(names(img),"eid")) img[[n]] <- numeric_clean(img[[n]])
# Absent fields remain planned tests with failed status; none silently dropped.
for(n in setdiff(D$column,names(img))) img[[n]] <- NA_real_
if(!"p26521_i2" %in% names(img)) img$p26521_i2 <- NA_real_

# Match participants by eid; use exactly the primary-cohort covs.all.
joined <- dat[dat$eid %in% img$eid,c("eid",exposures,covs),drop=FALSE]
ii <- match(joined$eid,img$eid)
if(!nrow(joined)) stop("No cohort / MRI overlap")
base <- data.frame(eid=joined$eid)
cn <- paste0("cov",seq_along(covs))
for(j in seq_along(covs)) base[[cn[j]]] <- joined[[covs[j]]]
base$icv <- img$p26521_i2[ii]
base$icv[base$icv<=0] <- NA_real_
write_tab(data.frame(variable=c(covs,"ICV (volume/WMH only)"),
                    source=c(rep("FLT primary cohort",length(covs)),"p26521_i2")),"00_covariates_used.csv")
chronology <- "not verified: FLT/MRI dates not assessed"
qc <- data.frame(field_id=D$field_id,family=D$family,available=D$available,
                 N_overlap=nrow(joined),N_observed=vapply(D$column,function(n)sum(is.finite(img[[n]][ii])),integer(1)))
write_tab(qc,"00_observed_MRI_counts.csv")
# Participants with no observed outcome cannot enter any of these models.
# This is NOT a requirement to have all modalities: one observed field suffices.
has_any <- Reduce(`|`,lapply(D$column,function(n)is.finite(img[[n]][ii])))
write_tab(data.frame(overlap=nrow(joined),any_observed=sum(has_any),none_observed=sum(!has_any)),"00_analysis_eligibility.csv")
joined <- joined[has_any,,drop=FALSE];base <- base[has_any,,drop=FALSE]
img <- img[ii[has_any],unique(c("eid",D$column)),drop=FALSE];ii <- seq_len(nrow(img))
rm(dat);invisible(gc())

run_pair <- function(e,k) {
    if(!D$available[k]) stop("Field absent from raw table")
    terms <- cn
    bcols <- c(cn,if(D$adjust_icv[k]) "icv")
    d <- base[bcols];d$x <- joined[[e]];d$y <- img[[D$column[k]]][ii]
    d$y[d$y<0] <- NA_real_
    if(D$family[k]=="ThicknessDK") d$y[d$y==0] <- NA_real_
    if(D$family[k]=="FA27") d$y[d$y>1] <- NA_real_
    # All MRI outcomes, including WMH: direct scale(), no log transformation.
    d <- droplevels(d[complete(d),,drop=FALSE])
    if(nrow(d)<MIN_N) stop("Insufficient observed complete cases: ",nrow(d))
    xsd <- sd(d$x);ysd <- sd(d$y)
    d$xz <- standardize(d$x);d$yz <- standardize(d$y)
    constant <- cn[vapply(d[cn],function(x)length(unique(x))<2,logical(1))]
    if(length(constant)) stop("Constant inherited covariate: ",paste(covs[match(constant,cn)],collapse=", "))
    if(D$adjust_icv[k]) {d$icv <- standardize(d$icv);terms <- c(terms,"icv")}
    fit <- lm(reformulate(c("xz",terms),response="yz"),data=d,na.action=na.fail,x=TRUE)
    if(fit$rank<length(coef(fit))) stop("Rank-deficient model")
    beta <- unname(coef(fit)["xz"])
    df <- fit$df.residual
    bread <- summary(fit)$cov.unscaled
    influence <- drop(fit$x[,rownames(bread),drop=FALSE] %*% bread[,"xz"])
    se <- sqrt(nrow(d)/df * sum((residuals(fit)*influence)^2))
    if(!is.finite(beta)||!is.finite(se)||se<=0||df<=0) stop("Invalid coefficient or SE")
    ci <- beta+c(-1,1)*qt(.975,df)*se
    c(N=nrow(d),beta=beta,SE=se,LCI=ci[1],HCI=ci[2],P=2*pt(-abs(beta/se),df),
      exposure_SD=xsd,MRI_SD=ysd,beta_per_score_point=beta/xsd)
}
rows <- list();n <- 0L
for(e in exposures) {
    message("START MRI ",e,"; planned=",nrow(D))
    for(k in seq_len(nrow(D))) {
        err <- "";warn <- character()
        z <- tryCatch(withCallingHandlers(run_pair(e,k),warning=function(w){warn <<- c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(x){err <<- conditionMessage(x);NULL})
        if(is.null(z)) z <- setNames(rep(NA_real_,9),c("N","beta","SE","LCI","HCI","P","exposure_SD","MRI_SD","beta_per_score_point"))
        n <- n+1L;rows[[n]] <- cbind(data.frame(exposure=e,food_liking=names(liking_selected)[match(e,liking_selected)],exposure_type="FLT",model=MODEL),D[k,],as.data.frame(as.list(z)),data.frame(status=if(nzchar(err))"failed" else "ok",error=err,warnings=paste(unique(warn),collapse=" | ")))
    }
    message("DONE MRI ",e)
    # Checkpoint only completed exposure blocks, with explicit partial status.
    saveRDS(list(completed_exposure=e,partial_results=do.call(rbind,rows)),file.path(OUT,"checkpoint_partial.rds"))
}
r <- do.call(rbind,rows);rownames(r)<-NULL
r$FDR_family <- r$FDR_exposure_all <- r$FDR_FLTs_family <- NA_real_
for(e in exposures) {
    i <- which(r$exposure==e)
    r$FDR_exposure_all[i] <- p.adjust(r$P[i],"BH",n=nrow(D))
    for(f in unique(D$family)) {
        j <- i[r$family[i]==f]
        r$FDR_family[j] <- p.adjust(r$P[j],"BH",n=sum(D$family==f))
    }
}
for(f in unique(D$family)) {
    j <- which(r$family==f & r$exposure_type=="FLT")
    r$FDR_FLTs_family[j] <- p.adjust(r$P[j],"BH",n=sum(D$family==f)*length(exposures))
}
r$FDR_all_tests <- p.adjust(r$P,"BH",n=nrow(D)*length(exposures))
# Primary beta retains the s1_3 interpretation: per one point higher FLT.
# MRI outcomes are z-scaled within the complete-case sample for each model.
r$beta_standardized <- r$beta
r$SE_standardized <- r$SE
r$LCI_standardized <- r$LCI
r$HCI_standardized <- r$HCI
r$beta <- r$beta_standardized/r$exposure_SD
r$SE <- r$SE_standardized/r$exposure_SD
r$LCI <- r$LCI_standardized/r$exposure_SD
r$HCI <- r$HCI_standardized/r$exposure_SD
# Also retain s1_3's within-phenotype correction (3 FLTs), as a secondary scope.
r$FDR_trait <- NA_real_
for (id in D$field_id) {
    j <- which(r$field_id==id)
    r$FDR_trait[j] <- p.adjust(r$P[j], "BH", n=length(exposures))
}
r$significant <- !is.na(r$FDR_family) & r$FDR_family<0.05
write_tab(r,"mri_lm_results.csv")
write_tab(r[r$status!="ok" | nzchar(r$warnings),],"05_errors_warnings.csv")
counts <- do.call(rbind,lapply(split(r,interaction(r$exposure,r$family,drop=TRUE)),function(q)
    data.frame(exposure=q$exposure[1],family=q$family[1],planned=nrow(q),completed=sum(q$status=="ok"),P05=sum(q$P<.05,na.rm=TRUE),FDR05=sum(q$significant))))
write_tab(counts,"06_summary_counts.csv")
saveRDS(list(results=r,dictionary=D,cohort_config=cfg,liking_selected=liking_selected,
    covariates=covs.all,standard_error="HC1",
    primary_FDR="BH within FLT and MRI family",chronology=chronology),file.path(OUT,"mri_lm_results.rds"))
writeLines(c(
 "Three FLTs fitted separately; primary cohort from s0, no disease-specific eligibility filter.",
 "Observed instance-2 MRI outcomes/ICV only; no outcome or exposure imputation in this script.",
 "260 planned MRI measures: 66 thickness, 27 FA, 27 MD, 139 volumes and 1 WMH.",
 "All outcomes directly standardized per complete-case model, including WMH (no log transform).",
 "beta/SE/LCI/HCI: outcome SD per one FLT point; *_standardized: per one FLT SD.",
 "HC1 robust SE and residual-df two-sided t inference; minimum N=100.",
 paste("Adjustment covariates:",paste(covs.all,collapse=", ")),
 "Other covariates match s1_3; volumes/WMH additionally adjust observed ICV.",
 "Primary FDR_family: within FLT and MRI family; FDR_trait: within MRI phenotype (3 FLTs).",
 "FDR_FLTs_family pools 3 FLTs per modality; FDR_exposure_all pools 260 MRIs per FLT; FDR_all_tests pools all 780 tests.",
 "Failed/missing tests remain in planned multiplicity counts.",
 paste("Chronology:",chronology),
 "Associations only; no causal interpretation. Inspect QC and residual diagnostics before publication."
),file.path(OUT,"README_analysis.txt"))
capture.output(sessionInfo(),file=file.path(OUT,"sessionInfo.txt"))
if(file.exists(file.path(OUT,"checkpoint_partial.rds"))) unlink(file.path(OUT,"checkpoint_partial.rds"))
print(counts)
if(all(r$status=="failed")) stop("All models failed; inspect error table")
message("DONE output: ",OUT)
