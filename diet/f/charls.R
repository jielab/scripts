#!/usr/bin/env Rscript
# CHARLS has no source-specific food adapter in the supplied repository.
# This module scores declared, measured components. It never equates a broad
# frequency question with a quantitative nutrient/serving measurement.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L || !args[1] %in% c("disco","phe","assoc"))
    stop("Usage: Rscript charls.R disco|phe|assoc")
root <- Sys.getenv("DIET_ROOT")
source(file.path(root,"f/common.R"))
out <- diet_env("DIET_OUTDIR",required=TRUE)
input <- diet_env("DIET_CHARLS_INPUT",file.path(diet_env("DIET_CHARLS"),"food_components.tsv"))
phenofile <- diet_env("DIET_CHARLS_PHE",file.path(diet_env("DIET_CHARLS"),"phenotypes.tsv"))
config <- diet_env("DIET_CHARLS_CONFIG",file.path(root,"config/charls.json"))
if(!requireNamespace("jsonlite",quietly=TRUE)) stop("Install jsonlite (see install.sh)")
cfg <- jsonlite::fromJSON(config,simplifyVector=FALSE)
if(args[1]=="disco") {
    raw <- diet_read(input)
    raw$id <- diet_id(raw,"charls")
    res <- data.frame(id=raw$id,stringsAsFactors=FALSE)
    available <- list()
    columns <- function(score, required, units) {
        spec <- cfg$scores[[score]]
        if (is.null(spec) || !isTRUE(spec$enabled)) return(NULL)
        if(!identical(spec$unit,units)) stop(score,": config unit must be ",units)
        map <- spec$columns
        missing <- required[!required %in% names(map)]
        if(!length(missing)) missing <- required[!unlist(map[required]) %in% names(raw)]
        if(length(missing)) {
            available[[length(available)+1L]] <<- data.frame(score=score,status="unavailable",
                reason=paste("Missing measured components:",paste(missing,collapse=", ")))
            return(NULL)
        }
        x <- as.data.frame(lapply(required,function(k) {
            v <- raw[[map[[k]]]]
            if(is.factor(v)) v <- as.character(v)
            suppressWarnings(as.numeric(v))
        }))
        names(x) <- required
        if(any(vapply(x,function(v) any(v<0 | is.infinite(v),na.rm=TRUE),logical(1))))
            stop(score,": negative or infinite components; recode source missing codes before scoring")
        available[[length(available)+1L]] <<- data.frame(score=score,status="available",
            reason=paste("Measured components; definition",spec$definition))
        x
    }
    # Fixed manuscript recipe; no fitting to disease outcomes.
    foods <- c("Apples_Pears","Dried_fruit","Legumes","Olive_oil","Tomatoes","Eggs",
               "Green_leafy_vegetables","Whole_grains","Wine","Processed_meats","Sweetened_beverages")
    d <- columns("DOSE", foods,"servings/day")
    if(!is.null(d)) {
        d$eid <- raw$id; attr(d,"intake_unit") <- "servings/day"
        source(file.path(root,"f/dose/scoring.R"))
        z <- calculate_dose(d)
        res$diet.dose9half.sum <- z$dose9half
    } else res$diet.dose9half.sum <- NA_real_
    # MAHA domains are the measured, energy-adjusted inputs to the UKB rule.
    # They must be constructed by a declared CHARLS measurement mapping.
    d <- columns("MAHA",c("protein_gkg","dairy_srv","veg_srv","fruit_srv","wholegrain_srv",
                         "fat_srv","upf_srv","alcohol_srv","sodium_mg"),"declared_domain_units")
    if(!is.null(d)) {
        up <- function(v,t) pmin(pmax(v/t,0),1)*10
        down <- function(v,good,bad) pmin(pmax((bad-v)/(bad-good),0),1)*10
        pts <- data.frame(protein=up(d$protein_gkg,1.2),dairy=up(d$dairy_srv,3),veg=up(d$veg_srv,3),
            fruit=up(d$fruit_srv,2),wholegrain=up(d$wholegrain_srv,2),fat=up(d$fat_srv,1),
            upf=down(d$upf_srv,0,4),alcohol=down(d$alcohol_srv,0,2.5),sodium=down(d$sodium_mg,2300,5000))
        res$diet.maha.sum <- rowSums(pts,na.rm=FALSE)
    } else res$diet.maha.sum <- NA_real_
    # Declared component-point routes allow validated local DASH/MEDI variants
    # without pretending that CHARLS measured every original nutrient/food.
    for(score in c("DASH","MEDI")) {
        spec <- cfg$scores[[score]]
        nm <- paste0("diet.",tolower(score),".sum")
        if(!isTRUE(spec$enabled)) {res[[nm]] <- NA_real_; next}
        required <- unlist(spec$components)
        if(!length(required)) {res[[nm]] <- NA_real_; next}
        d <- columns(score,required,"component_points")
        if(is.null(d)) {res[[nm]] <- NA_real_; next}
        bounds <- spec$bounds
        for (k in required) if(any(d[[k]]<bounds[[1]] | d[[k]]>bounds[[2]],na.rm=TRUE))
            stop(score," component ",k," outside declared scoring bounds")
        res[[nm]] <- rowSums(d,na.rm=FALSE)
    }
    attr(res,"score_columns") <- grep("^diet\\.",names(res),value=TRUE)
    attr(res,"measurement_mapping") <- cfg$scores
    attr(res,"availability") <- available
    attr(res,"profile") <- "CHARLS measured-component adapter; see definitions in config"
    diet_save(res,file.path(out,".cache/disco/charls.charls.rds"))
    if(length(available)) diet_tsv(do.call(rbind,available),file.path(out,"disco/charls.availability.tsv"))
} else if(args[1]=="phe") {
    scores <- readRDS(file.path(out,"disco/charls.disco.rds"))
    raw <- diet_read(phenofile); raw$id <- diet_id(raw,"charls")
    if(!all(c("start_date","end_date") %in% names(raw))) stop("CHARLS phenotypes need start_date and end_date")
    raw$start_date <- as.Date(raw$start_date); raw$end_date <- as.Date(raw$end_date)
    data <- merge(scores,raw,by="id",all=FALSE,suffixes=c("",".phe"))
    if(!nrow(data)) stop("CHARLS score/phenotype IDs have no overlap")
    outcomes <- unlist(cfg$outcomes)
    long <- list()
    for(y in outcomes) {
        event_name <- paste0(y,"_date")
        if(!event_name %in% names(data)) next
        status_name <- paste0(y,"_status")
        if(!status_name %in% names(data)) stop("CHARLS needs ",status_name,
            ": 0=confirmed absent through end_date, 1=known event, NA=unknown")
        status <- suppressWarnings(as.numeric(as.character(data[[status_name]])))
        if(any(!is.na(status) & !status %in% c(0,1))) stop(status_name, " must be 0, 1, or NA")
        event_date <- as.Date(data[[event_name]])
        if(any(status==0 & !is.na(event_date),na.rm=TRUE)) stop(y,": status=0 conflicts with an observed event date")
        observed <- !is.na(status) & (status==0 | !is.na(event_date))
        end <- data$end_date
        if(y!="death" && "death_date" %in% names(data)) {
            death <- as.Date(data$death_date)
            if(!"death_status" %in% names(data)) stop("CHARLS death_date requires death_status (0/1/NA)")
            ds <- suppressWarnings(as.numeric(as.character(data$death_status)))
            if(any(!is.na(ds) & !ds %in% c(0,1))) stop("death_status must be 0, 1, or NA")
            if(any(ds==0 & !is.na(death),na.rm=TRUE)) stop("death_status=0 conflicts with death_date")
            observed <- observed & !is.na(ds) & (ds==0 | !is.na(death))
            end <- pmin(end,death,na.rm=TRUE)
        }
        prevalent <- !is.na(event_date) & event_date<=data$start_date
        keep <- observed & !prevalent & !is.na(data$start_date) & !is.na(data$end_date) & !is.na(end) & end>data$start_date
        keep[is.na(keep)] <- FALSE
        d <- data[keep,,drop=FALSE]
        e <- event_date[keep]; stop_date <- end[keep]
        d$event <- as.integer(!is.na(e)&e<=stop_date)
        d$time <- as.numeric(pmin(e,stop_date,na.rm=TRUE)-d$start_date)/365.25
        d$outcome <- y; long[[y]] <- d
    }
    if(!length(long)) stop("No observed outcome dates; expected *_date columns from config")
    diet_save(list(data=data,outcomes=long,id="id",profile="charls"),file.path(out,"phe/charls.charls.rds"))
} else {
    if(!requireNamespace("survival",quietly=TRUE)) stop("Install survival")
    x <- readRDS(file.path(out,"phe/charls.charls.rds"))
    results <- list(); failed <- list()
    scores <- diet_score_columns(x$data); covs <- unlist(cfg$covariates)
    miss <- setdiff(covs,names(x$data))
    if(length(miss)) stop("Missing configured covariates: ",paste(miss,collapse=", "))
    for(y in names(x$outcomes)) for(sc in scores) {
        d <- x$outcomes[[y]]; vars <- c(sc,covs,"event","time")
        d <- d[complete.cases(d[vars]) & is.finite(d$time) & d$time>0,,drop=FALSE]
        sd0 <- sd(d[[sc]])
        if(nrow(d)<20 || sum(d$event)<5 || !is.finite(sd0) || sd0<=0) next
        d$exposure <- (d[[sc]]-mean(d[[sc]]))/sd0
        formula <- as.formula(paste("survival::Surv(time,event) ~",paste(c("exposure",covs),collapse=" + ")))
        fit <- tryCatch(survival::coxph(formula,data=d,ties="efron"),error=function(e)e)
        if(inherits(fit,"error")) {failed[[length(failed)+1L]]<-data.frame(outcome=y,score=sc,reason=conditionMessage(fit)); next}
        s <- summary(fit); b <- coef(fit)[["exposure"]]; se <- sqrt(vcov(fit)["exposure","exposure"])
        results[[length(results)+1L]] <- data.frame(study="charls",cohort="charls",outcome=y,diet=sc,
          model="configured_covariates",scale="SD",effect_type="HR",beta=b,se=se,
          estimate=exp(b),lower=exp(b-1.959964*se),upper=exp(b+1.959964*se),
          p=2*pnorm(-abs(b/se)),n=nrow(d),events=sum(d$event),adjustment=paste(covs,collapse=";"),
          status="ok",stringsAsFactors=FALSE)
    }
    if(!length(results)) stop("No estimable CHARLS score/outcome model; inspect component coverage and events")
    ans <- do.call(rbind,results); ans$fdr <- p.adjust(ans$p,"BH")
    diet_tsv(ans,file.path(out,"assoc/charls.charls.tsv"))
    if(length(failed)) diet_tsv(do.call(rbind,failed),file.path(out,"assoc/charls.failures.tsv"))
}
