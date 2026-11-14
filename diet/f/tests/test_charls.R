#!/usr/bin/env Rscript
# Meaningful synthetic checks for the declared CHARLS component adapter.
# No participant data or source-specific CHARLS variable mappings are used.
argv <- commandArgs(trailingOnly=TRUE)
this <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1L])
root <- if(length(argv)) normalizePath(argv[1]) else normalizePath(file.path(dirname(this),"../.."))
if(!all(vapply(c("jsonlite","survival"),requireNamespace,logical(1),quietly=TRUE)))
    stop("Tests require jsonlite and survival")
rscript <- Sys.getenv("DIET_TEST_RSCRIPT",file.path(R.home("bin"),"Rscript"))
tmp <- tempfile("diet-charls-tests-",tmpdir=Sys.getenv("DIET_TEST_TMPDIR",tempdir()))
dir.create(tmp,recursive=TRUE)
Sys.setenv(DIET_ROOT=root)
source(file.path(root,"f/common.R"))
source(file.path(root,"f/dose/scoring.R"))
cfg <- jsonlite::fromJSON(file.path(root,"config/charls.json"),simplifyVector=FALSE)
results <- list()
check <- function(name,code) {
    e <- tryCatch({force(code);NULL},error=function(e)e)
    ok <- is.null(e)
    results[[length(results)+1L]] <<- data.frame(test=name,passed=ok,
        detail=if(ok)"" else conditionMessage(e),stringsAsFactors=FALSE)
    cat(if(ok)"PASS" else "FAIL",name,if(!ok)paste0(": ",conditionMessage(e))else "","\n")
    invisible(ok)
}
expect_error <- function(code,pattern=NULL) {
    e<-tryCatch({force(code);NULL},error=function(e)e)
    if(is.null(e))stop("Expected an error")
    if(!is.null(pattern)&&!grepl(pattern,conditionMessage(e),ignore.case=TRUE))stop(conditionMessage(e))
}
run <- function(stage,input,phe=NULL,out,config=cfg,success=TRUE,pattern=NULL) {
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
    conf<-file.path(out,"config.json");jsonlite::write_json(config,conf,auto_unbox=TRUE,pretty=TRUE)
    Sys.setenv(DIET_CHARLS_INPUT=input,DIET_OUTDIR=out,DIET_CHARLS_CONFIG=conf)
    if(!is.null(phe))Sys.setenv(DIET_CHARLS_PHE=phe)
    log<-suppressWarnings(system2(rscript,c(shQuote(file.path(root,"f/charls.R")),stage),stdout=TRUE,stderr=TRUE))
    status<-attr(log,"status");if(is.null(status))status<-0L
    writeLines(log,file.path(out,paste0(stage,".log")))
    if(success && status!=0L)stop(paste(log,collapse="\n"))
    if(!success && status==0L)stop("Expected worker to reject invalid input")
    if(!is.null(pattern)&&!any(grepl(pattern,log,ignore.case=TRUE)))stop(paste(log,collapse="\n"))
    invisible(log)
}
save_input<-function(x,name){p<-file.path(tmp,paste0(name,".rds"));saveRDS(x,p);p}
foods<-dose_food_rules("dose11")$Food
favorable<-data.frame(eid="A",matrix(c(1,1,1,1,1,.5,1.56,3,2,1,0),nrow=1),check.names=FALSE)
names(favorable)<-c("eid",foods);attr(favorable,"intake_unit")<-"servings/day"

check("DOSE frozen boundaries and half-domain missingness",{
    z<-calculate_dose(favorable);stopifnot(z$dose11==11,z$dose9==9,z$dose9half==9)
    d<-favorable;d$Whole_grains<-1;stopifnot(calculate_dose(d)$dose9half==8)
    d$Whole_grains<-1+1e-8;stopifnot(calculate_dose(d)$dose9half==9)
    d$Whole_grains<-3+1e-8;stopifnot(calculate_dose(d)$dose9half==8)
    d<-favorable;d$Eggs<-.5+1e-8;stopifnot(calculate_dose(d)$dose9half==8)
    d<-favorable;d$Dried_fruit<-NA_real_;z<-calculate_dose(d)
    stopifnot(is.na(z$dose11),is.na(z$dose9half),z$dose9==9)
    d<-favorable;d[foods]<-0;stopifnot(calculate_dose(d)$dose9half==1)
})
check("DOSE rejects invalid intake, duplicate IDs, and wrong units",{
    d<-favorable;d$Eggs<- -1;expect_error(calculate_dose(d),"Negative")
    d<-favorable;d$Eggs<-Inf;expect_error(calculate_dose(d),"infinite")
    d<-rbind(favorable,favorable);expect_error(calculate_dose(d),"unique")
    d<-favorable;attr(d,"intake_unit")<-"g/day";expect_error(calculate_dose(d),"servings/day")
})

make_components<-function(n=240L){
    set.seed(19061)
    x<-data.frame(id=sprintf("C%04d",seq_len(n)),stringsAsFactors=FALSE)
    for(f in foods)x[[f]]<-runif(n,0,2)
    x$Whole_grains<-runif(n,0,4);x$Sweetened_beverages<-sample(c(0,0,.5,1),n,TRUE)
    # Measured, energy-adjusted domain inputs, not invented CHARLS raw mappings.
    x$protein_gkg<-runif(n,.6,1.8);x$dairy_srv<-runif(n,0,4)
    x$veg_srv<-runif(n,0,5);x$fruit_srv<-runif(n,0,3);x$wholegrain_srv<-runif(n,0,3)
    x$fat_srv<-runif(n,0,2);x$upf_srv<-runif(n,0,6)
    x$alcohol_srv<-runif(n,0,3);x$sodium_mg<-runif(n,1200,6500)
    for(v in unlist(cfg$scores$DASH$columns))x[[v]]<-sample(1:5,n,TRUE)
    for(v in unlist(cfg$scores$MEDI$columns))x[[v]]<-sample(0:1,n,TRUE)
    x
}
baseline<-make_components()
check("Missing complete component columns leave whole scores unavailable",{
    x<-baseline;x[c("Wine","sodium_mg","DASH_sodium","MEDI_alcohol")]<-NULL
    out<-file.path(tmp,"missing_columns");run("disco",save_input(x,"missing_columns"),out=out)
    z<-readRDS(file.path(out,".cache/disco/charls.charls.rds"))
    stopifnot(nrow(z)==nrow(x),all(is.na(as.matrix(z[grep("^diet\\.",names(z))]))))
})
check("Missing nutrition and nonnumeric factor values never become zero or level codes",{
    x<-baseline
    x[1,setdiff(names(x),"id")]<-NA
    favorable_maha<-c(protein_gkg=1.2,dairy_srv=3,veg_srv=3,fruit_srv=2,wholegrain_srv=2,fat_srv=1,upf_srv=0,alcohol_srv=0,sodium_mg=2300)
    for(v in names(favorable_maha))x[3,v]<-favorable_maha[[v]]
    p<-as.character(x$protein_gkg);p[2]<-"Unknown";x$protein_gkg<-factor(p)
    out<-file.path(tmp,"missing_values");run("disco",save_input(x,"missing_values"),out=out)
    z<-readRDS(file.path(out,".cache/disco/charls.charls.rds"))
    stopifnot(all(is.na(z[1,grep("^diet\\.",names(z))])),is.na(z$diet.maha.sum[2]),abs(z$diet.maha.sum[3]-90)<1e-10)
})
check("Repeated participant IDs and negative inputs are rejected",{
    x<-baseline;x$id[2]<-x$id[1]
    run("disco",save_input(x,"duplicate_ids"),out=file.path(tmp,"duplicate_ids"),success=FALSE,pattern="unique")
    x<-baseline;x$id[2]<-paste0(" ",x$id[1]," ")
    run("disco",save_input(x,"duplicate_trim"),out=file.path(tmp,"duplicate_trim"),success=FALSE,pattern="unique")
    x<-baseline;x$Eggs[1]<- -9
    run("disco",save_input(x,"negative"),out=file.path(tmp,"negative"),success=FALSE,pattern="negative")
})

# Complete pipeline, including outcome-specific status and dated censoring.
main_out<-file.path(tmp,"pipeline")
input<-save_input(baseline,"pipeline_components")
phe<-data.frame(id=baseline$id,start_date=as.Date("2010-01-01"),end_date=as.Date("2010-01-01")+2000,
                age=40+seq_len(nrow(baseline))%%35,sex=seq_len(nrow(baseline))%%2)
for(j in seq_along(cfg$outcomes)){
    y<-cfg$outcomes[[j]];ix<-seq_len(nrow(phe))
    status<-as.integer((ix+j)%%3==0)
    date<-rep(as.Date(NA),nrow(phe));date[status==1]<-phe$start_date[status==1]+300+(ix[status==1]*13+j*101)%%1400
    phe[[paste0(y,"_status")]]<-status;phe[[paste0(y,"_date")]]<-date
}
phe$cad_status[1:7]<-c(1,NA,1,0,1,1,1)
phe$cad_date[1:7]<-phe$start_date[1:7]+c(-1,NA,NA,NA,500,300,2500)
phe$death_status[1:7]<-c(0,0,0,0,1,1,0)
phe$death_date[1:7]<-phe$start_date[1:7]+c(NA,NA,NA,NA,300,500,NA)
phe$death_status[9]<-1;phe$death_date[9]<-as.Date(NA)
phe$start_date[10]<-as.Date(NA);phe$end_date[11]<-phe$start_date[11]
phe$cad_status[12]<-1;phe$cad_date[12]<-phe$start_date[12]
phe$end_date[15]<-as.Date(NA);phe$death_status[15]<-1;phe$death_date[15]<-phe$start_date[15]+1000
phe$age[20]<-NA_real_
phenofile<-save_input(phe,"pipeline_phenotypes")
pipeline_ok<-check("CHARLS disco and canonical publication run end to end",{
    run("disco",input,out=main_out)
    z<-diet_publish_disco("charls",main_out,c(charls=file.path(main_out,".cache/disco/charls.charls.rds")))
    stopifnot(nrow(z)==nrow(baseline),!anyDuplicated(z$id),identical(z$id,z$ID),length(diet_score_columns(z))==4L)
    run("phe",input,phenofile,out=main_out)
})
if(pipeline_ok){
    check("Unknown or undated events are excluded; prevalent cases never enter risk sets",{
        x<-readRDS(file.path(main_out,"phe/charls.charls.rds"));cad<-x$outcomes$cad
        excluded<-baseline$id[c(1,2,3,9,10,11,12,15)]
        stopifnot(!any(excluded %in% cad$id))
    })
    check("Events and death censoring produce correct positive follow-up",{
        x<-readRDS(file.path(main_out,"phe/charls.charls.rds"));cad<-x$outcomes$cad
        rows<-cad[match(baseline$id[c(4,5,6,7)],cad$id),]
        stopifnot(!anyNA(rows$id),identical(rows$event,c(0L,0L,1L,0L)))
        stopifnot(isTRUE(all.equal(rows$time,c(2000,300,300,2000)/365.25,tolerance=1e-12)))
        for(y in names(x$outcomes))stopifnot(all(is.finite(x$outcomes[[y]]$time)&x$outcomes[[y]]$time>0))
    })
    check("Association sample counts and BH family match all estimable configured models",{
        run("assoc",input,phenofile,out=main_out)
        a<-read.delim(file.path(main_out,"assoc/charls.charls.tsv"),check.names=FALSE)
        stopifnot(nrow(a)==16L,all(a$status=="ok"),!anyNA(a$p))
        stopifnot(isTRUE(all.equal(a$fdr,p.adjust(a$p,"BH"),tolerance=1e-12)))
        x<-readRDS(file.path(main_out,"phe/charls.charls.rds"))
        for(i in seq_len(nrow(a))){
            d<-x$outcomes[[a$outcome[i]]];v<-c(a$diet[i],unlist(cfg$covariates),"event","time")
            d<-d[complete.cases(d[v])&is.finite(d$time)&d$time>0,,drop=FALSE]
            stopifnot(a$n[i]==nrow(d),a$events[i]==sum(d$event))
        }
    })
}
summary<-do.call(rbind,results)
if(any(!summary$passed)){
    print(summary[!summary$passed,],row.names=FALSE)
    cat("Synthetic fixtures retained at",tmp,"\n")
    quit(save="no",status=1)
}
unlink(tmp,recursive=TRUE)
cat(nrow(summary),"CHARLS/DOSE checks passed. These are synthetic contracts, not manuscript replication.\n")
