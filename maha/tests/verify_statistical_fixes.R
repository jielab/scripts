suppressPackageStartupMessages({library(dplyr); library(tidyr); library(survival)})
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)[1]
root <- dirname(dirname(normalizePath(sub("^--file=","",arg))))
chns <- parse(file.path(root,"f","MAHA_chns.R"),encoding="UTF-8")
nhanes <- parse(file.path(root,"f","MAHA_nhanes.R"),encoding="UTF-8")
wanted <- c("num","recode_chns_school_years","zstd","good_covariates",
            "fit_time_updated_score","pick","colnum","first_num","cache_identical_inputs",
            "prep_risk_plot_dat")
extract <- function(x) {
  if(missing(x)) return(invisible(NULL))
  if(is.call(x) && identical(x[[1]],as.name("<-")) && is.symbol(x[[2]]) &&
     as.character(x[[2]]) %in% wanted) eval(x,.GlobalEnv)
  if(is.call(x)||is.expression(x)||is.pairlist(x)) for(y in as.list(x)) extract(y)
}
extract(chns)
x <- c(0,11:16,21:29,31:36,10,17,30,37,99,NA)
expected <- c(0,1:6,7:12,10:12,13:18,rep(NA_real_,6))
stopifnot(identical(recode_chns_school_years(x),as.numeric(expected)))
lines <- readLines(file.path(root,"f","MAHA_chns.R"),encoding="UTF-8")
start <- which(lines=="# Education"); end <- which(lines=="# Household income")
baseline_wave <- 2011
for(field in c("A11","EDUYRS")) {
  educ <- tibble(IDIND=1:3,WAVE=2011)
  educ[[field]] <- c(11,21,34)
  eval(parse(text=lines[start:(end-1)]))
  stopifnot(identical(educb$education_years,if(field=="A11") c(1,7,16) else c(11,21,34)))
}
set.seed(42)
d <- tibble(IDIND=rep(1:200,each=2),tstart=rep(c(0,1),200),tstop=rep(c(1,2),200),
            surv_event=rep(c(0,1),200),score_cum=rnorm(400))
build_time_updated <- function(d,mortality_definition) d
res <- fit_time_updated_score(d,"score",character(),"end_wave")
d$score_z <- zstd(d$score_cum)
ref <- broom::tidy(coxph(Surv(tstart,tstop,surv_event)~score_z+cluster(IDIND),data=d,ties="efron"))
stopifnot(isTRUE(all.equal(res$se,ref$robust.se)),
          isTRUE(all.equal(res$CI_low,exp(res$beta-1.96*res$se))),
          isTRUE(all.equal(res$CI_high,exp(res$beta+1.96*res$se))),
          isTRUE(all.equal(res$P,2*pnorm(-abs(res$beta/res$se)))))
count <- 0L
fn <- cache_identical_inputs(function(d,mortality_definition) {
  count <<- count+1L
  transform(d,result=nchar(mortality_definition))
})
z <- data.frame(value=1:3)
a <- fn(z,"end_wave"); b <- fn(z,"end_wave")
stopifnot(identical(a,b),count==1L)
b$value[1] <- 99
stopifnot(identical(fn(z,"end_wave"),a),count==1L)
invisible(fn(z,"last_contact")); invisible(fn(data.frame(value=2:4),"end_wave"))
stopifnot(count==3L)
extract(nhanes)
risk <- tibble(risk=c(.2,.3),lower=c(.1,.25),upper=c(.4,.6),n=c(10,20),stratum="low")
r <- prep_risk_plot_dat(risk,"stratum")
stopifnot(!"sd" %in% names(r),identical(r$mean,risk$risk*100),
          identical(r$lower_plot,risk$lower*100),identical(r$upper_plot,risk$upper*100))

# Scoped publication assembly must not require another cohort's temporary files.
pub <- parse(file.path(root,"f","MAHA_publication.R"),encoding="UTF-8")
for(e in pub) if(is.call(e)&&identical(e[[1]],as.name("<-"))&&is.symbol(e[[2]])&&
  as.character(e[[2]]) %in% c("supplement_renumbering","finalize_supplement_renumbering")) eval(e)
scratch <- tempfile("maha-publication-test-"); dir.create(scratch)
pub_path <- function(f) file.path(scratch,f)
selected <- supplement_renumbering %>% filter(grepl("chns",tmp_stem))
for(stem in selected$tmp_stem) for(suffix in c(".png",".out.xlsx"))
  writeBin(as.raw(rep(65L,1100L)),pub_path(paste0(stem,suffix)))
sentinel <- pub_path("FigS8.ukb.construct_stress_tests.png")
writeLines("unchanged",sentinel)
finalize_supplement_renumbering("chns")
for(stem in selected$final_stem) for(suffix in c(".png",".out.xlsx"))
  stopifnot(file.exists(pub_path(paste0(stem,suffix))))
stopifnot(identical(readLines(sentinel),"unchanged"))
cat("PASS: schooling mapping, EDUYRS passthrough, robust SE/CI/P, exact-input cache, original risk CI, scoped publication assembly\n")
