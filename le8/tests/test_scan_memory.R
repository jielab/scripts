# Run with Rscript tests/test_scan_memory.R from the repository root.
# Extract only functions: no participant data, pipeline startup or installation.
suppressPackageStartupMessages({library(dplyr); library(tibble); library(survival)})
wanted <- c("parallel_map", "cox_scan", "cox_scan_delayed_entry", "bt")
for (expr in parse("f/comm.f.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      as.character(expr[[2]]) %in% wanted) eval(expr)
}
set.seed(29)
n <- 1200L
d <- data.frame(t=rexp(n), event=rbinom(n,1,.4), age=rnorm(n,55,8),
                a=rnorm(n), b=rnorm(n), constant=1)
d$entry <- runif(n,40,60); d$exit <- d$entry+d$t
d$a[seq(1,n,17)] <- NA_real_
d$age[seq(1,n,23)] <- NA_real_
scan <- function() list(
  cox_scan(d,c("a","b","constant"),"age",time_var="t",event_var="event"),
  cox_scan_delayed_entry(d,c("a","b","constant"),"age",
    entry_var="entry",exit_var="exit",event_var="event"))
optimized_map <- parallel_map
parallel_map <- function(x,fun) lapply(x,fun)
reference <- scan()
parallel_map <- optimized_map
for (cores in c(1L,2L)) {
  N_CORES <- cores
  stopifnot(isTRUE(all.equal(reference,scan(),tolerance=1e-12)))
  stopifnot(identical(parallel_map(integer(),identity),list()))
  stopifnot(identical(unname(unlist(parallel_map(1:5,function(x)x*x))),
                      c(1L,4L,9L,16L,25L)))
}
cat("PASS: serial/fork scans preserve Cox estimates, missing-data handling and task order\n")

N_CORES <- 2L
stopifnot(identical(parallel_map(1:2,function(x)NULL),list(NULL,NULL)))
failed <- suppressWarnings(try(parallel_map(1:2,function(x)stop("test worker failure")),silent=TRUE))
stopifnot(inherits(failed,"try-error"))
if (.Platform$OS.type != "windows") {
  failed <- suppressWarnings(try(parallel_map(1:2,function(x) {
    if(x==1L) tools::pskill(Sys.getpid(),signal=9L)
    x
  }),silent=TRUE))
  stopifnot(inherits(failed,"try-error"),grepl("refusing partial results",as.character(failed)))
}
cat("PASS: worker failure/SIGKILL cannot silently produce partial scans; intentional NULL is preserved\n")

# Verify the narrow mediation sample preserves participant selection and values.
for (expr in parse("f/c4_connect.R")) {
  if (is.call(expr) && identical(expr[[1]],as.name("<-")) &&
      identical(expr[[2]],as.name("stratified_sample"))) eval(expr)
}
set.seed(30)
wide <- as.data.frame(matrix(rnorm(6000L*100L),nrow=6000L))
wide$eid <- seq_len(nrow(wide)); wide$event <- rep(c(0,1,NA),length.out=nrow(wide))
cols <- c("eid","event","V1","V2")
set.seed(31); full_sample <- stratified_sample(wide,"event",1000L)
set.seed(31); narrow_sample <- stratified_sample(wide[,cols,drop=FALSE],"event",1000L)
stopifnot(identical(full_sample[,cols],narrow_sample),
          as.numeric(object.size(narrow_sample)) < as.numeric(object.size(full_sample))/10)
cat("PASS: narrow C4 sampling preserves participants and values; fixture table size reduced >90%\n")
