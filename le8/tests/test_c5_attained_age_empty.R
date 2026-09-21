# Rscript tests/test_c5_attained_age_empty.R
suppressPackageStartupMessages({
  library(dplyr);library(purrr);library(survival);library(ggplot2)
})
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))[1]
fdir<-file.path(dirname(normalizePath(script)),"..","f")
# Load only function definitions; never start the pipeline or read cohort inputs.
for(expr in parse(file.path(fdir,"c5_legacy.R"))){
  if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&
     as.character(expr[[2]])%in%c("attained_age_score_sensitivity","plot_attained_age_score_sensitivity"))eval(expr)
}
bt<-function(x)paste0("`",x,"`")
set.seed(27)
n<-400L
d<-tibble(eid=seq_len(n),age=runif(n,40,70),clinical=rnorm(n),
  time=runif(n,1,10),event=rep(0:1,n/2))|>
  mutate(.attained_entry=age,.attained_exit=age+time)
scores<-tibble(eid=d$eid,method="Reference",split="validation",score_z=rnorm(n))
run<-function(dat=d,rows=scores,methods="Reference",covars=c("age","clinical"))
  attained_age_score_sensitivity(rows,dat,methods,covars,"time","event")
empty<-run(methods="Absent")
stopifnot(nrow(empty)==0L,identical(names(empty),
  c("method","time_scale","beta","se","N","events","lo","hi")))
check_empty<-function(z)stopifnot(identical(z,empty))
# The reported failure: a method exists, but both complete-case fits are empty.
for(usable in c(0L,8L,95L)){
  missing<-d;missing$clinical[seq_len(n)>usable]<-NA_real_
  check_empty(run(missing))
}
few_events<-d;few_events$event<-0L;few_events$event[1:19]<-1L
check_empty(run(few_events))
check_empty(run(rows=mutate(scores,split="training")))
check_empty(run(rows=mutate(scores,score_z=NA_real_)))
# Model errors also preserve the output schema (one observed factor level).
bad_factor<-d;bad_factor$clinical<-factor(rep("only",n))
check_empty(run(bad_factor))
# Empty output is renderable as the existing unavailable panel.
blank_plot<-function(title,label)ggplot()+labs(title=title,subtitle=label)
p<-plot_attained_age_score_sensitivity(empty)
grDevices::pdf(NULL)
tryCatch(stopifnot(grepl("No attained-age score model was estimable",p$labels$subtitle),
  inherits(ggplotGrob(p),"gtable")),finally=grDevices::dev.off())
# Compare successful fits against direct Cox models, including confidence limits.
z<-run()
direct<-list(
  coxph(Surv(time,event)~score_z+age+clinical,data=inner_join(d,scores,by="eid")),
  coxph(Surv(.attained_entry,.attained_exit,event)~score_z+clinical,data=inner_join(d,scores,by="eid")))
stopifnot(nrow(z)==2L,all(z$N==n),all(z$events==n/2))
for(i in seq_along(direct)){
  sm<-coef(summary(direct[[i]]))
  stopifnot(isTRUE(all.equal(z$beta[i],unname(sm["score_z","coef"]))),
    isTRUE(all.equal(z$se[i],unname(sm["score_z","se(coef)"]))))
}
stopifnot(isTRUE(all.equal(z$lo,z$beta-1.96*z$se)),
  isTRUE(all.equal(z$hi,z$beta+1.96*z$se)))
# Missing attained-age fields must not discard a valid baseline-time analysis.
missing_age<-d;missing_age$.attained_entry<-NA_real_;missing_age$.attained_exit<-NA_real_
stopifnot(isTRUE(all.equal(run(missing_age),z[1,])))
# Age adjustment is not needed for the attained-age fit; baseline time is not needed either.
missing_baseline<-d;missing_baseline$age<-NA_real_;missing_baseline$time<-NA_real_
stopifnot(isTRUE(all.equal(run(missing_baseline),z[2,])))
# One unestimable method cannot discard another method's valid results.
mixed<-bind_rows(scores,mutate(scores[1:8,],method="Too small"))
stopifnot(isTRUE(all.equal(run(rows=mixed,methods=c("Too small","Reference")),z)))
# Reordered cached score rows are still joined by participant ID.
stopifnot(isTRUE(all.equal(run(rows=scores[n:1,]),z,tolerance=1e-10)))
cat("C5 attained-age empty and estimable sensitivity regression tests passed.\n")
