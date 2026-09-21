# Run from the repository root: Rscript tests/test_c4_focus_empty.R
suppressPackageStartupMessages({
  library(dplyr);library(tibble);library(purrr);library(ggplot2);library(patchwork)
})
# Load definitions without initializing the analysis or reading participant data.
for(expr in parse("f/c4_focus.R")) {
  if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&
     as.character(expr[[2]])%in%c("focus_contrasts","focus_plot"))eval(expr)
}
metrics<-tibble(model=c("Clinical","NS_5"),stratum="All",landmark=0,
  horizon=10,AUC=c(.6,.65),budget=c(0,5),paradigm=c("Clinical","NS"))
empty<-focus_contrasts(metrics,tibble(),5)
stopifnot(nrow(empty)==0,all(c("stratum","landmark","delta_AUC")%in%names(empty)),
  identical(focus_contrasts(tibble(),tibble(),5),empty))
paired<-bind_rows(metrics,mutate(metrics[2,],model="YS_Yin_5",AUC=.7))
z<-focus_contrasts(paired,tibble(),5)
stopifnot(nrow(z)==1,z$model=="YS_Yin_5",z$reference=="NS_5",
  abs(z$delta_AUC-.05)<1e-12)
boot<-bind_rows(tibble(model="NS_5",replicate=1:25,AUC=.65),
  tibble(model="YS_Yin_5",replicate=1:25,AUC=.7))|>
  mutate(stratum="All",landmark=0,horizon=10)
zb<-focus_contrasts(paired,boot,5)
stopifnot(zb$valid==25,abs(zb$delta_lo-.05)<1e-12,abs(zb$delta_hi-.05)<1e-12)
theme_5c<-function(base_size)theme_classic(base_size)
le8_queue_figure<-function(p,...) {
  # Building the full patchwork also evaluates each panel's aesthetic mappings.
  grob<-patchwork::patchworkGrob(p)
  stopifnot(inherits(grob,"gtable"))
}
le8_flush_figures<-function(...)invisible(NULL)
proxy<-tibble(component="bmi",model="NS_5",delta_R2=.01)
pillars<-tibble(component="bmi",n=0,cohort="Yin")
for(contrasts in list(empty,tibble(),z,zb,mutate(z,stratum="Other")))
  focus_plot(paired,contrasts,proxy,pillars,tempdir())
cat("C4 focus empty and paired comparison regression tests passed.\n")
