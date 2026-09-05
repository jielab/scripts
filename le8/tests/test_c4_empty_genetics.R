# Run without sourcing comm.f.R or loading participant data.
# Usage: Rscript test_c4_empty_genetics.R /path/to/f/c4_connect.R
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(patchwork); library(forcats)
})
for (expr in parse(args[[1]])) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      as.character(expr[[2]]) %in% c("make_genetic_edges", "plot_c4")) eval(expr)
}
FDR_CUT <- .05
safe_log <- log
disease <- tibble(term=c("A","B","C"), beta=c(.1,.2,.3), p.value=c(.01,.02,.03))
membership <- tibble(feature=c("A","B","C"),set="YS",primary_component="diet")
disc <- tibble(feature=c("A","B","C"),component="cad_prs",component_var="cad_prs",
  r=c(.2,.3,.4),z=c(4,5,6),p.value=.001,FDR_component=.01)
rep_scan <- disc |> mutate(r=c(.2,-.3,.4),FDR_component=c(.01,.01,.2))
normal <- make_genetic_edges(disc,rep_scan,disease,membership)
stopifnot(identical(normal |> filter(replicated) |> pull(feature), "A"))
stopifnot(normal$same_direction[match("B",normal$feature)]==FALSE)
cat("PASS: populated scans retain direction/FDR replication criteria\n")
empty_cases <- list(list(tibble(),tibble()),list(disc,tibble()),list(tibble(),rep_scan),
  list(disc[0,],rep_scan[0,]))
for (inputs in empty_cases) {
  edges <- make_genetic_edges(inputs[[1]],inputs[[2]],disease,membership)
  stopifnot(nrow(edges)==0L,identical(names(edges),names(normal)),
    identical(vapply(edges,typeof,character(1)),vapply(normal,typeof,character(1))),
    identical(edges |> filter(replicated) |> pull(feature) |> unique(), character()))
}
cat("PASS: absent and one-sided scans preserve schema and export empty anchors\n")
# Exercise the actual six-panel plot function, including ggplot rendering.
names.le8 <- "diet"; cols_le8 <- c(diet="#2166AC"); Y <- "cvd_cad"
cap <- function(x,m) pmax(-m,pmin(m,x))
theme_5c <- function(size=10) theme_minimal(base_size=size)
forest_theme <- theme_5c
blank_plot <- function(title,subtitle=NULL) ggplot()+theme_void()+labs(title=title,subtitle=subtitle)
scan <- tibble(feature="A",component="diet",split="discovery",z=4)
primary <- tibble(feature="A",strict_YS=TRUE,YS_model=TRUE,FDR_disc=.01,z_disc=4)
sets <- list(primary=primary,
  membership=tibble(feature="A",set="YS",group="Other",primary_component="diet",disease_beta=.1,disease_p=.01),
  YS_edges=tibble(feature="A",component="diet",FDR_disc=.01))
med <- tibble(indirect_beta=double(),indirect_p=double(),feature=character(),component=character())
plots_built <- character()
save_plot <- function(p,filename,...,outdir) {
  if(inherits(p,"patchwork")) patchwork::patchworkGrob(p) else ggplotGrob(p)
  plots_built <<- c(plots_built,filename)
}
genetics <- make_genetic_edges(tibble(),tibble(),disease,membership)
for(layer in c("protein","metabolite")) plot_c4(scan,sets,med,genetics,layer,tempdir())
stopifnot(length(plots_built)==12L)
cat("PASS: all six C4 panels render for both layers without genetic edges\n")
