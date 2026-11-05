# Execute native analysis functions on sparse synthetic data, without UKB inputs.
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[1]) else normalizePath(".")
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(purrr); library(survival);
    library(ggplot2); library(patchwork)})
env <- new.env(parent=globalenv())
wanted <- c("le8_empty_association", "le8_empty_mr", "le8_read_c1_optional", "cox_scan",
    "cox_scan_delayed_entry", "logistic_scan", "component_association", "le8_fit_budget_model",
    "le8_prepare_prediction_matrix", "focus_proxy_accuracy", "run_c4_nonlin", "c3_select_mr",
    "pgs_calibrate", "pgs_complete", "pgs_folds", "pgs_bind", "pgs_num", "make_proxy_sets", "focus_plot")
extract <- function(expr) {
    if (missing(expr) || !is.call(expr)) return(invisible(NULL))
    if (identical(expr[[1]], as.name("<-")) && length(expr)==3L && is.symbol(expr[[2]]) &&
        as.character(expr[[2]]) %in% wanted) {
        eval(expr, env)
        return(invisible(NULL))
    }
    if (identical(expr[[1]], as.name("function"))) return(invisible(NULL))
    for (part in as.list(expr)[-1L]) extract(part)
}
for (file in c("0.common.R", "c1.correlate.R", "c2.cause.R", "c3.coloc.R", "c4.connect.R"))
    for (expr in parse(file.path(root, "f", file))) extract(expr)
stopifnot(all(vapply(wanted, exists, logical(1), envir=env, inherits=FALSE)))
env$parallel_map <- function(x, f, ...) lapply(x, f)
env$bt <- function(x) paste0("`", x, "`")
env$`%||%` <- function(x,y) if (is.null(x)) y else x
env$Y <- "toy"
set.seed(57)
d <- data.frame(eid=as.character(seq_len(650)), age=rnorm(650), protein=rnorm(650),
    toy.t2e=c(rep(2,5),rep(13,645)), toy.Yt2e=c(rep(1,5),rep(0,645)), prevalent=rbinom(650,1,.3))
incident <- env$cox_scan(d, "protein", "age")
prevalent <- env$logistic_scan(d, "protein", "age", "prevalent")
stopifnot(nrow(incident)==1, is.na(incident$p.value), incident$status=="skipped_insufficient_sample_or_events",
    nrow(prevalent)==1, is.finite(prevalent$p.value))
# A sufficient layer runs; a sparse feature is skipped independently within it.
large <- d
large$toy.Yt2e <- rbinom(650,1,.3)
large$toy.t2e <- ifelse(large$toy.Yt2e, runif(650,1,9), 13)
large$missing <- large$protein
large$missing[which(large$toy.Yt2e==1)[-seq_len(5)]] <- NA_real_
scan <- env$cox_scan(large, c("protein","missing"), "age")
stopifnot(is.finite(scan$p.value[scan$term=="protein"]), is.na(scan$p.value[scan$term=="missing"]))
stopifnot(nrow(env$cox_scan(d, character(), "age"))==0)

# Low incident support does not remove prevalent decomposition fits.
d$.prevalent <- d$prevalent
inc <- env$component_association(d, "protein", "age", "toy.t2e", "toy.Yt2e")
prev <- env$component_association(d, "protein", "age", "toy.t2e", "toy.Yt2e", TRUE)
stopifnot(is.na(inc["p"]), is.finite(prev["p"]))

# An unestimable PGS calibration fold keeps its schema for C2 decomposition.
small <- data.frame(eid=letters[1:4], .m=1:4, .g=4:1, .prev=rep(0,4))
cal <- env$pgs_calibrate(small, character())
stopifnot(nrow(cal$folds)==0, "slope" %in% names(cal$folds), all(is.na(cal$data$.M)))

# Zero disease events cannot fit a risk model; LE8 reconstruction still works.
d$toy.Yt2e <- 0
d$p1 <- d$protein + rnorm(650, sd=.2)
d$p2 <- rnorm(650)
env$le8_hash_object <- function(...) "synthetic-hash"
fit <- env$le8_fit_budget_model(d[1:500,], d[501:650,], "age", "protein", "toy.t2e", "toy.Yt2e")
proxy <- env$focus_proxy_accuracy(d[1:500,], d[501:650,], "protein", "p1", "age", "YS_1")
stopifnot(grepl("insufficient",fit$status), nrow(proxy)==1, proxy$R2_basic_omics>.8)
# A Cox fit does not require censored controls when event times are observed.
all_events <- d; all_events$toy.Yt2e <- 1L; all_events$toy.t2e <- runif(nrow(d),1,9)
fit_all <- env$le8_fit_budget_model(all_events[1:500,],all_events[501:650,],"age","protein","toy.t2e","toy.Yt2e")
stopifnot(fit_all$status=="ok")

# Replicated LE8 proxies remain available when disease ranking is unestimable.
env$FDR_CUT <- .05; env$YSP_PLUS_N <- 1L; env$NS_N <- 1L
env$truthy <- function(x) toupper(x) %in% c("TRUE", "1")
replicated <- tibble(feature="A", component="bmi", component_var="bmi.pts",
    r=.5, z=5, p.value=.001, FDR_component=.001)
sets <- env$make_proxy_sets(replicated, replicated, env$le8_empty_association(c("A","B")),
    tibble(feature=c("A","B")))
stopifnot(identical(sets$YS,"A"), !length(sets$NS), !length(sets$YSP_plus))

# Empty disease metrics still generate a figure with the independent proxy panels.
env$theme_5c <- function(...) theme_minimal()
env$blank_plot <- function(title, subtitle="") ggplot()+labs(title=title,subtitle=subtitle)
env$le8_queue_figure <- function(plot, ...) stopifnot(inherits(plot,"patchwork"))
env$le8_flush_figures <- function(...) invisible(NULL)
env$focus_plot(tibble(model=character(),stratum=character(),landmark=numeric(),AUC=numeric()),
    tibble(),proxy,tibble(component="p1",cohort="Yin",n=1L),tempdir())

scratch <- tempfile("le8-sparse-native-"); dir.create(scratch)
env$out.prot <- env$out.met <- scratch
env$LE8_JOB <- "c4_connect"
env$LE8_REUSE_RESULTS <- FALSE
env$vars.basic <- "age"
env$vars.le8 <- c("p1","p2")
env$read_all <- function(...) d
env$filter_analysis_cohort <- identity
env$make_outcome <- function(x, ...) x
env$setwd2 <- function(...) invisible(NULL)
env$le8_job_dir <- function(out, module) file.path(out, module)
env$c4_cache_current <- function(...) FALSE
env$le8_stage_fingerprint <- function(...) "synthetic-fingerprint"
env$le8_participant_groups <- function(x) x$eid
env$write_raw_csv <- function(x, name, dir) {
    dir.create(dir,recursive=TRUE,showWarnings=FALSE)
    write.csv(x, file.path(dir,name),row.names=FALSE)
}
env$plot_c4_nonlin <- function(tests, curves, nadir, ...) {
    stopifnot(nrow(tests)==0, nrow(curves)==0, all(c("variable","nadir_x","nadir_HR") %in% names(nadir)))
}
env$plot_c4_penalty <- function(...) invisible(NULL)
env$finalize_outputs <- function(...) invisible(NULL)
env$run_c4_nonlin("prot")
stopifnot(file.exists(file.path(scratch,"c4_connect/c4.nonlin.res.rds")),
    readRDS(file.path(scratch,"c4_connect/c4.penalty.res.rds"))$status$status=="unavailable")

# Missing upstream associations preserve assay identities, without inventing estimates.
fallback <- env$le8_read_c1_optional(scratch,"protein",c("A","B"))
stopifnot(identical(fallback$association$term,c("A","B")), all(is.na(fallback$association$p.value)))
selected <- env$c3_select_mr(env$le8_empty_mr(),"protein",.05,200L)
stopifnot(nrow(selected)==0, "instrument_snps" %in% names(selected))
unlink(scratch, recursive=TRUE)
cat("Sparse native C1-C4 regression checks passed\n")
