# Exercise the actual C1 paired Cox and C2 cross-fitted decomposition with
# bare and suffixed score names. Only fixture I/O and figure rendering are stubbed.
suppressPackageStartupMessages({library(dplyr); library(purrr); library(data.table); library(survival)})
script <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
root <- dirname(dirname(normalizePath(script)))
load_functions <- function(file) {
  walk <- function(x) {
    if (missing(x) || (!is.call(x) && !is.expression(x))) return(invisible(NULL))
    if (is.call(x) && identical(x[[1]], as.name('<-')) && length(x) == 3L &&
        is.symbol(x[[2]]) && is.call(x[[3]]) && identical(x[[3]][[1]], as.name('function'))) {
      eval(x, envir=.GlobalEnv); return(invisible(NULL))
    }
    for (item in as.list(x)) walk(item)
  }
  walk(parse(file))
}
for (file in c('0.common.R', 'c1.correlate.R', 'c2.cause.R')) load_functions(file.path(root, 'f', file))
if (exists('split', envir=.GlobalEnv, inherits=FALSE)) rm(split, envir=.GlobalEnv)
scratch <- tempfile('le8-pgs-join-', tmpdir='/tmp'); dir.create(scratch)
.le8_identity_hashes <- new.env(parent=emptyenv())
SEED <- 2026L; N_CORES <- 1L; Y <- 'cvd_cad'; vars.basic <- 'age'; C2_FIXED_TOP <- c('PCSK9', 'LPA')
out.prot <- out.met <- scratch
Sys.unsetenv(c('LE8_GROUP_FILE', 'LE8_GROUP_COLUMN', 'PGS_GROUP_COLUMN', 'LE8_OUTER_ROSTER'))
Sys.setenv(C1_DIRECTION_ANCHORS='PCSK9,LPA', C2_LE4_COVARS='', C2_TREATMENT_VARS='',
           C2_STATE_FEATURES='PCSK9,LPA', C2_RUN_STATE_PROJECTION='TRUE', C2_DECOMP_BOOT='3')
set.seed(SEED); n <- 2400L
scores <- data.frame(eid=as.character(seq_len(n)), PCSK9=rnorm(n), LPA=rnorm(n))
ph <- data.frame(eid=scores$eid, age=rnorm(n, 60, 7), .group=rep(seq_len(n/2), each=2))
biom <- data.frame(eid=scores$eid, PCSK9=.6*scores$PCSK9+rnorm(n), LPA=.4*scores$LPA+rnorm(n))
ph$.prevalent <- rbinom(n, 1, plogis(-1+.4*biom$PCSK9))
ph$cvd_cad.t2e <- pmin(rexp(n, exp(.3*biom$PCSK9)/12), 16)
ph$cvd_cad.Yt2e <- as.integer(ph$cvd_cad.t2e < 16)
dat <- inner_join(ph, biom, by='eid')
# Pre-existing internal-looking columns and shuffled string IDs cannot overwrite data.
dat$.le8_pgs_1 <- 17
j <- le8_join_pgs(dat, scores[n:1,], c(PCSK9='PCSK9', LPA='LPA'))
stopifnot(identical(j$data$PCSK9, dat$PCSK9), all(j$data$.le8_pgs_1 == 17),
          identical(j$data[[j$score_map[['PCSK9']]]], scores$PCSK9),
          !anyDuplicated(names(j$data)), identical(le8_join_pgs(dat,scores,character())$data,dat))
bad <- tryCatch(le8_join_pgs(dat, scores[c(1,1),], c(PCSK9='PCSK9')), error=identity)
stopifnot(inherits(bad,'error'))
cat('PASS PGS join preserves distinct measured values, scores, IDs and existing columns\n')
find_c1_pgs_file <- find_c2_score_file <- function(layer) 'fixture'
read_c1_pgs <- read_c2_scores <- function(file) scores
read_prot <- function() biom
read_all <- function(...) ph
filter_analysis_cohort <- function(d) d
make_outcome <- function(d, ...) d
make_prevalent_status <- function(d, ...) d$.prevalent
le8_job_dir <- function(...) scratch
write_raw_csv <- function(d, name, dir) data.table::fwrite(d,file.path(dir,name))
le8_plot_c1_temporal <- function(...) invisible(NULL)
read_c2_heritability <- function(...) list(data=tibble(feature=character(),snp_h2=numeric(),snp_h2_se=numeric()),status=tibble(status='unavailable'))
bare1 <- le8_c1_additions(dat, 'protein', 'age', 'cvd_cad.t2e', 'cvd_cad.Yt2e')
bare2 <- run_individual_genetic_decomposition('protein', scratch, C2_FIXED_TOP)
stopifnot(nrow(bare1$paired_associations) == 8L, all(is.finite(bare1$paired_associations$beta)),
          all(bare2$summary$status == 'ok'), all(bare2$summary$score_column == C2_FIXED_TOP),
          all(is.finite(bare2$summary$incident_beta_genetic)), bare2$state_projection$status$status == 'ok',
          setequal(fread(file.path(scratch,'c2.state_PGS_calibration.csv'))$pgs, C2_FIXED_TOP))
names(scores)[-1] <- paste0(names(scores)[-1], '.pgs')
suffixed1 <- le8_c1_additions(dat, 'protein', 'age', 'cvd_cad.t2e', 'cvd_cad.Yt2e')
suffixed2 <- run_individual_genetic_decomposition('protein', scratch, C2_FIXED_TOP)
stopifnot(isTRUE(all.equal(bare1, suffixed1)),
          isTRUE(all.equal(select(bare2$summary, -score_column), select(suffixed2$summary, -score_column))),
          isTRUE(all.equal(bare2$folds, suffixed2$folds)),
          isTRUE(all.equal(bare2$trajectory, suffixed2$trajectory)),
          isTRUE(all.equal(bare2$state_projection, suffixed2$state_projection)),
          all(suffixed2$summary$score_column == paste0(C2_FIXED_TOP,'.pgs')))
cat('PASS C1 separate/joint Cox fits identical for bare and suffixed PGS columns\n')
cat('PASS C2 decomposition, temporal scans and frozen state projection identical; source names retained\n')
unlink(scratch,recursive=TRUE)
