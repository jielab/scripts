# Synthetic export regression: no participant data and no network.
suppressPackageStartupMessages(library(dplyr))
entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(entry), ".."))
source(file.path(root, "../0f/phenotype.R"))
source(file.path(root, "f/gwas_export.R"))
Sys.setenv(ADULT_ENABLE = "FALSE")
D <- as.Date
dat <- data.frame(eid = letters[1:9], ethnic.c = c(rep("White", 8), "Other"),
    center = "c1", gen_batch = c(-11, 1, NA, 10, 20, 30, 40, 50, 60),
    birth_date = D("1960-01-01"), date_attend = D(c(rep("2010-01-01", 8), "2025-01-01")),
    date_death = D(c(rep(NA, 5), "2015-01-01", NA, NA, NA)),
    date_lost = D(c(rep(NA, 6), "2012-01-01", NA, NA)))
for (nm in c("tdi", "edu.sco", "age", "sex", "bmi", "height", "drug.lipid", "drug.dm", "drug.htn", "bald")) dat[[nm]] <- 1
for (j in 1:40) dat[[paste0("PC", j)]] <- seq_len(nrow(dat)) + j / 100
for (trait in ukb_gwas_traits) {
    dat[[paste0("fod_icd10_", trait)]] <- D(c("2011-01-01", "2009-01-01", NA, NA,
        "2024-01-01", "2016-01-01", "2012-01-01", "1900-01-01", NA))
    dat[[paste0("cnt_icd10_", trait)]] <- c(1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 0L)
}
out <- ukb_gwas_events(dat, date_follow_end)
expected <- c(1, NA, 0, NA, 0, 0, 1, NA, NA)
for (trait in ukb_gwas_traits) {
    stopifnot(identical(out[[paste0(trait, ".Yt2e")]], expected),
        out[[paste0(trait, ".t2e")]][1] == 365 / 365.25,
        out[[paste0(trait, ".t2e")]][6] == as.numeric(D("2015-01-01") - D("2010-01-01")) / 365.25)
}
stopifnot(identical(as.character(out$batch), c("b-11", "b1", NA, "b10", "b20", "b30", "b40", "b50", "b60")))
# Disease evidence from an unrelated trait must not exclude a control.
dat$cnt_icd10_cvd_cad[3] <- 1L
out <- ukb_gwas_events(dat, date_follow_end)
stopifnot(is.na(out$cvd_cad.Yt2e[3]), out$t2dm.Yt2e[3] == 0, out$ad.Yt2e[3] == 0)

# Run the actual publication stage, including LDL alignment and staged writes.
indir <- tempfile("ukb-gwas-test-"); dir.create(indir)
dir.create(file.path(indir, "Rdata")); outdir <- indir
saveRDS(dat, file.path(indir, "Rdata/all.rds"))
saveRDS(data.frame(eid = rev(dat$eid), bb_LDL = 9:1), file.path(indir, "Rdata/bbc.rds"))
read_rds <- function(name) readRDS(file.path(indir, "Rdata", paste0(name, ".rds")))
exprs <- parse(file.path(root, "f/phe.R"))
stage <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("run_if")) &&
    identical(x[[2]], "phe4gwas"), as.list(exprs))[[1]][[3]]
eval(stage)
exported <- read.table(file.path(indir, "ukb.phe"), header = TRUE, check.names = FALSE)
stopifnot(identical(head(names(exported), 5), c("FID", "IID", "ethnic.c", "center", "batch")),
    identical(grep("^PC[0-9]+$", names(exported), value = TRUE), paste0("PC", 1:10)),
    !any(grepl("(r2e|bi2e)$", names(exported))), nrow(exported) == 8L,
    identical(exported$ldl, 1:8), is.na(exported$batch[3]),
    nrow(read_rds("all")) == 9L,
    all(c(outer(ukb_gwas_traits, c(".t2e", ".Yt2e"), paste0)) %in% names(exported)))
# Re-running on all.rds must not multiply columns or rows.
before <- tools::md5sum(file.path(indir, "ukb.phe"))
eval(stage)
stopifnot(identical(before, tools::md5sum(file.path(indir, "ukb.phe"))))
bad <- read_rds("all"); bad$PC10 <- NULL; saveRDS(bad, file.path(indir, "Rdata/all.rds"))
error <- tryCatch({ eval(stage); NULL }, error = identity)
stopifnot(inherits(error, "error"), identical(before, tools::md5sum(file.path(indir, "ukb.phe"))))
unlink(indir, recursive = TRUE)
cat("PASS: 16 disease outcomes, censoring, batch/PC schema, safe repeat/export failure\n")
