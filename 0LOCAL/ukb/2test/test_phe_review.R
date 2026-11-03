# Synthetic regressions for the 2026-11-02 phenotype review.
# Run from any directory: Rscript /path/to/scripts/ukb/2test/test_phe_review.R
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
stopifnot(requireNamespace("impute", quietly = TRUE),
          requireNamespace("R.utils", quietly = TRUE))
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), "..", ".."))
source(file.path(root, "0f/phenotype.R"))
source(file.path(root, "ukb/f/date_utils.R"))
D <- as.Date
checks <- character()
ok <- function(s) { checks <<- c(checks, s); cat("PASS:", s, "\n") }
err <- function(expr) inherits(tryCatch({ force(expr); NULL }, error = identity), "error")
exprs <- parse(file.path(root, "ukb/f/phe.R"))
stage <- function(name) {
  hit <- which(vapply(exprs, function(x) is.call(x) &&
    identical(x[[1L]], as.name("run_if")) && identical(x[[2L]], name), logical(1)))
  stopifnot(length(hit) == 1L)
  exprs[[hit]][[3L]]
}
for (x in exprs) {
  if (is.call(x) && identical(x[[1L]], as.name("<-")) &&
      as.character(x[[2L]])[1L] %in% c("phe_num", "derive_ukb_uacr")) eval(x)
}

stopifnot(length(date_follow_end) == 1L, date_follow_end == D("2023-04-01"),
          identical(ukb_date_parse(c("2010-02-03", "03/02/2010", "2010/02/03")),
                    rep(D("2010-02-03"), 3)),
          all(is.na(ukb_date_valid(c("1900/01/01", "2037-07-07")))))
ok("fixed administrative cutoff and equivalent date formats")

code <- data.frame(eid = c("B", "A"), p20002_i0 = c("200|201", "100|101"),
                   p20002_i1 = c("202", "102"))
dates <- data.frame(eid = c("A", "B"), p20008_i1_a0 = c(2010, 2012),
                    p20008_i0_a1 = c(2001, 2003), p20008_i0_a0 = c(2000, 2002))
z <- ukb_expand_codes(code, dates, "srd")
stopifnot(identical(z$eid, c("A", "B")), identical(z$i1, c("102", "202")),
          identical(z$i2, c("101", "201")), identical(z$i3, c("100", "200")),
          err(ukb_expand_codes(rbind(code, code[1, ]), dates, "srd")))
ok("code/date pairing by eid, assessment instance and array")

deaths <- data.frame(eid = c("A", "B"), p40000_i0 = c("2024-01-05", "1900-01-01"),
  p40001_i0 = c("C50", "I21"), p40000_i1 = c("2024-01-07", "2019-03-04"), p40001_i1 = c("I21", "I25"))
dd <- ukb_death_dates(deaths)
stopifnot(identical(dd$date_death, D(c("2024-01-05", "2019-03-04"))),
          identical(dd$date_cvd_death, D(c("2024-01-07", "2019-03-04")))); ok("true death dates and corresponding causes retained")

cc <- data.frame(eid = as.character(1:4), i1 = c("7410", "0410", "4109", "4140"))
dc <- data.frame(eid = as.character(1:4), d1 = rep(D("2010-01-01"), 4))
ans <- get_icd_fast("icd9", cc, dc,
  data.frame(pattern = c("41[0-4]", "410"), trait = c("cvd_cad", "mi")), verbose = FALSE)
stopifnot(identical(ans$cnt_icd9_mi, c(0L, 0L, 1L, 0L)),
          identical(ans$cnt_icd9_cvd_cad, c(0L, 0L, 1L, 1L)), all(is.na(ans$fod_icd9_mi[1:2])))
cc$i1 <- c("1074", "1075", "11074", "10740")
ans <- get_icd_fast("srd", cc, dc, data.frame(pattern = "1074|1075", trait = "cad"), verbose = FALSE)
stopifnot(identical(ans$cnt_srd_cad, c(1L, 1L, 0L, 0L)))
ok("ICD prefixes and exact SRD codes; false substring matches excluded")

previous_cap <- Sys.getenv("SRD_YEAR_CAP", unset = NA_character_)
Sys.unsetenv("SRD_YEAR_CAP")
stopifnot(ukb_srd_year_cap() == 2026L, ukb_srd_year_cap("2028") == 2028L,
          err(ukb_srd_year_cap("2026.5")), err(ukb_srd_date(2027.1, 1960, TRUE, 2026, TRUE)),
          identical(ukb_srd_date(c(-1, -3, 2000.5, 2024.5), rep(1960, 4), TRUE, 2026, TRUE),
                    D(c(NA, NA, "2000-07-02", "2024-07-02"))))
ok("SRD input-year limit independent of follow-up; excessive years fail explicitly")

dat <- data.frame(eid = as.character(1:7), birth_date = rep(D("1960-01-15"), 7),
  date_attend = rep(D("2000-01-01"), 7), date_lost = rep(D(NA), 7), date_death = rep(D("2027-01-01"), 7),
  sex = rep(0:1, length.out = 7), fod_srd_cad = D(c("1999-01-01", "2000-01-01", "2010-01-01",
    "2023-04-01", "2024-07-02", NA, NA)), cnt_srd_cad = c(1L, 1L, 1L, 1L, 1L, 1L, 0L))
timed <- t2e(dat, NA, "fod_srd_cad", "birth_date", "date_attend", "date_lost", "date_death", date_follow_end, "srd", "day")
ob <- adu_observations(dat, "fod_srd_cad", date_follow_end, "cnt_srd_cad")
stopifnot(identical(timed$srd.Yt2e, c(NA_real_, NA_real_, 1, 1, 0, NA_real_, 0)),
          identical(ob$status, c(1L, 1L, 1L, 1L, 0L, NA_integer_, 0L)),
          timed$srd.t2e[5] == as.numeric(date_follow_end - dat$date_attend[5]),
          identical(timed$date_death, dat$date_death))
ok("t2e/ADuLT distinguish prevalence, incident disease, later diagnosis and undated evidence")

tmp <- tempfile("phe-review-"); dir.create(tmp); dir.create(file.path(tmp, "rap")); dir.create(file.path(tmp, "Rdata")); dir.create(file.path(tmp, "common"))
saveRDS(data.frame(eid = c("B", "A"), birth_year = c(1960, 1950)), file.path(tmp, "Rdata/phe.rds"))
input <- data.frame(eid = c("A", "B"), p20008_i0_a0 = c(2024.5, 2025.6))
fwrite(input, file.path(tmp, "rap/srd.year.tab.gz"), sep = "\t", compress = "gzip")
e <- new.env(parent = globalenv()); e$indir <- tmp
eval(stage("srd"), e)
meta_file <- file.path(tmp, "Rdata/srd.date.meta.rds"); output_file <- file.path(tmp, "rap/srd.date.tab.gz")
meta <- readRDS(meta_file)
stopifnot(meta$year_cap == 2026L, identical(meta$administrative_end, date_follow_end),
          all(!is.na(fread(output_file)$p20008_i0_a0)))
before <- tools::md5sum(c(meta_file, output_file))
input$p20008_i0_a0[2] <- 2027.1
fwrite(input, file.path(tmp, "rap/srd.year.tab.gz"), sep = "\t", compress = "gzip")
stopifnot(err(eval(stage("srd"), e)), identical(before, tools::md5sum(c(meta_file, output_file))))
ok("actual SRD stage writes metadata and preserves previous files on rejected input")
if (is.na(previous_cap)) Sys.unsetenv("SRD_YEAR_CAP") else Sys.setenv(SRD_YEAR_CAP = previous_cap)

one <- clean_biom(data.frame(eid = as.character(1:12), single = seq_len(12)))
none <- clean_biom(data.frame(eid = as.character(1:12), absent = rep(NA_real_, 12)))
stopifnot(identical(one$eid, as.character(1:12)), all(one$single == seq_len(12)),
          identical(names(none), "eid"), identical(none$eid, as.character(1:12)))
ok("one/zero retained biomarkers preserve dimensions and participant IDs")

signed <- data.frame(eid = as.character(1:20), NPX_A = seq(-3, 3, length.out = 20), NPX_B = 20:1)
signed_clean <- clean_biom(signed, imp = FALSE)
stopifnot(identical(signed_clean$eid, signed$eid),
          isTRUE(all.equal(signed_clean$NPX_A, signed$NPX_A)),
          isTRUE(all.equal(signed_clean$NPX_B, signed$NPX_B)))
u <- derive_ukb_uacr(data.frame(eid = c("a", "b", "c"), bb_UALB = c(10, NA, -1), bb_UCRE = c(8840, 8840, 8840)))
stopifnot(identical(u$bb_UACR, c(10, NA_real_, NA_real_)), identical(u$bb_UACR_assessed, c(TRUE, FALSE, FALSE)))
ok("signed NPX and raw numeric UACR behavior preserved")

# Execute the production Met block, including real KNN and final invalid-value mask.
parts <- as.list(stage("biom"))[-1L]
texts <- vapply(parts, function(x) paste(deparse(x), collapse = " "), character(1))
start <- which(grepl("/common/met.lst", texts, fixed = TRUE)); finish <- which(grepl('saveRDS(met,', texts, fixed = TRUE))
stopifnot(length(start) == 1L, length(finish) == 1L)
raw <- data.frame(eid = as.character(1:40), p23460_i0 = 2 + seq_len(40), p23461_i0 = 5 + seq_len(40))
for (j in 1:30) raw[[paste0("p", 24000 + j, "_i0")]] <- seq_len(40) + j / 2
mapping <- c("data_field\tmet_name\tfull_name", "p23460\tAla\tAlanine", "p23461\tGln\tGlutamine",
  paste0("p", 24001:24030, "\tM", 1:30, "\tMeasurement"), "p23460/p23461\tAla2Gln\tRatio")
writeLines(mapping, file.path(tmp, "common/met.lst"))
run_met <- function(values) {
  fwrite(values, file.path(tmp, "rap/met.tab.gz"), sep = "\t", compress = "gzip")
  env <- new.env(parent = globalenv()); env$indir <- tmp
  for (x in parts[start:finish]) eval(x, env)
  readRDS(file.path(tmp, "Rdata/met.rds"))
}
normal <- run_met(raw)
stopifnot(identical(normal$eid, raw$eid), isTRUE(all.equal(normal$Ala, as.numeric(raw$p23460_i0))),
          isTRUE(all.equal(normal$Ala2Gln, raw$p23460_i0 / raw$p23461_i0)))
ok("production Met block preserves ordinary observed values and valid ratios")
bad <- raw; bad$p23460_i0[1] <- -8; bad$p23461_i0[1] <- -2
bad$p23461_i0[2] <- 0; bad$p23460_i0[3] <- 0; bad$p23460_i0[4] <- Inf
bad$p23460_i0[5] <- NA_real_
clean <- run_met(bad)
stopifnot(nrow(clean) == 40L, all(is.na(clean$Ala2Gln[c(1, 2, 4)])), clean$Ala2Gln[3] == 0,
          is.na(clean$Ala[1]), is.na(clean$Gln[1]), !is.na(clean$Ala2Gln[5]))
ok("invalid ratio operands remain NA after KNN; valid zero and ordinary missing data retain their policies")

for (x in parts) if (is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name("prepare_ukb_imaging"))) eval(x)
fwrite(data.frame(eid = c("a", "b"), p26517_i2 = c(1, 2), p26517_i3 = c(3, 4)),
  file.path(tmp, "rap/img.tab.gz"), sep = "\t", compress = "gzip")
writeLines("p26517\tvolume\tall", file.path(tmp, "common/img.lst"))
extra <- data.frame(eid = c("a", "b"), p53_i2 = rep("2020-01-01", 2))
fwrite(extra, file.path(tmp, "rap/img.extra.tab.gz"), sep = "\t", compress = "gzip")
img <- prepare_ukb_imaging(tmp)
stopifnot(identical(img$eid, c("a", "b")), identical(names(img), c("eid", "img_p26517")),
          identical(img$img_p26517, 1:2))
img_paths <- file.path(tmp, "Rdata", c("img.raw.rds", "img.visits.rds", "img.dictionary.csv"))
img_before <- tools::md5sum(img_paths)
fwrite(rbind(extra, extra[1, ]), file.path(tmp, "rap/img.extra.tab.gz"), sep = "\t", compress = "gzip")
stopifnot(err(prepare_ukb_imaging(tmp)), identical(img_before, tools::md5sum(img_paths)))
ok("i2 imaging preserved; duplicate extra IDs fail before changing outputs")
unlink(tmp, recursive = TRUE)
cat("TOTAL_PASS:", length(checks), "\n")
