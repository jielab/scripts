suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(survival) })
source("f/ukb_followup.R")
source("f/nhanes_meat_sensitivity.R")
load_functions <- function(path, env = .GlobalEnv) {
  for (x in parse(path, encoding = "UTF-8")) {
    if (is.call(x) && identical(x[[1]], as.name("<-")) && is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function"))) eval(x, env)
  }
}
# Boundary cases: before/on/after landmark, event at censoring, loss and death.
d <- data.frame(eid = as.character(1:8), diet_landmark = as.Date(rep("2010-01-01", 8)),
  date_lost = as.Date(c(NA, NA, NA, "2011-01-01", "2011-01-01", NA, NA, "2010-01-01")),
  date_death = as.Date(c(NA, NA, NA, NA, NA, "2011-01-01", "2009-12-31", NA)),
  fod_icd10_test = as.Date(c("2009-12-31", "2010-01-01", "2010-01-02", "2011-01-02", "2011-01-01", NA, NA, NA)),
  icd10Ct_cvd = 1L)
x <- maha_outcome_followup(d, "test", "2012-01-01")
stopifnot(identical(x$test.Yt2e, c(NA_integer_, NA_integer_, 1L, 0L, 1L, 0L, NA_integer_, NA_integer_)),
  all(x$test.t2e[!is.na(x$test.t2e)] > 0), abs(x$test.t2e[3] - 1/365.25) < 1e-12)
x <- maha_outcome_followup(d, "death", "2012-01-01")
stopifnot(x$death.Yt2e[6] == 1, is.na(x$death.Yt2e[7]))
d$cnt_icd10_test <- c(rep(0L, 5), 1L, 0L, 0L)
stopifnot(is.na(maha_outcome_followup(d, "test", "2012-01-01")$test.Yt2e[6]))
raw <- data.frame(eid = 1:2)
for (k in 0:4) raw[[paste0("p105010_i", k)]] <- as.Date(rep(NA, 2))
raw$p105010_i0 <- as.Date(c("2009-01-01", "2009-01-01"))
raw$p105010_i4 <- as.Date(c("2012-01-01", NA))
raw$p20082_i4 <- c(1, 30) # A late recall rejected by duration QC still bounds the window.
w <- maha_webq_window(raw)
stopifnot(w$webq_window_end[1] == as.Date("2012-01-01"), w$webq_undated[2])
z <- maha_set_landmark(data.frame(eid=1:2, date_attend=as.Date(rep("2008-01-01",2)), age=50, score=5), w, "score")
stopifnot(z$diet_landmark[1] == as.Date("2012-01-01"), is.na(z$diet_landmark[2]), z$webq_timing_excluded[2])
# P-value underflow is flagged, not confused with a measured P of 1e-308.
r <- data.frame(snp = rep(c("A", "B"), each = 3), p = c(1e-80, .01, .1, 0, .03, .7))
rp <- maha_phewas_p(r)
stopifnot(identical(rp$p_raw, r$p), identical(rp$p, r$p), rp$logp_display[1] == 50,
  rp$p_underflow[4], isTRUE(all.equal(rp$FDR[1:3], p.adjust(r$p[1:3], "BH"))),
  isTRUE(all.equal(rp$FDR[4:6], p.adjust(r$p[4:6], "BH"))))
stopifnot(identical(maha_processed_meat_flag(c("beef steak", "pork chop", "hamburger", "ham, sliced", "sausage", "corned beef", "lamb")),
                    c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE)))
# Exercise the actual food/score functions; original score must be unchanged.
load_functions("f/MAHA_nhanes.R")
old <- new.env(parent = .GlobalEnv)
load_functions("reviewer_revision/code_before_20260906/MAHA_nhanes.R", old)
food_vars <- c("fruit","berry","vegetable","green_leafy","allium","legumes","nuts","dairy","lowfat_dairy","whole_grain","refined_grain","fish","poultry","red_processed_meat","processed_meat_proxy","sweets_pastries","ssb","fried_fast","coffee_tea")
nutrient_vars <- c("kcal","protein_g","carb_g","sugar_g","fiber_g","fat_g","sfat_g","mufa_g","pufa_g","chol_mg","sodium_mg","alcohol_g")
diet.lst <- c(maha="MAHA", dash="DASH", mind="MIND", medi="MEDI")
group_pct_th <- .4
set.seed(606)
dd <- as.data.frame(matrix(runif(1000 * length(c(food_vars, nutrient_vars)), 1, 100), nrow=1000))
names(dd) <- c(food_vars,nutrient_vars)
dd$processed_meat_proxy <- pmin(dd$red_processed_meat, dd$processed_meat_proxy)
dd$sex <- sample(1:2, 1000, TRUE); dd$weight_kg <- runif(1000, 45, 100)
newscore <- construct_diet_scores(dd); oldscore <- old$construct_diet_scores(dd)
for (nm in grep("^diet[.]", names(oldscore), value = TRUE)) stopifnot(isTRUE(all.equal(newscore[[nm]], oldscore[[nm]])))
for (nm in c("maha_processed", "maha_nomeat")) stopifnot(any(abs(newscore[[paste0("diet.",nm,".sum")]] - newscore$diet.maha.sum) > 0))
newscore$age <- runif(1000,40,80); newscore$death <- rbinom(1000,1,.3)
newscore$death_time_y <- runif(1000,.1,15); newscore$wt <- runif(1000,.5,2)
newscore$strata <- factor(rep(1:10,each=100)); newscore$psu <- factor(rep(1:20,each=50))
sens <- nhanes_meat_sensitivity(newscore, "age")
stopifnot(nrow(sens$mortality_common_sample) == 6, length(unique(sens$mortality_common_sample$N)) == 1,
  all(is.finite(sens$mortality_common_sample$HR_per_SD)))
cat("Timing boundaries, missing dates, raw P/FDR, original score invariance, and weighted sensitivity tests PASSED.\n")
# Tiny real PheWAS call validates schema and resumable cache (no UKB participant data).
ph <- data.frame(id=as.character(1:1000), check.names=FALSE)
ph[["250.2"]] <- rbinom(1000,1,.2); ph[["401.1"]] <- rbinom(1000,1,.4)
geno <- data.frame(eid=ph$id, exposure=rnorm(1000), age=runif(1000,40,80))
td <- tempfile(); dir.create(td)
saveRDS(ph,file.path(td,"phenotypes.rds"))
pa <- maha_raw_phewas(geno, "exposure", "age", file.path(td,"phenotypes.rds"), file.path(td,"cache"))
pb <- maha_raw_phewas(geno, "exposure", "age", file.path(td,"phenotypes.rds"), file.path(td,"cache"))
stopifnot(nrow(pa$res) == 2, identical(pa$res, pb$res), identical(pa$signature, pb$signature), identical(pa$res$p, pa$res$p_raw))
geno$exposure[1] <- 100
pc <- maha_raw_phewas(geno, "exposure", "age", file.path(td,"phenotypes.rds"), file.path(td,"cache"))
stopifnot(!identical(pa$signature,pc$signature))
cat("Actual PheWAS call, unchanged-cache reuse and input-change invalidation PASSED.\n")
