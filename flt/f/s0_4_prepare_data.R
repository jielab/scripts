source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
library(dplyr)
library(readxl)

project <- Sys.getenv("FLT_OUT")
indir <- Sys.getenv("FLT_UKB")
outdir <- file.path(project, "00data")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Cohort and covariates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
covs.all <- c("age", "sex", "edu", "inc", "smoke_status", "alcohol_status", "bmi", "tdi", "pa_met_tot") # "ethnic.c3"

dat0 <- as.data.frame(readRDS(file.path(indir, "Rdata/all.rds")))
flt <- as.data.frame(readRDS(file.path(Sys.getenv("FLT_PREP"), "flt.rds")))
brain <- as.data.frame(readRDS(file.path(Sys.getenv("FLT_PREP"), "fod.icd10.brain.rds")))

dat0$eid <- as.character(dat0$eid)
flt$eid <- as.character(flt$eid)
brain$eid <- as.character(brain$eid)
flt_vars <- read_excel(Sys.getenv("FLT_DICTIONARY"), sheet = "item182")$name # 187L
flt <- flt[c("eid", flt_vars, "p20750")]
flt[flt_vars] <- lapply(flt[flt_vars], function(x) replace(x, !is.finite(x), NA_real_))

brain_def <- data.table::fread(file.path(indir, "common/icd10_brain.lst"), sep = "\t", header = FALSE, fill = TRUE, select = 1:3, col.names = c("pattern", "name", "label"), data.table = FALSE)
Ys <- brain_def$name
brain_def$label[is.na(brain_def$label) | !nzchar(brain_def$label)] <- Ys[is.na(brain_def$label) | !nzchar(brain_def$label)]
outcomes <- data.frame(outcome = Ys, label = brain_def$label, pattern = brain_def$pattern, date_field = paste0("fod_icd10_", Ys), count_field = paste0("cnt_icd10_", Ys))
Ys <- outcomes$outcome
brain <- brain[c("eid", outcomes$date_field, outcomes$count_field)]
for (nm in c("dat0", "flt", "brain")) {
  ids <- get(nm)$eid
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop(nm, ": eid 缺失或重复")
}

dat <- flt %>%
  inner_join(dat0 %>% select(eid, date_attend, date_lost, date_death,
    age, sex, ethnic.c, edu, inc, smoke_status, alcohol_status, bmi, tdi, pa_met_tot), by="eid", relationship="one-to-one", na_matches="never") %>%
  left_join(brain, by="eid", relationship="one-to-one", na_matches="never")
cat("FLT 队列人数：", nrow(dat), "\n") # FLT 队列人数： 181835  
rm(dat0, flt, brain) 
date_cols <- c("date_attend", "date_lost", "date_death", outcomes$date_field)
dat[date_cols] <- lapply(dat[date_cols], as.Date)

# Covariate
has_code <- function(x, codes) {grepl(paste0("(^|\\|)(", paste(codes, collapse = "|"), ")(\\||$)"), x)}
current_status <- function(x) {factor(ifelse(x %in% c(0, 1), "No", ifelse(x == 2, "Yes", NA)), c("No", "Yes"))}
dat <- dat %>% mutate(
  sex = factor(sex, 0:1, c("Female", "Male")),
  ethnic.c3 = factor(ifelse(is.na(ethnic.c), NA, ifelse(ethnic.c == "White", "White", "Others")), c("White", "Others")),
  smoke_status = current_status(smoke_status),
  alcohol_status = current_status(alcohol_status),
  edu = factor(case_when(has_code(edu, 1) ~ "high", has_code(edu, c(2, 5, 6)) ~ "middle", has_code(edu, c(3, 4, 7, -7)) ~ "low", TRUE ~ NA_character_), c("low", "middle", "high")),
  inc = factor(case_when(inc == 1 ~ "low", inc %in% c(2, 3) ~ "middle", inc %in% c(4, 5) ~ "high", TRUE ~ NA_character_), levels = c("low", "middle", "high")),
  pa_met_tot = ifelse(is.finite(pa_met_tot) & pa_met_tot >= 0, pa_met_tot, NA_real_)
)
for (v in covs.all) {if (is.numeric(dat[[v]])) dat[[v]][!is.finite(dat[[v]])] <- NA_real_}
dat$complete_covariates <- complete.cases(dat[covs.all])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Baseline diseases, follow-up and outcomes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 从入组日期起算；分结局排除既往病，以及有诊断但无日期等无法确定状态者。
# 与 0phe.f.R 的随访截止日期保持一致。
source(file.path(Sys.getenv("FLT_SHARED"), "phenotype.R"))
dat$start_date <- dat$date_attend
dat$end0 <- pmin(dat$date_lost, dat$date_death, date_follow_end, na.rm = TRUE)
dat$valid_followup <- !is.na(dat$start_date) & dat$end0 > dat$start_date
tab_event <- list()
for (i in seq_along(Ys)) {
  Y <- Ys[i]
  d <- dat[[outcomes$date_field[i]]]
  count <- dat[[outcomes$count_field[i]]]
  # 排除计数/日期矛盾及无效计数，避免将未知状态当作无病。
  consistent <- is.finite(count) & count >= 0 & count == floor(count) & ((count == 0 & is.na(d)) | (count > 0 & !is.na(d)))
  eligible <- consistent & dat$valid_followup & !is.na(count) & (count == 0 | !is.na(d)) & (is.na(d) | d > dat$start_date)
  eligible[is.na(eligible)] <- FALSE
  event <- !is.na(d) & d > dat$start_date & d <= dat$end0
  time <- as.numeric(pmin(d, dat$end0, na.rm = TRUE) - dat$start_date) / 365.25
  dat[[paste0(Y, ".eligible")]] <- eligible
  dat[[paste0(Y, ".Yt2e")]] <- ifelse(eligible, as.integer(event), NA_integer_)
  dat[[paste0(Y, ".t2e")]] <- ifelse(eligible, time, NA_real_)
  followup <- time[eligible]
  tab_event[[i]] <- data.frame(outcome = Y, excluded_inconsistent = sum(!consistent), total = sum(eligible), cases = sum(event[eligible]),
    person_years = sum(followup), mean_fu = if (length(followup)) mean(followup) else NA_real_,
    median_fu = if (length(followup)) median(followup) else NA_real_,
    complete_covariates = sum(dat$complete_covariates[eligible]))
}
tab_event <- bind_rows(tab_event)
dat.base <- dat[rowSums(dat[paste0(Ys, ".eligible")]) > 0, , drop = FALSE]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Covariates and missing values
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stopifnot(all(covs.all %in% names(dat.base)))
covs_miss <- bind_rows(lapply(covs.all, function(v) {
  x <- dat.base[[v]]
  tibble(variable=v, class=class(x)[1], n_total=length(x),
    n_missing=sum(is.na(x)), pct_missing=round(mean(is.na(x))*100, 2),
    n_nonmissing=sum(!is.na(x)),
    n_level=if (is.factor(x) || is.character(x)) n_distinct(x, na.rm=TRUE) else NA_integer_)
})) %>% arrange(desc(pct_missing), variable)
print(covs_miss, n=Inf)

# 连续变量填均值，分类变量填众数；食物偏好和结局不插补。
dat.base.imp <- dat.base
for (v in covs.all) {
  x <- dat.base[[v]]
  if (all(is.na(x))) stop("协变量全部缺失，无法插补：", v)
  value <- if (is.numeric(x)) mean(x, na.rm = TRUE) else names(which.max(table(x)))
  dat.base.imp[[v]][is.na(x)] <- value
}
dat.base.imp$complete_covariates <- complete.cases(dat.base.imp[covs.all])
stopifnot(all(dat.base.imp$complete_covariates))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat.base 保留未插补协变量；完整病例分析还须筛选 complete_covariates。
analysis_config <- list(Ys = Ys, flt_vars = flt_vars, covs.all = covs.all, outcomes = outcomes, date_follow_end = date_follow_end,
  primary_data_file = "dat.base.imp.rds", complete_case_data_file = "dat.base.rds")
saveRDS(dat.base, file.path(outdir, "dat.base.rds"))
saveRDS(dat.base.imp, file.path(outdir, "dat.base.imp.rds"))
saveRDS(analysis_config, file.path(outdir, "analysis_config.rds"))
write.csv(tab_event, file.path(outdir, "tab_event.csv"), row.names = FALSE)
print(tab_event)
