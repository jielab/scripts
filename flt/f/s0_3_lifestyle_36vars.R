source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(data.table)
  library(tidyverse)
  library(ggrepel)
  library(lubridate)
  library(broom)
  library(purrr)
  library(survival)
  library(ggplot2)
  library(forcats)
  library(stringr)
})

indir = Sys.getenv("FLT_UKB")
source(file.path(Sys.getenv("FLT_SHARED"), "phenotype.R"))


## ------------------------------------------------------------------
## 提取 19项人体测量＋17项生活方式/社会因素
## ------------------------------------------------------------------
fn	<- fread(paste0(indir, "/common/cmd_anthropometry_lifestyle.lst"), sep = "\t", select = 1:3, fill = TRUE, header = FALSE) %>% as.data.frame()
fn$V1[duplicated(fn$V1)]; fn$V2[duplicated(fn$V2)]
phe0 <- fread(flt_tab("cmd_anthropometry_lifestyle"), sep = "\t", header = TRUE) %>% as.data.frame()
names(phe0)	<- names(phe0) %>% gsub("_i0$|_i0_a0$", "", .)
setdiff(fn$V1, names(phe0)); setdiff(names(phe0), fn$V1)
names(phe0)	<- sapply(names(phe0), replace_code, sep = "_", mapping_df = fn) %>% as.character()

phe <- phe0
phe$eid <- as.character(phe$eid)
for (v in setdiff(fn$V2, c("eid", "edu"))) {
  x <- suppressWarnings(as.numeric(phe[[v]]))
  if (any(!is.na(phe[[v]]) & is.na(x))) stop("存在非数值内容：", v)
  x[!is.finite(x)] <- NA_real_
  if (v != "tdi") x[!is.na(x) & x < 0] <- NA_real_ # Townsend负值有效
  phe[[v]] <- x
}
for (v in c("tv_time", "computer_time")) phe[[v]][which(phe0[[v]] == -10)] <- 0.5 # -10表示不足1小时，按0.5小时处理是分析约定；缺失值不补0
keep_codes <- function(x, codes) { x[!x %in% codes] <- NA_real_; x }
for (v in c("smoke_status", "alcohol_status")) phe[[v]] <- keep_codes(phe[[v]], 0:2)
phe$alcohol_freq <- 6 - keep_codes(phe$alcohol_freq, 1:6) # 0=从不，5=最频繁
phe$inc <- keep_codes(phe$inc, 1:5) # 收入、手机和休闲步行保留类别，不换算金额/分钟。
for (v in c("liking_gym", "liking_stairs", "liking_sweat")) phe[[v]] <- keep_codes(phe[[v]], 1:9)
for (v in c("days_pa_mod", "days_pa_vig")) phe[[v]] <- keep_codes(phe[[v]], 0:7)
for (v in c("sleep_duration", "tv_time", "computer_time")) phe[[v]][which(phe[[v]] > 24)] <- NA_real_
for (v in c("walking_duration", "dura_pa_mod", "dura_pa_vig")) phe[[v]][which(phe[[v]] > 1440)] <- NA_real_
for (v in grep("fat_pct", names(phe), value = TRUE)) phe[[v]][which(phe[[v]] > 100)] <- NA_real_
# MVPA为分钟/周，未做IPAQ截断或MET加权；0天记0，天数>0但时长缺失仍为NA。
weekly_minutes <- function(days, minutes) ifelse(!is.na(days) & days == 0, 0, days * minutes)
phe$summed_mvpa <- weekly_minutes(phe$days_pa_mod, phe$dura_pa_mod) + weekly_minutes(phe$days_pa_vig, phe$dura_pa_vig)
# 教育多选：出现1记大学学历；出现2~6或-7/7记非大学；其余记缺失。
edu_cols <- grep("^edu($|_i0_a[0-9]+$)", names(phe0), value = TRUE)
edu_text <- do.call(paste, c(lapply(phe0[edu_cols], as.character), sep = "|"))
edu_values <- unique(edu_text)
edu_lookup <- vapply(strsplit(edu_values, "[|,;\\[\\] ]+", perl = TRUE), function(z) {
  z <- suppressWarnings(as.numeric(z))
  if (any(z == 1, na.rm = TRUE)) return(1)
  if (any(z %in% c(2:6, -7, 7))) return(0)
  NA_real_
}, numeric(1))
phe$education_degree <- unname(edu_lookup[match(edu_text, edu_values)])
body_vars <- fn$V2[match(paste0("p", c(23099, 23100, 21001, 23111:23126)), fn$V1)]
lifestyle_vars <- c("smoke_status", "alcohol_status", "alcohol_freq", "sleep_duration",
  "phone_call_time", "tv_time", "computer_time", "pa_met_tot", "summed_mvpa", "walking_duration",
  "walk_time", "liking_gym", "liking_stairs", "liking_sweat", "education_degree", "tdi", "inc")
analysis_vars <- c(body_vars, lifestyle_vars)
dat_lifestyle <- phe[c("eid", analysis_vars)]
saveRDS(dat_lifestyle, file.path(Sys.getenv("FLT_PREP"), "lifestyle_36vars.rds"))
