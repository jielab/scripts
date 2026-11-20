pacman::p_load(readxl, dplyr, data.table, tidyverse, ggrepel, lubridate, broom, purrr, survival, ggplot2, forcats, stringr)

dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:/")
indir = paste0(dir0, "/data/phe/ukb/phe")
invisible(lapply(c("phe.f.R", "assoc.f.R", "plot.f.R"), function(f) source(file.path(dir0, "scripts", "f", f))))


## ------------------------------------------------------------------
## 提取 icd10_brain.lst
## ------------------------------------------------------------------
# 11 mental disorders【doi:10.1038/s44220-024-00226-0】
# 13 neurological and psychiatric disorders【10.1038/s41467-026-72529-y】
code0 <- readRDS(paste0(indir, "/Rdata/icd10.code0.rds"))
date0 <- readRDS(paste0(indir, "/Rdata/icd10.date0.rds"))
fn <- fread(paste0(indir, "/common/icd10_brain.lst"), sep = "\t", header = FALSE, fill = TRUE) %>% as.data.frame()
colnames(fn)[1:3] <- c("pattern", "name", "label")
dat.other <- data.frame(eid = code0$eid, stringsAsFactors = FALSE)
for(i in 1:nrow(fn)) {print(i); dat.other <- get_icd("icd10", code0, date0, dat.other, fn$pattern[i], fn$name[i])}
saveRDS(dat.other, paste0(indir, "/Rdata/fod.icd10.brain.rds"))


## ------------------------------------------------------------------
## 提取 behavioral.lst
## behavioral phenotypes：5 cognitive, 7 emotional tests
## ------------------------------------------------------------------
## bash ❗
dir0=/mnt/d
phedir=$dir0/data/phe/ukb/phe
for dat in behavioral; do # 需要在 *.lst 添加第三列【all或i0】
	echo RUN: $dat
	awk 'BEGIN{print "eid"}{print $1}' $phedir/common/$dat.lst > $dat.id
	dos2unix $dat.id 
	sort $dat.id | uniq -d
	cat $dat.id | sed 's/^p//; 1d' | sort -nc # 查看排序
	grep -wf $dat.id pheno.id.uniq | wc -l
	grep -vwf pheno.id.uniq $dat.id | head # 显示UKB里面没有的
	grep -nwf $dat.id pheno.id.2col | sed 's/:/ /; s/ /\t/g' > $dat.id.tmp1
	python3 $dir0/scripts/f/0join_file.py -i "$dat.id.tmp1,TAB,1 $phedir/common/$dat.lst,TAB,0" -o $dat.id.tmp2
	cols=$(awk '{if ($3 == "eid" || $3 ~ /^p[0-9]+$/ || $3 ~ /_i0(_a[0-9]+)?$/) printf "%s,", $1}' "$dat.id.tmp2" | sed 's/,$//')
	zcat raw/pheno.tab.gz | cut -f $cols | gzip -f > $dat.tab.gz
done


## QC
# p20507 无后缀; p399_i0_a1 p3
fn <- fread(file.path(indir, "common/behavioral.lst"), header = FALSE, sep = "\t", select = 1:3, fill = TRUE)
phe0 <- as.data.frame(fread(file.path(indir, "rap/behavioral.tab.gz")))
field <- sub("_i0(_a[0-9]+)?$", "", names(phe0))
array <- ifelse(grepl("_a[0-9]+$", names(phe0)), sub(".*(_a[0-9]+)$", "\\1", names(phe0)), ""); array[array == "_a0"] <- ""
idx <- match(field, fn$V1)
keep <- !is.na(idx)
names(phe0)[keep] <- paste0(fn$V2[idx[keep]], array[keep])

phe <- phe0
items <- phe[, setdiff(names(phe), "eid"), drop = FALSE]
for (nm in names(items)) { # 负值和非有限值设为 NA
  x <- items[[nm]]
  if (!is.numeric(x) && !all(is.na(x))) stop("题目不是数值型: ", nm)
  x <- as.numeric(x)
  x[!is.finite(x) | x < 0] <- NA_real_
  items[[nm]] <- x
}
items$pairs_matching <- items$pairs_matching_a2  # 第二轮：6 对卡片

## 2. 共用函数：检查编码、完整题目求和、秩-logit 标准化
valid <- function(x, allowed) {
  bad <- !is.na(x) & !(x %in% allowed)
  if (any(bad)) warning("Unexpected item codes replaced by NA: ", paste(unique(x[bad]), collapse = ", "))
  x[bad] <- NA_real_
  x
}
score <- function(nms, allowed, offset = 0) {
  x <- lapply(items[nms], function(x) valid(x, allowed) - offset)
  rowSums(as.data.frame(x), na.rm = FALSE)
}
logit_scale <- function(x) {
  ans <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n <- sum(ok)
  if (n < 2L || length(unique(x[ok])) < 2L) return(ans)
  z <- qlogis((rank(x[ok], ties.method = "average") - 0.5) / n)
  ans[ok] <- as.numeric(scale(z))
  ans
}

## 3. 生活满意度：未就业且其余四题完整者，按四题均值 × 5 折算
life <- c("job_satisfy", "health_satisfy", "family_satisfy", "friend_satisfy", "finance_satisfy")
not_employed <- !is.na(items$job_satisfy) & items$job_satisfy == 7
items$job_satisfy[not_employed] <- NA_real_
life_items <- as.data.frame(lapply(items[life], valid, allowed = 1:6))
other_complete <- complete.cases(life_items[-1])
eligible <- not_employed & other_complete
life_total <- rowSums(life_items, na.rm = FALSE)
# 仅“未就业且其余四题完整”者折算；其他缺题仍为 NA。
life_total[eligible] <- rowSums(life_items[eligible, -1, drop = FALSE]) * 5 / 4

## 4. 连续评分：越高越差；除上述未就业情况外，缺任一题即为 NA
phq4 <- c("depress_freq", "disinterest_freq", "restless_freq", "tired_freq")
phq9 <- paste0("phq9_", c("inadequacy", "concentration", "depression", "appetite",
                         "self_harm", "disinterest", "sleep", "psychomotor", "tiredness"))
gad7 <- paste0("gad7_", c("irritability", "anxiety", "uncontrolled_worry", "foreboding",
                         "relaxing", "restlessness", "excessive_worry"))
cidi <- c(paste0("cidi_", c("concentration", "death_thoughts", "disinterest", "depression",
                           "tiredness", "worthlessness", "sleep")), "weight_change")
items$meaning <- 6 - valid(items$happy_life, 1:5)
items$weight_change <- as.numeric(valid(items$cidi_weight, 0:3) > 0)
raw <- list(
  fluid_intelligence = -items$fluid_intelligence,
  reaction_time = items$reaction_time,
  pairs_matching = items$pairs_matching,
  numeric_memory = -items$numeric_memory,
  life_satisfaction = life_total,
  subjective_wellbeing = score(c("happy_general", "happy_health", "meaning"), 1:6),
  phq4 = score(phq4, 1:4, 1),  # S20 四题组合，并非标准临床 PHQ-4
  phq9 = score(phq9, 1:4, 1),
  cidi = score(cidi, 0:1),     # 任何体重变化均记 1
  gad7 = score(gad7, 1:4, 1)
)
continuous <- names(raw)
for (nm in names(raw)) {
  phe[[paste0(nm, ".raw")]] <- raw[[nm]]
  phe[[paste0(nm, ".z")]] <- logit_scale(raw[[nm]])
}

## 5. 二分类评分：1 表示较差
phe$prospective_memory.b <- as.integer(valid(items$prospective_memory, 0:2) != 1)
gp <- valid(items$mental_gp, 0:1)
psy <- valid(items$mental_psychiatrist, 0:1)
phe$general_mental_health.b <- as.integer(gp == 1 | psy == 1)
# 任一 Yes→1；两题均 No→0；其余→NA。

## 6. QC 与保存
score_names <- c(paste0(continuous, ".raw"), paste0(continuous, ".z"), "prospective_memory.b", "general_mental_health.b")
qc <- rbindlist(lapply(score_names, function(nm) {
  x <- phe[[nm]]
  ok <- is.finite(x)
  data.frame(variable = nm, n = sum(ok), missing = sum(!ok), min = if (any(ok)) min(x[ok]) else NA_real_, max = if (any(ok)) max(x[ok]) else NA_real_)
}))
# 折算、CIDI 体重二值化、秩-logit 公式均为实现假设，未经作者原始代码确认。
# 在线问卷不等于招募基线；随访须复用同一参考分布和标准化参数。
attr(phe, "scoring_notes") <- c(
  "Pairs matching: i0_a2; array choice not specified in S20",
  "Complete-item sums; no imputation for skip patterns",
  "Life satisfaction: five-item sum; job=7 and other four complete -> four-item mean*5; implementation assumption, not verified author code",
  "CIDI weight: any change=1; no change=0",
  "Logit: qlogis((average_rank-0.5)/n), then z-score; implementation assumption")
print(qc)

saveRDS(phe, paste0(indir, "/Rdata/behavioral.rds"))


## ------------------------------------------------------------------
## gwas下载/清洗
## ------------------------------------------------------------------
# https://www.nature.com/articles/s41467-026-72529-y
# out <- "D:\data\gwas\brain\raw"