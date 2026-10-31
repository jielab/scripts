#!/usr/bin/env Rscript
# foodwas: preprocess, cohort, patterns, items, groups, trend, rcs, plot, plot_myth
# Launch with dose.sh; source this file to register definitions for reuse.
if (!exists("dose_step_code", mode="function")) {
    source(file.path(Sys.getenv("DOSE_SCRIPT_DIR", if (.Platform$OS.type == "windows")
        "D:/scripts/dose/f" else "/mnt/d/scripts/dose/f"), "0.common.R"), local=TRUE)
}

# ===========================================================================
# foodwas / preprocess
# ===========================================================================
dose_register_step("foodwas", "preprocess", expression({
pacman::p_load(readxl, tidyverse, data.table, survival)

dir0 <- dose_drive()
indir <- dose_ukb_dir()
invisible(lapply(c("phenotype.R", "association.R", "plotting.R"), \(f) source(file.path(dose_helper_dir(), f))))


phase_rds <- function(names) file.path(indir, "Rdata", paste0(names, ".rds"))
pattern_files <- phase_rds(paste0("diet.", c("medi24", "dash", "mind", "hpdi", "ahei", "phdi", "modern", "digm", "maha", "hlcd", "hlfd", "redii", "fds")))
clinical_files <- phase_rds(c("phe", "bbc", "move", "ses", "death", "fod.icd10", "fod.icd9", "fod.opcs4", "fod.srd", "fod.gp", "fod.ref", "pca", "gen", "prs"))
foods_file <- file.path(dose_files_dir(), "foods.xlsx")
upf_columns <- names(read_excel(foods_file, sheet="Sheet2", n_max=0))
upf_files <- phase_rds(c("diet.upf", if (all(c("UPF_2", "Energy_KJ/100g") %in% upf_columns)) "diet.upf_2"))

dose_cached_phase("qc", inputs=c(file.path(indir, "rap/diet.tab.gz"), dose_definition_file("diet")), outputs=phase_rds("diet0"), code={
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 数据QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fn <- fread(dose_definition_file("diet"), sep = "\t", select = 1:5, header = FALSE, fill = TRUE) %>% as.data.frame()
	# Field 20088 is labelled types_sauce in the current dictionary.
	fn$V2[fn$V1 == "p20088"] <- "type_sauce"
	fn$V4 <- as.integer(fn$V4)
	stopifnot(!anyDuplicated(fn$V1), !anyDuplicated(fn$V2))
diet0 <- fread(paste0(indir, "/rap/diet.tab.gz"), sep = "\t", header = TRUE) %>% as.data.frame()
	names(diet0) <- sapply(names(diet0), replace_code, sep = "_", mapping_df = fn[, 1:2]) %>% as.character()
	stopifnot(!anyDuplicated(names(diet0)), !anyDuplicated(diet0$eid))
base_name <- \(x) sub("_i\\d+.*$", "", x)
get_col <- \(x) names(diet0)[base_name(names(diet0)) %in% x]
code_qc <- function(dat, cols, codes) {if (!length(cols)) return(tibble(variable = character(), code = numeric(), n = integer())); map_dfr(cols, \(v) tibble(variable = v, code = codes, n = map_int(codes, \(z) sum(dat[[v]] == z, na.rm = TRUE)))) %>% filter(n > 0)}
date.var <- fn$V2[fn$V1 %in% c("p105010", "p20078", "p105030")]
diet0 <- diet0 %>% mutate(across(all_of(get_col(date.var)), as.Date))
# Touchscreen饮食变量：-10表示“少于1”，按0.5处理，而不是设为NA
less1.var <- c("touch.veg_cooked", "touch.veg_raw", "touch.fruit_fresh", "touch.fruit_dried", "touch.bread", "touch.cereal", "touch.tea", "touch.coffee", "touch.water")
	less1.col <- get_col(less1.var)
	less1.qc <- code_qc(diet0, less1.col, -10)
	diet0 <- diet0 %>% mutate(across(all_of(less1.col), \(x) replace(x, !is.na(x) & x == -10, 0.5)))
# 缺失编码
num.col <- names(diet0)[vapply(diet0, is.numeric, logical(1))]
	neg.code <- c(-1, -2, -3, -7, -121, -818)
	neg.qc <- code_qc(diet0, num.col, neg.code)
	diet0 <- neg2NA(diet0, neg.code)
# 原始WebQ回答中的份量特殊编码；排除连续营养素、重量、性能、体力活动和废弃营养素
cat.no.recode <- c(100117, 100118, 100114, 100108, 100098)
	recode.var <- fn$V2[!is.na(fn$V4) & !fn$V4 %in% cat.no.recode]
	recode.col <- get_col(recode.var); recode.col <- recode.col[vapply(diet0[recode.col], is.numeric, logical(1))]
special.code <- c(200, 300, 400, 500, 600, 444, 555)
	special.qc <- code_qc(diet0, recode.col, special.code)
	diet0 <- diet0 %>% mutate(across(all_of(recode.col), \(x) case_when(x == 200 ~ 2, x == 300 ~ 3, x == 400 ~ 4, x == 500 ~ 5, x == 600 ~ 6, x == 444 ~ 0.25, x == 555 ~ 0.5, TRUE ~ x)))
# “Varied”按Oxford WebQ营养计算规则处理
sugar.var <- c("p100370", "p100490", "p100900")       # sugar：1 tsp
sweetener.var <- c("p100380", "p100500", "p100910")   # sweetener：类推为1 tsp/tablet
half.var <- c("p100260", "p100280", "p100320", "p100350", "p100360", "p100460", "p100470", "p100480", "p100880", "p100890") # milk/decaf/dried fruit：一半有、一半无
meat.var <- c("p103120", "p103130")                   # fat/skin removed
stopifnot(all(c(sugar.var, sweetener.var, half.var, meat.var) %in% base_name(names(diet0))))
sugar.col <- get_col(sugar.var); sugar111.qc <- code_qc(diet0, sugar.col, 111)
sweetener.col <- get_col(sweetener.var); sweetener111.qc <- code_qc(diet0, sweetener.col, 111)
half.col <- get_col(half.var); half111.qc <- code_qc(diet0, half.col, 111)
meat.col <- get_col(meat.var); meat.code.qc <- code_qc(diet0, meat.col, c(111, 222, 313))
diet0 <- diet0 %>% mutate(
	across(all_of(c(sugar.col, sweetener.col)), \(x) replace(x, !is.na(x) & x == 111, 1)),
	across(all_of(half.col), \(x) replace(x, !is.na(x) & x == 111, 0.5)),
	across(all_of(meat.col), \(x) case_when(x %in% c(111, 222) ~ 0.5, x == 313 ~ NA_real_, TRUE ~ x))
)
dose_save_rds(diet0, paste0(indir, "/Rdata/diet0.rds"))
# QC后检查
num.col <- names(diet0)[vapply(diet0, is.numeric, logical(1))]
	qc.after <- bind_rows(code_qc(diet0, less1.col, -10), code_qc(diet0, num.col, neg.code), code_qc(diet0, recode.col, special.code), code_qc(diet0, c(sugar.col, sweetener.col, half.col), 111), code_qc(diet0, meat.col, c(111, 222, 313)))
	if (nrow(qc.after) > 0) print(qc.after, n = Inf); stopifnot(nrow(qc.after) == 0)
code111.review <- code_qc(diet0, setdiff(recode.col, c(sugar.col, sweetener.col, half.col, meat.col)), 111) # 其他原始WebQ变量中不应再有111
	if (nrow(code111.review) > 0) print(code111.review, n = Inf); stopifnot(nrow(code111.review) == 0)
diet.qc <- bind_rows(less1.qc %>% mutate(action = "-10 -> 0.5"), neg.qc %>% mutate(action = "set NA"), special.qc %>% mutate(action = "recode"), sugar111.qc %>% mutate(action = "111 -> 1 tsp"), sweetener111.qc %>% mutate(action = "111 -> 1 tsp/tablet (assumption)"), half111.qc %>% mutate(action = "111 -> 0.5"), meat.code.qc %>% mutate(action = case_when(code %in% c(111, 222) ~ "code -> 0.5", code == 313 ~ "set NA"))) %>% select(action, everything())
# fwrite(diet.qc, "diet.qc.tsv", sep = "\t")
})

dose_cached_phase("medito", inputs=phase_rds("diet0"), outputs=phase_rds("diet.medito"), code={
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 MEDI 改良版地中海（湘雅PMID 40285703), 基线触屏问卷
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# type_milk[1418], type_spread [1428], type_bread [1448], type_cereal [1468]
diet0 <- readRDS(file.path(indir, "Rdata/diet0.rds")) %>% dplyr::select(matches("^(eid$|touch\\.)")) %>% rename_with(~ gsub("touch\\.|_i0$", "", .x))
invisible(gc()) # Release the full QC table before scoring the touchscreen subset.
diet <- diet0 %>% mutate(
	fruit.medito = ifelse((fruit_fresh + fruit_dried/5) >= 3, 1, 0), # ≥3 servings/day
	veg.medito = ifelse((veg_cooked + veg_raw) >= 3, 1, 0), # ≥3 servings/day
	fish.medito = ifelse((fish_oily + fish_other) >= 2, 1, 0), # ≥2 times/week
	meat_processed.medito = ifelse(meat_processed <= 1, 1, 0), # ≤ once a week
	meat_red.medito = ifelse((beef + lamb + pork) <= 2, 1, 0), # ≤2 times/week
	dairy.medito = ifelse(cheese <= 4 & type_milk %in% c(2, 3), 1, 0), # cheese ≤4 times/week
	# oil.medito = ifelse((bread > 2 * 2 * 7 & type_spread %in% c(1,3)), 0, 1),
	oil.medito = ifelse( (((type_spread %in% c(0,2)) | (type_spread == 3 & type_butter %in% c(2, 6, 7, 8))) & bread <= 2 * 2 * 7), 1, 0), # ≤2 servings/day	
	bread.b = ifelse(type_bread == 3, bread, 0), bread.w = ifelse(type_bread %in% c(1, 2, 4), bread, 0),
	cereal.b = ifelse(type_cereal %in% c(1, 3, 4), cereal, 0), cereal.w = ifelse(type_cereal %in% c(2, 5), cereal, 0),
	bread.b.medito = ifelse((bread.b + cereal.b >= 3 * 7), 1, 0), # ≥3 servings/day
	bread.w.medito = ifelse((bread.w + cereal.w <= 2 * 7), 1, 0), # ≤2 servings/day
	sugar.medito = ifelse(grepl("4", never_eat), 1, 0)
) %>% mutate(
	diet.medito.sum = rowSums1(across(matches("\\.medito$"))),
	diet.medito.pts = cut(diet.medito.sum, breaks = c(-Inf, 2, 4, 6, 8, Inf), labels = c(0, 25, 50, 80, 100), right = FALSE) %>% as.character() %>% as.numeric()
) %>% dplyr::select(eid, diet.medito.sum, diet.medito.pts, matches("\\.medito")) 
table(diet$diet.medito.pts, useNA = "always")
dose_save_rds(diet, paste0(indir, "/Rdata/diet.medito.rds"))
rm(diet0, diet)
invisible(gc())
})

dose_cached_phase("webq", inputs=c(phase_rds("diet0"), foods_file), outputs=phase_rds("diet.fudan"), code={
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 复旦（PMID 40603580) 修订版, 24小时饮食回顾（Oxford WebQ）数据处理 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1）纳入人群
inst <- 0:4 # *_i0 到 *_i4
field <- read_excel(file.path(dose_files_dir(), "foods.xlsx"), sheet = "Sheet2")
diet0 <- readRDS(paste0(indir, "/Rdata/diet0.rds")) %>% rename_with(~ sub("^touch\\.vari", "vari", .x), matches("^touch\\.vari")) %>%
	select(!matches("^touch\\.")) %>% filter(!is.na(p20077)) # 👫 原始有Oxford WebQ记录的人群 210739
# 填写时长：仅做QC描述，不作为参与者或单次recall的排除条件
dura.col <- paste0("p20082_i", inst) # 填写时间低于第5百分位或超过24小时的单次问卷
	dura <- diet0 %>% select(eid, all_of(dura.col))
	q <- quantile(as.matrix(dura[, -1]), na.rm = TRUE, probs = c(0.05, 0.95)); cat("WebQ填写时长P5/P95（分钟）:",paste(round(q, 1), collapse = " / "), "\n")
# dura <- dura %>% mutate(across(all_of(dura.col), ~ ifelse(.x >= q[1] & .x <= 60*24, .x, NA_real_)), enroll = rowSums(!is.na(pick(all_of(dura.col))))); table(dura$enroll)
# 排除 Field 1548 回答 Often（编码3）的参与者
# vari <- diet0 %>% select(eid, matches("^vari")) %>% mutate(vari = as.integer(rowSums(across(-eid, ~ .x == 3), na.rm = TRUE) > 0)); table(vari$vari, useNA = "ifany") 
# keep <- dura %>% inner_join(vari %>% select(eid, vari), by = "eid") %>% filter(enroll >= 1, vari == 0) 
# diet <- diet0 %>% filter(eid %in% keep$eid) %>% select(-all_of(dura.col)) %>% inner_join(keep %>% select(eid, enroll, vari, all_of(dura.col)), by = "eid") 

vari <- diet0 %>% # 🏮仅使用基线 Field 1548:vari_i0
	transmute(eid, vari = as.integer(!is.na(vari_i0) & vari_i0 == 3)); table(vari$vari, useNA = "ifany")
keep <- vari %>% filter(vari == 0)
diet <- diet0 %>% inner_join(keep, by = "eid")
rm(diet0, dura, vari, keep)
invisible(gc())
cat("不排除填写时长、仅排除基线饮食经常变化者：", nrow(diet), "\n") # 👫194434
	
# 有效recall = 填写时长有效 + 能量合理；无效recall整体作废
energy.col <- paste0("p26002_i", inst); date.col <- paste0("p105010_i", inst); valid.col <- paste0("valid_i", inst)
diet <- diet %>% mutate(across(all_of(date.col), as.Date))
for (j in seq_along(inst)) {
	energy <- diet[[energy.col[j]]]
	ok <- !is.na(diet[[dura.col[j]]]) & !is.na(energy) & ((diet$sex == 0 & energy >= 500*4.184 & energy <= 3500*4.184) | (diet$sex == 1 & energy >= 800*4.184 & energy <= 4200*4.184))
	ok[is.na(ok)] <- FALSE
	diet[[valid.col[j]]] <- ok; diet[[energy.col[j]]][!ok] <- NA_real_; diet[[date.col[j]]][!ok] <- as.Date(NA)
}
diet <- diet %>% mutate( # 根据最终有效recall重新计算次数、平均能量和起始日期
	enroll = rowSums(pick(all_of(valid.col))),  # number of dietary assessments
	energy.kJ = rowMeans1(pick(all_of(energy.col)), per = 0.8),
	start_date = do.call(pmax, c(pick(all_of(date.col)), list(na.rm = TRUE))) # 最后一次有效recall日期，作为平均饮食暴露的随访起点
)
diet$start_date[rowSums(!is.na(diet[date.col])) == 0] <- as.Date(NA)
diet <- diet %>% filter(enroll >= 1, !is.na(energy.kJ), !is.na(start_date))
cat("最终有效WebQ饮食人群：", nrow(diet), "\n")  # 👫187372 
table(diet$enroll)

# 2）食物处理
# 🍞 面包、🥥油、🥣汤和🥫Sauce分类
breads <- c("bread", "baguette", "bap", "roll")
soups <- c("canned", "homemade"); soup_code <- c("1"="fish", "2"="meat", "3"="pasta", "4"="pulse", "5"="vegetable", "6"="other")
sauce_code <- list(sauce.Jamhoney=332, sauce.cream=333, sauce.pnutbutter=334, sauce.yeast=335, sauce.hummus=336, sauce.guacamole=337, sauce.chutney=338, sauce.ketchup=339,sauce.brown=340, sauce.mayonnaise=c(341,342), sauce.saladdressing=343, sauce.olives=344, sauce.pesto=345, sauce.tomato=346, sauce.cheese=347, sauce.white=348, sauce.gravy=349, sauce.other=350)
for (k in inst) {
	# 🍞 面包分类：多个类型等比例分配,代码1、5为w,代码2–4为b
	for (f in breads) {
		type_col <- sprintf("type_%s_i%d", f, k); intake_col <- sprintf("%s_i%d", f, k); type_x <- as.character(diet[[type_col]]); intake_x <- as.numeric(diet[[intake_col]])
		valid <- !is.na(type_x) & str_detect(type_x, "^[1-5](\\|[1-5])*$")
		n_all <- ifelse(valid, str_count(type_x, fixed("|")) + 1, NA_real_); n_w <- ifelse(valid, str_count(type_x, "[15]"), NA_real_); n_b <- ifelse(valid, str_count(type_x, "[2-4]"), NA_real_)
		diet[[sprintf("%s.w_i%d", f, k)]] <- ifelse(is.na(intake_x), NA_real_, ifelse(valid, intake_x*n_w/n_all, NA_real_)); diet[[sprintf("%s.b_i%d", f, k)]] <- ifelse(is.na(intake_x), NA_real_, ifelse(valid, intake_x*n_b/n_all, NA_real_)); diet[[sprintf("%s.unanswered_i%d", f, k)]] <- ifelse(is.na(intake_x), NA_real_, ifelse(!valid, intake_x, 0))
	}
	# 🥥 油分为好坏：选择354–357中的任一种，记为10 g
	raw_col <- sprintf("type_oil_i%d", k); oil_x <- as.character(diet[[raw_col]])
		diet[[sprintf("oil_i%d", k)]] <- ifelse(is.na(oil_x), NA_real_, ifelse(str_detect(oil_x, "(^|\\|)35[4-7](\\||$)"), 10, 0))
	# 🥣 汤：多个配料等比例分配
	for (f in soups) {
		type_col <- sprintf("type_%s_i%d", f, k); intake_col <- sprintf("%s_i%d", f, k); type_x <- as.character(diet[[type_col]]); intake_x <- as.numeric(diet[[intake_col]])
		valid <- !is.na(type_x) & str_detect(type_x, "^[1-6](\\|[1-6])*$")
		n_all <- ifelse(valid, str_count(type_x, fixed("|")) + 1, NA_real_)
		for (cc in names(soup_code)) {
			selected <- str_detect(type_x, paste0("(^|\\|)", cc, "(\\||$)"))
			diet[[sprintf("%s.%s_i%d", f, soup_code[[cc]], k)]] <- ifelse(is.na(intake_x), NA_real_, ifelse(valid & selected, intake_x/n_all, 0))
		}
		diet[[sprintf("%s.unanswered_i%d", f, k)]] <- ifelse(is.na(intake_x), NA_real_, ifelse(!valid, intake_x, 0))
	}
	# 🥫 Sauce：每个选中的类型记为1
	raw_col <- sprintf("type_sauce_i%d", k); sauce_x <- as.character(diet[[raw_col]])
	for (nm in names(sauce_code)) {
		codes <- sauce_code[[nm]]
		pat <- paste0("(^|\\|)(", paste(codes, collapse="|"), ")(\\||$)")
		diet[[sprintf("%s_i%d", nm, k)]] <- ifelse(is.na(sauce_x), NA_real_, as.numeric(str_detect(sauce_x, pat)))
	}
}
# diet <- diet %>% select(-matches(paste0("^(", paste(c(breads), collapse = "|"), ")_i")), -matches("^type"))
food.col <- unlist(lapply(field$name, \(x) paste0(x, "_i", inst)), use.names=FALSE)
dat1 <- diet %>% select(eid, all_of(food.col))
for (k in inst) { # 有效recall中未回答食物设为0；无效recall全部设为NA
	valid_k <- diet[[paste0("valid_i", k)]]; idx_k <- grep(paste0("_i", k, "$"), names(dat1))
	dat1[idx_k] <- lapply(dat1[idx_k], \(x) {x <- suppressWarnings(as.numeric(as.character(x))); ifelse(!valid_k, NA_real_, ifelse(is.na(x), 0, x))})
}
dat.food <- dat1 %>% select(eid) # 多次recall取均值
	for (i in field$name) dat.food[[i]] <- rowMeans1(dat1[paste0(i, "_i", inst)], per=0.8)
wt_vars <- strsplit("p26150 p26131 p26133 p26096 p26102 p26154 p26110 oil p26062 p26063 p26111 p26112 p26064 p26153 p26152 p26151 p26067 p26138", " ")[[1]] # 原始单位为重量的food items转换为portion数
	portion <- field %>% filter(name %in% wt_vars) %>% transmute(name, Portion=as.numeric(Portion)) %>% filter(!is.na(Portion), Portion > 0) %>% tibble::deframe()
dat.food <- dat.food %>% mutate(across(all_of(names(portion)), ~ .x/portion[cur_column()])); 
length(names(dat.food))-1 # 222

# 🚩 Oxford WebQ食物重量和营养素
field.w <- read_excel(file.path(dose_files_dir(), "foods.xlsx"), sheet="Sheet3")
dat1 <- diet %>% select(eid, all_of(unlist(lapply(field.w$name, \(x) paste0(x, "_i", inst)), use.names=FALSE)))
weight_names <- field.w %>% filter(Category != "nutrient") %>% pull(name) %>% unique() # 93
nutrient_names <- field.w %>% filter(Category == "nutrient") %>% pull(name) %>% unique() # 63
for (k in inst) {
	valid_k <- diet[[paste0("valid_i", k)]]; idx_k <- grep(paste0("_i", k, "$"), names(dat1)); base_k <- sub("_i\\d+$", "", names(dat1)[idx_k])
	# 食物组重量：有效recall中NA -> 0；无效recall全部 -> NA
	idx_wk <- idx_k[base_k %in% weight_names]
	if (length(idx_wk)) dat1[idx_wk] <- lapply(dat1[idx_wk], \(x) {x <- suppressWarnings(as.numeric(as.character(x))); ifelse(!valid_k, NA_real_, ifelse(is.na(x), 0, x))})
	# 营养素：无效recall -> NA；有效recall中的原始NA保持NA
	idx_nk <- idx_k[base_k %in% nutrient_names]
	if (length(idx_nk)) dat1[idx_nk] <- lapply(dat1[idx_nk], \(x) {x <- suppressWarnings(as.numeric(as.character(x))); ifelse(!valid_k, NA_real_, x)})
}
dat.weight <- dat1 %>% select(eid)  #🏮取均值
for (i in field.w$name)
	dat.weight[[paste0("w.", i)]] <- rowMeans1(dat1[paste0(i, "_i", inst)], per=0.8)
stopifnot(identical(dat.weight$eid, diet$eid))
length(names(dat.weight))-1 # ✔ 156
dat.food <- dat.food %>% left_join(dat.weight, by="eid") %>% left_join(diet %>% transmute(eid, sex, weight=weight_i0, energy.kJ, start_date, enroll, vari), by="eid")
dose_save_rds(dat.food, paste0(indir, "/Rdata/diet.fudan.rds"))
rm(diet, dat1, dat.weight)
invisible(gc())
message("Food/weight recall averages saved; scoring dietary patterns")
})

dose_cached_phase("patterns", inputs=phase_rds("diet.fudan"), outputs=pattern_files, code={
dat.food <- readRDS(file.path(indir, "Rdata/diet.fudan.rds"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 🧠MODERN 用了复旦文章的代码
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field <- field %>% filter(MODERN > 0)
cols <- 'Berries_and_Citrus Green_leafy_vegetables Olive_oil Potatoes Eggs Poultry Sweetened_beverages' %>% strsplit(' ') %>% unlist() 
modern_map <- list(  # 7个组
	Berries_and_Citrus = c("p104470", "p104480", "p104490", "p104530", "p104540"),
	Green_leafy_vegetables = c("p104140", "p104160", "p104240", "p104300", "p104310", "p104370"),
	Olive_oil = c("p26110"),
	Potatoes = c("p104030", "p104050"),
	Eggs = c("p102940", "p102950", "p102960", "p102970", "p102980"),
	Poultry = c("p103060"),
	Sweetened_beverages = c("p100170", "p100180", "p100220", "p100230", "p100550", "p100160", "p100540")
)
dat1 <- dat.food %>% select(eid, all_of(unique(unlist(modern_map))))
for (g in names(modern_map)) {dat1[[g]] <- rowSums1(dat1[, modern_map[[g]], drop = FALSE])}
dat.modern <- dat1 %>% select(eid, all_of(cols)) %>% transmute(eid,
	Berries_and_Citrus.modern	  = ifelse(Berries_and_Citrus > 0 & Berries_and_Citrus <= 2, 1, 0),
	Green_leafy_vegetables.modern = ifelse(Green_leafy_vegetables > 0.25, 1, 0),
	Olive_oil.modern			  = ifelse(Olive_oil > 0, 1, 0),
	Potatoes.modern 			  = ifelse(Potatoes > 0 & Potatoes <= 0.75, 1, 0),
	Eggs.modern					  = ifelse(Eggs > 0 & Eggs <= 1, 1, 0),
	Poultry.modern 				  = ifelse(Poultry > 0 & Poultry <= 0.5, 1, 0),
	Sweetened_beverages.modern 	  = ifelse(Sweetened_beverages == 0, 1, 0)
) %>% mutate(diet.modern.sum = rowSums1(across(ends_with(".modern")))) %>% select(eid, diet.modern.sum, ends_with(".modern"))
dose_save_rds(dat.modern, paste0(indir, "/Rdata/diet.modern.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 🧠MIND (参考文献 PMID: 41245532)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field <- field %>% filter (MIND > 0)
cols <- 'Green_leafy_vegetables Other_vegetables Berries Olive_oil Red_and_processed_meat Poultry Butter_and_margarine Cheese_not_low_fat Wine Legumes Nuts Seafood_not_fried Sweets_desserts_and_sugary_drinks Whole_grain Fast_fried_foods' %>% strsplit(' ') %>% unlist()
mind_map <- list( # 15个组 0–15
	Green_leafy_vegetables = c("p104080", "p104160", "p104240", "p104300", "p104370"), # p104090 *0.5
	Other_vegetables = c("p104060", "p104070", "p104100", "p104120", "p104130", "p104140", "p104150", "p104170", "p104180", "p104190", "p104200", "p104210", "p104220", "p104230", "p104250", "p104260", "p104270", "p104280", "p104290", "p104310", "p104320", "p104330", "p104340", "p104350", "p104360", "p104380", "homemade.vegetable", "canned.vegetable", "sauce.guacamole", "sauce.tomato"), # p104090 *0.5
	Berries = c("p104470"),
	Olive_oil = c("p26110"),
	Red_and_processed_meat = c("p103030", "p103020", "p103040", "p103010", "p103070", "p103080", "p103090"),
	Poultry = c("p103060"),
	Butter_and_margarine = c("p26062", "p26063", "p26111", "p26112", "p26154"),
	Cheese_not_low_fat = c("p102820", "p102830", "p102840", "p102860", "p102880", "p102890", "p102900", "p102910"),
	Wine = c("p26153", "p26152", "p26151"),
	Legumes = c("p104110", "p104000", "p104010", "p103270", "sauce.hummus", "canned.pulse", "homemade.pulse"),
	Nuts = c("p102440", "p102430", "p102420", "p102410", "sauce.pnutbutter", "p102450"),
	Seafood_not_fried = c("p103160", "p103190", "p103150", "p103230", "p103200", "p103210", "p103220", "canned.fish", "homemade.fish"),
	Sweets_desserts_and_sugary_drinks = c("p102350", "p102340", "p102260", "p102310", "p102300", "p102290", "p102280", "p102270", "p102360", "p102330", "p102320", "p102380", "p102190", "p102220", "p102200", "p102180", "p102060", "p102210", "p102140", "p102150", "p102230", "p102170", "p102120", "p102370", "p101970", "p101980", "p101990", "p102010", "p102020", "p102030", "p102050", "p102070", "p100530", "p100550", "p100170", "p100180"),
	Whole_grain = c("p100770", "p100800", "p100810", "p100840", "p100850", "bread.b", "baguette.b", "bap.b", "roll.b", "p101250", "p101260", "p102720", "p102740", "p102770", "p102780"),
	Fast_fried_foods = c("p103050", "p103170", "p103180", "p104020", "p102460")
)
dat1 <- dat.food %>% select(eid, all_of(unique(c(unlist(mind_map, use.names = FALSE), "p104090"))))
for (g in names(mind_map)) {dat1[[g]] <- rowSums1(dat1[, mind_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% mutate(
	Green_leafy_vegetables = Green_leafy_vegetables + 0.5 * p104090, 
	Other_vegetables = Other_vegetables + 0.5 * p104090
) %>% select(eid, all_of(cols))
weekly_cols <- 'Green_leafy_vegetables Other_vegetables Berries Red_and_processed_meat Poultry Cheese_not_low_fat Wine Legumes Nuts Seafood_not_fried Sweets_desserts_and_sugary_drinks Fast_fried_foods' %>% strsplit(' ') %>% unlist()
dat1 <- dat1 %>% mutate(across(all_of(weekly_cols), ~ .x * 7))
dat.mind <- dat1 %>% transmute(eid,
	Green_leafy_vegetables.mind = case_when(Green_leafy_vegetables <= 2 ~ 0, Green_leafy_vegetables >= 6 ~ 1, TRUE ~ 0.5),
	Other_vegetables.mind = case_when(Other_vegetables < 5 ~ 0, Other_vegetables >= 7 ~ 1, TRUE ~ 0.5),
	Berries.mind = case_when(Berries < 1 ~ 0, Berries >= 2 ~ 1, TRUE ~ 0.5),
	Olive_oil.mind = ifelse(Olive_oil > 0, 1, 0),
	Red_and_processed_meat.mind = case_when(Red_and_processed_meat > 6 ~ 0, Red_and_processed_meat < 4 ~ 1, TRUE ~ 0.5),
	Poultry.mind = case_when(Poultry < 1 ~ 0, Poultry > 2 ~ 1, TRUE ~ 0.5),
	Butter_and_margarine.mind = case_when(Butter_and_margarine > 2 ~ 0, Butter_and_margarine < 1 ~ 1, TRUE ~ 0.5),
	Cheese_not_low_fat.mind = case_when(Cheese_not_low_fat >= 6 ~ 0, Cheese_not_low_fat < 1 ~ 1, TRUE ~ 0.5),
	Wine.mind = case_when(Wine < 0.25 | Wine > 7 ~ 0, dplyr::near(Wine, 7) ~ 1, TRUE ~ 0.5),
	Legumes.mind = case_when(Legumes < 1 ~ 0, Legumes > 3 ~ 1, TRUE ~ 0.5),
	Nuts.mind = case_when(Nuts < 0.5 ~ 0, Nuts > 4 ~ 1, TRUE ~ 0.5),
	Seafood_not_fried.mind = case_when(Seafood_not_fried < 0.25 ~ 0, Seafood_not_fried >= 1 ~ 1, TRUE ~ 0.5),
	Sweets_desserts_and_sugary_drinks.mind = case_when(Sweets_desserts_and_sugary_drinks >= 7 ~ 0, Sweets_desserts_and_sugary_drinks < 5 ~ 1, TRUE ~ 0.5),
	Whole_grain.mind = case_when(Whole_grain < 1 ~ 0, Whole_grain > 2 ~ 1, TRUE ~ 0.5),
	Fast_fried_foods.mind = case_when(Fast_fried_foods > 3 ~ 0, Fast_fried_foods < 1 ~ 1, TRUE ~ 0.5)
) %>% mutate(diet.mind.sum = rowSums(across(ends_with(".mind")))) %>% select(eid, diet.mind.sum, ends_with(".mind"))
dose_save_rds(dat.mind, paste0(indir, "/Rdata/diet.mind.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Medi Pattern score (参考文献 PMID: 41245532)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field <- field %>% filter (MedDiet > 0) 
score_up <- function(x, high) {x <- as.numeric(x); pmin(pmax(x / high, 0), 1)} # 正向组分：0摄入得0分，达到推荐摄入量得1分，中间线性计分
score_down <- function(x, low, high) {x <- as.numeric(x); pmin(pmax((high - x) / (high - low), 0), 1)} # 反向组分：不超过low得1分，达到或超过high得0分，中间线性计分
cols <- 'Vegetables Fruit_including_fruit_juice Olive_oil Red_and_processed_meat White_meat Butter_and_margarine Sugary_drinks Wine Legumes Nuts Seafood Sweets_and_desserts Sofrito' %>% strsplit(' ') %>% unlist()
medi_map <- list(# 13个组 0–13
	Vegetables = c("p104060", "p104070", "p104080", "p104090", "p104100", "p104120", "p104130", "p104140", "p104150", "p104160", "p104170", "p104180", "p104190", "p104200", "p104210", "p104220", "p104230", "p104240", "p104250", "p104260", "p104270", "p104280", "p104290", "p104300", "p104310", "p104320", "p104330", "p104340", "p104350", "p104360", "p104370", "p104380", "canned.vegetable", "homemade.vegetable", "sauce.guacamole"),
	Fruit_including_fruit_juice = c("p104410", "p104420", "p104430", "p104440", "p104450", "p104460", "p104470", "p104480", "p104490", "p104500", "p104510", "p104520", "p104530", "p104540", "p104550", "p104560", "p104570", "p104580", "p104590", "p102490", "p100200", "p100190", "p100210"),
	Olive_oil = c("p26110"),
	Red_and_processed_meat = c("p103030", "p103020", "p103040", "p103010", "p103070", "p103080", "p103090"),
	White_meat = c("p103060", "p103050"),
	Butter_and_margarine = c("p26062", "p26063", "p26111", "p26112", "p26154"),
	Sugary_drinks = c("p100170", "p100160", "p100180"),
	Wine = c("p26153", "p26152"),
	Legumes = c("p104110", "p104000", "p104010", "p103270", "sauce.hummus", "canned.pulse", "homemade.pulse"),
	Nuts = c("p102410", "p102430", "p102420", "p102440", "p102450", "sauce.pnutbutter"),
	Seafood = c("p103170", "p103180", "p103160", "p103190", "p103150", "p103230", "p103200", "p103210", "p103220", "canned.fish", "homemade.fish"),
	Sweets_and_desserts = c("p102350", "p102340", "p102260", "p102310", "p102300", "p102290", "p102280", "p102270", "p102360", "p102330", "p102320", "p102380", "p102190", "p102220", "p102200", "p102180", "p102060", "p102210", "p102140", "p102150", "p102230", "p102170", "p102120", "p102370", "p101970", "p101980", "p101990", "p102010", "p102020", "p102030", "p102050", "p102070"),
	Sofrito = c("sauce.tomato")
)
dat1 <- dat.food %>% select(eid, all_of(unique(unlist(medi_map, use.names = FALSE))))
for (g in names(medi_map)) {dat1[[g]] <- rowSums1(dat1[, medi_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% select(eid, all_of(cols))
weekly_cols <- 'Wine Legumes Nuts Seafood Sweets_and_desserts Sofrito' %>% strsplit(' ') %>% unlist()
dat1 <- dat1 %>% mutate(across(all_of(weekly_cols), ~ .x * 7))
dat.medi24 <- dat1 %>% transmute(eid,
	Vegetables.medi = score_up(Vegetables, high = 2),
	Fruit_including_fruit_juice.medi = score_up(Fruit_including_fruit_juice, high = 3),
	Olive_oil.medi = as.numeric(Olive_oil > 0),
	Red_and_processed_meat.medi = score_down(Red_and_processed_meat, low = 1, high = 2),
	White_meat.medi = as.numeric(White_meat > Red_and_processed_meat),
	Butter_and_margarine.medi = score_down(Butter_and_margarine, low = 1, high = 2),
	Sugary_drinks.medi = score_down(Sugary_drinks, low = 1, high = 2),
	Wine.medi = score_up(Wine, high = 7),
	Legumes.medi = score_up(Legumes, high = 3),
	Nuts.medi = score_up(Nuts, high = 3),
	Seafood.medi = score_up(Seafood, high = 3),
	Sweets_and_desserts.medi = score_down(Sweets_and_desserts, low = 2, high = 4),
	Sofrito.medi = score_up(Sofrito, high = 2)
) %>% mutate(diet.medi24.sum = rowSums1(across(ends_with(".medi")))) %>% select(eid, diet.medi24.sum, ends_with(".medi"))
dose_save_rds(dat.medi24, paste0(indir, "/Rdata/diet.medi24.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 DASH (同上，参考文献 PMID: 41245532)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field <- field %>% filter (DASH > 0)
score_q5_zero <- function(x, reverse=FALSE){x <- as.numeric(x); out <- rep(NA_integer_, length(x)); i0 <- !is.na(x)&x==0; i1 <- !is.na(x)&x>0; out[i0] <- 1L; out[i1] <- dplyr::ntile(x[i1], 4)+1L; if(reverse) out <- 6L-out; out} # 对大量 0 的变量，用“0 单独一档 + 其余做分位”
cols <- 'Vegetables Fruit_including_fruit_juice Red_and_processed_meat Nuts_and_legumes Whole_grains Low_fat_dairy Sugary_drinks Sodium' %>% strsplit(' ') %>% unlist()
dash_map <- list( # 8组 8–39
	Vegetables = c("p104060", "p104070", "p104080", "p104090", "p104100", "p104120", "p104130", "p104140", "p104150", "p104160", "p104170", "p104180", "p104190", "p104200", "p104210", "p104220", "p104230", "p104240", "p104250", "p104260", "p104270", "p104280", "p104290", "p104300", "p104310", "p104320", "p104330", "p104340", "p104350", "p104360", "p104370", "p104380", "canned.vegetable", "homemade.vegetable", "sauce.guacamole", "sauce.tomato"),
	Fruit_including_fruit_juice = c("p104410", "p104420", "p104430", "p104440", "p104450", "p104460", "p104470", "p104480", "p104490", "p104500", "p104510", "p104520", "p104530", "p104540", "p104550", "p104560", "p104570", "p104580", "p104590", "p102490", "p100200", "p100190", "p100210"),
	Red_and_processed_meat = c("p103030", "p103020", "p103040", "p103010", "p103070", "p103080", "p103090"),
	Nuts_and_legumes = c("p104110", "p104000", "p104010", "p103270", "sauce.hummus", "canned.pulse", "homemade.pulse", "p102410", "p102430", "p102420", "p102440", "p102450", "sauce.pnutbutter"),
	Whole_grains = c("p100770", "p100800", "p100810", "p100840", "p100850", "bread.b", "baguette.b", "bap.b", "roll.b", "p101250", "p101260", "p102720", "p102740", "p102770", "p102780"),
	Low_fat_dairy = c("p102850", "p102810", "p102870", "p26131", "p26133", "p26102"),
	Sugary_drinks = c("p100170", "p100160", "p100180"),
	Sodium = c("w.p26052")
)
dat1 <- dat.food %>% select(eid, all_of(unique(unlist(dash_map, use.names = FALSE))))
for (g in names(dash_map)) {dat1[[g]] <- rowSums1(dat1[, dash_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% select(eid, all_of(cols))
dat.dash <- dat1 %>% mutate(
	across(c(Vegetables, Fruit_including_fruit_juice, Nuts_and_legumes, Whole_grains, Low_fat_dairy), ~ score_q5_zero(.x, reverse = FALSE), .names = "{.col}.dash"),
	across(c(Red_and_processed_meat, Sugary_drinks, Sodium), ~ score_q5_zero(.x, reverse = TRUE), .names = "{.col}.dash")
) %>% mutate(diet.dash.sum = rowSums1(across(ends_with(".dash")))) %>% select(eid, diet.dash.sum, ends_with(".dash"))
dose_save_rds(dat.dash, paste0(indir, "/Rdata/diet.dash.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 AHEI-2010 (同上，参考文献 PMID: 41245532)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field <- field %>% filter (AHEI > 0)
score_up10 <- function(x, high) {x <- as.numeric(x); pmin(pmax(x / high, 0), 1) * 10} # 0摄入为0分，达到或超过推荐摄入量为10分，中间线性计分
score_down10 <- function(x, high) {x <- as.numeric(x); pmin(pmax(1 - x / high, 0), 1) * 10} # 0摄入为10分，达到或超过上限为0分，中间线性计分
score_mid10 <- function(x, lo, hi, high0) {x <- as.numeric(x); ifelse(is.na(x), NA_real_, ifelse(x >= high0, 0, ifelse(x >= lo & x <= hi, 10, ifelse(x < lo, pmin(pmax(x / lo, 0), 1) * 10, pmin(pmax((high0 - x) / (high0 - hi), 0), 1) * 10))))} # 中等摄入最佳：lo-hi得10分；0摄入及达到high0时得0分；两侧线性计分
cols <- 'Vegetables Fruit_not_fruit_juice Red_and_processed_meat Nuts_and_legumes Sugary_drinks_and_fruit_juices Whole_grains Trans_fat Long_chain_n3_fats Polyunsaturated_fatty_acids Sodium Alcohol' %>% strsplit(' ') %>% unlist()
ahei_map <- list( # 11组 0–110
	Vegetables = c("p104060", "p104070", "p104080", "p104090", "p104100", "p104120", "p104130", "p104140", "p104150", "p104160", "p104170", "p104180", "p104190", "p104200", "p104210", "p104220", "p104230", "p104240", "p104250", "p104260", "p104270", "p104280", "p104290", "p104300", "p104310", "p104320", "p104330", "p104340", "p104350", "p104360", "p104370", "p104380", "canned.vegetable", "homemade.vegetable", "sauce.guacamole", "sauce.tomato"),
	Fruit_not_fruit_juice = c("p104410", "p104420", "p104430", "p104440", "p104450", "p104460", "p104470", "p104480", "p104490", "p104500", "p104510", "p104520", "p104530", "p104540", "p104550", "p104560", "p104570", "p104580", "p104590", "p102490"),
	Red_and_processed_meat = c("p103030", "p103020", "p103040", "p103010", "p103070", "p103080", "p103090"),
	Nuts_and_legumes = c("p104110", "p104000", "p104010", "p103270", "sauce.hummus", "canned.pulse", "homemade.pulse", "p102410", "p102430", "p102420", "p102440", "p102450", "sauce.pnutbutter"),
	Sugary_drinks_and_fruit_juices = c("p100190", "p100200", "p100210", "p100220", "p100170", "p100180", "p100530", "p100550"),
	Whole_grains = c("w.p26076", "w.p26077", "w.p26074", "w.p26114"), # w.p26105 × 0.5
	Trans_fat = c("w.p26155"),
	Long_chain_n3_fats = c("p103160"),
	Polyunsaturated_fatty_acids = c("w.p26015", "w.p26016"),
	Sodium = c("w.p26052"),
	Alcohol = c("p26153", "p26152", "p26151", "p26067", "p26138")
)
dat1 <- dat.food %>% select(eid, sex, energy.kJ, all_of(unique(c(unlist(ahei_map, use.names = FALSE), "w.p26105"))))
for (g in names(ahei_map)) {dat1[[g]] <- rowSums1(dat1[, ahei_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% mutate(
	Whole_grains = Whole_grains + 0.5 * w.p26105,
	energy_kcal = energy.kJ * 0.239,
	Trans_fat = ifelse(energy_kcal > 0, Trans_fat * 9 / energy_kcal * 100, NA_real_),
	Polyunsaturated_fatty_acids = ifelse(energy_kcal > 0, Polyunsaturated_fatty_acids * 9 / energy_kcal * 100, NA_real_)
) %>% select(eid, all_of(cols), sex)
dat.ahei <- dat1 %>% transmute(eid,
	Vegetables.ahei = score_up10(Vegetables, high = 5),
	Fruit_not_fruit_juice.ahei = score_up10(Fruit_not_fruit_juice, high = 4),
	Red_and_processed_meat.ahei = score_down10(Red_and_processed_meat, high = 1.5),
	Nuts_and_legumes.ahei = score_up10(Nuts_and_legumes, high = 1),
	Sugary_drinks_and_fruit_juices.ahei = score_down10(Sugary_drinks_and_fruit_juices, high = 1),
	Whole_grains.ahei = ifelse(sex == 0, score_up10(Whole_grains, high = 75), score_up10(Whole_grains, high = 90)),
	Trans_fat.ahei = pmin(pmax((4 - Trans_fat) / (4 - 0.5), 0), 1) * 10,
	Long_chain_n3_fats.ahei = score_up10(Long_chain_n3_fats, high = 2/7),
	Polyunsaturated_fatty_acids.ahei = pmin(pmax((Polyunsaturated_fatty_acids - 2) / (10 - 2), 0), 1) * 10,
	Sodium.ahei = {dec <- dplyr::ntile(Sodium, 10); 10 * (10 - dec) / 9},
	Alcohol.ahei = ifelse(sex == 0, score_mid10(Alcohol, lo = 0.5, hi = 1.5, high0 = 2.5), score_mid10(Alcohol, lo = 0.5, hi = 2.0, high0 = 3.5))
) %>% mutate(diet.ahei.sum = rowSums1(across(ends_with(".ahei")))) %>% select(eid, diet.ahei.sum, ends_with(".ahei"))
dose_save_rds(dat.ahei, paste0(indir, "/Rdata/diet.ahei.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 healthful Plant-Based Index (hPDI) (同上，参考文献 PMID: 41245532)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field <- field %>% filter (hPDI > 0)
cols <- 'Vegetables Fruit_not_fruit_juice Legumes_and_vegetarian_protein_alternatives Nuts Whole_grains Tea_and_coffee Refined_grains Potatoes Sugary_drinks Fruit_juices Sweets_and_desserts Animal_fat Dairy Eggs Seafood Meat Miscellaneous_animal_based_foods' %>% strsplit(' ') %>% unlist()
hpdi_map <- list( # 17组，17–85分
	Vegetables = c("p104060", "p104070", "p104080", "p104090", "p104100", "p104120", "p104130", "p104140", "p104150", "p104160", "p104170", "p104180", "p104190", "p104200", "p104210", "p104220", "p104230", "p104240", "p104250", "p104260", "p104270", "p104280", "p104290", "p104300", "p104310", "p104320", "p104330", "p104340", "p104350", "p104360", "p104370", "p104380", "canned.vegetable", "homemade.vegetable", "sauce.guacamole", "sauce.tomato"),
	Fruit_not_fruit_juice = c("p104410", "p104420", "p104430", "p104440", "p104450", "p104460", "p104470", "p104480", "p104490", "p104500", "p104510", "p104520", "p104530", "p104540", "p104550", "p104560", "p104570", "p104580", "p104590", "p102490"),
	Legumes_and_vegetarian_protein_alternatives = c("p104110", "p104000", "p104010", "p103270", "sauce.hummus", "canned.pulse", "homemade.pulse", "p103260", "p103280", "p103290"),
	Nuts = c("p102440", "p102430", "p102420", "p102410", "sauce.pnutbutter", "p102450"),
	Whole_grains = c("p100770", "p100800", "p100810", "p100840", "p100850", "bread.b", "baguette.b", "bap.b", "roll.b", "p101250", "p101260", "p102720", "p102740", "p102770", "p102780"),
	Tea_and_coffee = c("p100250", "p100270", "p100290", "p100300", "p100310", "p100330", "p100400", "p100410", "p100420", "p100430", "p100440"),
	Refined_grains = c("p100820", "p100830", "bread.w", "baguette.w", "bap.w", "roll.w", "p101230", "p101240", "p102710", "p102730", "p102010", "p102020", "p102050", "p102070", "p102470", "p102500", "p102760"),
	Potatoes = c("p104020", "p104030", "p104050"),
	Sugary_drinks = c("p100170", "p100160", "p100180"),
	Fruit_juices = c("p100190", "p100200", "p100210"),
	Sweets_and_desserts = c("p102350", "p102340", "p102260", "p102310", "p102300", "p102290", "p102280", "p102270", "p102360", "p102330", "p102320", "p102380", "p102190", "p102180", "p102200", "p102060", "p102210", "p102230", "p102170", "p102370", "p101970", "p101980", "p101990", "p102020", "p102030"),
	Animal_fat = c("p26062", "p26063"),
	Dairy = c("p100230", "p100530", "p26096", "p26102", "p102120", "p102810", "p102850", "p102870", "p102820", "p102830", "p102840", "p102860", "p102880", "p102890", "p102900", "p102910", "p102220", "p102140", "p102150"),
	Eggs = c("p102940", "p102950", "p102960", "p102970", "p102980"),
	Seafood = c("p103160", "p103190", "p103150", "p103230", "p103200", "p103210", "p103220", "canned.fish", "homemade.fish", "p103170", "p103180"),
	Meat = c("p103010", "p103020", "p103030", "p103040", "p103070", "p103080", "p103050", "p103060", "p103090", "p103100", "canned.meat", "homemade.meat"),
	Miscellaneous_animal_based_foods = c("p102000", "p102040")
)
dat1 <- dat.food %>% select(eid, all_of(unique(unlist(hpdi_map, use.names = FALSE))))
for (g in names(hpdi_map)) {dat1[[g]] <- rowSums1(dat1[, hpdi_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% select(eid, all_of(cols))
dat.hpdi <- dat1 %>% mutate(
	across(c(Vegetables, Fruit_not_fruit_juice, Legumes_and_vegetarian_protein_alternatives, Nuts, Whole_grains, Tea_and_coffee), ~ score_q5_zero(.x, reverse = FALSE), .names = "{.col}.hpdi"),
	across(c(Refined_grains, Potatoes, Sugary_drinks, Fruit_juices, Sweets_and_desserts, Animal_fat, Dairy, Eggs, Seafood, Meat, Miscellaneous_animal_based_foods), ~ score_q5_zero(.x, reverse = TRUE), .names = "{.col}.hpdi")
) %>% mutate(diet.hpdi.sum = rowSums1(across(ends_with(".hpdi")))) %>% select(eid, diet.hpdi.sum, ends_with(".hpdi"))
dose_save_rds(dat.hpdi, paste0(indir, "/Rdata/diet.hpdi.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 🌏 PHDI (同上，参考文献 PMID: 41245532)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field.w <- field.w %>% filter (PHDI > 0) %>% mutate(name = paste0("w.", name))
score_up10 <- function(x, thr) {x <- as.numeric(x); pmin(pmax(x / thr, 0), 1) * 10}
score_down10 <- function(x, thr) {x <- as.numeric(x); pmin(pmax(1 - x / thr, 0), 1) * 10}
score_down10_range <- function(x, low, high) {x <- as.numeric(x); out <- ifelse(is.na(x), NA_real_, ifelse(x <= low, 10, ifelse(x >= high, 0, 10 * (high - x) / (high - low)))); pmin(pmax(out, 0), 10)}
cols <- 'Whole_grains Starchy_vegetables Vegetables Whole_fruits Dairy_food Red_and_processed_meat Poultry Eggs Fish_and_shellfish Nuts_and_seeds Legumes Soybean Added_fat_saturated Added_fat_unsaturated Added_sugar_and_fruit_juices' %>% strsplit(' ') %>% unlist()
phdi_map <- list( # 15组；Legumes和Soybean权重各为0.5；总分0–140
	Whole_grains = c("w.p26076", "w.p26077", "w.p26074", "w.p26114"), # w.p26105 × 0.5
	Starchy_vegetables = c("w.p26115", "w.p26119", "w.p26120", "w.p26118"),
	Vegetables = c("w.p26098", "w.p26146", "w.p26143", "w.p26125", "w.p26147", "w.p26123", "w.p26144"),
	Whole_fruits = c("w.p26089", "w.p26090", "w.p26091", "w.p26093"),
	Dairy_food = c("w.p26084", "w.p26087", "w.p26096", "w.p26102", "w.p26099", "w.p26103", "w.p26131", "w.p26133", "w.p26150", "w.p26154"),
	Red_and_processed_meat = c("w.p26066", "w.p26100", "w.p26117", "w.p26122", "w.p26104"),
	Poultry = c("w.p26069", "w.p26121"),
	Eggs = c("w.p26088"),
	Fish_and_shellfish = c("w.p26109", "w.p26149", "w.p26132", "w.p26070"),
	Nuts_and_seeds = c("w.p26106", "w.p26108", "w.p26107"),
	Legumes = c("w.p26101"),
	Soybean = c("w.p26136", "w.p26086", "w.p26137"),
	Added_fat_saturated = c("w.p26014", "w.p26155"),
	Added_fat_unsaturated = c("w.p26032", "w.p26015", "w.p26016"),
	Added_sugar_and_fruit_juices = c("w.p26012")
)
dat1 <- dat.food %>% select(eid, sex, energy.kJ, all_of(unique(c(unlist(phdi_map, use.names = FALSE), "w.p26105"))))
for (g in names(phdi_map)) {dat1[[g]] <- rowSums1(dat1[, phdi_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% mutate(
	Whole_grains = Whole_grains + 0.5 * w.p26105,
	energy_kcal = energy.kJ * 0.239,
	Added_fat_saturated = ifelse(energy_kcal > 0, Added_fat_saturated * 9 / energy_kcal * 100, NA_real_),
	Added_fat_unsaturated = ifelse(energy_kcal > 0, Added_fat_unsaturated * 9 / energy_kcal * 100, NA_real_),
	Added_sugar_and_fruit_juices = ifelse(energy_kcal > 0, Added_sugar_and_fruit_juices * 4 / energy_kcal * 100, NA_real_)
) %>% select(eid, all_of(cols), sex)
dat.phdi <- dat1 %>% transmute(eid,
	Whole_grains.phdi = case_when(sex == 0 ~ score_up10(Whole_grains, 75), sex == 1 ~ score_up10(Whole_grains, 90), TRUE ~ NA_real_),
	Starchy_vegetables.phdi = score_down10_range(Starchy_vegetables, 50, 200),
	Vegetables.phdi = score_up10(Vegetables, 300),
	Whole_fruits.phdi = score_up10(Whole_fruits, 200),
	Dairy_food.phdi = score_down10_range(Dairy_food, 250, 1000),
	Red_and_processed_meat.phdi = score_down10_range(Red_and_processed_meat, 14, 100),
	Poultry.phdi = score_down10_range(Poultry, 29, 100),
	Eggs.phdi = score_down10_range(Eggs, 13, 120),
	Fish_and_shellfish.phdi = score_up10(Fish_and_shellfish, 28),
	Nuts_and_seeds.phdi = score_up10(Nuts_and_seeds, 50),
	Legumes.phdi = score_up10(Legumes, 100),
	Soybean.phdi = score_up10(Soybean, 50),
	Added_fat_saturated.phdi = score_down10(Added_fat_saturated, 10),
	Added_fat_unsaturated.phdi = pmin(pmax((Added_fat_unsaturated - 3.5) / (21 - 3.5), 0), 1) * 10,
	Added_sugar_and_fruit_juices.phdi = score_down10_range(Added_sugar_and_fruit_juices, 5, 25)
) %>% mutate(diet.phdi.sum = rowSums1(select(., ends_with(".phdi"), -Legumes.phdi, -Soybean.phdi)) + 0.5 * rowSums1(select(., Legumes.phdi, Soybean.phdi))
) %>% select(eid, diet.phdi.sum, ends_with(".phdi"))
dose_save_rds(dat.phdi, paste0(indir, "/Rdata/diet.phdi.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DI-GM (同上，参考文献 PMID: 41245532)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field.w <- field.w %>% filter (DIGM > 0) %>% mutate(name = paste0("w.", name))
get_med_by_sex <- function(x, sex){tapply(x, sex, median, na.rm = TRUE)} # 性别特异性中位数评分
score_benefit <- function(x, sex, med_by_sex){med <- med_by_sex[as.character(sex)]; ifelse(is.na(x) | is.na(sex), NA_real_, ifelse(is.na(med), NA_real_, ifelse(med == 0, x > 0, x >= med)))}
score_harm <- function(x, sex, med_by_sex){med <- med_by_sex[as.character(sex)]; ifelse(is.na(x) | is.na(sex), NA_real_, ifelse(is.na(med), NA_real_, ifelse(med == 0, x == 0, x < med)))}
cols <- 'Avocados Broccoli Coffee Fermented_dairy Green_tea Fiber Soybean Whole_grains Refined_grains Processed_meat Red_meat High_fat_diet' %>% strsplit(' ') %>% unlist()
digm_map <- list( # 12组，0–12分
	Coffee = c("w.p26081", "w.p26082"),
	Fermented_dairy = c("w.p26096", "w.p26102", "w.p26099", "w.p26103"),
	Fiber = c("w.p26017"),
	Soybean = c("w.p26136", "w.p26086", "w.p26137"),
	Whole_grains = c("w.p26076", "w.p26077", "w.p26074", "w.p26114"), # w.p26105 × 0.5
	Refined_grains = c("w.p26075", "w.p26078", "w.p26073", "w.p26072", "w.p26113", "w.p26083", "w.p26128", "w.p26097"),
	Processed_meat = c("w.p26122"),
	Red_meat = c("w.p26066", "w.p26100", "w.p26117", "w.p26104"),
	High_fat_diet = c("w.p26008")
)
dat1 <- dat.food %>% select(eid, sex, energy.kJ, p104100, p104140, p100420, all_of(unique(c(unlist(digm_map, use.names = FALSE), "w.p26105"))))
for (g in names(digm_map)) {dat1[[g]] <- rowSums1(dat1[, digm_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% mutate(
	Avocados = p104100, Broccoli = p104140, Green_tea = p100420,
	Whole_grains = Whole_grains + 0.5 * w.p26105,
	energy_kcal = energy.kJ * 0.239,
	High_fat_diet = ifelse(energy_kcal > 0, High_fat_diet * 9 / energy_kcal * 100, NA_real_)
) %>% select(eid, all_of(cols), sex)
benefit_vars <- c("Avocados", "Broccoli", "Coffee", "Fermented_dairy", "Green_tea", "Fiber", "Soybean", "Whole_grains")
harm_vars <- c("Refined_grains", "Processed_meat", "Red_meat")
meds <- c(lapply(dat1[benefit_vars], get_med_by_sex, sex = dat1$sex), lapply(dat1[harm_vars], get_med_by_sex, sex = dat1$sex))
dat.digm <- dat1 %>% mutate(
	across(all_of(benefit_vars), ~ score_benefit(.x, sex, meds[[cur_column()]]), .names = "{.col}.digm"),
	across(all_of(harm_vars), ~ score_harm(.x, sex, meds[[cur_column()]]), .names = "{.col}.digm"),
	High_fat_diet.digm = ifelse(is.na(High_fat_diet), NA_real_, as.numeric(High_fat_diet < 40))
) %>% mutate(diet.digm.sum = rowSums1(across(ends_with(".digm")))) %>% select(eid, diet.digm.sum, ends_with(".digm"))
dose_save_rds(dat.digm, paste0(indir, "/Rdata/diet.digm.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 HLCD: Healthy Low-Carbohydrate Diet (同上，参考文献 PMID: 41245532)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field.w <- field.w %>% filter (HLCD > 0) %>% mutate(name = paste0("w.", name))
score_11grp <- function(x, sex, reverse = FALSE){x <- as.numeric(x); sex <- as.numeric(sex); g <- ave(x, sex, FUN = function(v) if(all(is.na(v))) rep(NA_real_, length(v)) else dplyr::ntile(v, 11) - 1); if(reverse) g <- 10 - g; g}
cols <- 'Low_quality_carbohydrate Vegetable_protein Unsaturated_fat' %>% strsplit(' ') %>% unlist()
hlcd_map <- list( # 3组，0–30分
	Low_quality_carbohydrate = c("w.p26075", "w.p26073", "w.p26072", "w.p26113", "w.p26083", "w.p26128", "w.p26097", "w.p26079", "w.p26134", "w.p26064", "w.p26095", "w.p26115", "w.p26119", "w.p26120", "w.p26118"),
	Vegetable_protein = c("w.p26006"),
	Unsaturated_fat = c("w.p26015", "w.p26016", "w.p26032")
)
dat1 <- dat.food %>% select(eid, sex, energy.kJ, all_of(unique(unlist(hlcd_map, use.names = FALSE))))
for (g in names(hlcd_map)) {dat1[[g]] <- rowSums1(dat1[, hlcd_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% mutate(
	energy_kcal = energy.kJ * 0.239,
	Low_quality_carbohydrate = ifelse(energy_kcal > 0, Low_quality_carbohydrate * 4 / energy_kcal * 100, NA_real_),
	Vegetable_protein = ifelse(energy_kcal > 0, Vegetable_protein * 4 / energy_kcal * 100, NA_real_),
	Unsaturated_fat = ifelse(energy_kcal > 0, Unsaturated_fat * 9 / energy_kcal * 100, NA_real_)
) %>% select(eid, all_of(cols), sex)
dat.hlcd <- dat1 %>% transmute(eid,
	Low_quality_carbohydrate.hlcd = score_11grp(Low_quality_carbohydrate, sex, reverse = TRUE),
	Vegetable_protein.hlcd = score_11grp(Vegetable_protein, sex),
	Unsaturated_fat.hlcd = score_11grp(Unsaturated_fat, sex)
) %>% mutate(diet.hlcd.sum = rowSums1(across(ends_with(".hlcd")))) %>% select(eid, diet.hlcd.sum, ends_with(".hlcd"))
dose_save_rds(dat.hlcd, paste0(indir, "/Rdata/diet.hlcd.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 HLFD: Healthy Low-Fat Diet (同上，参考文献 PMID: 41245532)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field.w <- field.w %>% filter(HLFD > 0) %>% mutate(name = paste0("w.", name))
cols <- 'High_quality_carbohydrate Vegetable_protein Saturated_fat' %>% strsplit(' ') %>% unlist()
hlfd_map <- list( # 3组，0–30分
	High_quality_carbohydrate = c("w.p26076", "w.p26077", "w.p26074", "w.p26114", "w.p26089", "w.p26090", "w.p26091", "w.p26093", "w.p26101", "w.p26098", "w.p26146", "w.p26143", "w.p26125", "w.p26147", "w.p26123", "w.p26144"), # w.p26105 × 0.5
	Vegetable_protein = c("w.p26006"),
	Saturated_fat = c("w.p26014")
)
dat1 <- dat.food %>% select(eid, sex, energy.kJ, all_of(unique(c(unlist(hlfd_map, use.names = FALSE), "w.p26105"))))
for (g in names(hlfd_map)) {dat1[[g]] <- rowSums1(dat1[, hlfd_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% mutate(
	High_quality_carbohydrate = High_quality_carbohydrate + 0.5 * w.p26105,
	energy_kcal = energy.kJ * 0.239,
	High_quality_carbohydrate = ifelse(energy_kcal > 0, High_quality_carbohydrate * 4 / energy_kcal * 100, NA_real_),
	Vegetable_protein = ifelse(energy_kcal > 0, Vegetable_protein * 4 / energy_kcal * 100, NA_real_),
	Saturated_fat = ifelse(energy_kcal > 0, Saturated_fat * 9 / energy_kcal * 100, NA_real_)
) %>% select(eid, all_of(cols), sex)
dat.hlfd <- dat1 %>% transmute(eid,
	High_quality_carbohydrate.hlfd = score_11grp(High_quality_carbohydrate, sex),
	Vegetable_protein.hlfd = score_11grp(Vegetable_protein, sex),
	Saturated_fat.hlfd = score_11grp(Saturated_fat, sex, reverse = TRUE)
) %>% mutate(diet.hlfd.sum = rowSums1(across(ends_with(".hlfd")))) %>% select(eid, diet.hlfd.sum)
dose_save_rds(dat.hlfd, paste0(indir, "/Rdata/diet.hlfd.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩  rEDII (reversed energy‑adjusted diet inflammatory index) (同上，参考文献 PMID: 41245532,37480061)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field.w <- field.w %>% filter(rEDII != 0) %>% mutate(name = paste0("w.", name))
get_z  <- function(x) {x <- as.numeric(x); s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(NA_real_, length(x)) else (x - mean(x, na.rm = TRUE)) / s}
get_cp <- function(z) {(pnorm(as.numeric(z)) * 2) - 1}
cols <- paste(
	'Alcohol Carbohydrate Total_protein Total_fat Polyunsaturated_fat Saturated_fat Monounsaturated_fat',
	'Trans_fat Omega3_fatty_acids Omega6_fatty_acids Dietary_fiber Cholesterol Beta_carotene Vitamin_B12',
	'Vitamin_B6 Niacin Selenium Thiamin Riboflavin Vitamin_A Vitamin_C',
	'Vitamin_D Vitamin_E Folate Iron Magnesium Zinc Tea'
) %>% strsplit(' ') %>% unlist()
redii_map <- list( # 28组；Energy仅作为密度法分母
	Alcohol = c("w.p26030"), Carbohydrate = c("w.p26013"), Total_protein = c("w.p26005"), Total_fat = c("w.p26008"), Polyunsaturated_fat = c("w.p26015", "w.p26016"), Saturated_fat = c("w.p26014"), Monounsaturated_fat = c("w.p26032"),
	Trans_fat = c("w.p26155"), Omega3_fatty_acids = c("w.p26015"), Omega6_fatty_acids = c("w.p26016"), Dietary_fiber = c("w.p26017"), Cholesterol = c("w.p26037"), Beta_carotene = c("w.p26039"), Vitamin_B12 = c("w.p26021"),
	Vitamin_B6 = c("w.p26020"), Niacin = c("w.p26054"), Selenium = c("w.p26058"), Thiamin = c("w.p26034"), Riboflavin = c("w.p26035"), Vitamin_A = c("w.p26061"), Vitamin_C = c("w.p26023"),
	Vitamin_D = c("w.p26029"), Vitamin_E = c("w.p26028"), Folate = c("w.p26022"), Iron = c("w.p26019"), Magnesium = c("w.p26025"), Zinc = c("w.p26033"), Tea = c("w.p26141", "w.p26142")
)
redii_beta <- c(
	Alcohol = -0.278, Carbohydrate = 0.097, Total_protein = 0.021, Total_fat = 0.298, Polyunsaturated_fat = -0.337, Saturated_fat = 0.373, Monounsaturated_fat = -0.009,
	Trans_fat = 0.229, Omega3_fatty_acids = -0.436, Omega6_fatty_acids = -0.159, Dietary_fiber = -0.663, Cholesterol = 0.110, Beta_carotene = -0.584, Vitamin_B12 = 0.106,
	Vitamin_B6 = -0.365, Niacin = -0.246, Selenium = -0.191, Thiamin = -0.098, Riboflavin = -0.068, Vitamin_A = -0.401, Vitamin_C = -0.424,
	Vitamin_D = -0.446, Vitamin_E = -0.419, Folate = -0.190, Iron = 0.032, Magnesium = -0.484, Zinc = -0.313, Tea = -0.536
)
stopifnot(identical(cols, names(redii_map)), identical(cols, names(redii_beta)))
dat1 <- dat.food %>% select(eid, energy.kJ, all_of(unique(unlist(redii_map, use.names = FALSE))))
for (g in names(redii_map)) {dat1[[g]] <- rowSums1(dat1[, redii_map[[g]], drop = FALSE])}
dat1 <- dat1 %>% mutate(
	energy_kcal = energy.kJ * 0.239,
	across(all_of(cols), ~ ifelse(energy_kcal > 0, .x / energy_kcal * 1000, NA_real_))
) %>% select(eid, all_of(cols))
dat.redii <- dat1 %>% mutate(across(all_of(cols), ~ get_cp(get_z(.x)) * redii_beta[[cur_column()]], .names = "{.col}.edii")
) %>% mutate(diet.edii.sum = rowSums1(across(ends_with(".edii"))), diet.redii.sum = -diet.edii.sum) %>% select(eid, diet.redii.sum)
dose_save_rds(dat.redii, paste0(indir, "/Rdata/diet.redii.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 🍟美国MAHA（根据MAHA指南PDF文件）
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
score_up10 <- function(x, target) {x <- as.numeric(x); pmin(pmax(x / target, 0), 1) * 10}	# 越多越好 (Up): 达到 target 得10分
score_down10 <- function(x, limit_good, limit_bad) {x <- as.numeric(x); pmin(pmax((limit_bad - x) / (limit_bad - limit_good), 0), 1) * 10} # 越少越好 (Down): < limit_good 得10分，> limit_bad 得0分
score_mid10 <- function(x, lo, hi) {x <- as.numeric(x); ifelse(x >= lo, 10, (x / lo) * 10)}	# 区间最好 (Mid): lo~hi 区间得10分 (MAHA 蛋白质: <1.2分低, 1.2-1.6满分)[ 逻辑: >= lo 且 <= hi 满分; < lo 线性下降; > hi 暂不扣分(Prioritize)]
score_min10 <- function(x, target) {x <- as.numeric(x); pmin(pmax(x / target, 0), 1) * 10} # 修改：新增最低达标得满分函数，使 MAHA 蛋白评分与 Methods 表述一致
dairy_vars 		<- 'p26150 p26096 p102820 p102830 p102840 p102860 p102880 p102890 p102900 p102910' %>% strsplit(' ') %>% unlist()
veg_vars 		<- 'p104090 p104240 p104370 p104140 p104160 p104180 p104300 p104310 p104220 p104230 p104260 p104340 p104350 p104150 p104170 p104290 p104330 p104130 p104190 p104270 p104360 p104060 p104070 p104200 p104210 p104250 p104320 p104380' %>% strsplit(' ') %>% unlist()
fruit_vars 		<- 'p104450 p104560 p104470 p104480 p104490 p104530 p104540 p104430 p104420 p104550 p104410 p104580 p104500 p104460 p104510 p104570 p104520 p104440 p104590' %>% strsplit(' ') %>% unlist()
fat_vars 		<- 'p26110 oil p26062 p26063 p102440 p102450 p103160 p104100' %>% strsplit(' ') %>% unlist()
wholegrain_vars <- 'p100770 p100810 p100800 p100840 p100850 bread.b baguette.b bap.b roll.b p102720 p102740 p102780' %>% strsplit(' ') %>% unlist()
refined_vars 	<- 'p100820 p101230 p101240 p101250 p101260 p101270 p102710 p102730 p102750 p102760 p102770 bread.w baguette.w bap.w roll.w p100830 p100860' %>% strsplit(' ') %>% unlist()
ssb_vars 		<- 'p100170 p100180 p100550 p100160 p100540 p100530 p100230 p100190 p100200 p100210 p100220' %>% strsplit(' ') %>% unlist()
sweet_vars 		<- 'p26064 p102120 p102140 p102150 p102170 p102180 p102190 p102200 p102210 p102220 p102230 p102260 p102270 p102280 p102290 p102300 p102310 p102320 p102330 p102340 p102350 p102360 p102370 p102380 p102060 p102010 p102020 p102050 p102070 p101970 p101980 p101990 p102030' %>% strsplit(' ') %>% unlist()
hp_vars 		<- 'p102460 p102470 p102480 p102500 p102040 p102000 p104020 p103050 p103170 p103180 p103010 p103070 p103080' %>% strsplit(' ') %>% unlist()
alcohol_vars 	<- 'p26067 p26138 p26151 p26152 p26153' %>% strsplit(' ') %>% unlist()
dat.maha <- dat.food %>% mutate(
	ef 				= if_else(!is.na(energy.kJ) & energy.kJ > 0, 8368 / energy.kJ, NA_real_),
	protein.maha 	= score_min10(w.p26005 / weight, 1.2), # 修改：原来用 score_mid10(lo=1.2, hi=1.6)，现改为达到最低推荐值即满分
	dairy.maha		= score_up10(rowSums1(across(all_of(dairy_vars))) * ef, 3),
	veg.maha		= score_up10(rowSums1(across(all_of(veg_vars))) * ef, 3), 
	fruit.maha 		= score_up10(rowSums1(across(all_of(fruit_vars))) * ef, 2),
	wholegrain.maha = score_up10(rowSums1(across(all_of(wholegrain_vars))) * ef, 2), 
	fat.maha 		= score_up10(rowSums1(across(all_of(fat_vars))) * ef, 1), # 修改：用 rowSums1 替代 rowSums(..., na.rm = TRUE)，并与其他 servings 成分一样做能量校正
	upf.maha		= score_down10(rowSums1(across(all_of(c(refined_vars, ssb_vars, sweet_vars, hp_vars)))) * ef, 0, 4), # 修改：用 rowSums1 替代 rowSums(..., na.rm = TRUE)
	alcohol.maha	= score_down10(rowSums1(across(all_of(alcohol_vars))), 0, 2.5),
	sodium.maha		= score_down10(w.p26052, 2300, 5000) # 低于指南要求的 2300mg 给予满分 10 分；随着摄入量增加分数线性下降，设定 5000mg 为 0 分
) %>% mutate(diet.maha.sum = rowSums1(across(ends_with(".maha")))) %>% select(eid, diet.maha.sum, ends_with(".maha")) 
dose_save_rds(dat.maha, paste0(indir, "/Rdata/diet.maha.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FDS: Flavodiet Score（servings/day） 【PMID: 40456886, 39292460, 38778045】
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat_field <- field %>% filter (FDS > 0)
cols <- 'tea red_wine apples berries grapes oranges grapefruit sweet_peppers onions dark_chocolate' %>% strsplit(' ') %>% unlist()
fds_map <- list( # 10组，单位均为servings/day
	tea = c("p100400", "p100420"),
	red_wine = c("p26152"),
	apples = c("p104450"), berries = c("p104470"), grapes = c("p104500"), oranges = c("p104530"), grapefruit = c("p104490"),
	sweet_peppers = c("p104290"), onions = c("p104260"),
	dark_chocolate = c("p102290")
)
dat1 <- dat.food %>% select(eid, all_of(unique(unlist(fds_map, use.names = FALSE))))
for (g in names(fds_map)) {dat1[[g]] <- rowSums1(dat1[, fds_map[[g]], drop = FALSE])}
dat.fds <- dat1 %>% select(eid, all_of(cols)) %>% rename_with(~ paste0(.x, ".fds"), all_of(cols)) %>%
	mutate(diet.fds.sum = rowSums1(across(ends_with(".fds")))) %>% select(eid, diet.fds.sum, ends_with(".fds"))
dose_save_rds(dat.fds, paste0(indir, "/Rdata/diet.fds.rds"))
# Flavonoid subclass intake（mg/day）❓还没完成, 需要食物成分数据库，不能由servings/day直接得到
})

dose_cached_phase("diet_merge", inputs=c(phase_rds(c("diet.medito", "diet.fudan")), pattern_files), outputs=phase_rds("diet"), code={
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 👨‍👩‍👧‍👦 merge all together
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lst <- c("medito", "medi24", "dash", "mind", "hpdi", "ahei", "phdi", "modern", "digm", "maha", "hlcd", "hlfd", "redii", "fds") # "upf"
dat0 <- dose_merge_participant_rds(file.path(indir, "Rdata", paste0("diet.", lst, ".rds")))
sum_cols <- names(dat0)[grepl("\\.sum$", names(dat0))]
lapply(dat0 %>% select(matches("\\.sum$")), summary)
lapply(dat0 %>% select(matches("\\.(pts|q5)$")), table, useNA = "always")
dat.fudan <- readRDS(file.path(indir, "Rdata/diet.fudan.rds")) %>% transmute(eid = as.character(eid))
dat0 <- dat0 %>% left_join(dat.fudan, by = "eid")
dose_save_rds(dat0, paste0(indir, "/Rdata/diet.rds"))
rm(dat0, dat.fudan)
invisible(gc())
})

dose_cached_phase("clinical_merge", inputs=c(clinical_files, phase_rds("diet")), outputs=phase_rds("all"), code={
sum_cols <- grep("\\.sum$", names(readRDS(file.path(indir, "Rdata/diet.rds"))), value=TRUE)
invisible(gc())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 👨‍👩‍👧‍👦 与其他数据合并
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(tidyverse, data.table, survival)

dir0 <- dose_drive()
indir <- dose_ukb_dir()
invisible(lapply(c("phenotype.R", "association.R", "plotting.R"), \(f) source(file.path(dose_helper_dir(), f))))
setwd2(file.path(dose_output_dir(), "ukb/00data"))

pacman::p_load(matrixStats)
lst <- c("phe", "bbc", "diet", "move", "ses", "death", "fod.icd10", "fod.icd9", "fod.opcs4", "fod.srd", "fod.gp", "fod.ref", "pca", "gen", "prs")
dat0 <- dose_merge_participant_rds(file.path(indir, "Rdata", paste0(lst, ".rds")))
dat0 <- dat0 %>% filter(!grepl("^-", eid))
invisible(gc())
message("Participant merge complete; deriving covariates")
dat0 <- dat0 %>% mutate(
	eid = as.character(eid),
	fod_icd10_death = date_death, 
	drug.lipid = ifelse(matrixStats::rowSums2(as.matrix(across(starts_with("drug.big3"), ~ grepl("1", as.character(.))))) > 0, 1, 0),
	drug.dm = ifelse(matrixStats::rowSums2(as.matrix(across(starts_with("drug.big3"), ~ grepl("3", as.character(.))))) > 0, 1, 0),
	drug.htn = ifelse(matrixStats::rowSums2(as.matrix(across(starts_with("drug.big3"), ~ grepl("2", as.character(.))))) > 0, 1, 0),
	date_mi = do.call(pmin, c(across(any_of(c("fod_icd10_cvd_cad", "fod_icd9_cvd_cad", "fod_srd_cvd_cad"))), list(na.rm = TRUE))),
	date_stroke_i = do.call(pmin, c(across(matches("^fod_icd.*_cvd_stroke_i$")), list(na.rm = TRUE))),
	date_stroke_h = do.call(pmin, c(across(matches("^fod_icd.*_cvd_stroke_(i|s)h$")), list(na.rm = TRUE))),
	date_stroke = pmin(date_stroke_i, date_stroke_h, na.rm = TRUE),
	date_mace3 = pmin(date_mi, date_stroke, date_cvd_death, na.rm = TRUE)
) %>% mutate( # ⓸➕⓸ = LE8
	# 饮食 🥢
	diet.sum = diet.medito.sum, diet.pts = diet.medito.pts,
	# 运动 🏃‍ 
	across(c(days_pa_mod, dura_pa_mod, days_pa_vig, dura_pa_vig), ~ ifelse(.x < 0, NA, .x)),
	pa_mod_w = days_pa_mod * pmin(dura_pa_mod, 180),
	pa_vig_w = days_pa_vig * pmin(dura_pa_vig, 180),
	pa_mins = rowSums1(data.frame(pa_mod_w, 2 * pa_vig_w), per = 0.5),
	pa.pts = cut(pa_mins, breaks = c(0, 1, 30, 60, 90, 120, 150, Inf), labels = c(0, 20, 40, 60, 80, 90, 100), right = FALSE) %>% as.character() %>% as.numeric(),
	# 吸烟 🚬
	smoke.c = factor(smoke_status, levels = 0:2, labels = c("never", "previous", "current")),
	smoke_quit_years = ifelse(smoke_quit_age > 0 & age >= smoke_quit_age, age - smoke_quit_age, NA_real_),
	smoke_base = case_when(smoke_status == 0 ~ 100, smoke_status == 2 ~ 0, smoke_status == 1 & smoke_quit_years >= 5 ~ 75, smoke_status == 1 & smoke_quit_years >= 1 ~ 50, smoke_status == 1 ~ 25, TRUE ~ NA_real_),
	smoke_shs_home = if_else(smoke_2nd_home > 0, 1L, 0L),
	smoke.pts = pmax(smoke_base - if_else(smoke_shs_home == 1L & smoke_base > 0, 20, 0), 0),	
	# 睡眠 🛏
	sleep.pts = cut(sleep_duration, breaks = c(0, 4, 5, 6, 7, 9, 10, 24), labels = c(0, 20, 40, 70, 100, 90, 40), right = FALSE) %>% as.character() %>% as.numeric(),
	# BMI 🧖
	bmi.c = cut(bmi, breaks = c(10, 18.5, 25, 30, 100), labels = c("lean", "healthy", "overweight", "obese")),
	bmi.c = factor(bmi.c, levels = c("healthy", "lean", "overweight", "obese")),
	bmi.pts = cut(bmi, breaks = c(0, 25, 30, 35, 40, Inf), labels = c(100, 70, 30, 15, 0), right = FALSE) %>% as.character() %>% as.numeric(),
	# 血脂 🍟
	nonhdl = pmax((bb_TC - bb_HDL) * 38.67, 0), # Units of measurement are mmol/L, 转换为 mg/dL
	nonhdl.c = cut(nonhdl, breaks = c(0, 130, 160, 190, 220, Inf), labels = c(100, 60, 40, 20, 0), right = FALSE) %>% as.character() %>% as.numeric(),
	nonhdl.pts = ifelse(drug.lipid == 1 & nonhdl.c > 0, pmax(nonhdl.c - 20, 0), nonhdl.c),
	# 血糖 🍬
	hba1c_ngsp = bb_HBA1C * 0.0915 + 2.15, # NGSP%
	dm.yes = if_else(!is.na(fod_srd_t2dm) | !is.na(fod_icd10_t2dm) | !is.na(fod_ref_t2dm) | drug.dm == 1 | (!is.na(hba1c_ngsp) & hba1c_ngsp >= 6.5), 1, 0),
	hba1c.pts = case_when(is.na(hba1c_ngsp) ~ NA_real_, dm.yes == 0 & hba1c_ngsp < 5.7 ~ 100, dm.yes == 0 & hba1c_ngsp < 6.5 ~ 60, hba1c_ngsp < 7 ~ 40, hba1c_ngsp < 8 ~ 30, hba1c_ngsp < 9 ~ 20, hba1c_ngsp < 10 ~ 10, TRUE ~ 0),
	# 血压 🦆 sbp_man[93], dbp_man[94], dbp_auto[4079], sbp_auto[4080]
	sbp = rowMeans(across(matches("^sbp_.*_i0.*")), na.rm = TRUE),
	dbp = rowMeans(across(matches("^dbp_.*_i0.*")), na.rm = TRUE),
	htn.yes = if_else(!is.na(fod_srd_cvd_htn) | !is.na(fod_icd10_cvd_htn)  | !is.na(fod_ref_cvd_htn) | drug.htn == 1, 1L, 0L),
	bp_base = case_when(sbp >= 160 | dbp >= 100 ~ 0, (sbp >= 140 & sbp < 160) | (dbp >= 90 & dbp < 100) ~ 25, (sbp >= 130 & sbp < 140) | (dbp >= 80 & dbp < 90) ~ 50, (sbp >= 120 & sbp < 130) & (dbp > 0 & dbp < 80) ~ 75, (sbp > 0 & sbp < 120) & (dbp > 0 & dbp < 80) ~ 100, TRUE ~ NA_real_),
	bp.pts = if_else(coalesce(drug.htn, 0) == 1 & bp_base > 0, pmax(bp_base - 20, 0), bp_base),
	# bp1.pts = case_when(is.na(sbp) | is.na(dbp) ~ NA_real_, TRUE ~ pmax(pmin(100 - 2 * pmax(sbp - 115, 0) - 2 * pmax(dbp - 75, 0), 100), 0)),
	# bp2.pts = case_when(is.na(htn.yes) & is.na(drug.htn) ~ NA_real_, TRUE ~ pmax(100 - 30 * coalesce(htn.yes, 0) - 10 * coalesce(drug.htn, 0), 0))	
	# 喝酒 🍷
	alcohol.c = factor(alcohol_status, levels = 0:2, labels = c("never", "previous", "current")),
	alcohol_days = case_match(alcohol_freq, 1~7, 2~3.5, 3~1.5, 4~0.5, 5~0.1, 6~0, .default=NA),
	alcohol_units = case_match(alcohol_amount, 1~1.5, 2~3.5, 3~5.5, 4~8, 5~12, .default=NA),
	# alcohol.pts = case_when(alcohol_status %in% c(0, 1) ~ 100, sex == 0 ~ as.numeric(as.character(cut(alcohol_days * alcohol_units * 8 / 7, c(0, 7, 14, 21, 28, Inf), c(100, 80, 50, 20, 0), right = FALSE))), sex == 1 ~ as.numeric(as.character(cut(alcohol_days * alcohol_units * 8 / 7, c(0, 14, 28, 42, 56, Inf), c(100, 80, 50, 20, 0), right = FALSE))), TRUE ~ NA_real_),
	alcohol.pts = case_when(alcohol_status %in% c(0, 1) ~ 100, is.na(alcohol_status) ~ NA_real_, sex == 0 ~ as.numeric(as.character(cut(coalesce(alcohol_days * alcohol_units * 8 / 7, 14.1), c(0, 7, 14, 21, 28, Inf), c(100, 80, 50, 20, 0), right = FALSE))), sex == 1 ~ as.numeric(as.character(cut(coalesce(alcohol_days * alcohol_units * 8 / 7, 28.1), c(0, 14, 28, 42, 56, Inf), c(100, 80, 50, 20, 0), right = FALSE))), TRUE ~ NA_real_)
) %>% mutate(
	le8.sco = rowMeans1(across(all_of(vars.le8)), 0.7), # 📍允许30%的NA值
	le8r.sco = -le8.sco, le8r.3c = f3c(le8r.sco),	
#	cad.prs = cad.6m.score_sum, cad.prs.3c = f3c(cad.prs), # 🏮
	cvh.3c = cut(le8.sco, breaks = c(0, 50, 80, Inf), right = FALSE, labels = c("unhealthy", "average", "healthy")),
	cvh.3c = factor(cvh.3c, levels = c("healthy", "average", "unhealthy"))
)
dose_save_rds(dat0, paste0(indir, "/Rdata/all.rds"))
message("Combined participant data saved: all.rds")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std_10_90 <- function(x){  # 统一标准化到同一尺度（10th→90th 的增量）
	x <- as.numeric(x); q <- quantile(x, c(0.10, 0.90), na.rm = TRUE, names = FALSE); if(!all(is.finite(q)) || q[2] == q[1]) rep(NA_real_, length(x)) else (x - q[1]) / (q[2] - q[1])
}
score_q5 <- function(x){ # 真正五分位 -> 0/25/50/75/100
	x <- as.numeric(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x); out[ok] <- c(0, 25, 50, 75, 100)[dplyr::ntile(x[ok], 5)]; out
}
score_0_100 <- function(x){
	x <- as.numeric(x); ok <- is.finite(x); if(sum(ok) < 2) rep(NA_real_, length(x)) else {mn <- min(x[ok]); mx <- max(x[ok]); if(mx == mn) rep(NA_real_, length(x)) else {out <- rep(NA_real_, length(x)); out[ok] <- (x[ok] - mn) / (mx - mn) * 100; out}}
}
score_0_100_q5 <- function(x){ # 先缩放到 0-100，再按固定宽度切成 5 档
	x <- as.numeric(x); ok <- is.finite(x); if(sum(ok) < 2) rep(NA_real_, length(x)) else {mn <- min(x[ok]); mx <- max(x[ok]); if(mx == mn) rep(NA_real_, length(x)) else {out <- rep(NA_real_, length(x)); out[ok] <- (x[ok] - mn) / (mx - mn) * 100; as.numeric(as.character(cut(out, breaks = c(-Inf, 20, 40, 60, 80, Inf), labels = c(0, 25, 50, 75, 100), right = FALSE)))}}
}
dat <- dat0 %>% select(eid, all_of(sum_cols), matches("\\.(pts|q5)$"))
rm(dat0)
invisible(gc()) # Score QC needs no clinical/genetic columns from the full join.
dat <- dat %>% mutate(
	across(all_of(sum_cols), std_10_90, .names = "{sub('\\\\.sum$', '', .col)}.std"),
    across(all_of(sum_cols), score_0_100, .names = "{sub('\\\\.sum$', '', .col)}.s100"),
	across(all_of(sum_cols), score_0_100_q5, .names = "{sub('\\\\.sum$', '', .col)}.pts"),
	across(all_of(sum_cols), score_q5, .names = "{sub('\\\\.sum$', '', .col)}.q5"),
	across(all_of(sum_cols), f3c, .names = "{sub('\\\\.sum$', '', .col)}.3c")
)
lapply(dat %>% select(matches("\\.sum$")), summary)
lapply(dat %>% select(matches("\\.(pts|q5)$")), table, useNA = "always")
rm(dat)
invisible(gc())






#❓❓❓❓❓❓❓❓❓❓❓❓❓❓❓需要检查
})

dose_cached_phase("upf", inputs=c(phase_rds("diet.fudan"), foods_file), outputs=upf_files, code={
dat.food <- readRDS(file.path(indir, "Rdata/diet.fudan.rds"))
field <- read_excel(file.path(dose_files_dir(), "foods.xlsx"), sheet="Sheet2")
field.w <- read_excel(file.path(dose_files_dir(), "foods.xlsx"), sheet="Sheet3")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 upf: ultra-processed foods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (同上，参考文献 PMID: 41245532)
# dat_field.w <- field.w %>% filter (upf > 0) %>% mutate(name = paste0("w.", name))
cols <- 'Beverages Dairy_products Fruits_vegetables Meat_fish_eggs Sauces_soups Starchy_foods Sugary_snacks Savory_snacks' %>% strsplit(' ') %>% unlist()
	cols <- paste0("weight_", cols)
upf_map <- list( # 8组，共35个UPF食物
	weight_Beverages = c("w.p26126", "w.p26127", "w.p26138"),
	weight_Dairy_products = c("w.p26086", "w.p26084", "w.p26087", "w.p26102", "w.p26124"),
	weight_Fruits_vegetables = c("w.p26064", "w.p26090", "w.p26144"),
	weight_Meat_fish_eggs = c("w.p26122", "w.p26137", "w.p26145", "w.p26149"),
	weight_Sauces_soups = c("w.p26106", "w.p26129", "w.p26130", "w.p26110"),
	weight_Starchy_foods = c("w.p26075", "w.p26068", "w.p26076", "w.p26119", "w.p26105", "w.p26078", "w.p26116", "w.p26073", "w.p26079", "w.p26097"),
	weight_Sugary_snacks = c("w.p26080", "w.p26085", "w.p26140"),
	weight_Savory_snacks = c("w.p26108", "w.p26083", "w.p26134")
)
# 所有食物和饮料重量变量，用于计算总摄入重量
weight_vars <- field.w %>% filter(startsWith(Food.Group, "w.")) %>% pull(name)#  p26000
	weight_vars <- unique(weight_vars[!is.na(weight_vars) & weight_vars != ""])
	weight_vars <- setdiff(weight_vars, "p26000") # 若包含已经汇总的总重量字段，则排除，避免重复计入
	weight_vars <- paste0("w.", weight_vars)
stopifnot(length(unlist(upf_map, use.names = FALSE)) == 35, !anyDuplicated(unlist(upf_map, use.names = FALSE))) # 检查UPF映射：应为35个食物，且不能重复
dat1 <- dat.food %>% select(eid, all_of(unique(c(unlist(upf_map, use.names = FALSE), weight_vars))))
for (g in names(upf_map)) {dat1[[paste0(g, ".upf")]] <- rowSums1(dat1[, upf_map[[g]], drop = FALSE])}
dat.upf <- dat1 %>% mutate(
	diet.total.w = rowSums1(across(all_of(weight_vars))), # field.w中全部食物和饮料重量之和，g/day
	diet.upf.sum = rowSums1(across(all_of(paste0(cols, ".upf")))), # 35个UPF食物摄入重量之和，g/day
	diet.upf.prop = ifelse(diet.total.w > 0, diet.upf.sum / diet.total.w * 100, NA_real_) # UPF重量占全部饮食重量的百分比
) %>% select(eid, all_of(paste0(cols, ".upf")), diet.upf.prop)
dose_save_rds(dat.upf, paste0(indir, "/Rdata/diet.upf.rds"))

# 参考文献 PMID: 37348154
if (all(c("UPF_2", "Energy_KJ/100g") %in% names(field))) {
calc_grams <- function(dat, info) {x <- as.matrix(dat[, info$name, drop = FALSE]); x[is.na(x)] <- 0; sweep(x, 2, info$Portion, "*")}
calc_energy <- function(dat, info) {g <- calc_grams(dat, info); sweep(g, 2, info$`Energy_KJ/100g` / 100, "*")}
calc_serv <- function(dat, info) {x <- as.matrix(dat[, info$name, drop = FALSE]); x[is.na(x)] <- 0; x}
add_upf_metric <- function(dat, info, prefix) {
	serv <- rowSums1(calc_serv(dat, info)); gram <- rowSums1(calc_grams(dat, info)); ener <- rowSums1(calc_energy(dat, info))
	tibble::tibble(
		!!paste0(prefix, "_servings_day")  := serv, # 每日UPF摄入份数
		!!paste0(prefix, "_grams_day")     := gram, # 每日UPF摄入克数
		!!paste0(prefix, "_grams_pct")     := ifelse(dat$food_total_weight > 0, gram / dat$food_total_weight * 100, NA_real_), # UPF克数占总食物重量比例
		!!paste0(prefix, "_energy_kj_day") := ener, # 每日来自UPF的能量
		!!paste0(prefix, "_energy_pct")    := ifelse(dat$energy.kJ > 0, ener / dat$energy.kJ * 100, NA_real_) # UPF能量占总能量比例
	)
}
dat_field <- field %>% filter(UPF_2 > 0) %>% filter(!str_detect(name, "^bap\\."), !str_detect(name, "^roll\\."))
add_bread_rows <- tibble::tribble(
	~Category, ~Food.Group, ~Food.Item, ~FiledID, ~Portion, ~name, ~`Energy_KJ/100g`, ~UPF_2,
	"Cereal_Bread", "Refined_grains", "Bap intake",        101090, 90, "bap.total",  1065, 1,
	"Cereal_Bread", "Refined_grains", "Bread roll intake", 101160, 60, "roll.total", 1084, 1
)
dat_field <- bind_rows(dat_field, add_bread_rows) %>% distinct(name, .keep_all = TRUE)
coffee_vars <- c("p100250", "p100270", "p100290", "p100300", "p100310", "p100330")
sweetbev_vars <- c("p100160", "p100170", "p100180", "p100530")
alcohol_vars <- c("p26151", "p26138")
sweetener_vars <- c("p100380", "p100500")
dat_field <- dat_field %>% mutate(
	upf_subgroup = case_when(
		name %in% coffee_vars ~ "coffee",
		name %in% sweetbev_vars ~ "sweetened_beverages",
		name %in% alcohol_vars ~ "alcoholic_beverages",
		name %in% sweetener_vars ~ "artificial_sweeteners",
		Category == "Dairy_dairy_free_products" ~ "dairy_products",
		Category %in% c("Vegetables", "Legumes_Nuts") ~ "fruits_vegetables",
		Category %in% c("Meats", "Meat_substitutes") ~ "meat_fish_eggs",
		name %in% c("p102530") ~ "sauces_soups",
		Category == "Cereal_Bread" | name %in% c("p102000", "bap.total", "roll.total") ~ "starchy_foods",
		Food.Group == "Snacks_and_pastries" ~ "sugary_snacks",
		name %in% c("p102040", "p102460", "p102470", "p102480", "p102500") ~ "savory_snacks",
		TRUE ~ "other_upf"
	),
	upf_beverage = upf_subgroup %in% c("coffee", "sweetened_beverages", "alcoholic_beverages", "artificial_sweeteners")
)
info_total <- dat_field
info_bev <- dat_field %>% filter(upf_beverage)
info_nonbev <- dat_field %>% filter(!upf_beverage)
dat.upf <- bind_cols(
	dat.food %>% select(eid, food_total_weight, energy.kJ),
	add_upf_metric(dat.food, info_total, "diet.upf"),
	add_upf_metric(dat.food, info_bev, "diet.upf_beverage"),
	add_upf_metric(dat.food, info_nonbev, "diet.upf_without_beverage")
)
for(g in sort(unique(dat_field$upf_subgroup))) {
	info_g <- dat_field %>% filter(upf_subgroup == g)
	dat.upf <- bind_cols(dat.upf, add_upf_metric(dat.food, info_g, paste0("diet.upf_", g)))
}
dat.upf <- dat.upf %>% select(eid, starts_with("diet.upf_"))
dose_save_rds(dat.upf, paste0(indir, "/Rdata/diet.upf_2.rds"))
} else {
    message("Skipping optional UPF_2 appendix: foods.xlsx lacks UPF_2 / Energy_KJ/100g. Main DOSE inputs are complete.")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Dietary Diversity Score (DDS) 【PMID: 39952880】
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 主 DDS：DDS18_total

# 健康多样性：DDS_healthy

# 非健康/加工多样性：DDS_unhealthy 或 DDS_processing

# Shannon diversity index
})
}))

# ===========================================================================
# foodwas / cohort
# ===========================================================================
dose_register_step("foodwas", "cohort", expression({
pacman::p_load(tidyverse, survival, writexl)

dir0 <- dose_drive()
indir <- dose_ukb_dir()
source(file.path(dose_helper_dir(), "phenotype.R"))
DIR_PROJECT <- dose_output_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data"); dir.create(DIR_DATA, recursive=TRUE, showWarnings=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Cohort and covariates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS(file.path(indir, "Rdata/all.rds")) %>% mutate(eid=as.character(eid))
diet <- readRDS(file.path(indir, "Rdata/diet.fudan.rds")) %>% mutate(eid=as.character(eid))
# The WebQ file is authoritative for these fields; all.rds may already contain
# older copies from another analysis, which would otherwise acquire .x/.y suffixes.
diet_cohort <- diet %>% select(eid, energy.kJ, enroll, start_date)
dat <- dat0 %>% select(-any_of(setdiff(names(diet_cohort), "eid"))) %>%
    inner_join(diet_cohort, by="eid") %>%
	mutate(start_date=as.Date(start_date), date_lost=as.Date(date_lost), date_death=as.Date(date_death))
rm(dat0, diet, diet_cohort)
invisible(gc())
stopifnot(all(dat$enroll >= 1), all(!is.na(dat$energy.kJ)), all(!is.na(dat$start_date)))
cat("有效 WebQ 人群：", nrow(dat), "\n") # 187372 

has_code <- function(x, code) !is.na(x) & grepl(paste0("(^|\\|)", code, "(\\||$)"), as.character(x))
region_levels <- c(
	"London", "North West", "North East", "Yorkshire", "West Midlands",
	"East Midlands", "South East", "South West", "Scotland", "Wales"
)

dat <- dat %>% mutate(
	sex.b=factor(sex, 0:1, c("Female","Male")),
	ethnic.c3=factor(ifelse(is.na(ethnic.c), NA, ifelse(ethnic.c=="White","White","Others")), c("White","Others")),
	edu.c=factor(case_when(grepl("(^|\\|)1(\\||$)",edu)~"high", grepl("(^|\\|)(2|5|6)(\\||$)",edu)~"middle", grepl("(^|\\|)(3|4|7)(\\||$)",edu)~"low"), c("low","middle","high")),
	center=as.character(center),
	# Assessment-centre region definition retained from the original 01_data.R.
	# region/region_code are metadata for geographic 10-fold ML splitting only;
	# they are NOT added to M1/M2/M3 and are NOT used as strata in the current Cox models.
	region=case_when(
		center %in% c("c11012","c11020","c11018") ~ "London",
		center %in% c("c11001","c11016","c10003","c11024","c11025","c11008") ~ "North West",
		center %in% c("c11017","c11009","c11027") ~ "North East",
		center %in% c("c11010","c11014") ~ "Yorkshire",
		center %in% c("c11021","c11006") ~ "West Midlands",
		center %in% c("c11013") ~ "East Midlands",
		center %in% c("c11002","c11007","c11026") ~ "South East",
		center %in% c("c11011","c11028") ~ "South West",
		center %in% c("c11005","c11004") ~ "Scotland",
		center %in% c("c11003","c11022","c11023") ~ "Wales",
		TRUE ~ NA_character_
	),
	region=factor(region, levels=region_levels),
	region_code=as.integer(region)-1L,
	smoke.c=factor(smoke_status, 0:2, c("never","previous","current")),
	alcohol.c=factor(alcohol_status, 0:2, c("never","previous","current")),
	days_pa_mod.c=factor(ifelse(days_pa_mod %in% 0:7, days_pa_mod, NA_integer_), 0:7), # Field 884
	age.b=factor(ifelse(age<60,"<60","≥60"), c("<60","≥60")),
	bmi.b=factor(ifelse(bmi<25,"<25",">=25"), c("<25",">=25"))
)
# First-degree family history of stroke: Fields 20107/20110/20111, Coding 1010 Group 1.
dat <- dose_add_family_history(dat)
fh_cols <- c("fhist_father","fhist_mother","fhist_siblings"); stopifnot(all(fh_cols %in% names(dat)))
fh_yes <- Reduce(`|`, lapply(dat[fh_cols], has_code, code=2))
fh_info <- Reduce(`|`, lapply(dat[fh_cols], \(x) !is.na(x) & grepl("(^|\\|)(-17|1|2|6|8|9|10)(\\||$)", as.character(x))))
dat$fhist_stroke <- factor(ifelse(fh_yes,"Yes",ifelse(fh_info,"No",NA)), c("No","Yes"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Baseline diseases, follow-up and outcomes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dxs <- c("primary_ahtn","t2dm","cvd"); drugs <- c("drug.lipid","drug.htn","drug.dm")
date.cvd <- paste0("fod_icd10_", c("cvd_cad","cvd_hfail","cvd_afib","cvd_pad")); stopifnot(all(date.cvd %in% names(dat)), all(drugs %in% names(dat)))
dat$fod_icd10_cvd <- do.call(pmin, c(dat[date.cvd], list(na.rm=TRUE)))
dat$fod_icd10_cvd[!is.finite(as.numeric(dat$fod_icd10_cvd))] <- as.Date(NA)
for(v in dxs) dat[[v]] <- !is.na(dat[[paste0("fod_icd10_",v)]]) & dat[[paste0("fod_icd10_",v)]] <= dat$start_date

dat <- dat %>% mutate(
	prior_qi_stroke=!is.na(fod_icd10_qi_stroke) & fod_icd10_qi_stroke<=start_date,
	across(all_of(dxs), ~factor(ifelse(.x,"Yes","No"), c("No","Yes"))),
	across(all_of(drugs), ~factor(.x, 0:1, c("No","Yes"))),
	end0=pmin(date_lost,date_death,as.Date(date_follow_end),na.rm=TRUE)
) %>% filter(!prior_qi_stroke, end0>start_date)

Ys <- c("qi_stroke","cvd_stroke_i","cvd_stroke_ih","cvd_stroke_sh")
out_map <- c(qi_stroke="Overall stroke", cvd_stroke_i="Ischemic stroke", cvd_stroke_ih="Intracerebral hemorrhage", cvd_stroke_sh="Subarachnoid hemorrhage")
# Outcome hierarchy locked for DOSE development on 2026-09-04.
# Do not switch the primary outcome in response to post hoc P-value or AUC variation.
primary_development_outcome <- "cvd_stroke_i"
key_secondary_outcome <- "qi_stroke"
secondary_subtype_outcomes <- c("cvd_stroke_ih", "cvd_stroke_sh")
for(Y in Ys){
	old <- grep(paste0("^",Y,"\\.(Yt2e|Yr2e|t2e|r2e|b2e)$"), names(dat), value=TRUE); if(length(old)) dat[old] <- NULL
	dat <- t2e(dat, NA, paste0("fod_icd10_",Y), "birth_date", "start_date", "date_lost", "date_death", date_follow_end, Y, "year")
}
dat.base <- dat %>% mutate(across(ends_with("Yt2e"), as.integer))
rm(dat)
invisible(gc())
cat("主分析人数 =", nrow(dat.base), "\n") # 186345

# Geographic fold QC for machine learning.
# Expected coding:
# 0 London; 1 North West; 2 North East; 3 Yorkshire; 4 West Midlands;
# 5 East Midlands; 6 South East; 7 South West; 8 Scotland; 9 Wales.
if(anyNA(dat.base$region_code)){
	bad_center <- dat.base %>% filter(is.na(region_code)) %>% count(center, sort=TRUE)
	print(tibble::as_tibble(bad_center), n=Inf)
	stop("Missing region_code found; check unmapped assessment centre codes above.")
}
stopifnot(identical(sort(unique(dat.base$region_code)), 0:9))
tab_region <- dat.base %>% dplyr::count(region, region_code, name="n") %>% arrange(region_code)
print(tibble::as_tibble(tab_region), n=Inf)
write.csv(tab_region, file.path(DIR_DATA, "region_code_summary.csv"), row.names=FALSE)

# Event summary
event_stat <- function(Y){
	y <- paste0(Y,".Yt2e"); t <- paste0(Y,".t2e"); case <- dat.base[[y]]==1; ncase <- sum(case,na.rm=TRUE); n <- sum(!is.na(dat.base[[y]])); py <- sum(dat.base[[t]],na.rm=TRUE)
	tibble(outcome=out_map[Y], total=n, case=ncase, case_pct=round(ncase/n*100,1), person_years=round(py,1), ir_1e5=round(ncase/py*1e5,1),
		mean_ageonset=round(mean(dat.base$age[case]+dat.base[[t]][case],na.rm=TRUE),1), mean_fu=round(mean(dat.base[[t]],na.rm=TRUE),1),
		sd_fu=round(sd(dat.base[[t]],na.rm=TRUE),1), median_fu=round(median(dat.base[[t]],na.rm=TRUE),1))
}
tab_event <- map_dfr(Ys,event_stat); print(tab_event); write.csv(tab_event,file.path(DIR_DATA,"tab_event.csv"),row.names=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Models and missing covariates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# M1: age + sex + ethnicity + education + TDI + energy + smoking + PA + BMI (+ alcohol for patterns without alcohol)
# M2: M1 + family history + HTN + T2DM ; primary/full model
# M3: M2 + lipid/HTN/DM medication; sensitivity model
# FoodWAS does not adjust for alcohol.
dxs <- c("primary_ahtn","t2dm") # 不校正"cvd"
covs.m1 <- c("age","sex.b","ethnic.c3","edu.c","tdi","energy.kJ","smoke.c","days_pa_mod.c","bmi")
covs.alcohol <- "alcohol.c"; covs.m2.add <- c("fhist_stroke",dxs); covs.m3.add <- drugs
covs.all <- unique(c(covs.m1,covs.alcohol,covs.m2.add,covs.m3.add)); stopifnot(all(covs.all %in% names(dat.base)))

covs_miss <- tibble(variable=covs.all) %>% mutate(
	class=map_chr(variable,~class(dat.base[[.x]])[1]), n_total=nrow(dat.base), n_missing=map_int(variable,~sum(is.na(dat.base[[.x]]))),
	pct_missing=round(n_missing/n_total*100,2), n_nonmissing=n_total-n_missing,
	n_level=map_int(variable,~if(is.factor(dat.base[[.x]])||is.character(dat.base[[.x]])) n_distinct(dat.base[[.x]],na.rm=TRUE) else NA_integer_)
) %>% arrange(desc(pct_missing),variable)
print(covs_miss,n=Inf); write.csv(covs_miss,file.path(DIR_DATA,"covariate_missingness.csv"),row.names=FALSE)
write_xlsx(covs_miss %>% filter(n_missing>0) %>% transmute(Covariates=variable,`No. of sample with missing covariates`=n_missing,`Missing (%)`=pct_missing), file.path(DIR_DATA,"Table_S_Missing_Covariates.xlsx"))

mode_imp <- function(x){x[is.na(x)] <- names(which.max(table(x,useNA="no"))); x}
dat.base.imp <- dat.base %>% mutate(
	across(c(days_pa_mod.c,fhist_stroke), ~fct_na_value_to_level(.x,"Unknown")),
	across(c(edu.c,ethnic.c3,smoke.c,alcohol.c), mode_imp),
	across(c(bmi,tdi), ~replace(.x,is.na(.x),mean(.x,na.rm=TRUE)))
)
stopifnot(sum(is.na(dat.base.imp[covs.all]))==0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
analysis_config <- list(
	Ys=Ys, out_map=out_map,
	primary_development_outcome=primary_development_outcome,
	key_secondary_outcome=key_secondary_outcome,
	secondary_subtype_outcomes=secondary_subtype_outcomes,
	covs.m1=covs.m1, covs.alcohol=covs.alcohol,
	covs.m2.add=covs.m2.add, covs.m3.add=covs.m3.add, covs.all=covs.all,
	dxs=dxs, drugs=drugs,
	primary_model="M2", sensitivity_model="M3",
	region_var="region", fold_var="region_code",
	date_follow_end=date_follow_end
)
dose_save_rds(dat.base,file.path(DIR_DATA,"dat.base.rds")); dose_save_rds(dat.base.imp,file.path(DIR_DATA,"dat.base.imp.rds")); dose_save_rds(analysis_config,file.path(DIR_DATA,"analysis_config.rds"))
cat("Saved:",DIR_DATA,"\n")
}))

# ===========================================================================
# foodwas / patterns
# ===========================================================================
dose_register_step("foodwas", "patterns", expression({
pacman::p_load(tidyverse, survival, broom)
DIR_PROJECT <- dose_output_dir()
DIR_DATA <- file.path(DIR_PROJECT,"ukb/00data"); DIR_OUT <- file.path(DIR_PROJECT,"ukb/01pattern"); dir.create(DIR_OUT,recursive=TRUE,showWarnings=FALSE)
# all.rds already contains diet.rds, therefore dat.base.imp already contains diet.*.sum scores.
cfg <- readRDS(file.path(DIR_DATA,"analysis_config.rds")); dat <- readRDS(file.path(DIR_DATA,"dat.base.imp.rds"))
spec <- tribble(
	~diet,~alcohol_in_score,~include,
	"medi24",TRUE,TRUE,"mind",TRUE,TRUE,"ahei",TRUE,TRUE,"redii",TRUE,TRUE,"fds",TRUE,TRUE,
	"dash",FALSE,TRUE,"hpdi",FALSE,TRUE,"phdi",FALSE,TRUE,"digm",FALSE,TRUE,"hlcd",FALSE,TRUE,"hlfd",FALSE,TRUE,"modern",FALSE,TRUE,
	"medito",FALSE,FALSE,"maha",TRUE,FALSE
) %>% mutate(variable=paste0("diet.",diet,".sum"),available=variable %in% names(dat))
use <- spec %>% filter(include,available); write.csv(spec,file.path(DIR_OUT,"pattern_spec.csv"),row.names=FALSE)
if(!nrow(use)) stop("No requested diet.*.sum variables were found in dat.base.imp.rds")
needed <- unique(c("eid", use$variable, cfg$covs.all,
    unlist(lapply(cfg$Ys, function(y) paste0(y, c(".t2e", ".Yt2e"))))))
dat <- dat %>% select(all_of(needed))
invisible(gc())

model_covs <- function(d,m){
	m1 <- c(cfg$covs.m1,if(!use$alcohol_in_score[match(d,use$diet)]) cfg$covs.alcohol)
	switch(m,M1=m1,M2=c(m1,cfg$covs.m2.add),M3=c(m1,cfg$covs.m2.add,cfg$covs.m3.add)) %>% unique()
}

# Pattern exposure: P10-P90, per 1-SD, and quartiles (Q1-Q4; MYTH-like quantile cutpoints).
res.cont <- list(); res.cat <- list(); k <- 1L
for(d in use$diet){
	x0 <- paste0("diet.",d,".sum"); xp <- paste0("diet.",d,".p1090"); xz <- paste0("diet.",d,".z"); xq <- paste0("diet.",d,".q")
	dat[[xp]] <- std_10_90(dat[[x0]]); dat[[xz]] <- std_z(dat[[x0]]); dat[[xq]] <- qcat(dat[[x0]])
	for(Y in cfg$Ys) for(m in c("M1","M2","M3")){
		co <- model_covs(d,m)
		res.cont[[k]] <- bind_rows(
			cox_cont(dat,Y,xp,co) %>% mutate(scale="P10-P90"),
			cox_cont(dat,Y,xz,co) %>% mutate(scale="Per 1-SD")
		) %>% mutate(diet=d,outcome=Y,model=m,.before=1)
		r <- cox_cat(dat,Y,xq,x0,co,FALSE)
		if(!is.null(r)) res.cat[[k]] <- r$contrast %>%
			mutate(diet=d,outcome=Y,model=m,P_global=r$global$P_global[1],P_trend=r$global$P_trend[1],.before=1)
		k <- k+1L
	}
}
write.csv(bind_rows(res.cont),file.path(DIR_OUT,"pattern_continuous.csv"),row.names=FALSE)
write.csv(bind_rows(res.cat),file.path(DIR_OUT,"pattern_categories.csv"),row.names=FALSE)
cat("Pattern complete; primary model = M2; patterns =",paste(use$diet,collapse=", "),"\n")
}))

# ===========================================================================
# foodwas / items
# ===========================================================================
dose_register_step("foodwas", "items", expression({
pacman::p_load(readxl, tidyverse, survival, broom)

# ======================================================================
# foodwas.R
# Item-level FoodWAS for incident ischemic stroke
#
# Association analysis follows the original 01_data.R / dementia-style code:
#   A. item-level FoodWAS: non-intake vs intake; M1 + M2; FDR across items
#
# M2 is primary. No M3, no alcohol adjustment, no strata(region).
# This script is independent of foodwas.R; output names are unchanged.
# The two post-ML combined food groups are intentionally NOT created here.
# ======================================================================

dir0 <- dose_drive()
indir <- dose_ukb_dir()
DIR_PROJECT <- dose_output_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_OUT <- file.path(DIR_PROJECT, "ukb/02foodwas")
DIR_ITEM <- file.path(DIR_OUT, "item")
DIR_GROUP <- file.path(DIR_OUT, "group")
DIR_ML <- file.path(DIR_PROJECT, "ukb/03ml")
dir.create(DIR_ITEM, recursive=TRUE, showWarnings=FALSE)
dir.create(DIR_GROUP, recursive=TRUE, showWarnings=FALSE)
dir.create(DIR_ML, recursive=TRUE, showWarnings=FALSE)

source(file.path(dose_helper_dir(), "phenotype.R"))

# ----------------------------------------------------------------------
# 1. Cohort, detailed food items, outcome and covariates
# ----------------------------------------------------------------------
cfg <- readRDS(file.path(DIR_DATA, "analysis_config.rds"))
base <- readRDS(file.path(DIR_DATA, "dat.base.imp.rds")) %>% mutate(eid=as.character(eid))
food <- readRDS(file.path(indir, "Rdata/diet.fudan.rds")) %>% mutate(eid=as.character(eid)) %>% semi_join(base %>% select(eid), by="eid")

stopifnot(
	all(c("region", "region_code") %in% names(base)),
	!anyNA(base$region_code),
	identical(sort(unique(base$region_code)), 0:9),
	!anyDuplicated(base$eid),
	!anyDuplicated(food$eid),
	nrow(food)==nrow(base),
	setequal(food$eid, base$eid)
)

field <- read_excel(file.path(dose_files_dir(), "foods.xlsx"), sheet="Sheet4") # %>% filter(doseFood>0)
# field$Food.Group[field$name=="p103100"] <- "Organ_meats"
stopifnot(
	!anyDuplicated(field$name),
	all(field$name %in% names(food)),
	all(!is.na(field$Food.Group))
)

# Prespecified primary development outcome; do not switch based on post hoc results.
Y <- "cvd_stroke_i"
y <- paste0(Y, ".Yt2e")
t <- paste0(Y, ".t2e")
cat("Outcome:", Y)

# FoodWAS does not adjust for alcohol.
covs.m1 <- cfg$covs.m1
covs.m2 <- unique(c(covs.m1, cfg$covs.m2.add))
stopifnot(!"alcohol.c" %in% covs.m1, !"alcohol.c" %in% covs.m2)

keep <- unique(c("eid", "region", "region_code", t, y, covs.m2))

dat.item <- base %>%
	select(all_of(keep)) %>%
	left_join(food %>% select(eid, all_of(field$name)), by="eid")

stopifnot(
	nrow(dat.item)==nrow(base),
	all(c(t, y, covs.m2, field$name) %in% names(dat.item))
)
rm(base, food)
invisible(gc())

write.csv(field, file.path(DIR_OUT, "food_group_map.csv"), row.names=FALSE)
cat("FoodWAS cohort =", nrow(dat.item), "; items =", nrow(field), "\n")

# ----------------------------------------------------------------------
# 2. Item-level FoodWAS: 0 = non-intake, 1 = intake
#    This follows the original dementia-style FoodWAS logic.
# ----------------------------------------------------------------------
run_item_cox <- function(data, exposure, Y, covs){
	y <- paste0(Y, ".Yt2e")
	t <- paste0(Y, ".t2e")

	d <- data %>%
		select(all_of(c(t, y, exposure, covs))) %>%
		rename(pheno=all_of(exposure)) %>%
		drop_na()

	d$pheno <- ifelse(d$pheno==0, 0L, 1L)
	if(length(unique(d$pheno))<2) return(NULL)

	st <- d %>%
		group_by(pheno) %>%
		summarise(total=n(), case=sum(.data[[y]], na.rm=TRUE), .groups="drop")

	fit <- coxph(
		reformulate(c("pheno", covs), response=paste0("Surv(", t, ",", y, ")")),
		data=d
	)

	td <- broom::tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>%
		filter(term=="pheno")

	zph_p <- tryCatch(
		cox.zph(fit)$table["pheno", "p"],
		error=function(e) NA_real_
	)

	tibble(
		exposure=exposure,
		outcome=Y,
		class=1:2,
		level=c("Non-intake", "Intake"),
		total=st$total[match(c(0,1), st$pheno)],
		case=st$case[match(c(0,1), st$pheno)],
		HR=c(1, td$estimate),
		LCI=c(1, td$conf.low),
		HCI=c(1, td$conf.high),
		P=c(1, td$p.value),
		res.zph=zph_p
	)
}

run_item_model <- function(model, covs){
	res <- lapply(field$name, function(v) run_item_cox(dat.item, v, Y, covs)) %>%
		bind_rows() %>%
		mutate(model=model, .after=outcome) %>%
		arrange(exposure, class)

	# Equivalent to the original FDR(p, 2): adjust the actual intake-vs-no-intake
	# coefficient P values across food items; reference rows remain FDR=1.
	res$FDR <- 1
	idx <- which(res$class==2 & is.finite(res$P))
	if(length(idx)) res$FDR[idx] <- p.adjust(res$P[idx], method="fdr")
	res
}

item_m1 <- run_item_model("M1", covs.m1)
item_m2 <- run_item_model("M2", covs.m2)
item_all <- bind_rows(item_m1, item_m2)

write.csv(item_m1, file.path(DIR_ITEM, paste0("Food_", Y, "_M1.csv")), row.names=FALSE)
write.csv(item_m2, file.path(DIR_ITEM, paste0("Food_", Y, ".csv")), row.names=FALSE)
write.csv(item_all, file.path(DIR_ITEM, paste0("Food_", Y, "_all_models.csv")), row.names=FALSE)

item_sig_m2 <- item_m2 %>%
	filter(class==2, FDR<0.05) %>%
	arrange(FDR, P)
write.csv(item_sig_m2, file.path(DIR_ITEM, "Need_Item_FDR.csv"), row.names=FALSE)

cat("Item-level M2 FDR<0.05 =", nrow(item_sig_m2), "\n")


cat("Item-level FoodWAS complete\n")
cat(" - cohort:", nrow(dat.item), "\n")
cat(" - detailed food items:", nrow(field), "\n")
cat(" - item M1 FDR<0.05:", sum(item_m1$class==2 & item_m1$FDR<0.05, na.rm=TRUE), "\n")
cat(" - item M2 FDR<0.05:", nrow(item_sig_m2), "\n")
cat(" - output folder:", DIR_ITEM, "\n")
}))

# ===========================================================================
# foodwas / groups
# ===========================================================================
dose_register_step("foodwas", "groups", expression({
# ======================================================================
# foodwas.R
# Food-group FoodWAS + candidate preparation for incident ischemic stroke
#
# Association analysis follows the original 01_data.R / dementia-style code:
#   B. food-group FoodWAS: quintile-style grouping with fallback; all-reference
#      Cox; M1 + M2; raw coefficient P values retained (no group-level FDR)
#
# Candidate preparation for downstream ML is a hybrid of the original
# dementia-style FoodWAS screen and the MYTH manuscript ML preparation:
#   M2 group Cox P < 0.05 in any non-reference comparison
#       -> candidate food groups
#       -> Spearman correlation / hierarchical clustering
#       -> remove highly dependent groups at |rho| >= 0.60
#       -> one group per cluster
#       -> residualize each selected food on all M2 covariates
#       -> z-standardize residuals
#       -> Need_Food_lst.txt + Food_resid_z.csv + ML_input.csv for ML
#
# M2 is primary. No M3, no alcohol adjustment, no strata(region).
# This script is independent of foodwas.R; group/ML output names are unchanged.
# The two post-ML combined food groups are intentionally NOT created here.
# ======================================================================
pacman::p_load(readxl, tidyverse, survival, broom)

dir0 <- dose_drive()
indir <- dose_ukb_dir()
DIR_PROJECT <- dose_output_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_OUT <- file.path(DIR_PROJECT, "ukb/02foodwas")
	DIR_ITEM <- file.path(DIR_OUT, "item")
	DIR_GROUP <- file.path(DIR_OUT, "group")
DIR_ML <- file.path(DIR_PROJECT, "ukb/03ml")
	dir.create(DIR_ITEM, recursive=TRUE, showWarnings=FALSE)
	dir.create(DIR_GROUP, recursive=TRUE, showWarnings=FALSE)
	dir.create(DIR_ML, recursive=TRUE, showWarnings=FALSE)

source(file.path(dose_helper_dir(), "phenotype.R"))

# ----------------------------------------------------------------------
# Cohort, detailed food items, outcome and covariates
# ----------------------------------------------------------------------
cfg <- readRDS(file.path(DIR_DATA, "analysis_config.rds"))
base <- readRDS(file.path(DIR_DATA, "dat.base.imp.rds")) %>% mutate(eid=as.character(eid))
food <- readRDS(file.path(indir, "Rdata/diet.fudan.rds")) %>% mutate(eid=as.character(eid)) %>% semi_join(base %>% select(eid), by="eid")
field <- read_excel(file.path(dose_files_dir(), "foods.xlsx"), sheet="Sheet4") # %>% filter(doseFood>0)
#	field$Food.Group[field$name=="p103100"] <- "Organ_meats"

Y <- "cvd_stroke_i"
cat("Outcome:", Y)
y <- paste0(Y, ".Yt2e"); t <- paste0(Y, ".t2e")

# FoodWAS does not adjust for alcohol.
covs.m1 <- cfg$covs.m1; covs.m2 <- unique(c(covs.m1, cfg$covs.m2.add))
	stopifnot(!"alcohol.c" %in% covs.m1, !"alcohol.c" %in% covs.m2)
keep <- unique(c("eid", "region", "region_code", t, y, covs.m2))
write.csv(field, file.path(DIR_OUT, "food_group_map.csv"), row.names=FALSE)
cat("Group FoodWAS cohort =", nrow(base), "; items used to build groups =", nrow(field), "\n")

# ----------------------------------------------------------------------
# Build ORIGINAL food groups
# ----------------------------------------------------------------------
grp.link <- split(field$name, field$Food.Group)
grp <- food %>% select(eid)
for(g in names(grp.link)) grp[[g]] <- rowSums1(food[grp.link[[g]]], per=.8)
# grp$Citrus_ApplesPears_Driedfruit <- grp$Citrus + grp$Apples_Pears + grp$Dried_fruit
group_vars <- setdiff(names(grp), "eid")
dat.group <- base %>% select(all_of(keep)) %>% left_join(grp, by="eid")
base <- base %>% select(eid) # Retain only IDs for the final ML alignment checks.
rm(food, grp)
invisible(gc())
dose_save_rds(dat.group, file.path(DIR_OUT, "dat.group.rds"))
cat("Food groups =", length(group_vars), "\n")

# ----------------------------------------------------------------------
# quintile-style food-group categorization
#    - ordinary quintiles when cut points are unique
#    - fewer categories when quantile cut points collapse
#    - binary 0 vs >0 when only two unique cut points remain
# ----------------------------------------------------------------------
func_quin <- function(x){
	x <- as.numeric(x)
	ok <- is.finite(x)
	out <- rep(NA_integer_, length(x))
	if(sum(ok)<2) return(factor(out))

	br0 <- quantile(x[ok], probs=seq(0,1,1/5), na.rm=TRUE, names=FALSE)
	bru <- unique(br0)

	if(length(bru)==2){
		out[ok] <- ifelse(x[ok]>0, 2L, 1L)
	} else {
		br <- if(length(bru)<6) c(min(x[ok])-1, bru) else bru
		out[ok] <- as.integer(cut(x[ok], breaks=br, include.lowest=TRUE, labels=FALSE))
	}

	lv <- sort(unique(out[!is.na(out)]))
	out <- match(out, lv)
	factor(out, levels=seq_along(lv), labels=paste0("C",seq_along(lv)), ordered=FALSE)
}

group_q <- dat.group %>%
	select(eid, all_of(c(t, y, covs.m2, group_vars)))

for(g in group_vars) group_q[[g]] <- func_quin(group_q[[g]])

cat_meta <- tibble(
	exposure=group_vars,
	n_cat=sapply(group_vars, function(g) nlevels(group_q[[g]])),
	zero_pct=sapply(group_vars, function(g) mean(dat.group[[g]]==0, na.rm=TRUE))
)
write.csv(cat_meta, file.path(DIR_GROUP, "Group_categories.csv"), row.names=FALSE)

# ----------------------------------------------------------------------
# 5. Statistics: median / range / n / case for each intake category
#    Preserve the original stat() + for + rbind() organization.
# ----------------------------------------------------------------------
stat <- function(x, data_q, data_num, Y){
	y <- paste0(Y, ".Yt2e")

	temp <- data_q %>%
		select(eid, all_of(y), quantile=all_of(x)) %>%
		left_join(data_num %>% select(eid, Food=all_of(x)), by="eid") %>%
		mutate(quantile=as.character(quantile))

	temp %>%
		group_by(quantile) %>%
		summarise(
			Median=round(median(Food, na.rm=TRUE), 2),
			Range=paste(round(range(Food, na.rm=TRUE), 2), collapse="~"),
			n=n(),
			case=sum(.data[[y]], na.rm=TRUE),
			.groups="drop"
		) %>%
		mutate(pheno=x, .before="quantile") %>%
		mutate(`Consumption level`=paste0(Median, " (", Range, ")"), .after="quantile")
}

stat_data <- data.frame()
for(i in group_vars){
	stat0 <- stat(i, group_q, dat.group, Y)
	stat_data <- rbind(stat_data, stat0)
}
write.csv(stat_data, file.path(DIR_GROUP, "Group_stat.csv"), row.names=FALSE)

# ----------------------------------------------------------------------
# 6. Food-group all-reference Cox, M1 + M2
#    Keep the original all-reference purpose but only run observed references
#    once (avoids duplicate ref=4/5 runs in food groups with <4/5 categories).
# ----------------------------------------------------------------------
run_group_cox <- function(data_q, group_vars, Y, covs, model){
	y <- paste0(Y, ".Yt2e"); t <- paste0(Y, ".t2e")
	out <- list()
	z <- 1L

	for(component in group_vars){
		d0 <- data_q %>% select(all_of(c(t, y, component, covs))) %>% rename(pheno=all_of(component)) %>% drop_na()

		d0$pheno <- droplevels(factor(d0$pheno, ordered=FALSE))
		lv0 <- levels(d0$pheno)
		if(length(lv0)<2) next

		# PH test is reference-invariant; compute once per food/model.
		fit.ph <- coxph(
			reformulate(c("pheno", covs), response=paste0("Surv(", t, ",", y, ")")),
			data=d0
		)
		zph_p <- tryCatch(cox.zph(fit.ph)$table["pheno", "p"], error=function(e) NA_real_)

		for(ref in lv0){
			d <- d0
			d$pheno <- relevel(d$pheno, ref=ref)

			fit <- coxph(
				reformulate(c("pheno", covs), response=paste0("Surv(", t, ",", y, ")")),
				data=d
			)

			td <- broom::tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>%
				filter(grepl("^pheno", term))

			lv <- levels(d$pheno)
			res <- tibble(
				Ref=ref,
				exposure=component,
				outcome=Y,
				model=model,
				class=as.character(lv),
				HR=1,
				LCI=1,
				HCI=1,
				P=1,
				res.zph=zph_p
			)

			if(nrow(td)){
				idx <- match(sub("^pheno", "", td$term), res$class)
				res$HR[idx] <- td$estimate
				res$LCI[idx] <- td$conf.low
				res$HCI[idx] <- td$conf.high
				res$P[idx] <- td$p.value
			}

			out[[z]] <- res
			z <- z + 1L
		}
	}

	bind_rows(out)
}

group_m1 <- run_group_cox(group_q, group_vars, Y, covs.m1, "M1")
	write.csv(group_m1, file.path(DIR_GROUP, paste0("Group_", Y, "_M1.csv")), row.names=FALSE)
group_m2 <- run_group_cox(group_q, group_vars, Y, covs.m2, "M2")
	write.csv(group_m2, file.path(DIR_GROUP, paste0("Group_", Y, ".csv")), row.names=FALSE)
group_all <- bind_rows(group_m1, group_m2)
	write.csv(group_all, file.path(DIR_GROUP, paste0("Group_", Y, "_all_models.csv")), row.names=FALSE)

# ----------------------------------------------------------------------
# 7. Permissive M2 candidate screen for ML: any non-reference P < 0.05
# ----------------------------------------------------------------------
sig_m2 <- group_m2 %>% filter(Ref!=class, P<0.05) %>% arrange(P)
need_lst <- sig_m2 %>% distinct(exposure) %>% pull(exposure); need_lst
write.csv(sig_m2, file.path(DIR_GROUP, "Need_Food_P005.csv"), row.names=FALSE)
writeLines(need_lst, file.path(DIR_GROUP, "Need_Food_before_cluster.txt"))
cat("M2 group FoodWAS P<0.05 before clustering =", length(need_lst), "\n")

# ----------------------------------------------------------------------
# 8. MYTH-like ML preparation
#    Manuscript: exclude food groups with high mutual dependence
#    (|Spearman rho| > 0.60), then residualize selected food intake on all
#    covariates and standardize the residuals for ML.
# ----------------------------------------------------------------------
cor_threshold <- 0.60
cut_height <- 1 - cor_threshold

if(length(need_lst)==1){
	cor_mat <- matrix(1, nrow=1, ncol=1, dimnames=list(need_lst, need_lst))
	clusters <- setNames(1L, need_lst)
	need_lst1 <- need_lst
	cor_pairs <- tibble(food1=character(), food2=character(), rho=numeric(), abs_rho=numeric())
} else {
	cor_mat <- cor(dat.group %>% select(all_of(need_lst)), method="spearman", use="pairwise.complete.obs")

	# dist = 1-|rho|.  Therefore |rho| >= 0.60 corresponds to distance <= 0.40.
	dist_mat <- as.dist(1 - abs(cor_mat))
	hc <- hclust(dist_mat, method="average")
	clusters <- cutree(hc, h=cut_height)

	pdf(file.path(DIR_GROUP, "Group_clustering.pdf"), width=10, height=7)
	plot(hc, main=paste0("Hierarchical clustering of candidate food groups (|rho| >= ", cor_threshold, ")"), xlab="", sub="", cex=.8)
	abline(h=cut_height, col="red", lty=2)
	bounded_dev_off()

	# Keep one food per cluster, matching the MYTH public-code organization.
	need_lst1 <- tapply(names(clusters), clusters, function(x) x[1]) %>% unname()

	# Pairwise correlation table without factor-comparison warnings.
	idx <- which(upper.tri(cor_mat), arr.ind=TRUE)
	cor_pairs <- tibble(
		food1=rownames(cor_mat)[idx[,1]], food2=colnames(cor_mat)[idx[,2]], rho=cor_mat[idx]
	) %>% mutate(abs_rho=abs(rho)) %>% arrange(desc(abs_rho))
}

write.csv(cor_mat, file.path(DIR_GROUP, "Group_correlation.csv"), row.names=TRUE)
write.csv(cor_pairs, file.path(DIR_GROUP, "Group_correlation_pairs.csv"), row.names=FALSE)
write.csv(
	tibble(exposure=names(clusters), cluster=as.integer(clusters)),
	file.path(DIR_GROUP, "Group_clusters.csv"),
	row.names=FALSE
)
writeLines(need_lst1, file.path(DIR_GROUP, "Need_Food_lst.txt"))

# ----------------------------------------------------------------------
# 9. Residual method + z-standardization for downstream ML
#    For each selected food group: food ~ all M2 covariates.
#    Standardized residuals are used as ML features.
# ----------------------------------------------------------------------
residualize_food <- function(x, data, covs){
	d <- data %>%
		select(all_of(c(x, covs))) %>%
		rename(food=all_of(x))

	fit <- lm(
		reformulate(covs, response="food"),
		data=d,
		na.action=na.exclude
	)

	r <- residuals(fit)
	s <- sd(r, na.rm=TRUE)
	z <- if(is.finite(s) && s>0) (r - mean(r, na.rm=TRUE))/s else rep(NA_real_, length(r))

	list(
		z=as.numeric(z),
		info=tibble(
			exposure=x,
			n=sum(is.finite(r)),
			R2=summary(fit)$r.squared,
			resid_mean=mean(r, na.rm=TRUE),
			resid_sd=sd(r, na.rm=TRUE)
		)
	)
}

food_resid_z <- dat.group %>% select(eid)
resid_info <- list()
for(g in need_lst1){
	rz <- residualize_food(g, dat.group, covs.m2)
	food_resid_z[[g]] <- rz$z
	resid_info[[g]] <- rz$info
}
resid_info <- bind_rows(resid_info)

write.csv(food_resid_z, file.path(DIR_GROUP, "Food_resid_z.csv"), row.names=FALSE)
write.csv(resid_info, file.path(DIR_GROUP, "Food_residualization_QC.csv"), row.names=FALSE)

# ML-ready dataset: geographic fold + outcome/time + residualized standardized food features.
# This is the direct input for the s2 Python module; no separate s2_0 preparation step is needed.
dat.ml <- dat.group %>% select(eid, region, region_code, all_of(c(t, y))) %>% left_join(food_resid_z, by="eid")

stopifnot(
	nrow(dat.ml)==nrow(base),
	setequal(dat.ml$eid, base$eid),
	!anyNA(dat.ml$region_code),
	identical(sort(unique(dat.ml$region_code)), 0:9),
	!anyNA(dat.ml[, c(t, y, need_lst1), drop=FALSE])
)

dose_save_rds(dat.ml, file.path(DIR_GROUP, "dat.foodwas.ml.rds"))
write.csv(dat.ml, file.path(DIR_ML, "ML_input.csv"), row.names=FALSE)

fold_stat <- dat.ml %>% group_by(region, region_code) %>% summarise(n=n(), case=sum(.data[[y]]==1), case_pct=round(mean(.data[[y]]==1)*100, 3), .groups="drop") %>% arrange(region_code)
write.csv(fold_stat, file.path(DIR_ML, "Region_code_fold_summary.csv"), row.names=FALSE)

# Keep raw selected groups for audit/reproducibility.
dat.ml.raw <- dat.group %>% select(eid, region, region_code, all_of(c(t, y, covs.m2, need_lst1)))
dose_save_rds(dat.ml.raw, file.path(DIR_GROUP, "dat.foodwas.ml.raw.rds"))

cat("Group-level FoodWAS + ML preparation complete\n")
cat(" - cohort:", nrow(dat.group), "\n")
cat(" - detailed food items used to build groups:", nrow(field), "\n")
cat(" - original food groups:", length(group_vars), "\n")
cat(" - M2 group P<0.05 before clustering:", length(need_lst), "\n")
cat(" - clustering correlation threshold |rho| >=", cor_threshold, "\n")
cat(" - groups after clustering:", length(need_lst1), "\n")
if(nrow(cor_pairs)>0){
	cat(" - max |rho| among candidates:", round(max(cor_pairs$abs_rho, na.rm=TRUE), 3), "\n")
	cat(" - pairs with |rho| >=", cor_threshold, ":", sum(cor_pairs$abs_rho>=cor_threshold, na.rm=TRUE), "\n")
}
cat(" - final ML list:", paste(need_lst1, collapse=", "), "\n")
cat(" - residualized ML matrix:", file.path(DIR_GROUP, "Food_resid_z.csv"), "\n")
cat(" - direct s2 ML input:", file.path(DIR_ML, "ML_input.csv"), "\n")
cat(" - Region_code folds: 0-9; summary:", file.path(DIR_ML, "Region_code_fold_summary.csv"), "\n")
}))

# ===========================================================================
# foodwas / trend
# ===========================================================================
dose_register_step("foodwas", "trend", expression({
pacman::p_load(tidyverse, survival, broom)

# ======================================================================
# foodwas.R
# P for trend across food-group intake categories.
#
# MYTH converted ordered intake categories to numeric values and fitted a
# Cox model.  Here each participant is assigned the MEDIAN raw intake of
# their category, which preserves the actual dose spacing and is standard
# for nutritional P-trend analyses.
#
# M1 + M2 are reported; M2 is primary. No alcohol, no M3, no strata(region).
# ======================================================================

DIR_PROJECT <- dose_output_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_IN <- file.path(DIR_PROJECT, "ukb/02foodwas")
DIR_GROUP <- file.path(DIR_IN, "group")

cfg <- readRDS(file.path(DIR_DATA, "analysis_config.rds"))
dat <- readRDS(file.path(DIR_IN, "dat.group.rds"))
group_map <- read.csv(file.path(DIR_IN, "food_group_map.csv"), check.names=FALSE)
group_vars <- unique(group_map$Food.Group)

# Prespecified primary development outcome; do not switch based on post hoc results.
Y <- "cvd_stroke_i"
y <- paste0(Y, ".Yt2e")
t <- paste0(Y, ".t2e")
cat("Outcome:", Y)

covs.m1 <- cfg$covs.m1
covs.m2 <- unique(c(covs.m1, cfg$covs.m2.add))
stopifnot(!"alcohol.c" %in% covs.m1, !"alcohol.c" %in% covs.m2)

func_quin <- function(x){
	x <- as.numeric(x); ok <- is.finite(x); out <- rep(NA_integer_, length(x))
	if(sum(ok)<2) return(factor(out))
	br0 <- quantile(x[ok], probs=seq(0,1,1/5), na.rm=TRUE, names=FALSE)
	bru <- unique(br0)
	if(length(bru)==2){
		out[ok] <- ifelse(x[ok]>0, 2L, 1L)
	} else {
		br <- if(length(bru)<6) c(min(x[ok])-1, bru) else bru
		out[ok] <- as.integer(cut(x[ok], breaks=br, include.lowest=TRUE, labels=FALSE))
	}
	lv <- sort(unique(out[!is.na(out)])); out <- match(out, lv)
	factor(out, levels=seq_along(lv), labels=paste0("C",seq_along(lv)), ordered=FALSE)
}

trend_one <- function(food, covs, model){
	vars <- unique(c(t, y, food, covs))
	d <- dat[, vars, drop=FALSE]
	d <- d[stats::complete.cases(d), , drop=FALSE]
	names(d)[names(d)==food] <- "intake"

	d$category <- func_quin(d$intake)
	if(nlevels(d$category)<2) return(NULL)

	med <- tapply(d$intake, d$category, median, na.rm=TRUE)
	d$trend <- unname(med[as.character(d$category)])

	fit <- coxph(
		reformulate(c("trend", covs), response=paste0("Surv(", t, ",", y, ")")),
		data=d
	)

	td <- broom::tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>% dplyr::filter(term=="trend")

	tibble(
		exposure=food,
		outcome=Y,
		model=model,
		total=fit$n,
		case=fit$nevent,
		n_cat=nlevels(d$category),
		category_medians=paste(names(med), signif(as.numeric(med),4), sep="=", collapse="; "),
		HR=td$estimate,
		LCI=td$conf.low,
		HCI=td$conf.high,
		P_trend=td$p.value,
		Linear=ifelse(td$p.value<0.05, "*", "")
	)
}

res <- dplyr::bind_rows(
	purrr::map_dfr(group_vars, ~trend_one(.x, covs.m1, "M1")),
	purrr::map_dfr(group_vars, ~trend_one(.x, covs.m2, "M2"))
) %>%
	dplyr::arrange(exposure, model)

write.csv(res, file.path(DIR_GROUP, paste0("Group_trend_", Y, ".csv")), row.names=FALSE)

cat("trend_P complete\n")
cat(" - food groups:", length(group_vars), "\n")
cat(" - M2 P-trend <0.05:", sum(res$model=="M2" & res$P_trend<0.05, na.rm=TRUE), "\n")
cat(" - output:", file.path(DIR_GROUP, paste0("Group_trend_", Y, ".csv")), "\n")
}))

# ===========================================================================
# foodwas / rcs
# ===========================================================================
dose_register_step("foodwas", "rcs", expression({
# ======================================================================
# foodwas.R
# Restricted cubic spline analysis for DOSE food-group intake
#
# v24 changes
# - fixes the all-39 failure caused by calling the non-exported rms::anova.rms;
#   uses stats::anova(fit), which dispatches to the rms method.
# - robust 3-knot rule: full-sample P10/P50/P90; if duplicated because of
#   zero inflation, use consumer (>0) P10/P50/P90, while fitting ALL people.
# - creates TWO RCS figure styles from the SAME fitted model:
#     1) MYTH / rcssci-like: grey intake distribution + red curve/ribbon
#     2) 01_plot.R style: blue histogram + red curve + dashed 95% CI
# - M2 primary model; no alcohol; no M3; no strata(region).
# ======================================================================

pacman::p_load(dplyr, tibble, ggplot2, rms, survival, patchwork)

DIR_PROJECT <- dose_output_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_IN <- file.path(DIR_PROJECT, "ukb/02foodwas")
DIR_GROUP <- file.path(DIR_IN, "group")
DIR_RCS <- file.path(DIR_GROUP, "RCS")
DIR_MYTH <- file.path(DIR_RCS, "myth_style")
DIR_OLD <- file.path(DIR_RCS, "01plot_style")
dir.create(DIR_RCS, recursive=TRUE, showWarnings=FALSE)
dir.create(DIR_MYTH, recursive=TRUE, showWarnings=FALSE)
dir.create(DIR_OLD, recursive=TRUE, showWarnings=FALSE)

RCS_SCOPE <- "all"   # "all", "candidates", or "selected"

cfg <- readRDS(file.path(DIR_DATA, "analysis_config.rds"))
dat <- readRDS(file.path(DIR_IN, "dat.group.rds"))
group_map <- read.csv(file.path(DIR_IN, "food_group_map.csv"), check.names=FALSE, stringsAsFactors=FALSE)
all_groups <- unique(group_map$Food.Group)

# Prespecified primary development outcome; do not switch based on post hoc results.
Y <- "cvd_stroke_i"
y <- paste0(Y, ".Yt2e")
t <- paste0(Y, ".t2e")
cat("Outcome:", Y)
covs <- unique(c(cfg$covs.m1, cfg$covs.m2.add))
stopifnot(!"alcohol.c" %in% covs)

read_list <- function(f){
    if(!file.exists(f)) return(character())
    x <- trimws(readLines(f, warn=FALSE))
    x[nzchar(x)]
}

group_vars <- switch(
    RCS_SCOPE,
    all=all_groups,
    candidates=read_list(file.path(DIR_GROUP, "Need_Food_before_cluster.txt")),
    selected=read_list(file.path(DIR_GROUP, "Need_Food_lst.txt")),
    stop("RCS_SCOPE must be 'all', 'candidates', or 'selected'.")
)
group_vars <- intersect(group_vars, names(dat))
if(!length(group_vars)) stop("No food groups available for RCS_SCOPE = ", RCS_SCOPE)

fmt_p <- function(p){
    ifelse(is.na(p), "NA", ifelse(p < .001, "< 0.001", sprintf("= %.3f", p)))
}

# P10/P50/P90; if zero inflation collapses knots, calculate knots in consumers.
get_k3 <- function(x){
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if(length(x) < 100 || length(unique(x)) < 5) return(NULL)

    k <- as.numeric(stats::quantile(x, c(.10,.50,.90), na.rm=TRUE, names=FALSE, type=2))
    source <- "all_P10_P50_P90"

    if(length(unique(k)) < 3){
        xp <- x[x > 0]
        if(length(xp) >= 20 && length(unique(xp)) >= 3){
            k <- as.numeric(stats::quantile(xp, c(.10,.50,.90), na.rm=TRUE, names=FALSE, type=2))
            source <- "consumer_P10_P50_P90"
        }
    }
    if(length(unique(k)) < 3) return(NULL)
    list(knots=k, source=source)
}

# Extract overall and nonlinear P from rms::anova output.
get_anova_p <- function(an, food){
    rn <- trimws(rownames(an))
    cn <- colnames(an)
    pcol <- which(cn == "P")
    if(!length(pcol)) pcol <- grep("^P$|Pr", cn)[1]
    if(!length(pcol) || is.na(pcol)) return(c(overall=NA_real_, nonlinear=NA_real_))

    # In rms output the exposure is the main row and Nonlinear is the child row.
    io <- which(rn == food)[1]
    if(is.na(io)) io <- grep(paste0("^", gsub("([\\W])", "\\\\\\1", food), "$"), rn)[1]
    if(is.na(io)) io <- which(!grepl("Nonlinear|TOTAL|Total|ERROR", rn, ignore.case=TRUE))[1]
    inl <- grep("Nonlinear", rn, ignore.case=TRUE)[1]

    c(
        overall=if(length(io) && !is.na(io)) as.numeric(an[io, pcol]) else NA_real_,
        nonlinear=if(length(inl) && !is.na(inl)) as.numeric(an[inl, pcol]) else NA_real_
    )
}

make_hist_df <- function(x, xmax, bins=20){
    z <- x[is.finite(x) & x <= xmax]
    h <- graphics::hist(z, breaks=bins, plot=FALSE)
    data.frame(
        x=h$mids,
        count=h$counts,
        pct=100*h$counts/sum(h$counts),
        width=diff(h$breaks)
    )
}

plot_myth <- function(pr, d, food, ref_value, pv, x995){
    hd <- make_hist_df(d[[food]], x995, 20)
    ymax <- max(pr$upper, na.rm=TRUE)
    if(!is.finite(ymax) || ymax <= 1) ymax <- 1.5
    ymax <- ymax * 1.08
    sf <- if(max(hd$pct, na.rm=TRUE) > 0) ymax/max(hd$pct, na.rm=TRUE) else 1

    ggplot2::ggplot() +
        ggplot2::geom_col(
            data=hd, ggplot2::aes(x=x, y=pct*sf), width=hd$width,
            fill="#f9f7f7", colour="grey75", linewidth=.2
        ) +
        ggplot2::geom_hline(yintercept=1, linetype=2, colour="grey55") +
        ggplot2::geom_ribbon(
            data=pr, ggplot2::aes(x=.data[[food]], ymin=lower, ymax=upper),
            fill="#e23e57", alpha=.10
        ) +
        ggplot2::geom_line(
            data=pr, ggplot2::aes(x=.data[[food]], y=yhat),
            colour="#e23e57", linewidth=1.0
        ) +
        ggplot2::geom_point(ggplot2::aes(x=ref_value, y=1), colour="#e23e57", size=2) +
        ggplot2::scale_y_continuous(
            "Hazard ratio (95% CI)",
            sec.axis=ggplot2::sec_axis(~./sf, name="Participants (%)")
        ) +
        ggplot2::annotate(
            "text", x=min(pr[[food]], na.rm=TRUE), y=ymax,
            label=paste0("P overall ", fmt_p(pv["overall"]), "\nP non-linear ", fmt_p(pv["nonlinear"])),
            hjust=0, vjust=1, size=3.3
        ) +
        ggplot2::labs(title=gsub("_", " ", food), x="Intake (servings/day)") +
        ggplot2::coord_cartesian(xlim=c(min(pr[[food]], na.rm=TRUE), x995), ylim=c(0, ymax)) +
        ggplot2::theme_bw(base_size=11) +
        ggplot2::theme(panel.grid=ggplot2::element_blank(), panel.border=ggplot2::element_blank(), axis.line=ggplot2::element_line())
}

# Style copied from the user's old 01_plot.R concept:
# blue histogram + red curve + dashed CI + secondary frequency axis.
plot_01 <- function(pr, d, food, pv, x995){
    dplot <- d[d[[food]] <= x995, , drop=FALSE]
    h <- graphics::hist(dplot[[food]], breaks=20, plot=FALSE)
    ym <- max(pr$upper, na.rm=TRUE)
    if(!is.finite(ym) || ym <= 1) ym <- 1.5
    ym <- ym * 1.08
    sc <- if(max(h$counts) > 0) ym/max(h$counts) else 1

    ggplot2::ggplot(pr, ggplot2::aes(x=.data[[food]])) +
        ggplot2::geom_histogram(
            data=dplot,
            ggplot2::aes(x=.data[[food]], y=ggplot2::after_stat(count)*sc),
            bins=20, fill="#A6CEE3", colour="white", alpha=.70,
            inherit.aes=FALSE
        ) +
        ggplot2::geom_hline(yintercept=1, linetype=2, colour="grey60") +
        ggplot2::geom_line(ggplot2::aes(y=yhat), colour="#E41A1C", linewidth=1.1) +
        ggplot2::geom_line(ggplot2::aes(y=lower), colour="#E41A1C", linetype=2, linewidth=.8) +
        ggplot2::geom_line(ggplot2::aes(y=upper), colour="#E41A1C", linetype=2, linewidth=.8) +
        ggplot2::scale_y_continuous("HR (95% CI)", sec.axis=ggplot2::sec_axis(~./sc, name="Frequency")) +
        ggplot2::annotate(
            "text",
            x=min(pr[[food]], na.rm=TRUE) + .02*diff(range(pr[[food]], na.rm=TRUE)), y=ym,
            label=paste0("P overall ", fmt_p(pv["overall"]), "\nP non-linear ", fmt_p(pv["nonlinear"])),
            hjust=0, vjust=1, size=3.5
        ) +
        ggplot2::labs(title=gsub("_", " ", food), x="Intake (servings/day)") +
        ggplot2::coord_cartesian(xlim=c(min(pr[[food]], na.rm=TRUE), x995)) +
        ggplot2::theme_bw(base_size=11) +
        ggplot2::theme(
            panel.grid=ggplot2::element_blank(),
            plot.title=ggplot2::element_text(hjust=.5, face="bold", size=11, margin=ggplot2::margin(b=10)),
            axis.title.y.left=ggplot2::element_text(size=9, margin=ggplot2::margin(r=5)),
            axis.title.y.right=ggplot2::element_text(size=9, margin=ggplot2::margin(l=5)),
            axis.text=ggplot2::element_text(size=8), plot.margin=ggplot2::margin(15,15,15,15)
        )
}

run_one <- function(food){
    cat("RCS:", food, "\n")

    tryCatch({
        vars <- unique(c(t, y, food, covs))
        miss <- setdiff(vars, names(dat))
        if(length(miss)) stop("Missing variables: ", paste(miss, collapse=", "))

        d <- dat[, vars, drop=FALSE]
        d <- d[stats::complete.cases(d), , drop=FALSE]
        d[[food]] <- as.numeric(d[[food]])
        if(nrow(d) < 100) stop("<100 complete observations")

        k0 <- get_k3(d[[food]])
        if(is.null(k0)){
            return(list(
                myth=NULL, old=NULL,
                stat=tibble::tibble(
                    exposure=food, status="skip_knots", n=nrow(d), case=sum(d[[y]]==1),
                    unique_n=length(unique(d[[food]])), zero_pct=mean(d[[food]]==0),
                    knot_source=NA_character_, knot1=NA_real_, knot2=NA_real_, knot3=NA_real_,
                    ref_value=NA_real_, P_Overall=NA_real_, P_Nonlinear=NA_real_, x_995=NA_real_,
                    message="<3 unique knots"
                )
            ))
        }
        k <- k0$knots

        # datadist must be discoverable by name when rms::Predict runs.
        dd_name <- paste0("dd_rcs_", gsub("[^A-Za-z0-9_]", "_", food))
        assign(dd_name, rms::datadist(d), envir=.GlobalEnv)
        old_dd <- getOption("datadist")
        options(datadist=dd_name)
        on.exit({
            options(datadist=old_dd)
            if(exists(dd_name, envir=.GlobalEnv, inherits=FALSE)) rm(list=dd_name, envir=.GlobalEnv)
        }, add=TRUE)

        bt <- function(z) paste0("`", z, "`")
        rcs_term <- paste0("rcs(", bt(food), ", c(", paste(format(k, digits=16, scientific=FALSE), collapse=","), "))")
        f <- stats::as.formula(paste0(
            "survival::Surv(", bt(t), ",", bt(y), ") ~ ", rcs_term,
            if(length(covs)) paste0(" + ", paste(bt(covs), collapse=" + ")) else ""
        ))

        fit <- rms::cph(f, data=d, x=TRUE, y=TRUE, surv=TRUE)

        # IMPORTANT: anova.rms is an S3 method and is not called as rms::anova.rms.
        # stats::anova dispatches to the correct rms method.
        an <- stats::anova(fit)
        pv <- get_anova_p(an, food)

        # Reference is the datadist adjustment value (normally median; often 0 for zero-inflated foods).
        dd_obj <- get(dd_name, envir=.GlobalEnv)
        ref_value <- tryCatch(as.numeric(dd_obj$limits["Adjust to", food]), error=function(e) NA_real_)
        if(!is.finite(ref_value)) ref_value <- stats::median(d[[food]], na.rm=TRUE)

        x995 <- as.numeric(stats::quantile(d[[food]], .995, na.rm=TRUE, names=FALSE))
        xmin <- min(d[[food]], na.rm=TRUE)
        if(!is.finite(x995) || x995 <= xmin) x995 <- max(d[[food]], na.rm=TRUE)
        xseq <- seq(xmin, x995, length.out=200)

        args <- c(
            list(object=fit),
            setNames(list(xseq), food),
            list(ref.zero=TRUE, fun=exp, conf.int=.95)
        )
        pr <- as.data.frame(do.call(rms::Predict, args))
        req <- c(food, "yhat", "lower", "upper")
        if(!all(req %in% names(pr))) stop("Predict output missing: ", paste(setdiff(req, names(pr)), collapse=", "))
        pr <- pr[is.finite(pr[[food]]) & is.finite(pr$yhat) & is.finite(pr$lower) & is.finite(pr$upper), , drop=FALSE]
        if(!nrow(pr)) stop("No finite RCS predictions")

        p_myth <- plot_myth(pr, d, food, ref_value, pv, x995)
        p_old <- plot_01(pr, d, food, pv, x995)
        bounded_ggsave(file.path(DIR_MYTH, paste0(food, "_", Y, ".pdf")), p_myth, width=5.2, height=4.2)
        bounded_ggsave(file.path(DIR_OLD, paste0(food, "_", Y, ".pdf")), p_old, width=5.2, height=4.2)

        list(
            myth=p_myth, old=p_old,
            stat=tibble::tibble(
                exposure=food, status="ok", n=nrow(d), case=sum(d[[y]]==1),
                unique_n=length(unique(d[[food]])), zero_pct=mean(d[[food]]==0),
                knot_source=k0$source, knot1=k[1], knot2=k[2], knot3=k[3],
                ref_value=ref_value, P_Overall=unname(pv["overall"]), P_Nonlinear=unname(pv["nonlinear"]),
                x_995=x995, message=""
            )
        )
    }, error=function(e){
        list(
            myth=NULL, old=NULL,
            stat=tibble::tibble(
                exposure=food, status="error", n=NA_integer_, case=NA_integer_,
                unique_n=NA_integer_, zero_pct=NA_real_, knot_source=NA_character_,
                knot1=NA_real_, knot2=NA_real_, knot3=NA_real_, ref_value=NA_real_,
                P_Overall=NA_real_, P_Nonlinear=NA_real_, x_995=NA_real_,
                message=conditionMessage(e)
            )
        )
    })
}

res <- setNames(vector("list", length(group_vars)), group_vars)
for(i in seq_along(group_vars)) res[[i]] <- run_one(group_vars[i])

rcs_summary <- dplyr::bind_rows(lapply(res, `[[`, "stat"))
write.csv(rcs_summary, file.path(DIR_RCS, "RCS_summary.csv"), row.names=FALSE)
write.csv(rcs_summary, file.path(DIR_RCS, paste0("RCS_", Y, ".csv")), row.names=FALSE)

ok <- rcs_summary$exposure[rcs_summary$status == "ok"]
myth_ok <- Filter(Negate(is.null), lapply(res[ok], `[[`, "myth"))
old_ok <- Filter(Negate(is.null), lapply(res[ok], `[[`, "old"))

save_grid <- function(pl, file, title){
    if(!length(pl)) return(invisible(NULL))
    p <- patchwork::wrap_plots(pl, ncol=4) +
        patchwork::plot_annotation(title=title, theme=ggplot2::theme(plot.title=ggplot2::element_text(size=16, face="bold", hjust=.5)))
    bounded_ggsave(file, p, width=18, height=max(4.5, ceiling(length(pl)/4)*4.0), limitsize=FALSE)
}

save_grid(myth_ok, file.path(DIR_RCS, paste0("fig_rcs_all_myth_", Y, ".pdf")),
          paste0("Restricted cubic spline associations with ", Y, " — MYTH/rcssci-like style"))
save_grid(old_ok, file.path(DIR_RCS, paste0("fig_rcs_all_01plot_", Y, ".pdf")),
          paste0("Restricted cubic spline associations with ", Y, " — 01_plot.R style"))

nl <- rcs_summary$exposure[rcs_summary$status=="ok" & !is.na(rcs_summary$P_Nonlinear) & rcs_summary$P_Nonlinear < .05]
save_grid(Filter(Negate(is.null), lapply(res[nl], `[[`, "myth")),
          file.path(DIR_RCS, paste0("fig_rcs_nonlinear_myth_", Y, ".pdf")),
          "Food groups with significant nonlinearity — MYTH/rcssci-like style")
save_grid(Filter(Negate(is.null), lapply(res[nl], `[[`, "old")),
          file.path(DIR_RCS, paste0("fig_rcs_nonlinear_01plot_", Y, ".pdf")),
          "Food groups with significant nonlinearity — 01_plot.R style")

cat("RCS complete\n")
cat(" - scope:", RCS_SCOPE, "\n")
cat(" - requested groups:", length(group_vars), "\n")
cat(" - successful:", sum(rcs_summary$status=="ok"), "\n")
cat(" - skipped:", sum(grepl("^skip", rcs_summary$status)), "\n")
cat(" - errors:", sum(rcs_summary$status=="error"), "\n")
cat(" - consumer-knot fallback:", sum(rcs_summary$knot_source=="consumer_P10_P50_P90", na.rm=TRUE), "\n")
cat(" - nonlinear P<0.05:", sum(rcs_summary$status=="ok" & rcs_summary$P_Nonlinear<.05, na.rm=TRUE), "\n")
cat(" - MYTH-style plots:", DIR_MYTH, "\n")
cat(" - 01_plot-style plots:", DIR_OLD, "\n")
cat(" - summary:", file.path(DIR_RCS, "RCS_summary.csv"), "\n")

# If anything still fails, print the first few actual messages immediately.
err <- rcs_summary[rcs_summary$status=="error", c("exposure","message"), drop=FALSE]
if(nrow(err)){
    cat("\nFirst RCS errors:\n")
    print(utils::head(err, 10), row.names=FALSE)
}
}))

# ===========================================================================
# foodwas / plot
# ===========================================================================
dose_register_step("foodwas", "plot", expression({
# ======================================================================
# foodwas.R
# v24: explicit dplyr/tidyr namespaces to avoid select()/filter() masking after rms/rcssci
# Visualization and result tables for DOSE FoodWAS
#
# Inputs from foodwas / groups:
#   ukb/02foodwas/food_group_map.csv
#   ukb/02foodwas/dat.group.rds
#   ukb/02foodwas/item/Food_cvd_stroke_i.csv
#   ukb/02foodwas/item/Food_cvd_stroke_i_M1.csv
#   ukb/02foodwas/group/Group_cvd_stroke_i.csv
#   ukb/02foodwas/group/Group_cvd_stroke_i_M1.csv
#   ukb/02foodwas/group/Group_stat.csv
#
# Main outputs:
#   group/Group_cvd_stroke_i_proper.csv
#   group/Group_proper_reference.csv
#   group/Group_asso_cvd_stroke_i.csv
#   item/Food_asso_cvd_stroke_i.csv
#   group/Group_correlation_pairs_all.csv
#   figures/food_group_forest.pdf
#   figures/food_group_circular_heatmap.pdf
#   figures/food_group_multiring.pdf
#   figures/food_item_manhattan.pdf
#   figures/food_group_mean.pdf
#   figures/food_group_corr.pdf
#   figures/food_group_corr_cluster_heatmap.pdf
#
# Notes:
# - M2 is the primary plotting model.
# - The two post-ML combined food groups are excluded if present.
# - RCS nonlinear-P ring is optional: it is added only when an RCS summary
#   file containing exposure/Exposure and P_Nonlinear is available.
# ======================================================================

pacman::p_load(
	dplyr, tidyr, ggplot2, ggrepel, stringr, scales,
	ggnewscale, cowplot, corrplot, pheatmap, tibble
)

DIR_PROJECT <- dose_output_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_OUT <- file.path(DIR_PROJECT, "ukb/02foodwas")
DIR_ITEM <- file.path(DIR_OUT, "item")
DIR_GROUP <- file.path(DIR_OUT, "group")
DIR_FIG <- file.path(DIR_GROUP, "figures")
dir.create(DIR_FIG, recursive=TRUE, showWarnings=FALSE)

cfg <- readRDS(file.path(DIR_DATA, "analysis_config.rds"))

# Prespecified primary development outcome; do not switch based on post hoc results.
Y <- "cvd_stroke_i"

out_label <- unname(cfg$out_map[Y])
excl <- "Citrus_ApplesPears_Driedfruit"

f_group_m2 <- file.path(DIR_GROUP, paste0("Group_", Y, ".csv"))
f_group_m1 <- file.path(DIR_GROUP, paste0("Group_", Y, "_M1.csv"))
f_group_stat <- file.path(DIR_GROUP, "Group_stat.csv")
f_item_m2 <- file.path(DIR_ITEM, paste0("Food_", Y, ".csv"))
f_item_m1 <- file.path(DIR_ITEM, paste0("Food_", Y, "_M1.csv"))
f_map <- file.path(DIR_OUT, "food_group_map.csv")
f_dat_group <- file.path(DIR_OUT, "dat.group.rds")

stopifnot(
	file.exists(f_group_m2), file.exists(f_group_m1), file.exists(f_group_stat),
	file.exists(f_item_m2), file.exists(f_item_m1), file.exists(f_map),
	file.exists(f_dat_group)
)

field <- read.csv(f_map, check.names=FALSE, stringsAsFactors=FALSE)
dat_group <- readRDS(f_dat_group)
group_m2 <- read.csv(f_group_m2, check.names=FALSE, stringsAsFactors=FALSE)
group_m1 <- read.csv(f_group_m1, check.names=FALSE, stringsAsFactors=FALSE)
group_stat <- read.csv(f_group_stat, check.names=FALSE, stringsAsFactors=FALSE)
item_m2 <- read.csv(f_item_m2, check.names=FALSE, stringsAsFactors=FALSE)
item_m1 <- read.csv(f_item_m1, check.names=FALSE, stringsAsFactors=FALSE)

norm_key <- function(x) stringr::str_to_lower(stringr::str_squish(as.character(x)))
fmt_p <- function(p){
	ifelse(is.na(p), NA_character_,
		ifelse(p == 1, "-", ifelse(p < .001, "<0.001", sprintf("%.3f", p))))
}

clean_category <- function(x){
	dplyr::case_when(
		x %in% c("Meats", "Fish_fish_dishes") ~ "Meats & Fish",
		x %in% c("Non_alcoholic_beverages", "Alcoholic_beverages") ~ "Water & beverage",
		x == "Fried_and_fast_foods" ~ "Fried & Fast food",
		x == "Legumes_Nuts" ~ "Legumes & Nuts",
		x == "Snack_Pastry" ~ "Snack & Pastry",
		x == "Cereal_Bread" ~ "Cereal & Bread",
		x == "Fats_Oil" ~ "Fats & Oil",
		x == "Dairy_dairy_free_products" ~ "Dairy",
		TRUE ~ stringr::str_replace_all(as.character(x), "_", " ")
	)
}

# Food-group category map: one category per food group.
if(!"Food.Item" %in% names(field)) field$Food.Item <- field$name

cat_map <- field %>%
	dplyr::filter(!is.na(Food.Group), Food.Group != "") %>%
	dplyr::mutate(Category=clean_category(Category)) %>%
	dplyr::group_by(Food.Group) %>%
	dplyr::summarise(
		Category={z <- Category[!is.na(Category)]; if(length(z)) z[1] else "Unknown"},
		.groups="drop"
	) %>%
	dplyr::mutate(Category=tidyr::replace_na(Category, "Unknown")) %>%
	dplyr::rename(exposure=Food.Group)

# Original food groups only.
group_vars <- intersect(unique(field$Food.Group), names(dat_group))
group_vars <- setdiff(group_vars, excl)

# ======================================================================
# 1. Choose the lowest-risk intake category as plotting reference
# ======================================================================
# Original plotting logic:
#   1) inspect the model in which C1 is the reference;
#   2) find the category with the lowest estimated HR;
#   3) use the all-reference block with that category as Ref.
# This makes the forest/circular plots easier to interpret.

pick_best_ref <- function(d){
	d <- d %>% dplyr::mutate(Ref=as.character(Ref), class=as.character(class))
	base_ref <- if("C1" %in% d$Ref) "C1" else d$Ref[1]
	x <- d %>% dplyr::filter(Ref == base_ref)
	if(!nrow(x)) return(base_ref)
	jj <- which.min(replace(x$HR, !is.finite(x$HR), Inf))
	if(!length(jj) || !is.finite(x$HR[jj])) return(base_ref)
	as.character(x$class[jj])
}

ref_tbl <- group_m2 %>%
	dplyr::filter(exposure %in% group_vars) %>%
	dplyr::group_by(exposure) %>%
	dplyr::group_modify(~tibble(best_ref=pick_best_ref(.x))) %>%
	dplyr::ungroup()

proper <- group_m2 %>%
	dplyr::mutate(Ref=as.character(Ref), class=as.character(class)) %>%
	dplyr::filter(exposure %in% group_vars) %>%
	dplyr::inner_join(ref_tbl, by="exposure") %>%
	dplyr::filter(Ref == best_ref) %>%
	dplyr::left_join(
		group_stat %>%
			dplyr::mutate(quantile=as.character(quantile)) %>%
			dplyr::select(pheno, quantile, `Consumption level`, Median, Range, n, case),
		by=c("exposure"="pheno", "class"="quantile")
	) %>%
	dplyr::left_join(cat_map, by="exposure") %>%
	dplyr::mutate(
		class_n=suppressWarnings(as.integer(sub("^C", "", class))),
		reference=class == best_ref,
		sig=!reference & P < .05,
		HR_CI=ifelse(reference, "Reference", sprintf("%.2f (%.2f-%.2f)", HR, LCI, HCI)),
		P_label=ifelse(reference, "-", fmt_p(P))
	) %>%
	dplyr::arrange(Category, exposure, class_n)

write.csv(ref_tbl, file.path(DIR_GROUP, "Group_proper_reference.csv"), row.names=FALSE)
write.csv(proper, file.path(DIR_GROUP, paste0("Group_", Y, "_proper.csv")), row.names=FALSE)

# ======================================================================
# 2. Food-group faceted forest plot
# ======================================================================

forest_dat <- proper %>%
	dplyr::mutate(
		exposure_lab=str_replace_all(exposure, "_", " "),
		class_lab=factor(class, levels=paste0("C", 1:5))
	)

p_forest <- ggplot(forest_dat, aes(x=HR, y=class_lab, xmin=LCI, xmax=HCI, colour=sig)) +
	geom_vline(xintercept=1, linetype="dashed", linewidth=.3, colour="grey55") +
	geom_errorbar(orientation="y", width=.20, linewidth=.4, na.rm=TRUE) +
	geom_point(size=1.5, na.rm=TRUE) +
	scale_colour_manual(values=c(`FALSE`="black", `TRUE`="#BE665E")) +
	facet_wrap(~ exposure_lab, ncol=4, scales="free_y") +
	labs(
		x="Hazard ratio (95% CI)", y="Intake category",
		colour="P < 0.05",
		title=paste0("Associations of food groups with incident ", out_label, " (M2)")
	) +
	theme_classic(base_size=8) +
	theme(
		strip.background=element_rect(fill="#c2e7ff", colour="black", linewidth=.1),
		strip.text=element_text(size=8), axis.text.y=element_text(size=7),
		legend.position="bottom",
		panel.border=element_rect(fill=NA, colour="black", linewidth=.15),
		panel.spacing=grid::unit(.8, "lines")
	)

bounded_ggsave(
	file.path(DIR_FIG, "food_group_forest.pdf"),
	p_forest, width=8, height=max(10, ceiling(n_distinct(forest_dat$exposure)/4)*2.0)
)

# ======================================================================
# 3. Food-item Manhattan-style plot (M2)
# ======================================================================

item_map <- field %>%
	dplyr::transmute(
		exposure=name,
		Category=clean_category(Category),
		Food.Group=Food.Group,
		Food.Item=Food.Item
	) %>%
	dplyr::mutate(Category=tidyr::replace_na(Category, "Unknown"))

item_plot <- item_m2 %>%
	dplyr::filter(class == 2) %>%
	dplyr::mutate(
		logP=-log10(pmax(P, 1e-300)),
		HR_type=ifelse(HR < 1, "HR < 1", "HR > 1")
	) %>%
	dplyr::left_join(item_map, by="exposure")

fdr_p <- if(any(item_plot$FDR < .05, na.rm=TRUE)) {
	max(item_plot$P[item_plot$FDR < .05], na.rm=TRUE)
} else NULL

p_lines <- c(.05, .01, .001)
item_grp <- item_plot %>%
	dplyr::group_by(Category, Food.Group) %>%
	dplyr::summarise(mP=min(P, na.rm=TRUE), sig=sum(P < .05, na.rm=TRUE), .groups="drop") %>%
	dplyr::arrange(Category, mP, desc(sig)) %>%
	dplyr::mutate(x=row_number())
item_ax <- item_grp %>%
	dplyr::group_by(Category) %>%
	dplyr::summarise(cx=mean(x), bnd=max(x)+.5, .groups="drop")

set.seed(12345)
item_df <- item_plot %>%
	dplyr::left_join(item_grp %>% dplyr::select(Food.Group, x), by="Food.Group") %>%
	dplyr::group_by(Food.Group) %>%
	dplyr::mutate(xj=x + runif(n(), -.30, .30)) %>%
	dplyr::ungroup()

p_item <- ggplot(item_df, aes(x=xj, y=logP)) +
	geom_vline(xintercept=head(item_ax$bnd, -1), linetype="dashed", colour="grey85", linewidth=.5) +
	geom_hline(yintercept=-log10(p_lines), linetype="dotted", colour="grey70", linewidth=.4) +
	annotate(
		"text", x=.5, y=-log10(p_lines)+.05,
		label=paste0("P==", p_lines), parse=TRUE,
		size=2.2, colour="grey50", hjust=0
	) +
	{if(!is.null(fdr_p)) geom_hline(yintercept=-log10(fdr_p), linetype="dashed", colour="#BE665E", linewidth=.6)} +
	{if(!is.null(fdr_p)) annotate("text", x=.5, y=-log10(fdr_p)+.10, label="FDR < 0.05", size=2.8, colour="#BE665E", hjust=0, fontface="bold")} +
	geom_point(aes(colour=Category, shape=HR_type, size=total), alpha=.72) +
	ggrepel::geom_text_repel(
		data=item_df %>% dplyr::filter(P < .05), aes(label=Food.Item),
		size=2.6, max.overlaps=50, box.padding=.30,
		segment.color="grey55", show.legend=FALSE, colour="black"
	) +
	scale_shape_manual(name="Effect direction", values=c("HR < 1"=16, "HR > 1"=17)) +
	scale_size_continuous(name="Sample size", range=c(2, 6)) +
	scale_x_continuous(
		breaks=item_grp$x,
		labels=str_replace_all(item_grp$Food.Group, "_", " "),
		expand=expansion(mult=c(.01, .01))
	) +
	labs(
		x="Food groups", y=expression(-log[10](italic(P)~value)),
		title=paste0("Associations between individual foods and incident ", out_label, " (M2)")
	) +
	theme_bw(base_size=11) +
	theme(
		panel.grid=element_blank(),
		axis.text.x=element_text(angle=60, hjust=1, size=7),
		legend.position="right"
	)

bounded_ggsave(file.path(DIR_FIG, "food_item_manhattan.pdf"), p_item, width=18, height=8)

# ======================================================================
# 4. Mean intake and food-group correlation plots
# ======================================================================

dat_grp <- dat_group[, group_vars, drop=FALSE]

dat_mean <- dat_grp %>%
	dplyr::summarise(across(everything(), ~mean(.x, na.rm=TRUE))) %>%
	tidyr::pivot_longer(everything(), names_to="Food.Group", values_to="mean_intake") %>%
	dplyr::left_join(cat_map %>% dplyr::rename(Food.Group=exposure), by="Food.Group") %>%
	dplyr::mutate(Food.Group=factor(Food.Group, levels=group_vars))

p_mean <- ggplot(dat_mean, aes(Food.Group, mean_intake, fill=Category)) +
	geom_col() +
	labs(x=NULL, y="Average consumption (servings/day)", title="Average consumption levels of food groups") +
	theme_bw(base_size=11) +
	theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1, size=7))
bounded_ggsave(file.path(DIR_FIG, "food_group_mean.pdf"), p_mean, width=20, height=5)

cor_mat <- cor(dat_grp, method="spearman", use="pairwise.complete.obs")

pdf(file.path(DIR_FIG, "food_group_corr.pdf"), width=15, height=15)
corrplot::corrplot(
	cor_mat, method="square", type="lower", diag=TRUE, order="hclust",
	tl.col="black", tl.srt=45, tl.cex=.8, number.cex=.65,
	cl.cex=.8, cl.ratio=.1, addCoef.col="black",
	col=colorRampPalette(c("#B2182B", "white", "#2166AC"))(200), col.lim=c(-1, 1)
)
bounded_dev_off()

# Keep clustering distance consistent with s1_2: 1 - |Spearman rho|.
cor_clust <- cor_mat
cor_clust[is.na(cor_clust)] <- 0
diag(cor_clust) <- 1
hc <- hclust(as.dist(1 - abs(cor_clust)), method="average")
pheatmap::pheatmap(
	cor_mat,
	cluster_rows=hc, cluster_cols=hc,
	color=colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
	breaks=seq(-1, 1, length.out=201), border_color=NA,
	fontsize_row=8, fontsize_col=8, angle_col=45,
	main="Spearman correlation and hierarchical clustering of food groups",
	filename=file.path(DIR_FIG, "food_group_corr_cluster_heatmap.pdf"),
	width=13, height=12
)

idx <- which(upper.tri(cor_mat), arr.ind=TRUE)
cor_long <- tibble(
	group1=rownames(cor_mat)[idx[,1]],
	group2=colnames(cor_mat)[idx[,2]],
	rho=cor_mat[idx]
) %>%
	dplyr::mutate(abs_rho=abs(rho)) %>%
	dplyr::arrange(desc(abs_rho))
write.csv(cor_long, file.path(DIR_GROUP, "Group_correlation_pairs_all.csv"), row.names=FALSE)

# ======================================================================
# 5. Circular heatmap: Category + C1-C5 log(HR)
# ======================================================================

circ <- proper %>%
	dplyr::mutate(
		logHR=log(HR),
		star=case_when(
			reference | is.na(P) ~ "",
			P < .001 ~ "***",
			P < .01 ~ "**",
			P < .05 ~ "*",
			TRUE ~ ""
		)
	) %>%
	dplyr::group_by(exposure, Category) %>%
	dplyr::mutate(
		score=min(ifelse(reference | is.na(P), 1, P), na.rm=TRUE),
		eff=ifelse(any(!reference), max(abs(logHR[!reference]), na.rm=TRUE), 0)
	) %>%
	dplyr::ungroup()

sig_groups <- group_m2 %>%
	dplyr::mutate(Ref=as.character(Ref), class=as.character(class)) %>%
	dplyr::filter(exposure %in% group_vars, Ref != class, P < .05) %>%
	dplyr::distinct(exposure) %>%
	dplyr::pull(exposure)

circ <- circ %>%
	dplyr::mutate(fw_sig=exposure %in% sig_groups) %>%
	dplyr::arrange(Category, desc(fw_sig), score, desc(eff), exposure) %>%
	dplyr::mutate(
		exposure=factor(exposure, levels=unique(exposure)),
		id=as.integer(exposure),
		lab=str_replace_all(as.character(exposure), "_", "\n")
	)

c_dat <- circ %>% dplyr::distinct(exposure, Category, id)
c_sep <- c_dat %>% dplyr::group_by(Category) %>% dplyr::summarise(x=max(id)+.5, .groups="drop") %>% head(-1)
l_dat <- circ %>% dplyr::distinct(exposure, lab, id) %>%
	dplyr::mutate(
		ang=90 - 360*id/n(),
		hj=ifelse(ang < -90, 1, 0),
		ang=ifelse(ang < -90, ang+180, ang)
	)
mx <- max(abs(circ$logHR[is.finite(circ$logHR)]), na.rm=TRUE)

p_circ <- ggplot(circ, aes(x=id)) +
	geom_tile(data=c_dat, aes(y=0, fill=Category), width=.98, height=.72, colour="white", linewidth=.3) +
	scale_fill_manual(
		values=setNames(grDevices::hcl.colors(n_distinct(c_dat$Category), "Set 3"), unique(c_dat$Category)),
		name="Category"
	) +
	ggnewscale::new_scale_fill() +
	geom_tile(aes(y=class_n, fill=logHR), width=.98, height=.98, colour="white", linewidth=.35) +
	geom_text(data=circ %>% dplyr::filter(star != ""), aes(y=class_n, label=star), size=2.2, fontface=2) +
	geom_point(data=circ %>% dplyr::filter(reference), aes(y=class_n), shape=21, size=1.5, stroke=.35, fill="white") +
	geom_vline(data=c_sep, aes(xintercept=x), linewidth=.25, colour="grey70") +
	geom_text(data=l_dat, aes(y=6.18, label=lab, angle=ang, hjust=hj), size=2.2, lineheight=.9) +
	coord_polar(clip="off") +
	scale_fill_gradient2(
		low="#3C78D8", mid="white", high="#D6604D",
		limits=c(-mx, mx), oob=scales::squish,
		name="log(HR)\nrelative to\nreference"
	) +
	scale_x_continuous(breaks=NULL) +
	scale_y_continuous(limits=c(-.65, 6.45), breaks=NULL) +
	labs(
		title=paste0("Associations of ", length(group_vars), " food groups with incident ", out_label),
		subtitle="Inner ring = food category; rings 1–5 = intake categories; hollow point = lowest-risk reference; * P < 0.05, ** P < 0.01, *** P < 0.001",
		x=NULL, y=NULL
	) +
	theme_void(base_size=10) +
	theme(
		legend.position="right",
		plot.title=element_text(hjust=.5, face="bold"),
		plot.subtitle=element_text(hjust=.5),
		plot.margin=margin(20, 40, 20, 40)
	)

bounded_ggsave(file.path(DIR_FIG, "food_group_circular_heatmap.pdf"), p_circ, width=11.5, height=11.5)

# ======================================================================
# 6. Compact multi-ring figure
#    Category + C1-C5 HR + optional RCS nonlinear P + mean intake.
# ======================================================================

# Optional RCS summary from either the old or future s1_5 output.
rcs_candidates <- c(
	file.path(DIR_GROUP, "RCS", "RCS_summary.csv"),
	file.path(DIR_GROUP, "RCS", paste0("RCS_", Y, ".csv")),
	file.path(DIR_GROUP, paste0("RCS_", Y, ".csv"))
)
rcs_file <- rcs_candidates[file.exists(rcs_candidates)][1]
if(length(rcs_file) && !is.na(rcs_file)){
	rcs_dat <- read.csv(rcs_file, check.names=FALSE, stringsAsFactors=FALSE)
	if("Exposure" %in% names(rcs_dat)) rcs_dat <- rcs_dat %>% dplyr::rename(exposure=Exposure)
	if(!"exposure" %in% names(rcs_dat)) rcs_dat <- tibble(exposure=character(), P_Nonlinear=numeric())
	if(!"P_Nonlinear" %in% names(rcs_dat)) rcs_dat$P_Nonlinear <- NA_real_
	rcs_dat <- rcs_dat %>% dplyr::select(exposure, P_Nonlinear) %>% dplyr::distinct(exposure, .keep_all=TRUE)
} else {
	rcs_dat <- tibble(exposure=group_vars, P_Nonlinear=NA_real_)
}

mean_map <- dat_mean %>% dplyr::transmute(exposure=as.character(Food.Group), mean_intake)

plot_dat <- proper %>%
	dplyr::left_join(mean_map, by="exposure") %>%
	dplyr::left_join(rcs_dat, by="exposure") %>%
	dplyr::mutate(
		logHR=log(HR),
		ref=reference,
		star=case_when(ref | is.na(P) ~ "", P<.001 ~ "***", P<.01 ~ "**", P<.05 ~ "*", TRUE ~ ""),
		nl_logp=ifelse(is.na(P_Nonlinear), NA_real_, -log10(pmax(P_Nonlinear, 1e-300))),
		nl_star=case_when(is.na(P_Nonlinear) ~ "", P_Nonlinear<.001 ~ "***", P_Nonlinear<.01 ~ "**", P_Nonlinear<.05 ~ "*", TRUE ~ "")
	) %>%
	dplyr::group_by(exposure, Category) %>%
	dplyr::mutate(
		score=min(ifelse(ref | is.na(P), 1, P), na.rm=TRUE),
		eff=ifelse(any(!ref), max(abs(logHR[!ref]), na.rm=TRUE), 0)
	) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(fw_sig=exposure %in% sig_groups) %>%
	dplyr::arrange(Category, desc(fw_sig), score, desc(eff), exposure) %>%
	dplyr::mutate(
		exposure=factor(exposure, levels=unique(exposure)),
		id=as.integer(exposure)
	)

n_id <- max(plot_dat$id)
gap_l <- -.35; gap_r <- n_id + .95; x_s <- gap_r-gap_l
y_c_min <- 1.02; y_c_max <- 1.18
q_t <- .38; q_p <- .41; y_q_b <- 1.35; y_q_t <- y_q_b + 4*q_p + q_t
has_rcs <- any(is.finite(plot_dat$nl_logp))
y_n_min <- y_q_t + .18; y_n_max <- y_n_min + .20
y_b_b <- if(has_rcs) y_n_max + .18 else y_q_t + .30

col_hr <- c("#5F8796", "#FAF7F0", "#BE665E")
col_nl <- c("#F2F5F5", "#6F8790")
col_mean <- c("#E3EEE5", "#6E9C7B")
pal_cat <- c(
	"Cereal & Bread"="#C5A56A", "Dairy"="#8FAF9A", "Eggs"="#D8BF73",
	"Fats & Oil"="#B4A37B", "Fried & Fast food"="#CF9273", "Fruits"="#D1A09A",
	"Legumes & Nuts"="#A1B36F", "Meats & Fish"="#91AFC1", "Snack & Pastry"="#C3899D",
	"Vegetables"="#7EA986", "Water & beverage"="#89A9BF", "Unknown"="#B8B8B8"
)

group_info <- plot_dat %>%
	dplyr::distinct(exposure, Category, id, fw_sig, mean_intake, nl_logp, nl_star) %>%
	dplyr::mutate(
		exp_key=norm_key(exposure),
		lab=str_replace_all(as.character(exposure), "_", " "),
		mean_h=scales::rescale(mean_intake, to=c(.20, .82)),
		bar_t=y_b_b + mean_h,
		ang0=(90 - 360*(id-gap_l)/x_s + 180)%%360 - 180,
		hj=ifelse(ang0 < -90 | ang0 > 90, 1, 0),
		ang=ifelse(ang0 < -90 | ang0 > 90, ang0+180, ang0),
		lab_y=max(y_b_b+mean_h, na.rm=TRUE)+.22,
		lab_col=ifelse(fw_sig, col_hr[3], "#2F3742")
	)

strip_dat <- plot_dat %>%
	dplyr::mutate(
		ymin=y_q_b + (class_n-1)*q_p,
		ymax=ymin+q_t,
		ymid=(ymin+ymax)/2,
		star_col=ifelse(logHR>0, col_hr[3], col_hr[1])
	)

c_sep2 <- group_info %>% dplyr::group_by(Category) %>% dplyr::summarise(x=max(id)+.5, .groups="drop") %>% head(-1)
mx_q <- max(abs(strip_dat$logHR[is.finite(strip_dat$logHR)]), na.rm=TRUE)
mx_nl <- if(has_rcs) max(group_info$nl_logp[is.finite(group_info$nl_logp)], na.rm=TRUE) else 0
mean_min <- min(group_info$mean_intake, na.rm=TRUE)
mean_max <- max(group_info$mean_intake, na.rm=TRUE)
cats <- sort(unique(group_info$Category))
cat_cols <- pal_cat[cats]; cat_cols[is.na(cat_cols)] <- "#B8B8B8"

b_y <- c(y_q_b-.08, y_q_t+.055, y_b_b-.055, y_b_b)
if(has_rcs) b_y <- c(y_q_b-.08, y_q_t+.055, y_n_min-.055, y_n_max+.055, y_b_b-.055, y_b_b)
b_ls <- data.frame(y=b_y)

p_main <- ggplot() +
	geom_segment(data=b_ls, aes(x=.5, xend=max(group_info$id)+.5, y=y, yend=y), colour="#E7E7E7", linewidth=.25) +
	geom_segment(data=c_sep2, aes(x=x, xend=x, y=y_c_min, yend=y_b_b+.80), colour="#ECECEC", linewidth=.22, linetype="dashed") +
	geom_rect(data=group_info, aes(xmin=id-.47, xmax=id+.47, ymin=y_c_min, ymax=y_c_max, fill=Category), colour="white", linewidth=.09) +
	scale_fill_manual(values=cat_cols, guide="none") +
	ggnewscale::new_scale_fill() +
	geom_rect(data=strip_dat, aes(xmin=id-.47, xmax=id+.47, ymin=ymin, ymax=ymax, fill=logHR), colour="white", linewidth=.07) +
	geom_point(data=strip_dat %>% dplyr::filter(ref), aes(x=id, y=ymid), shape=21, size=.48, stroke=.18, fill="white", colour="#1F2937") +
	geom_text(data=strip_dat %>% dplyr::filter(!ref, P<.05), aes(x=id, y=ymid, label=star, colour=star_col), size=2.55, fontface="bold") +
	scale_fill_gradient2(low=col_hr[1], mid=col_hr[2], high=col_hr[3], midpoint=0, limits=c(-mx_q,mx_q), oob=scales::squish, guide="none")

if(has_rcs){
	p_main <- p_main +
		ggnewscale::new_scale_fill() +
		geom_rect(data=group_info, aes(xmin=id-.44, xmax=id+.44, ymin=y_n_min, ymax=y_n_max, fill=nl_logp), colour="white", linewidth=.07, alpha=.90) +
		geom_text(data=group_info %>% dplyr::filter(nl_star != ""), aes(x=id, y=(y_n_min+y_n_max)/2, label=nl_star), size=2.45, fontface="bold", colour="#1F2937") +
		scale_fill_gradient(low=col_nl[1], high=col_nl[2], limits=c(0,max(2,mx_nl)), oob=scales::squish, guide="none")
}

p_main <- p_main +
	ggnewscale::new_scale_fill() +
	geom_rect(data=group_info, aes(xmin=id-.47, xmax=id+.47, ymin=y_b_b, ymax=bar_t, fill=mean_intake), colour="white", linewidth=.045, alpha=.93) +
	scale_fill_gradient(low=col_mean[1], high=col_mean[2], guide="none") +
	geom_text(data=group_info, aes(x=id, y=lab_y, label=lab, angle=ang, hjust=hj, colour=lab_col), size=3.0, lineheight=.78, fontface="bold") +
	scale_colour_identity() +
	coord_polar(clip="off") +
	scale_x_continuous(limits=c(gap_l,gap_r), breaks=NULL) +
	scale_y_continuous(limits=c(-.72, max(group_info$lab_y)+.88), breaks=NULL) +
	theme_void() + theme(plot.margin=margin(2,2,2,2))

# Compact legend panel.
mk_g <- function(y, cols, x0=.25, x1=.75, n=100){
	data.frame(x=seq(x0,x1,length.out=n), y=y, w=(x1-x0)/n*1.05, fill=colorRampPalette(cols)(n))
}
nr <- ceiling(length(cats)/2)
cat_df <- tibble(cat=cats, col=unname(cat_cols[cats]), i=seq_along(cats)) %>%
	dplyr::mutate(col_id=(i-1)%/%nr+1, row_id=(i-1)%%nr+1, x=ifelse(col_id==1,.12,.56), y=.88-(row_id-1)*.04)
b_hr <- mk_g(.365, col_hr)
b_mn <- mk_g(.085, col_mean)

p_leg <- ggplot() + theme_void() + coord_cartesian(xlim=0:1, ylim=0:1, clip="off") +
	annotate("text", .5, .97, label="Food Category", fontface="bold", size=3.0, colour="#1F2937") +
	geom_rect(data=cat_df, aes(xmin=x, xmax=x+.012, ymin=y-.014, ymax=y+.014, fill=col), colour=NA) +
	geom_text(data=cat_df, aes(x=x+.03, y=y, label=cat), hjust=0, size=1.90, colour="#4B5563") +
	annotate("text", .5, .425, label="C1-C5 log(HR)", size=2.15, fontface="bold", colour="#374151") +
	geom_tile(data=b_hr, aes(x,y,fill=fill), width=b_hr$w, height=.021) +
	annotate("text", c(.25,.50,.75), .323, label=c(sprintf("%.2f",-mx_q),"0",sprintf("%.2f",mx_q)), size=1.70, colour="#6B7280") +
	annotate("text", .5, .145, label="Average consumption level\n(servings/day)", size=1.95, fontface="bold", colour="#374151") +
	geom_tile(data=b_mn, aes(x,y,fill=fill), width=b_mn$w, height=.021) +
	annotate("text", c(.25,.75), .043, label=c(sprintf("%.2f",mean_min),sprintf("%.2f",mean_max)), size=1.70, colour="#6B7280") +
	scale_fill_identity()

if(has_rcs){
	b_nl <- mk_g(.225, col_nl)
	p_leg <- p_leg +
		annotate("text", .5, .285, label="RCS (-log10 nonlinear P)", size=2.05, fontface="bold", colour="#374151") +
		geom_tile(data=b_nl, aes(x,y,fill=fill), width=b_nl$w, height=.021) +
		annotate("text", c(.25,.75), .183, label=c("0",sprintf("%.1f",max(2,mx_nl))), size=1.70, colour="#6B7280")
}

p_multi <- cowplot::ggdraw() +
	cowplot::draw_plot(p_main, 0, 0, 1, 1) +
	cowplot::draw_plot(p_leg, .402, .402, .196, .196)
bounded_ggsave(file.path(DIR_FIG, "food_group_multiring.pdf"), p_multi, width=14, height=14)

# ======================================================================
# 7. Publication-ready association tables for M1 and M2
# ======================================================================

get_group_model <- function(f, model_name){
	read.csv(f, check.names=FALSE, stringsAsFactors=FALSE) %>%
		dplyr::mutate(Ref=as.character(Ref), class=as.character(class)) %>%
		dplyr::distinct(exposure, class, Ref, .keep_all=TRUE) %>%
		dplyr::mutate(
			HR_CI=ifelse(HR==1 & P==1, "Ref.", sprintf("%.2f (%.2f-%.2f)", HR, LCI, HCI)),
			P_val=fmt_p(P)
		) %>%
		dplyr::select(exposure, class, Ref, HR_CI, P_val) %>%
		tidyr::pivot_wider(names_from=Ref, values_from=c(HR_CI, P_val), names_sep="_Ref") %>%
		dplyr::rename_with(~paste0(model_name, "_", .x), -c(exposure, class))
}

tbl_group <- group_stat %>%
	dplyr::mutate(
		class=as.character(quantile),
		`Food group`=pheno,
		`Intake level (servings/day)`=`Consumption level`,
		Total=n, Event=case,
		Incidence=sprintf("%.2f%%", case/n*100)
	) %>%
	dplyr::select(`Food group`, class, `Intake level (servings/day)`, Total, Event, Incidence) %>%
	dplyr::left_join(get_group_model(f_group_m1, "Model1"), by=c("Food group"="exposure", "class"="class")) %>%
	dplyr::left_join(get_group_model(f_group_m2, "Model2"), by=c("Food group"="exposure", "class"="class"))
write.csv(tbl_group, file.path(DIR_GROUP, paste0("Group_asso_", Y, ".csv")), row.names=FALSE)

get_item_model <- function(f, model_name){
	read.csv(f, check.names=FALSE, stringsAsFactors=FALSE) %>%
		dplyr::filter(class==2) %>%
		dplyr::transmute(
		exposure,
		Total=total,
		Event=case,
		"HR_CI"=sprintf("%.2f (%.2f-%.2f)", HR, LCI, HCI),
		"P_value"=fmt_p(P),
		"FDR"=fmt_p(FDR)
		) %>%
		dplyr::rename_with(~paste0(model_name, "_", .x), c(HR_CI, P_value, FDR))
}

item_tbl <- item_m2 %>%
	dplyr::filter(class==2) %>%
	dplyr::select(exposure, Total=total, Event=case) %>%
	dplyr::left_join(get_item_model(f_item_m1, "Model1") %>% dplyr::select(-Total, -Event), by="exposure") %>%
	dplyr::left_join(get_item_model(f_item_m2, "Model2") %>% dplyr::select(-Total, -Event), by="exposure") %>%
	dplyr::left_join(item_map %>% dplyr::select(exposure, Food.Item, Food.Group, Category), by="exposure") %>%
	dplyr::select(Food.Item, Food.Group, Category, everything())
write.csv(item_tbl, file.path(DIR_ITEM, paste0("Food_asso_", Y, ".csv")), row.names=FALSE)

cat("food_group_plot complete\n")
cat(" - food groups:", length(group_vars), "\n")
cat(" - M2 groups with any P<0.05:", length(sig_groups), "\n")
cat(" - max |rho| across all groups:", round(max(cor_long$abs_rho, na.rm=TRUE), 3), "\n")
cat(" - proper-reference table:", file.path(DIR_GROUP, paste0("Group_", Y, "_proper.csv")), "\n")
cat(" - group association table:", file.path(DIR_GROUP, paste0("Group_asso_", Y, ".csv")), "\n")
cat(" - item association table:", file.path(DIR_ITEM, paste0("Food_asso_", Y, ".csv")), "\n")
cat(" - figures:", DIR_FIG, "\n")
if(has_rcs) cat(" - RCS nonlinear-P ring added from:", rcs_file, "\n") else cat(" - RCS nonlinear-P ring not added: no compatible RCS summary file found\n")
}))

# ===========================================================================
# foodwas / plot_myth
# ===========================================================================
dose_register_step("foodwas", "plot_myth", expression({
# ============================================================================
# foodwas.R
# MYTH-style circular FoodWAS result plot for the DOSE project
#
# Main visual encoding (after Fig. 2B of the MYTH Diet paper):
#   - one sector and one distinct color per food group
#   - radial bars: hazard ratios for intake categories
#   - grey radial lines/caps: 95% confidence intervals
#   - inner asterisk: Cox P for trend < 0.05 (M2)
#   - outer asterisk: each non-reference category-specific Cox P < 0.05
#   - outer text: each intake range is the outermost layer
#
# Outputs:
#   ukb/02foodwas/group/figures/food_group_MYTH_circular_ref_first_ordered.pdf
#   ukb/02foodwas/group/figures/food_group_MYTH_circular_ref_first_ordered.png
#   ukb/02foodwas/group/figures/food_group_MYTH_C1_sensitivity_ref_first_ordered.pdf
#   ukb/02foodwas/group/figures/food_group_MYTH_plot_data.csv
#   ukb/02foodwas/group/figures/food_group_MYTH_color_map.csv
#
# Statistical note:
#   Following the BHOD Fig. 2 display strategy, the main figure uses the
#   observed highest-risk intake category as reference for each food group.
#   A fixed-C1 sensitivity figure is also written. Data-driven reference
#   selection is exploratory and makes most non-reference HRs appear below 1.
#
# This script intentionally uses base R only, so it can run on a clean R install.
# ============================================================================

options(stringsAsFactors = FALSE)

DIR_PROJECT <- dose_output_dir()
DIR_FOODWAS <- file.path(DIR_PROJECT, "ukb/02foodwas")
DIR_GROUP <- file.path(DIR_FOODWAS, "group")
DIR_FIG <- file.path(DIR_GROUP, "figures")
dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# Use Arial consistently in both PDF and PNG outputs. Under WSL/Linux, Arial
# is loaded from the host Windows font directory through a temporary fontconfig
# file; on native Windows, Cairo resolves Arial directly.
PLOT_FONT <- "Arial"
if (.Platform$OS.type != "windows") {
  windows_font_dir <- "/mnt/c/Windows/Fonts"
  if (file.exists(file.path(windows_font_dir, "arial.ttf"))) {
    dose_fontconfig <- tempfile("dose_fontconfig_", fileext = ".conf")
    writeLines(
      c(
        '<?xml version="1.0"?>',
        '<!DOCTYPE fontconfig SYSTEM "urn:fontconfig:fonts.dtd">',
        '<fontconfig>',
        '  <include ignore_missing="yes">/etc/fonts/fonts.conf</include>',
        paste0("  <dir>", windows_font_dir, "</dir>"),
        '</fontconfig>'
      ),
      dose_fontconfig,
      useBytes = TRUE
    )
    Sys.setenv(FONTCONFIG_FILE = dose_fontconfig)
  } else {
    warning("Arial font files were not found; the graphics device may substitute another sans-serif font")
  }
}

Y <- "cvd_stroke_i"
EXCLUDE_GROUPS <- "Citrus_ApplesPears_Driedfruit"

f_model <- file.path(DIR_GROUP, paste0("Group_", Y, ".csv"))
f_stat <- file.path(DIR_GROUP, "Group_stat.csv")
f_map <- file.path(DIR_FOODWAS, "food_group_map.csv")
f_rcs <- file.path(DIR_GROUP, "RCS", "RCS_summary.csv")
f_trend <- file.path(DIR_GROUP, paste0("Group_trend_", Y, ".csv"))

required <- c(f_model, f_stat, f_map, f_trend)
if (any(!file.exists(required))) {
  stop("Missing required input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

model <- read.csv(f_model, check.names = FALSE)
stat <- read.csv(f_stat, check.names = FALSE)
food_map <- read.csv(f_map, check.names = FALSE)
trend <- read.csv(f_trend, check.names = FALSE)

needed_model <- c("Ref", "exposure", "class", "HR", "LCI", "HCI", "P")
needed_stat <- c("pheno", "quantile", "Consumption level")
if (!all(needed_model %in% names(model))) stop("Unexpected model columns in: ", f_model)
if (!all(needed_stat %in% names(stat))) stop("Unexpected category columns in: ", f_stat)
if (!all(c("exposure", "model", "P_trend") %in% names(trend))) {
  stop("Unexpected trend columns in: ", f_trend)
}

model$Ref <- as.character(model$Ref)
model$class <- as.character(model$class)
model$exposure <- as.character(model$exposure)
model <- model[!model$exposure %in% EXCLUDE_GROUPS, , drop = FALSE]
trend <- trend[trend$model == "M2", c("exposure", "P_trend"), drop = FALSE]
trend <- trend[!duplicated(trend$exposure), , drop = FALSE]

clean_category <- function(x) {
  x <- as.character(x)
  out <- x
  out[x %in% c("Meats", "Fish_fish_dishes")] <- "Meats & Fish"
  out[x %in% c("Non_alcoholic_beverages", "Alcoholic_beverages")] <- "Water & beverage"
  out[x == "Fried_and_fast_foods"] <- "Fried & Fast food"
  out[x == "Legumes_Nuts"] <- "Legumes & Nuts"
  out[x == "Snack_Pastry"] <- "Snack & Pastry"
  out[x == "Cereal_Bread"] <- "Cereal & Bread"
  out[x == "Fats_Oil"] <- "Fats & Oil"
  out[x == "Dairy_dairy_free_products"] <- "Dairy"
  out <- gsub("_", " ", out, fixed = TRUE)
  out[is.na(out) | out == ""] <- "Unknown"
  out
}

# One broad food category per food group.
map_ok <- !is.na(food_map$Food.Group) & food_map$Food.Group != ""
map0 <- food_map[map_ok, c("Food.Group", "Category"), drop = FALSE]
map0$Category <- clean_category(map0$Category)
map0 <- map0[!duplicated(map0$Food.Group), , drop = FALSE]
names(map0)[1] <- "exposure"

# RCS results are retained in the exported plotting data for reference, but
# asterisks are defined from Cox results below, not from RCS P values.
if (file.exists(f_rcs)) {
  rcs <- read.csv(f_rcs, check.names = FALSE)
  if ("Exposure" %in% names(rcs) && !"exposure" %in% names(rcs)) names(rcs)[names(rcs) == "Exposure"] <- "exposure"
  keep <- intersect(c("exposure", "P_Overall", "P_Nonlinear"), names(rcs))
  rcs <- rcs[, keep, drop = FALSE]
  if (!"P_Overall" %in% names(rcs)) rcs$P_Overall <- NA_real_
  if (!"P_Nonlinear" %in% names(rcs)) rcs$P_Nonlinear <- NA_real_
  rcs <- rcs[!duplicated(rcs$exposure), c("exposure", "P_Overall", "P_Nonlinear")]
} else {
  rcs <- data.frame(exposure = unique(model$exposure), P_Overall = NA_real_, P_Nonlinear = NA_real_)
}

pick_reference <- function(d, mode = c("highest_risk", "C1")) {
  mode <- match.arg(mode)
  if (mode == "C1") return(if ("C1" %in% d$Ref) "C1" else d$Ref[1])
  base_ref <- if ("C1" %in% d$Ref) "C1" else d$Ref[1]
  z <- d[d$Ref == base_ref & is.finite(d$HR), , drop = FALSE]
  if (!nrow(z)) return(base_ref)
  as.character(z$class[which.max(z$HR)])
}

prepare_plot_data <- function(reference_mode = c("highest_risk", "C1")) {
  reference_mode <- match.arg(reference_mode)
  spl <- split(model, model$exposure)
  best <- vapply(spl, pick_reference, character(1), mode = reference_mode)
  ans <- do.call(rbind, lapply(names(spl), function(nm) {
    d <- spl[[nm]]
    z <- d[d$Ref == best[[nm]], , drop = FALSE]
    z <- z[!duplicated(z$class), , drop = FALSE]
    z$plot_reference <- best[[nm]]
    z$reference_rule <- if (reference_mode == "C1") {
      "prespecified C1"
    } else {
      "observed highest-HR category"
    }
    z
  }))
  rownames(ans) <- NULL
  ans <- merge(ans, stat[, c("pheno", "quantile", "Consumption level")],
               by.x = c("exposure", "class"), by.y = c("pheno", "quantile"), all.x = TRUE)
  ans <- merge(ans, map0, by = "exposure", all.x = TRUE)
  ans <- merge(ans, rcs, by = "exposure", all.x = TRUE)
  ans <- merge(ans, trend, by = "exposure", all.x = TRUE)
  ans$Category[is.na(ans$Category)] <- "Unknown"
  ans$class_n <- suppressWarnings(as.integer(sub("^C", "", ans$class)))
  ans$is_reference <- ans$class == ans$plot_reference
  ans <- ans[is.finite(ans$class_n), , drop = FALSE]
  if (any(abs(ans$HR[ans$is_reference] - 1) > 1e-8, na.rm = TRUE)) {
    stop("At least one reference category does not have HR = 1")
  }

  # Follow the MYTH figure's visual progression around the circle so adjacent
  # sectors move smoothly through the palette instead of jumping between hues.
  # Cox significance remains category-specific. For example, if both C2 and
  # C3 differ from C1 at P < 0.05, both intake labels receive an outer star.
  ans$cox_sig <- !ans$is_reference & is.finite(ans$P) & ans$P < 0.05
  group_info <- unique(ans[, c("exposure", "Category", "P_trend")])
  group_info$order_key <- match(group_info$exposure, myth_food_order)
  group_info$order_key[is.na(group_info$order_key)] <-
    length(myth_food_order) + seq_len(sum(is.na(group_info$order_key)))
  group_info <- group_info[order(group_info$order_key), ]
  group_info$order_key <- NULL
  group_info$id <- seq_len(nrow(group_info))
  ans <- merge(ans, group_info[, c("exposure", "id")],
               by = "exposure", all.x = TRUE)
  # Within every food group, place the reference category first; retain the
  # natural C1-C5 order for all remaining categories.
  ans <- ans[order(ans$id, !ans$is_reference, ans$class_n), ]
  rownames(ans) <- NULL
  ans
}

lighten <- function(col, amount = 0) {
  rgb0 <- col2rgb(col) / 255
  rgb1 <- rgb0 + (1 - rgb0) * amount
  rgb(rgb1[1, ], rgb1[2, ], rgb1[3, ])
}

darken <- function(col, amount = 0) {
  rgb0 <- col2rgb(col) / 255
  rgb1 <- rgb0 * (1 - amount)
  rgb(rgb1[1, ], rgb1[2, ], rgb1[3, ])
}

# MYTH-style circular order and palette: teal fruit/vegetables, blue grains,
# purple legumes/egg/dairy/meat, beige/brown fats, orange/yellow beverages,
# and pink/red processed foods. Every DOSE food group retains a unique shade.
myth_food_order <- c(
  "Berries", "Citrus", "Apples_Pears", "Other_fruit", "Dried_fruit",
  "Green_leafy_vegetables", "Allium_vegetables", "Orange_vegetables",
  "Other_vegetables", "Tomatoes", "Potatoes",
  "Refined_grains", "Whole_grains", "Legumes", "Nuts_seeds", "Eggs",
  "High_fat_dairy", "Low_fat_dairy", "Processed_meats", "Red_meats",
  "Poultry", "Fish", "Organ_meats", "Butter", "Margarine",
  "Olive_oil", "Other_oil", "Beer_cider", "Coffee", "Fruit_juice",
  "Spirits", "Sweetened_beverages", "Tea", "Water", "Wine", "Pizza",
  "Snacks_and_pastries", "Fried_and_fast_foods"
)

myth_group_palette <- c(
  Berries = "#006B66", Citrus = "#238A7C", Apples_Pears = "#3FA087",
  Other_fruit = "#59B596", Dried_fruit = "#78C9A5",
  Green_leafy_vegetables = "#43A889", Allium_vegetables = "#72C6A0",
  Orange_vegetables = "#99D7AE", Other_vegetables = "#9CCFD1",
  Tomatoes = "#70BFA8", Potatoes = "#6DB5CF",
  Refined_grains = "#3E96BA", Whole_grains = "#0B6A9D",
  Legumes = "#084E79", Nuts_seeds = "#2D0E5C", Eggs = "#5D1987",
  High_fat_dairy = "#8E36C5", Low_fat_dairy = "#B26DDE",
  Processed_meats = "#D1ADEB", Red_meats = "#B47ED8",
  Poultry = "#EFD2B2", Fish = "#D7A978", Organ_meats = "#B87848",
  Butter = "#956039", Margarine = "#6A3317",
  Olive_oil = "#B55308", Other_oil = "#D96A00",
  Beer_cider = "#F18A00", Coffee = "#F8A700", Fruit_juice = "#FFC13A",
  Spirits = "#EFC061", Sweetened_beverages = "#F3D58B",
  Tea = "#F198B4", Water = "#F36F96", Wine = "#E44A7D",
  Pizza = "#C72D5C", Snacks_and_pastries = "#981042",
  Fried_and_fast_foods = "#650018"
)

get_group_colors <- function(exposures) {
  cols <- unname(myth_group_palette[exposures])
  missing <- is.na(cols)
  if (any(missing)) {
    n <- sum(missing)
    hues <- (seq_len(n) * 137.508 + 8) %% 360
    cols[missing] <- grDevices::hcl(h = hues, c = 72,
                                    l = rep(c(52, 68), length.out = n))
  }
  names(cols) <- exposures
  cols
}

polar_xy <- function(theta, radius) cbind(radius * cos(theta), radius * sin(theta))

annular_polygon <- function(theta0, theta1, r0, r1, col, border = NA) {
  tt <- seq(theta0, theta1, length.out = 24)
  xy <- rbind(polar_xy(tt, r0), polar_xy(rev(tt), r1))
  polygon(xy[, 1], xy[, 2], col = col, border = border)
}

radial_segment <- function(theta, r0, r1, col = "#AEB5BC", lwd = 0.7) {
  xy <- polar_xy(theta, c(r0, r1))
  segments(xy[1, 1], xy[1, 2], xy[2, 1], xy[2, 2], col = col, lwd = lwd)
}

draw_arc <- function(radius, theta0, theta1, col = "#D7DDE2", lwd = 0.5,
                     lty = 1) {
  tt <- seq(theta0, theta1, length.out = 720)
  xy <- polar_xy(tt, radius)
  lines(xy[, 1], xy[, 2], col = col, lwd = lwd, lty = lty)
}

draw_myth <- function(dat, output_pdf, output_png = NULL,
                      reference_mode = "highest_risk") {
  groups <- unique(dat[, c("id", "exposure", "Category", "P_trend")])
  groups <- groups[order(groups$id), ]
  n_group <- nrow(groups)
  if (!n_group) stop("No food groups to plot")

  group_cols <- get_group_colors(groups$exposure)

  # HR-to-radius transform: HR decreases from the inner edge to the outer edge.
  # The common HR=1 circle remains the bar baseline. With the highest-risk
  # category as reference, HRs below 1 therefore extend outwards from it.
  vals <- c(dat$LCI, dat$HCI)
  vals <- vals[is.finite(vals) & vals > 0]
  lower <- min(0.90, floor(min(vals) * 20) / 20)
  upper <- max(1.20, ceiling(max(vals) * 20) / 20)
  lower <- max(0.50, lower)
  upper <- min(2.00, upper)
  r_inner <- 0.43
  r_base <- 0.58
  r_outer <- 1.02
  hr_to_r <- function(x) {
    y <- ifelse(x <= 1,
                r_base + (1 - x) / (1 - lower) * (r_outer - r_base),
                r_inner + (upper - x) / (upper - 1) * (r_base - r_inner))
    pmin(r_outer, pmax(r_inner, y))
  }

  # A narrow 18-degree opening centred at 12 o'clock forms a C-shaped HR axis.
  # HR tick labels sit inside this opening rather than over food labels.
  opening <- 18 * pi / 180
  available_angle <- 2 * pi - opening
  theta_start <- pi / 2 - opening / 2
  gap <- 0.009
  group_width <- available_angle / n_group
  theta_end <- theta_start - available_angle

  render <- function(device = c("pdf", "png"), file) {
    device <- match.arg(device)
    tmp_file <- tempfile(
      pattern = paste0(".", tools::file_path_sans_ext(basename(file)), "_"),
      tmpdir = dirname(file),
      fileext = paste0(".", tools::file_ext(file))
    )
    device_open <- FALSE
    on.exit({
      if (device_open) try(bounded_dev_off(), silent = TRUE)
      if (file.exists(tmp_file)) unlink(tmp_file)
    }, add = TRUE)
    if (device == "pdf") {
      grDevices::cairo_pdf(tmp_file, width = 14, height = 14, family = PLOT_FONT)
    } else {
      bounded_png(tmp_file, width = 3000, height = 3000, res = 240, type = "cairo")
    }
    device_open <- TRUE
    par(mar = c(0.4, 0.4, 0.4, 4.0), xpd = NA, family = PLOT_FONT)
    plot.new()
    plot.window(xlim = c(-1.62, 2.05), ylim = c(-1.55, 1.55), asp = 1)

    # Alternating pale backgrounds stop at the HR panel. The outer intake-text
    # ring is deliberately white, with no sector fill under the labels.
    for (i in seq_len(n_group)) {
      th0 <- theta_start - i * group_width + gap / 2
      th1 <- theta_start - (i - 1) * group_width - gap / 2
      bg <- if (i %% 2) "#F3F6F8" else "#FCFDFD"
      annular_polygon(th0, th1, r_inner - 0.015, r_outer + 0.005, bg)
    }

    # C-shaped HR guide arcs and compact labels in the top opening.
    ticks <- sort(unique(c(lower, 1.00, 1.10, upper)))
    ticks <- ticks[ticks >= lower & ticks <= upper]
    for (z in ticks) {
      rr <- hr_to_r(z)
      draw_arc(rr, theta_end, theta_start,
               col = if (abs(z - 1) < 1e-8) "#7F8790" else "#D7DDE2",
               lwd = if (abs(z - 1) < 1e-8) 1.0 else 0.45)
      text(0, rr, labels = format(z, nsmall = 1), cex = 0.72,
           col = "#56606A", pos = 3, offset = 0.12)
    }

    # Inner food-group annulus and M2 Cox trend-test stars.
    for (i in seq_len(n_group)) {
      g <- groups[i, ]
      th0 <- theta_start - i * group_width + gap / 2
      th1 <- theta_start - (i - 1) * group_width - gap / 2
      base_col <- unname(group_cols[g$exposure])
      annular_polygon(th0, th1, 0.315, 0.410, base_col, border = "white")
      if (is.finite(g$P_trend) && g$P_trend < 0.05) {
        thm <- (th0 + th1) / 2
        xy <- polar_xy(thm, 0.365)
        text(xy[1], xy[2], "*", cex = 1.00, font = 2, col = "#17212B")
      }
    }
    draw_arc(0.302, theta_end, theta_start, col = "#D0D5DA", lty = 2)

    # Radial bars, point estimates, and 95% CIs.
    for (i in seq_len(n_group)) {
      g <- groups[i, ]
      d <- dat[dat$id == g$id, ]
      d <- d[order(!d$is_reference, d$class_n), ]
      k <- nrow(d)
      th_left <- theta_start - i * group_width + gap * 1.15
      th_right <- theta_start - (i - 1) * group_width - gap * 1.15
      # Start at the sector's leading edge and move clockwise. This preserves
      # the data order visually: reference first, then the remaining C groups.
      edges <- seq(th_right, th_left, length.out = k + 1)
      base_col <- unname(group_cols[g$exposure])
      nonref <- which(!d$is_reference)
      shades <- rep(0, k)
      if (length(nonref)) {
        shades[nonref] <- seq(0.34, 0.00,
                              length.out = max(length(nonref), 2))[seq_along(nonref)]
      }
      for (j in seq_len(k)) {
        th0 <- edges[j] - 0.0012
        th1 <- edges[j + 1] + 0.0012
        thm <- (th0 + th1) / 2
        # A reference category is exactly HR=1 and therefore has neither an
        # effect bar nor a confidence interval. Its intake label is retained.
        if (!isTRUE(d$is_reference[j])) {
          rr <- hr_to_r(d$HR[j])
          annular_polygon(th0, th1, r_base, rr, lighten(base_col, shades[j]))
        }
        if (!isTRUE(d$is_reference[j]) &&
            is.finite(d$LCI[j]) && is.finite(d$HCI[j])) {
          rl <- hr_to_r(d$LCI[j]); rh <- hr_to_r(d$HCI[j])
          radial_segment(thm, rl, rh, col = "#A5ADB5", lwd = 0.72)
          cap <- abs(th1 - th0) * 0.23
          radial_segment(thm - cap, rl, rl, col = NA)
          xy1 <- polar_xy(c(thm - cap, thm + cap), rl)
          xy2 <- polar_xy(c(thm - cap, thm + cap), rh)
          segments(xy1[1, 1], xy1[1, 2], xy1[2, 1], xy1[2, 2], col = "#A5ADB5", lwd = 0.55)
          segments(xy2[1, 1], xy2[1, 2], xy2[2, 1], xy2[2, 2], col = "#A5ADB5", lwd = 0.55)
        }
      }
    }

    # Each intake range follows its corresponding bar. Cox stars form the
    # inner layer and the range labels form the outermost circular layer.
    for (i in seq_len(n_group)) {
      g <- groups[i, ]
      d <- dat[dat$id == g$id, ]
      d <- d[order(!d$is_reference, d$class_n), ]
      k <- nrow(d)
      th_left <- theta_start - i * group_width + gap * 1.15
      th_right <- theta_start - (i - 1) * group_width - gap * 1.15
      # Use the same clockwise leading-edge order as the HR bars.
      edges <- seq(th_right, th_left, length.out = k + 1)
      intake <- as.character(d$`Consumption level`)
      intake[is.na(intake) | intake == ""] <- d$class[is.na(intake) | intake == ""]
      for (j in seq_len(k)) {
        th <- (edges[j] + edges[j + 1]) / 2
        deg <- th * 180 / pi
        rot <- deg
        hjust <- 0
        # Keep text upright on the left half while extending it outwards.
        if (cos(th) < 0) { rot <- rot + 180; hjust <- 1 }
        label_radius <- 1.082
        xy <- polar_xy(th, label_radius)
        text(xy[1], xy[2], intake[j], srt = rot,
             adj = c(hjust, 0.5), cex = 0.46,
             col = darken(unname(group_cols[g$exposure]), 0.12), font = 1)
        # The star belongs to this exact Cox contrast. It sits immediately
        # inside its intake label, so the label remains the outermost layer.
        # Multiple significant categories within one food group receive
        # multiple stars at their own category angles.
        if (isTRUE(d$cox_sig[j])) {
          xy_star <- polar_xy(th, 1.048)
          text(xy_star[1], xy_star[2], "*", cex = 0.90, font = 2,
               col = darken(unname(group_cols[g$exposure]), 0.12))
        }
      }
    }

    # Category separators.
    cat_change <- which(groups$Category[-1] != groups$Category[-n_group])
    if (length(cat_change)) {
      for (ii in cat_change) {
        th <- theta_start - ii * group_width
        xy <- polar_xy(th, c(0.305, r_outer + 0.005))
        segments(xy[1, 1], xy[1, 2], xy[2, 1], xy[2, 2], col = "white", lwd = 2.0)
      }
    }

    # Legend: food-group names inherit their broad-category colors.
    legend_x <- 1.49
    legend_y <- 1.38
    # Stack swatches edge-to-edge so the 38 food groups form one continuous
    # color strip while retaining enough room for the larger labels.
    line_h <- 0.052
    text(legend_x, legend_y + 0.08, "Food group", adj = 0, font = 2, cex = 0.92)
    for (i in seq_len(n_group)) {
      yy <- legend_y - (i - 1) * line_h
      col <- unname(group_cols[groups$exposure[i]])
      rect(legend_x, yy - line_h / 2, legend_x + 0.058, yy + line_h / 2,
           col = col, border = NA)
      text(legend_x + 0.065, yy, gsub("_", " ", groups$exposure[i], fixed = TRUE),
           adj = 0, cex = 0.60, col = "#17212B")
    }

    title_text <- if (reference_mode == "highest_risk") {
      paste0("Associations between ", n_group, " food groups and incident ischemic stroke")
    } else {
      paste0("Sensitivity analysis: C1-reference associations for ", n_group, " food groups")
    }
    text(0, 1.50, title_text, cex = 1.18, font = 2, col = "#17212B")
    text(0, -1.47,
         "HR decreases from inner to outer; reference = 1 (no bar); grey lines: 95% CI; inner *: P trend < 0.05; outer *: Cox P < 0.05",
         cex = 0.72, col = "#4D5863")
    if (reference_mode == "highest_risk") {
      text(0, -1.53,
           "Reference: observed highest-HR intake category for each food group",
           cex = 0.72, col = "#68737E")
    } else {
      text(0, -1.53, "Reference: prespecified lowest intake category (C1)",
           cex = 0.72, col = "#68737E")
    }

    bounded_dev_off()
    device_open <- FALSE
    if (!file.copy(tmp_file, file, overwrite = TRUE)) {
      stop("Unable to replace output file; close any PDF/image viewer using: ", file)
    }
    unlink(tmp_file)
  }

  render("pdf", output_pdf)
  if (!is.null(output_png)) render("png", output_png)
}

dat_myth <- prepare_plot_data("highest_risk")
dat_c1 <- prepare_plot_data("C1")

# Record both selections: displayed stars compare against the highest-risk
# reference; ML screens all contrasts and then removes correlated candidates.
# These different definitions need not produce identical sets on a new cohort.
f_ml_candidates <- file.path(DIR_GROUP, "Need_Food_lst.txt")
if (file.exists(f_ml_candidates)) {
  ml_candidates <- readLines(f_ml_candidates, warn = FALSE)
  starred_foods <- unique(dat_myth$exposure[dat_myth$cox_sig])
  audit_foods <- union(ml_candidates, unique(dat_myth$exposure))
  write.csv(data.frame(food=audit_foods, in_ML= audit_foods %in% ml_candidates,
                       starred_in_figure=audit_foods %in% starred_foods),
            file.path(DIR_FIG, "food_group_MYTH_ML_selection_audit.csv"), row.names=FALSE)
  if (!setequal(ml_candidates, starred_foods)) {
    message(
      "Figure/ML selection audit (different contrast and clustering rules). ML without displayed stars: ",
      paste(setdiff(ml_candidates, starred_foods), collapse = ", "),
      "; starred but not in ML: ",
      paste(setdiff(starred_foods, ml_candidates), collapse = ", ")
    )
  }
}

out_myth_pdf <- file.path(DIR_FIG, "food_group_MYTH_circular_ref_first_ordered.pdf")
out_myth_png <- file.path(DIR_FIG, "food_group_MYTH_circular_ref_first_ordered.png")
out_c1_pdf <- file.path(DIR_FIG, "food_group_MYTH_C1_sensitivity_ref_first_ordered.pdf")
out_data <- file.path(DIR_FIG, "food_group_MYTH_plot_data.csv")
out_colors <- file.path(DIR_FIG, "food_group_MYTH_color_map.csv")

draw_myth(dat_myth, out_myth_pdf, out_myth_png, "highest_risk")
draw_myth(dat_c1, out_c1_pdf, NULL, "C1")
write.csv(dat_myth, out_data, row.names = FALSE)
color_map <- unique(dat_myth[, c("id", "exposure", "Category")])
color_map <- color_map[order(color_map$id), ]
color_map$color <- unname(get_group_colors(color_map$exposure))
write.csv(color_map, out_colors, row.names = FALSE)

cat("MYTH-style FoodWAS figures completed.\n")
cat("  Main PDF:       ", out_myth_pdf, "\n", sep = "")
cat("  Main PNG:       ", out_myth_png, "\n", sep = "")
cat("  C1 sensitivity: ", out_c1_pdf, "\n", sep = "")
cat("  Plot data:      ", out_data, "\n", sep = "")
cat("  Color map:      ", out_colors, "\n", sep = "")
cat("  Food groups:    ", length(unique(dat_myth$exposure)), "\n", sep = "")
cat("  Inner trend P<.05: ", length(unique(dat_myth$exposure[is.finite(dat_myth$P_trend) & dat_myth$P_trend < .05])), "\n", sep = "")
cat("  Outer Cox P<.05 contrasts: ", sum(dat_myth$cox_sig, na.rm = TRUE), "\n", sep = "")
if (exists("ml_candidates")) {
  cat("  ML candidates: ", length(ml_candidates), "; starred in main figure: ",
      length(starred_foods), "; shared: ", length(intersect(ml_candidates, starred_foods)), "\n", sep = "")
}
}))

if (sys.nframe() == 0L) dose_run_module("foodwas")
