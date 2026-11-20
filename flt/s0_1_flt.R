pacman::p_load(tidyverse, data.table, survival, readxl)

dir0 <- "D:"
indir <- paste0(dir0, "/data/phe/ukb/phe")
invisible(lapply(c("phe.f.R", "mr.f.R", "assoc.f.R", "plot.f.R"), \(f) source(file.path(dir0, "scripts", "f", f))))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 数据QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fn	<- fread(paste0(indir, "/common/flt.lst"), sep = "\t", select = 1:3, fill = TRUE, header = FALSE) %>% as.data.frame()
fn$V1[duplicated(fn$V1)]; fn$V2[duplicated(fn$V2)]
phe0 <- fread(paste0(indir, "/rap/flt.tab.gz"), sep = "\t", header = TRUE) %>% as.data.frame()
names(phe0)	<- names(phe0) %>% gsub("_i0$|_i0_a0$", "", .)
setdiff(fn$V1, names(phe0)); setdiff(names(phe0), fn$V1)
names(phe0)	<- sapply(names(phe0), replace_code, sep = "_", mapping_df = fn) %>% as.character()
phe0 <- phe0 %>% mutate(across(all_of(unique(c(grep("date", names(phe0), value = TRUE), "p20750"))), as.Date))
phe0 <- neg2NA(phe0, c(-121,-818)) # -121＝从未尝试；-818＝不愿回答，均转为NA
saveRDS(phe0, paste0(indir, "/Rdata/phe.flt0.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 140个原始偏好 food-liking-traits＋ 4个 喜欢咖啡／茶”和“偏好加糖的饮法 【PMID: 35585065】
# 176 :36个派生偏好 【PMID: 40410146】
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 144个：140个原始偏好 + 4个 喜欢咖啡／茶”和“偏好加糖的饮法 【PMID: 35585065】
file <- file.path(dir0, "files/flt.xlsx") 
s1 <- read_excel(file, sheet = "Sheet1"); s2 <- read_excel(file, sheet = "Sheet2"); s4 <- read_excel(file, sheet = "Sheet4")
flt144 <- s1$name
group_map <- s2 %>% # 43 derived ： 从144 中最终进入 hierarchical factor tree 的 questionnaire traits生成
	pivot_longer(cols = starts_with("level_"), names_to = "level", values_to = "group") %>% 
	filter(!is.na(group), group != "") %>% distinct(group, name)
fatty_vars <- group_map %>% filter(group %in% c("F-Fatty_dairy", "F-Cheese")) %>% pull(name) %>% unique() # 添加F-Fatty：F-Fatty_dairy与F-Cheese覆盖的8个条目
group_map <- bind_rows(group_map, tibble(group = "F-Fatty", name = fatty_vars))
flt176_names <- s4$name
phe0 <- readRDS(file.path(indir, "Rdata/phe.flt0.rds"))
phe <- phe0 %>% mutate( # 咖啡／茶偏好：max至少一项有效；差值需要两项均有效,差值＝加糖评分－不加糖评分，正值表示更喜欢加糖
    coffee_max  = pmax(p20643, p20644, na.rm = TRUE),
    coffee_diff = p20643 - p20644,
    tea_max     = pmax(p20734, p20735, na.rm = TRUE),
    tea_diff    = p20734 - p20735
)
# 识别有Food-liking记录的参与者
has_fl <- rowSums(!is.na(phe[, flt144, drop = FALSE])) > 0
	if ("p20750" %in% names(phe)) {completed <- !is.na(phe$p20750) & nzchar(trimws(as.character(phe$p20750))); has_fl <- has_fl | completed}
	cat("总人数：", nrow(phe), "；无Food-liking记录：", sum(!has_fl), "；有Food-liking记录，QC前：", sum(has_fl), "\n")
# 保留问卷完成时间（随访起点）；不将其作为偏好条目参与评分或缺失率QC。
stopifnot("p20750" %in% names(phe))
flt <- phe[has_fl, c("eid", "p20750", flt144), drop = FALSE]
nrow(flt)# 181943
# 生成44个衍生组：144＋44＝188项：每组有效条目取均值，全部缺失则保留NA
for (g in unique(group_map$group)) {
    vars <- group_map$name[group_map$group == g]
    x <- rowMeans(flt[, vars, drop = FALSE], na.rm = TRUE)
    x[is.nan(x)] <- NA_real_
    flt[[g]] <- x
}

# 提取Sheet4的176项，并按这176项进行个体QC  【PMID: 40410146】
flt176 <- flt[, c("eid", "p20750", flt176_names), drop = FALSE]
keep <- rowMeans(is.na(flt176[, flt176_names, drop = FALSE])) <= 0.30  # 个体QC：176项中缺失超过30%者排除
cat("176项QC前：", nrow(flt176), "；因缺失超过30%排除：", sum(!keep), "；最终纳入，QC后：", sum(keep), "\n")
flt176 <- flt176[keep, , drop = FALSE]
nrow(flt176) # 181835
# saveRDS(flt176, file.path(indir, "Rdata/flt176.rds"))  # eid＋p20750＋176项

flt <- flt[keep, , drop = FALSE]
nrow(flt) # 181835
saveRDS(flt, file.path(indir, "Rdata/flt.rds")) # eid＋p20750＋188项

