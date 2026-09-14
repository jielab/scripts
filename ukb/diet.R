pacman::p_load(readxl, tidyverse, data.table, survival)
# UKB fields: touchscreen, 24-hour recall, nutrients, and Oxford WebQ.

dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
indir <- paste0(dir0, "/data/ukb/phe")
invisible(lapply(c("0phe.f.R", "mr.f.R", "assoc.f.R", "plot.f.R"), \(f) source(file.path(dir0, "scripts", "f", f))))


# 🚩 data QC
fn	<- fread(paste0(indir, "/common/diet.lst"), sep = "\t", select = 1:2, header = FALSE, fill = TRUE) %>% as.data.frame()
fn$V1[duplicated(fn$V1)]; fn$V2[duplicated(fn$V2)]
diet0 <- fread(paste0(indir, "/rap/diet.tab.gz"), sep = "\t", header = TRUE) %>% as.data.frame()
diet0 <- diet0 %>% mutate(across(grep("date", names(diet0), value = TRUE), as.Date))
names(diet0) <- sapply(names(diet0), replace_code, sep = "_", mapping_df = fn) %>% as.character()
diet0 <- neg2NA(diet0, c(-1, -2, -3, -7, -10, -121, -818))
diet0 <- diet0 %>% mutate(across(where(is.numeric), \(x) case_when(x==200~2, x==300~3, x==400~4, x==500~5, x==600~6, x==444~0.25, x==555~0.5, TRUE~x)))
diet0 %>% dplyr::select(where(~ any(stringr::str_detect(as.character(.x), stringr::fixed("|")), na.rm = TRUE))) %>% names()
saveRDS(diet0, paste0(indir, "/Rdata/diet0.rds")) 


# 🚩 MEDI(PMID 40285703) touch-screen
# type_milk[1418], type_spread [1428], type_bread [1448], type_cereal [1468]
diet0 <- readRDS(paste0(indir, "/Rdata/diet0.rds")) %>% as.data.frame() %>% dplyr::select(matches("^(eid$|touch\\.)")) %>% rename_with(~ gsub("touch\\.|_i0$", "", .x))
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
saveRDS(diet, paste0(indir, "/Rdata/diet.medito.rds"))


# 🚩 Oxford WebQ(Fudan, PMID 40603580)
inst <- 0:4 # *_i0 to *_i4
field <- read_excel(paste0(dir0, "/files/foods.xlsx"), sheet = "Sheet2")
diet0 <- readRDS(paste0(indir, "/Rdata/diet0.rds")) %>% dplyr::select(!matches("^touch\\.")) %>% filter(!is.na(p20077)) 
	dura <- diet0 %>% dplyr::select(eid, starts_with("p20082")) # finished too fast
	q <- quantile(as.matrix(dura[, -1]), na.rm = TRUE, probs = c(0.05, 0.95))
	dura <- dura %>% mutate(across(-eid, ~ ifelse(.x >= q[1] & .x <= 60 * 24, .x, NA_real_)), enroll = rowSums(!is.na(across(-eid)))); table(dura$enroll)
	vari <- diet0 %>% dplyr::select(eid, matches("vari")) %>% mutate(vari = as.integer(rowSums(across(-eid, ~ .x == 3, .names = NULL), na.rm = TRUE) > 0))
	keep <- dura %>% inner_join(vari %>% dplyr::select(eid, vari), by = "eid") %>% filter(enroll >= 1, vari == 0)
diet <- diet0 %>% filter(eid %in% keep$eid) %>% dplyr::select(-starts_with("p20082")) %>% inner_join(keep %>% dplyr::select(eid, enroll, starts_with("p20082")), by = "eid")
for (k in inst) {
	ok <- !is.na(diet[[paste0("p20082_i", k)]])
	diet[[paste0("p105010_i", k)]][!ok] <- as.Date(NA); diet[[paste0("p26002_i", k)]][!ok]  <- NA_real_
}
diet <- diet %>% mutate(
	across(starts_with("p105010_i"), as.Date), 
	start_date = do.call(pmax, c(across(starts_with("p105010_i")), list(na.rm = TRUE))),
	across(starts_with("p26002_i"), ~ ifelse((sex == 0) & (. < 500*4.184 | . > 3500*4.184) | (sex == 1) & (. < 800*4.184 | . > 4200*4.184), NA_real_, .)),
    energy.kJ = rowMeans1(across(starts_with("p26002_i")), per = 0.8)) %>% filter(!is.na(energy.kJ))
	nrow(diet) # 203463
breads <- c("bread", "baguette", "bap", "roll") # 
bread_code <- c("1" = "white", "2" = "mixed", "3" = "wholemeal", "4" = "seeded", "5" = "other")
for (k in inst) {
	for (f in breads) {
		type_col   <- sprintf("type_%s_i%d", f, k); intake_col <- sprintf("%s_i%d", f, k)
		valid_type <- !is.na(diet[[type_col]]) & str_detect(diet[[type_col]], "[1-5]"); total <- str_count(diet[[type_col]], "[1-5]")
		for (cc in names(bread_code)) {
			suf <- bread_code[[cc]]; pat <- paste0("(^|\\|)", cc, "(\\||$)"); cnt <- str_count(diet[[type_col]], pat)
			diet[[sprintf("%s.%s_i%d", f, suf, k)]] <- ifelse(!is.na(diet[[intake_col]]) & valid_type & total > 0 & cnt > 0, diet[[intake_col]] * cnt / total, 0)
		}
		diet[[sprintf("%s.unanswered_i%d", f, k)]] <- ifelse(is.na(diet[[intake_col]]), NA_real_, ifelse(!valid_type, diet[[intake_col]], 0))
	}
	raw_col <- sprintf("type_oil_i%d", k); new_col <- sprintf("oil_i%d", k) # 
	if (raw_col %in% names(diet)) { diet[[new_col]] <- ifelse(str_detect(diet[[raw_col]], "^35[4-7]$"), 10, 0)
	} else  {diet[[new_col]] <- NA_real_ }
}
soups <- c("canned", "homemade") # 
soup_code <- c("1" = "fish", "2" = "meat", "3" = "pasta", "4" = "pulse", "5" = "vegetable", "6" = "other")
for (k in inst) {
	for (f in soups) {    
		type_col   <- sprintf("type_%s_i%d", f, k); intake_col <- sprintf("%s_i%d", f, k)
		valid_type <- !is.na(diet[[type_col]]) & str_detect(diet[[type_col]], "[1-6]"); total <- str_count(diet[[type_col]], "[1-6]")
		for (cc in names(soup_code)) {
			suf <- soup_code[[cc]]; pat <- paste0("(^|\\|)", cc, "(\\||$)"); cnt <- str_count(diet[[type_col]], pat)
			diet[[sprintf("%s.%s_i%d", f, suf, k)]] <- ifelse(!is.na(diet[[intake_col]]) & valid_type & total > 0 & cnt > 0, diet[[intake_col]] * cnt / total, 0)
		}
    diet[[sprintf("%s.unanswered_i%d", f, k)]] <- ifelse(is.na(diet[[intake_col]]), NA_real_, ifelse(!valid_type, diet[[intake_col]], 0))
	}
}
sauce_code <- list(sauce.Jamhoney =c(332), sauce.cream =c(333), sauce.pnutbutter = c(334), sauce.yeast = c(335), sauce.hummus = c(336), sauce.guacamole = c(337), sauce.chutney = c(338), sauce.ketchup = c(339), sauce.brown = c(340), sauce.mayonnaise = c(341, 342), sauce.saladdressing = c(343), sauce.olives = c(344), sauce.pesto = c(345), sauce.tomato = c(346), sauce.cheese = c(347), sauce.white = c(348), sauce.gravy = c(349), sauce.other = c(350))
for (k in inst) {
	raw_col <- sprintf("types_sauce_i%d", k)
	if (raw_col %in% names(diet)) {
		for (nm in names(sauce_code)) {
			codes <- sauce_code[[nm]]; pat <- paste0("(^|\\|)(", paste(codes, collapse="|"), ")(\\||$)")
			diet[[sprintf("%s_i%d", nm, k)]] <- ifelse(is.na(diet[[raw_col]]), NA_real_, ifelse(str_detect(diet[[raw_col]], pat), 1, 0))
		}
	} else {for (nm in names(sauce_code)) diet[[sprintf("%s_i%d", nm, k)]] <- NA_real_}
}
diet <- diet %>% dplyr::select(-matches(paste0("^(", paste(c(breads), collapse = "|"), ")_i")), -matches("^type"))
dat1 <- diet %>% dplyr::select(eid)
	for (i in field$name) dat1 <- cbind(dat1, dplyr::select(diet, starts_with(i)))
	names(dat1)[duplicated(names(dat1))] 
for (k in inst) { 
	valid_k <- diet[[paste0("p26002_i", k)]]	
	idx  <- seq(2 + k, ncol(dat1), by = 5) 
	dat1[, idx] <- lapply(dat1[, idx, drop = FALSE], function(x) ifelse(is.na(x) & !is.na(valid_k), 0, x))
}
wt_vars <- "p26150 p26131 p26133 p26096 p26102 p26154 p26110 oil p26062 p26063 p26111 p26112 p26064 p26153 p26152 p26151 p26067 p26138" %>% strsplit(" ") %>% unlist() # 所有重量的都处理
portion <- field %>% filter(name %in% wt_vars) %>% dplyr::select(name, Portion) %>% mutate(Portion = as.numeric(Portion)) %>% tibble::deframe()
dat.food <- dat1 %>% dplyr::select(eid) #
	for (i in field$name) { dat.food[[i]] <- rowMeans1(dat1 %>% dplyr::select(starts_with(i)), per = 0.8)} 
	dat.food <- dat.food %>% mutate(across(all_of(names(portion)), ~ .x / portion[cur_column()])) # divide portion size  
length(names(dat.food)) #  235
stopifnot(identical(dat.food$eid, diet$eid))
dat.food <- cbind(dat.food, subset(diet, select=c("energy.kJ", "start_date", "enroll"))) 

# diet weight 
field.w <- read_excel(paste0(dir0, "/files/foods.xlsx"), sheet = "Sheet3")
dat1 <- diet %>% dplyr::select(eid)
	for (i in field.w$name) dat1 <- cbind(dat1, dplyr::select(diet, starts_with(i)))
	names(dat1)[duplicated(names(dat1))] 
weight_names   <- field.w %>% dplyr::filter(Category != "nutrient") %>% dplyr::pull(name) %>% unique()
nutrient_names <- field.w %>% dplyr::filter(Category == "nutrient") %>% dplyr::pull(name) %>% unique()
for (k in inst) {
	valid_k <- !is.na(diet[[paste0("p26002_i", k)]])  
	idx_k <- seq(2 + k, ncol(dat1), by = 5)
	idx_wk <- intersect(idx_k, which(sub("_i\\d+$", "", names(dat1)) %in% weight_names))
	if (length(idx_wk) > 0) {dat1[, idx_wk] <- lapply(dat1[, idx_wk, drop = FALSE], function(x) ifelse(is.na(x) & valid_k, 0, x))}
	idx_nk <- intersect(idx_k, which(sub("_i\\d+$", "", names(dat1)) %in% nutrient_names))
	if (length(idx_nk) > 0) {dat1[, idx_nk] <- lapply(dat1[, idx_nk, drop = FALSE], function(x) ifelse(!valid_k, NA_real_, suppressWarnings(as.numeric(as.character(x)))))}
}	
dat.weight <- dat1 %>% dplyr::select(eid) 
	for (i in field.w$name) { dat.weight[[paste0("w.", i)]] <- rowMeans1(dat1 %>% dplyr::select(starts_with(i)), per = 0.8)}
	length(names(dat.weight))  #  157
	stopifnot(identical(dat.weight$eid, diet$eid))	
dat.food <- dat.food %>% left_join(dat.weight, by = "eid") %>% left_join(diet %>% transmute(eid, sex, weight = weight_i0), by = "eid")
saveRDS(dat.food, paste0(indir, "/Rdata/diet.fudan.rds")) 


# 🚩 MODERN
cols <- 'Berries_and_Citrus Green_leafy_vegetables Olive_oil Potatoes Eggs Poultry Sweetened_beverages' %>% strsplit(' ') %>% unlist() 
dat_field <- field %>% filter(MODERN > 0)
dat1 <- dat.food
	grp.link <- split(dat_field$name, dat_field$Food.Group)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE], per = 1)}
dat1 <- dat1 %>% mutate(
	Berries_and_Citrus 	= Berries + Citrus,
	Sweetened_beverages = Sweetened_beverages + Fruit_juice
) %>% dplyr::select(eid, all_of(cols)) 
dat.modern <- dat1 %>% dplyr::select(eid, any_of(cols)) %>% transmute(eid,
	Berries_and_Citrus.modern	  = ifelse(Berries_and_Citrus > 0 & Berries_and_Citrus <= 2, 1, 0),
	Green_leafy_vegetables.modern = ifelse(Green_leafy_vegetables > 0.25, 1, 0),
	Olive_oil.modern			  = ifelse(Olive_oil > 0, 1, 0),
	Potatoes.modern 			  = ifelse(Potatoes > 0 & Potatoes <= 0.75, 1, 0),
	Eggs.modern					  = ifelse(Eggs > 0 & Eggs <= 1, 1, 0),
	Poultry.modern 				  = ifelse(Poultry > 0 & Poultry <= 0.5, 1, 0),
	Sweetened_beverages.modern 	  = ifelse(Sweetened_beverages == 0, 1, 0)
) %>% mutate(diet.modern.sum = rowSums1(across(ends_with(".modern")))) %>% dplyr::select(eid, diet.modern.sum) # ends_with(".modern")
saveRDS(dat.modern, paste0(indir, "/Rdata/diet.modern.rds"))


# 🚩 MIND
func_ter <- function(x){
	breaks = quantile(x, probs = seq(0, 1, 1/3), na.rm = TRUE) 
	if(length(unique(breaks)) == 2){a <- factor(ifelse(x>0, 1, 0))
	}else if(length(unique(breaks)) < 4){a <- cut(x,breaks = c(min(x, na.rm = TRUE)-1, unique(breaks)), include.lowest = T) 
	}else{a <- cut(x,breaks = unique(breaks),include.lowest = T)}
	return(a)
}
func_for1 <- function(x){a <- ifelse(x == 1,0,1); return(a)}	
func_for2 <- function(x){a <- ifelse(x == 1,0,ifelse(x == 2,0.5,1)); return(a)}
func_aginst1 <- function(x){a <- ifelse(x == 1,1,0); return(a)}
func_aginst2 <- function(x){a <- ifelse(x == 1,1,ifelse(x == 2,0.5,0)); return(a)}
func_wine <- function(x){a <- ifelse(x == 1,1,ifelse(x == 0 | x >1,0,0.5)); return(a)}
cols <- 'Berries Green_leafy_vegetables Other_vegetables Legumes Nuts_seeds Whole_grains Poultry Seafood Wine Olive_oil Red_processed_meat High_fat_diary Butter_and_Margarine Pastries_sweets Fried_and_fast_foods' %>% strsplit(' ') %>% unlist() 
dat_field <- field %>% filter (MIND > 0)
dat1 <- dat.food
	grp.link <- split(dat_field$name, dat_field$Food.Group)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE], per = 1)} 
dat1 <- dat1 %>% mutate(
	Green_leafy_vegetables = p104080 + p104090*0.5 + p104160 + p104240 + p104300 + p104370,
	Other_vegetables 	   = p104090*0.5 + p104140 + p104180 + p104310 + Allium_vegetables + Tomatoes + Orange_vegetables + Other_vegetables - p104080, 
	Red_processed_meat 	   = rowSums1(across(c(Red_meats, Processed_meats, Organ_meats))),
	Butter_and_Margarine   = rowSums1(across(c(Animal_fats, Plant_fats))),
	Seafood 			   = Fish + Fish_based_mixed_dishes, 
	Pastries_sweets		   = Snacks_and_pastries + Sweetened_beverages
) %>% dplyr::select(eid, all_of(cols)) 
dat.mind <- dat1 %>% mutate(across(-c(eid), func_ter)) %>% mutate(across(-c(eid), as.integer)) %>% mutate(
    across(all_of(c("Berries", "Nuts_seeds", "Olive_oil")), func_for1),
    across(all_of(c("Green_leafy_vegetables", "Other_vegetables", "Legumes", "Whole_grains","Poultry", "Seafood")), func_for2),
    across(all_of(c("Red_processed_meat", "High_fat_diary", "Butter_and_Margarine", "Pastries_sweets", "Fried_and_fast_foods")),  func_aginst2),
    Wine = func_wine(dat1$Wine)
) %>% rename_with(~paste0(., ".mind"), .cols = -eid
) %>% mutate(diet.mind.sum = rowSums1(across(ends_with(".mind")))) %>% dplyr::select(eid, diet.mind.sum)  # ends_with(".mind")
saveRDS(dat.mind, paste0(indir, "/Rdata/diet.mind.rds"))


# 🚩 Medi Pattern score (PMID: 41245532)
cols <- 'Vegetables Fruit_fruitjuice Olive_oil Red_processed_meat White_meat Butter_and_Margarine Sweetened_beverages Wine Legumes Nuts_seeds Seafood Sweets_desserts Sofrito' %>% strsplit(' ') %>% unlist() 
dat_field <- field %>% filter (MedDiet > 0) 
dat1 <- dat.food
	grp.link <- split(dat_field$name, dat_field$Food.Group)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE], per = 1)}
dat1 <- dat1 %>% mutate(
	Vegetables 		  	 = rowSums1(across(c(Green_leafy_vegetables, Allium_vegetables, Orange_vegetables, Other_vegetables, p104340, p104350))),
	Fruit_fruitjuice	 = rowSums1(across(c(Apples_Pears, Berries, Citrus, Dried_fruit, Other_fruit, Fruit_juice, p102490))),
	Olive_oil 		  	 = ifelse(Olive_oil > 0, 1, 0),
	Red_processed_meat 	 = rowSums1(across(c(Red_meats, Processed_meats, Organ_meats))),
	White_meat   	  	 = Poultry + p103050,
	Butter_and_Margarine = rowSums1(across(c(Animal_fats, Plant_fats))),
	Seafood 			 = rowSums1(across(c(Fish, Fish_based_mixed_dishes, p103170, p103180))),
	Sweets_desserts		 = Snacks_and_pastries - p102490,
	Sofrito 			 = sauce.tomato
	) %>% dplyr::select(eid, all_of(cols)) 
score_up <- function(x, thr_day) {x <- as.numeric(x); pmin(pmax(x / thr_day, 0), 1)}  # 0 at 0, 1 at >=thr (linear)
score_down <- function(x, low_day, high_day) {x <- as.numeric(x); out <- ifelse(is.na(x), NA_real_, ifelse(x <= low_day, 1, ifelse(x >= high_day, 0, 1 - (x - low_day) / (high_day - low_day)))); pmin(pmax(out, 0), 1)} # 1 at <=low, 0 at >=high (linear)
dat.medi24 <- dat1 %>% mutate( # 1) continuous(0-1)
    Vegetables.medi       = score_up(Vegetables, 2),
    Fruit_fruitjuice.medi = score_up(Fruit_fruitjuice, 3),
    Sofrito.medi          = score_up(Sofrito, 2/7),
    Wine.medi             = score_up(Wine, 1),
    Olive_oil.medi        = as.numeric(Olive_oil),
	across(c(Legumes, Nuts_seeds, Seafood), ~ score_up(.x, 3/7), .names = "{.col}.medi"),
	across(c(Red_processed_meat, Butter_and_Margarine, Sweetened_beverages), ~ score_down(.x, low_day = 1, high_day = 2), .names = "{.col}.medi"),
    Sweets_desserts.medi  = score_down(Sweets_desserts, low_day = 2/7, high_day = 4/7),
    White_meat.medi 	  = ifelse(is.na(White_meat) | is.na(Red_processed_meat), NA_real_, ifelse(White_meat >= Red_processed_meat, 1, 0))
) 
dat.medi24 <- dat.medi24 %>% mutate( # 2) integer(0/1)
    Vegetables.bin       = as.numeric(Vegetables >= 2),
    Fruit_fruitjuice.bin = as.numeric(Fruit_fruitjuice >= 3),
    Sofrito.bin          = as.numeric(Sofrito >= 2/7),
    Wine.bin             = as.numeric(Wine >= 1),
    Olive_oil.bin        = as.numeric(Olive_oil > 0),
    across(c(Legumes, Nuts_seeds, Seafood), ~ as.numeric(.x >= 3/7), .names = "{.col}.bin"),
    across(c(Red_processed_meat, Butter_and_Margarine, Sweetened_beverages), ~ as.numeric(.x <= 1), .names = "{.col}.bin"),
    Sweets_desserts.bin  = as.numeric(Sweets_desserts <= 2/7),
    White_meat.bin       = as.numeric(White_meat.medi)   
)
dat.medi24 <- dat.medi24 %>% mutate(
	diet.medi24.sum = rowSums1(across(ends_with(".medi"))),
	# diet.medi24_bin.sum = rowSums(across(ends_with(".bin"))),
) %>% dplyr::select(eid, diet.medi24.sum) #  ends_with(".medi")
saveRDS(dat.medi24, paste0(indir, "/Rdata/diet.medi24.rds"))


# 🚩 DASH (PMID: 41245532)
score_q5_zero <- function(x, reverse=FALSE){x <- as.numeric(x); out <- rep(NA_integer_, length(x)); i0 <- !is.na(x)&x==0; i1 <- !is.na(x)&x>0; out[i0] <- 1L; out[i1] <- ntile(x[i1], 4)+1L; if(reverse) out <- 6L-out; out} # 对大量 0 的变量，用“0 单独一档 + 其余做分位”
cols <- 'Vegetables Fruit_fruitjuice Red_processed_meat Legumes_Nuts Whole_grains Low_fat_diary Sweetened_beverages sodium' %>% strsplit(' ') %>% unlist() 
cols.new <- 'veg fruit red_meat nut whole_grain dairy sweet sodium' %>% strsplit(' ') %>% unlist() 

dat_field <- field %>% filter (DASH > 0)
dat1 <- dat.food
	grp.link <- split(dat_field$name, dat_field$Food.Group)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE], per = 1)} 
dat1 <- dat1 %>% mutate(
	Vegetables 		   = rowSums1(across(c(Green_leafy_vegetables, Allium_vegetables, Tomatoes, Orange_vegetables, Other_vegetables))),
	Fruit_fruitjuice   = rowSums1(across(c(Apples_Pears, Berries, Citrus, Dried_fruit, Other_fruit, Fruit_juice, p102490))),
	Red_processed_meat = rowSums1(across(c(Red_meats, Processed_meats, Organ_meats))),
	Legumes_Nuts 	   = Legumes + Nuts_seeds, 
	sodium 			   = w.p26052
	) %>% dplyr::select(eid, all_of(cols)) 
dat.dash <- dat1 %>% mutate(
    across(c(Vegetables, Fruit_fruitjuice, Legumes_Nuts, Whole_grains, Low_fat_diary), ~ score_q5_zero(.x, reverse = FALSE), .names = "{.col}.dash"),
	across(c(Red_processed_meat, Sweetened_beverages, sodium), ~ score_q5_zero(.x, reverse = TRUE), .names = "{.col}.dash")
) %>% mutate(diet.dash.sum = rowSums1(across(ends_with(".dash")))) %>% dplyr::select(eid, diet.dash.sum, ends_with(".dash")) %>%
	rename_with(~ paste0(cols.new[match(sub("\\.dash$", "", .x), cols)], ".dash"), all_of(paste0(cols, ".dash")))
saveRDS(dat.dash, paste0(indir, "/Rdata/diet.dash.rds"))


# 🚩 healthful Plant-Based Index (hPDI) (PMID: 41245532)
cols <- 'Vegetables Fruits Legumes_vegetarian Nuts_seeds Whole_grains Tea_coffee Refined_grains Potatoes Sugary_drinks Fruit_juice Sweets_desserts Animal_fats Dairy Eggs Seafood Meats Miscellaneous_animal_based_foods' %>% strsplit(' ') %>% unlist() 
dat_field <- field %>% filter (hPDI > 0)
dat1 <- dat.food
	grp.link <- split(dat_field$name, dat_field$Food.Group)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE], per = 1)}
dat1 <- dat1 %>% mutate(
	Vegetables  	   = rowSums1(across(c(Green_leafy_vegetables, Allium_vegetables, Tomatoes, Orange_vegetables, Other_vegetables))),
	Fruits			   = rowSums1(across(c(Apples_Pears, Berries, Citrus, Dried_fruit, Other_fruit, p102490))),
	Legumes_vegetarian = Legumes + Meat_substitutes,
	Tea_coffee 		   = Tea + Coffee, 
	Refined_grains 	   = Refined_grains + rowSums1(across(c(p102010, p102020, p102050, p102070, p102470, p102500))),
	Potatoes		   = Potatoes + p104020,
	Sugary_drinks      = Sweetened_beverages - rowSums1(across(c(p100230, p100530))),
	Sweets_desserts	   = Snacks_and_pastries - rowSums1(across(c(p102490, p102010, p102020, p102050, p102070, p102470, p102500, p102220, p102140, p102150, p102120))),
	Dairy 			   = rowSums1(across(c(Yogurt, High_fat_diary, Low_fat_diary, p100230, p100530, p102220, p102140, p102150, p102120))),
	Seafood 		   = rowSums1(across(c(Fish, Fish_based_mixed_dishes, p103170, p103180))),
	Meats			   = rowSums1(across(c(Poultry, Red_meats, Processed_meats, Organ_meats, Other_meats, p103050)))
) %>% dplyr::select(eid, all_of(cols)) 
dat.hpdi <- dat1 %>% mutate(
    across(all_of(c("Vegetables", "Fruits", "Legumes_vegetarian", "Nuts_seeds", "Whole_grains", "Tea_coffee")), ~ score_q5_zero(.x, reverse = FALSE), .names = "{.col}.hpdi"),
    across(all_of(c("Refined_grains", "Potatoes", "Sugary_drinks", "Fruit_juice", "Sweets_desserts", "Animal_fats", "Dairy", "Eggs", "Seafood", "Meats", "Miscellaneous_animal_based_foods")), ~ score_q5_zero(.x, reverse = TRUE),  .names = "{.col}.hpdi")
) %>% mutate(diet.hpdi.sum = rowSums1(across(ends_with(".hpdi")))) %>% dplyr::select(eid, diet.hpdi.sum)  # ends_with(".hpdi")
saveRDS(dat.hpdi, paste0(indir, "/Rdata/diet.hpdi.rds"))


# 🚩 AHEI-2010 (PMID: 41245532)
ahei_score_up10 <- function(x, max_x) {x <- as.numeric(x); pmin(pmax(x / max_x, 0), 1) * 10}
ahei_score_down10 <- function(x, max_x) {x <- as.numeric(x); pmin(pmax(1 - x / max_x, 0), 1) * 10}
ahei_score_mid10 <- function(x, lo, hi, high0){x <- as.numeric(x); ifelse(is.na(x), NA_real_, ifelse(x >= high0, 0, ifelse(x >= lo & x <= hi, 10, ifelse(x < lo, pmin(pmax(x/lo, 0), 1)*10, pmin(pmax((high0-x)/(high0-hi), 0), 1)*10))))} # 中等最好：lo~hi 得10；>high0 得0；<lo 线性上升
cols <- 'Vegetables Fruits Red_processed_meat Legumes_Nuts Sugary_drinks_fruitjuices weight_Whole_grains Trans_fat Oily_fish PUFA sodium Alcoholic_beverages' %>% strsplit(' ') %>% unlist() 
dat_field <- field %>% filter (AHEI > 0)
dat1 <- dat.food
	grp.link <- split(dat_field$name, dat_field$Food.Group)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE], per = 1)}
dat1 <- dat1 %>% mutate(
	Vegetables  	   = rowSums1(across(c(Green_leafy_vegetables, Allium_vegetables, Tomatoes, Orange_vegetables, Other_vegetables))),
	Fruits			   = rowSums1(across(c(Apples_Pears, Berries, Citrus, Dried_fruit, Other_fruit, p102490))),
	Red_processed_meat = rowSums1(across(c(Red_meats, Processed_meats, Organ_meats))), 
	Legumes_Nuts 	   = Legumes + Nuts_seeds, 
	Sugary_drinks_fruitjuices = Sweetened_beverages + Fruit_juice, 
	weight_Whole_grains = rowSums1(across(all_of(paste0("w.", c("p26076", "p26077", "p26074", "p26114"))))) + 0.5 * .data[["w.p26105"]], # g/day
	sodium 			   = w.p26052,
	Alcoholic_beverages= rowSums1(across(c(Wine, Beer_cider, Spirits))), 
	energy_kcal 	   = energy.kJ * 0.239, 
	Trans_fat 		   = (w.p26155 * 9 / energy_kcal) * 100,
	Oily_fish		   = p103160,
    PUFA 			   = ((w.p26015 + w.p26016) * 9 / energy_kcal) * 100
	) %>% dplyr::select(eid, all_of(cols), sex) 
dat.ahei <- dat1 %>% mutate(
    Vegetables.ahei 			 = ahei_score_up10(Vegetables, 5),
    Fruits.ahei 				 = ahei_score_up10(Fruits, 4),
    Legumes_Nuts.ahei 			 = ahei_score_up10(Legumes_Nuts, 1),
    Red_processed_meat.ahei		 = ahei_score_down10(Red_processed_meat, 1.5),
    Non_alcoholic_beverages.ahei = ahei_score_down10(Sugary_drinks_fruitjuices, 1),
    weight_Whole_grains.ahei 	 = ifelse(sex == 0, ahei_score_up10(weight_Whole_grains, 75), ahei_score_up10(weight_Whole_grains, 90)),
    Trans_fat.ahei 				 = pmin(pmax((4 - Trans_fat) / (4 - 0.5), 0), 1) * 10,
    Oily_fish.ahei 				 = ahei_score_up10(Oily_fish, 2/7),
    PUFA.ahei 					 = pmin(pmax((PUFA - 2) / (10 - 2), 0), 1) * 10,
	sodium.ahei 				 = { dec <- dplyr::ntile(sodium, 10); 10 * (10 - dec) / 9 },
    Alcoholic_beverages.ahei 	 = ifelse(sex == 0, ahei_score_mid10(Alcoholic_beverages, lo = 0.5, hi = 1.5, high0 = 2.5), ahei_score_mid10(Alcoholic_beverages, lo = 0.5, hi = 2.0, high0 = 3.5))
) %>% mutate(diet.ahei.sum = rowSums1(across(ends_with(".ahei")))) %>% dplyr::select(eid, diet.ahei.sum) # ends_with(".ahei")
saveRDS(dat.ahei, paste0(indir, "/Rdata/diet.ahei.rds"))


# 🚩 PHDI (PMID: 41245532)
phdi_score_up10 <- function(x, thr) {x <- as.numeric(x); pmin(pmax(x / thr, 0), 1) * 10}
phdi_score_down10 <- function(x, thr) {x <- as.numeric(x); pmin(pmax(1 - x / thr, 0), 1) * 10}
phdi_score_down10_range <- function(x, low, high) {x <- as.numeric(x); out <- ifelse(is.na(x), NA_real_, ifelse(x <= low, 10, ifelse(x >= high, 0, 10 * (high - x) / (high - low)))); pmin(pmax(out, 0), 10)}
cols <- 'weight_Whole_grains weight_Starchy_vegetables weight_Vegetables weight_fruits weight_Dairy weight_Red_processed_meat weight_Poultry weight_Eggs weight_Fish_shellfish weight_Nuts_Seeds weight_Legumes weight_Soybean n.Added_fat_Saturated n.Added_fat_unsaturated n.Added_sugar_and_fruit_juices' %>% strsplit(' ') %>% unlist() 
dat_field.w <- field.w %>% filter (PHDI > 0) %>% mutate(name = paste0("w.", name))
dat1 <- dat.food
	grp.link <- split(dat_field.w$name, dat_field.w$Category)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE], per = 1)}
dat1 <- dat1 %>% mutate(
	energy_kcal 				   = energy.kJ * 0.239,
	weight_Whole_grains 		   = rowSums1(across(all_of(paste0("w.", c("p26076", "p26077", "p26074", "p26114"))))) + 0.5 * .data[["w.p26105"]], # g/day	
	n.Added_fat_Saturated 		   = rowSums1(across(all_of(paste0("w.", c("p26014", "p26155"))))),
	n.Added_fat_unsaturated 	   = rowSums1(across(all_of(paste0("w.", c("p26032", "p26015", "p26016"))))),
	across(c(n.Added_fat_Saturated, n.Added_fat_unsaturated), ~ (.x * 9 / energy_kcal) * 100),
	n.Added_sugar_and_fruit_juices = (w.p26012 * 4 / energy_kcal) * 100
) %>% dplyr::select(eid, all_of(cols), sex) 
dat.phdi <- dat1 %>% mutate(
	Whole_grains.phdi 				  = case_when(sex == 0 ~ phdi_score_up10(weight_Whole_grains, 75), sex == 1 ~ phdi_score_up10(weight_Whole_grains, 90), TRUE ~ NA_real_),
	Starchy_vegetables.phdi 		  = phdi_score_down10_range(weight_Starchy_vegetables, 50, 200),
	Vegetables.phdi 				  = phdi_score_up10(weight_Vegetables, 300),
	Whole_fruits.phdi 				  = phdi_score_up10(weight_fruits, 200),
	Dairy.phdi						  = phdi_score_down10_range(weight_Dairy, 250, 1000),
	Red_processed_meat.phdi 		  = phdi_score_down10_range(weight_Red_processed_meat, 14, 100),
	Poultry.phdi 					  = phdi_score_down10_range(weight_Poultry, 29, 100),
	Eggs.phdi 						  = phdi_score_down10_range(weight_Eggs, 13, 120),
	Fish_shellfish.phdi 			  = phdi_score_up10(weight_Fish_shellfish, 28),
	Nuts_Seeds.phdi 				  = phdi_score_up10(weight_Nuts_Seeds, 50),
	Legumes.phdi 					  = phdi_score_up10(weight_Legumes, 100),
	Soybean.phdi 					  = phdi_score_up10(weight_Soybean, 50),
	Added_fat_Saturated.phdi 		  = phdi_score_down10(n.Added_fat_Saturated, 10),
	Added_fat_unsaturated.phdi 		  = pmin(pmax((n.Added_fat_unsaturated - 3.5) / (21 - 3.5), 0), 1) * 10,
	Added_sugar_and_fruit_juices.phdi = phdi_score_down10_range(n.Added_sugar_and_fruit_juices, 5, 25)
) %>% mutate(
	diet.phdi.sum = rowSums1(dplyr::select(., ends_with(".phdi"), -Legumes.phdi, -Soybean.phdi)) + 0.5 * rowSums1(dplyr::select(., Legumes.phdi, Soybean.phdi))
) %>% dplyr::select(eid, diet.phdi.sum) # ends_with(".phdi")
saveRDS(dat.phdi, paste0(indir, "/Rdata/diet.phdi.rds"))


# 🚩 DI-GM (PMID: 41245532)
get_med_by_sex <- function(x, sex){tapply(x, sex, median, na.rm = TRUE)}
score_benefit <- function(x, sex, med_by_sex){med <- med_by_sex[as.character(sex)]; ifelse(is.na(x) | is.na(sex), NA_real_, ifelse(is.na(med), NA_real_, ifelse(med == 0, x > 0, x >= med)))}
score_harm <- function(x, sex, med_by_sex){med <- med_by_sex[as.character(sex)]; ifelse(is.na(x) | is.na(sex), NA_real_, ifelse(is.na(med), NA_real_, ifelse(med == 0, x == 0, x < med)))}
dat_field.w <- field.w %>% filter (DIGM > 0) %>% mutate(name = paste0("w.", name))
cols <- 'Avocados Broccoli weight_Coffee weight_Dairy Green_tea Fiber weight_Soybean weight_Whole_grains weight_Refined_grains w.Processed_meat Red_meat High_fat_diet' %>% strsplit(' ') %>% unlist()
dat1 <- dat.food
	grp.link <- split(dat_field.w$name, dat_field.w$Category)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE], per = 1)}
dat1 <- dat1 %>% mutate(
	energy_kcal = energy.kJ * 0.239, 
	Avocados = p104100, 
	Broccoli = p104140, 
	weight_Coffee = weight_Non_alcoholic_beverages, 
	Green_tea = p100420, 
	Fiber = w.p26017, 
	weight_Whole_grains = rowSums1(across(all_of(paste0("w.", c("p26076", "p26077", "p26074", "p26114"))))) + 0.5 * .data[["w.p26105"]],
	w.Processed_meat = w.p26122, 
	Red_meat = rowSums1(across(all_of(paste0("w.", c("p26066", "p26100", "p26117", "p26104"))))),
	High_fat_diet = (w.p26008 * 9) / energy_kcal * 100 
) %>% dplyr::select(eid, all_of(cols), sex)
benefit_vars <- c("Avocados","Broccoli","weight_Coffee","weight_Dairy","Green_tea", "Fiber","weight_Soybean","weight_Whole_grains")
harm_vars    <- c("weight_Refined_grains","w.Processed_meat","Red_meat")
meds <- c(lapply(dat1[benefit_vars], get_med_by_sex, sex = dat1$sex), lapply(dat1[harm_vars], get_med_by_sex, sex = dat1$sex))
dat.digm <- dat1 %>% mutate(
	across(all_of(benefit_vars), ~ score_benefit(.x, sex, meds[[cur_column()]]), .names = "{.col}.digm"),
    across(all_of(harm_vars), ~ score_harm(.x, sex, meds[[cur_column()]]), .names = "{.col}.digm"),
    High_fat_diet.digm = ifelse(is.na(High_fat_diet), NA_real_, ifelse(High_fat_diet >= 40, 0, 1))
) %>% mutate(diet.digm.sum = rowSums1(dplyr::select(., ends_with(".digm")))) %>% dplyr::select(eid, diet.digm.sum) # ends_with(".digm")
saveRDS(dat.digm, paste0(indir, "/Rdata/diet.digm.rds"))


# 🚩 MAHA
maha_score_up10 <- function(x, target) {x <- as.numeric(x); pmin(pmax(x / target, 0), 1) * 10}
maha_score_down10 <- function(x, limit_good, limit_bad) {x <- as.numeric(x); pmin(pmax((limit_bad - x) / (limit_bad - limit_good), 0), 1) * 10}
maha_score_mid10 <- function(x, lo, hi, hi_bad = Inf) {
	x <- as.numeric(x)
	ifelse(is.na(x), NA_real_, ifelse(x < lo, pmin(pmax(x / lo, 0), 1) * 10, ifelse(x <= hi, 10, ifelse(is.finite(hi_bad), pmin(pmax((hi_bad - x) / (hi_bad - hi), 0), 1) * 10, 10))))
}
maha_score_min10 <- function(x, target) {x <- as.numeric(x); pmin(pmax(x / target, 0), 1) * 10}
rowsum_food <- function(dat, vars) rowSums1(dat %>% dplyr::select(any_of(vars)), per = 1)

dairy_vars 		<- 'p26150 p26096 p102820 p102830 p102840 p102860 p102880 p102890 p102900 p102910' %>% strsplit(' ') %>% unlist()
veg_vars 		<- 'p104090 p104240 p104370 p104140 p104160 p104180 p104300 p104310 p104220 p104230 p104260 p104340 p104350 p104150 p104170 p104290 p104330 p104130 p104190 p104270 p104360 p104060 p104070 p104200 p104210 p104250 p104320 p104380' %>% strsplit(' ') %>% unlist()
fruit_vars 		<- 'p104450 p104560 p104470 p104480 p104490 p104530 p104540 p104430 p104420 p104550 p104410 p104580 p104500 p104460 p104510 p104570 p104520 p104440 p104590' %>% strsplit(' ') %>% unlist()
fat_vars 		<- 'p26110 oil p26062 p26063 p102440 p102450 p103160 p104100' %>% strsplit(' ') %>% unlist()
wholegrain_vars <- 'p100770 p100810 p100800 p100840 p100850 bread.wholemeal baguette.wholemeal bap.wholemeal roll.wholemeal bread.seeded baguette.seeded bap.seeded roll.seeded p102720 p102740 p102780' %>% strsplit(' ') %>% unlist()
refined_vars 	<- 'p100820 p101230 p101240 p101250 p101260 p101270 p102710 p102730 p102750 p102760 p102770 bread.white baguette.white bap.white roll.white p100830 p100860' %>% strsplit(' ') %>% unlist()
ssb_vars 		<- 'p100170 p100180 p100550 p100160 p100540 p100530 p100230 p100190 p100200 p100210 p100220' %>% strsplit(' ') %>% unlist()
sweet_vars 		<- 'p26064 p102120 p102140 p102150 p102170 p102180 p102190 p102200 p102210 p102220 p102230 p102260 p102270 p102280 p102290 p102300 p102310 p102320 p102330 p102340 p102350 p102360 p102370 p102380 p102060 p102010 p102020 p102050 p102070 p101970 p101980 p101990 p102030' %>% strsplit(' ') %>% unlist()
hp_vars 		<- 'p102460 p102470 p102480 p102500 p102040 p102000 p104020 p103050 p103170 p103180 p103010 p103070 p103080' %>% strsplit(' ') %>% unlist()
alcohol_vars 	<- 'p26067 p26138 p26151 p26152 p26153' %>% strsplit(' ') %>% unlist()

base_maha <- dat.food %>% mutate(
	ef = if_else(!is.na(energy.kJ) & energy.kJ > 0, 8368 / energy.kJ, NA_real_),
	protein_gkg = w.p26005 / weight,
	dairy_srv = rowsum_food(., dairy_vars) * ef,
	veg_srv = rowsum_food(., veg_vars) * ef,
	fruit_srv = rowsum_food(., fruit_vars) * ef,
	wholegrain_srv = rowsum_food(., wholegrain_vars) * ef,
	fat_srv = rowsum_food(., fat_vars) * ef,
	upf_srv = rowsum_food(., c(refined_vars, ssb_vars, sweet_vars, hp_vars)) * ef,
	alcohol_srv = rowsum_food(., alcohol_vars),
	sodium_mg = w.p26052
)

# 
maha_primary <- base_maha %>% mutate(
	protein.maha = maha_score_min10(protein_gkg, 1.2),
	dairy.maha = maha_score_up10(dairy_srv, 3),
	veg.maha = maha_score_up10(veg_srv, 3),
	fruit.maha = maha_score_up10(fruit_srv, 2),
	wholegrain.maha = maha_score_up10(wholegrain_srv, 2),
	fat.maha = maha_score_up10(fat_srv, 1),
	upf.maha = maha_score_down10(upf_srv, 0, 4),
	alcohol.maha = maha_score_down10(alcohol_srv, 0, 2.5),
	sodium.maha = maha_score_down10(sodium_mg, 2300, 5000)
) %>% mutate(diet.maha.sum = rowSums1(across(ends_with('.maha')))) %>% dplyr::select(eid, diet.maha.sum, ends_with('.maha'))

# 
maha_bal <- base_maha %>% mutate(
	protein.maha_bal = maha_score_mid10(protein_gkg, 1.2, 1.6, 2.0),
	dairy.maha_bal = maha_score_mid10(dairy_srv, 1, 2, 4),
	veg.maha_bal = maha_score_up10(veg_srv, 3),
	fruit.maha_bal = maha_score_up10(fruit_srv, 2),
	wholegrain.maha_bal = maha_score_up10(wholegrain_srv, 2),
	fat.maha_bal = maha_score_mid10(fat_srv, 0.5, 1.0, 2.0),
	upf.maha_bal = maha_score_down10(upf_srv, 0, 3),
	alcohol.maha_bal = maha_score_down10(alcohol_srv, 0, 1.5),
	sodium.maha_bal = maha_score_down10(sodium_mg, 2000, 4000)
) %>% mutate(diet.maha_bal.sum = rowSums1(across(ends_with('.maha_bal')))) %>% dplyr::select(eid, diet.maha_bal.sum, ends_with('.maha_bal'))

# 
maha_strict <- base_maha %>% mutate(
	protein.maha_strict = maha_score_mid10(protein_gkg, 0.8, 1.2, 1.8),
	dairy.maha_strict = maha_score_mid10(dairy_srv, 0.5, 1.5, 3),
	veg.maha_strict = maha_score_up10(veg_srv, 4),
	fruit.maha_strict = maha_score_up10(fruit_srv, 2.5),
	wholegrain.maha_strict = maha_score_up10(wholegrain_srv, 3),
	fat.maha_strict = maha_score_mid10(fat_srv, 0.5, 1.0, 1.5),
	upf.maha_strict = maha_score_down10(upf_srv, 0, 2.5),
	alcohol.maha_strict = maha_score_down10(alcohol_srv, 0, 1),
	sodium.maha_strict = maha_score_down10(sodium_mg, 1500, 3000)
) %>% mutate(diet.maha_strict.sum = rowSums1(across(ends_with('.maha_strict')))) %>% dplyr::select(eid, diet.maha_strict.sum, ends_with('.maha_strict'))

# 
maha_nodairy <- base_maha %>% mutate(
	protein.maha_nodairy = maha_score_min10(protein_gkg, 1.2),
	veg.maha_nodairy = maha_score_up10(veg_srv, 3),
	fruit.maha_nodairy = maha_score_up10(fruit_srv, 2),
	wholegrain.maha_nodairy = maha_score_up10(wholegrain_srv, 2),
	fat.maha_nodairy = maha_score_up10(fat_srv, 1),
	upf.maha_nodairy = maha_score_down10(upf_srv, 0, 4),
	alcohol.maha_nodairy = maha_score_down10(alcohol_srv, 0, 2.5),
	sodium.maha_nodairy = maha_score_down10(sodium_mg, 2300, 5000)
) %>% mutate(diet.maha_nodairy.sum = rowSums1(across(ends_with('.maha_nodairy'))) * 9 / 8) %>% dplyr::select(eid, diet.maha_nodairy.sum, ends_with('.maha_nodairy'))

# 
maha_noprotein <- base_maha %>% mutate(
	dairy.maha_noprotein = maha_score_up10(dairy_srv, 3),
	veg.maha_noprotein = maha_score_up10(veg_srv, 3),
	fruit.maha_noprotein = maha_score_up10(fruit_srv, 2),
	wholegrain.maha_noprotein = maha_score_up10(wholegrain_srv, 2),
	fat.maha_noprotein = maha_score_up10(fat_srv, 1),
	upf.maha_noprotein = maha_score_down10(upf_srv, 0, 4),
	alcohol.maha_noprotein = maha_score_down10(alcohol_srv, 0, 2.5),
	sodium.maha_noprotein = maha_score_down10(sodium_mg, 2300, 5000)
) %>% mutate(diet.maha_noprotein.sum = rowSums1(across(ends_with('.maha_noprotein'))) * 9 / 8) %>% dplyr::select(eid, diet.maha_noprotein.sum, ends_with('.maha_noprotein'))

dat.maha <- maha_primary %>%
	left_join(maha_bal, by = 'eid') %>%
	left_join(maha_strict, by = 'eid') %>%
	left_join(maha_nodairy, by = 'eid') %>%
	left_join(maha_noprotein, by = 'eid')

saveRDS(dat.maha, paste0(indir, '/Rdata/diet.maha.rds'))


# 🚩 upf: ultra-processed foods
# cols <- 'Beverages Dairy_products Fruits_vegetables Meat_fish_eggs Sauces_soups Starchy_foods Sugary_snacks Savory_snacks' %>% strsplit(' ') %>% unlist()
cols <- 'bev dairy veg_fruit meat_fish_egg soup starch snack_sugar snack_savory' %>% strsplit(' ') %>% unlist(); 
bev 		<- 'p26126 p26127 p26138' %>% strsplit(' ') %>% unlist()
dairy 		<- 'p26086 p26084 p26087 p26102 p26124' %>% strsplit(' ') %>% unlist()
veg_fruit	<- 'p26064 p26090 p26144' %>% strsplit(' ') %>% unlist()
meat_fish_egg 		<- 'p26122 p26137 p26145 p26149' %>% strsplit(' ') %>% unlist()
soup 		<- 'p26106 p26129 p26130 p26110' %>% strsplit(' ') %>% unlist()
starch 		<- 'p26075 p26068 p26076 p26119 p26105 p26078 p26116 p26073 p26079 p26097' %>% strsplit(' ') %>% unlist()
snack_sugar 		<- 'p26080 p26085 p26140' %>% strsplit(' ') %>% unlist()
snack_savory 		<- 'p26108 p26083 p26134' %>% strsplit(' ') %>% unlist()
weight_vars <- field.w %>% filter(startsWith(Food.Group, "w.")) %>% pull(name) #  p26000
# dat_field.w <- field.w %>% filter (upf > 0) %>% mutate(name = paste0("w.", name))
dat.upf <- dat.food %>% mutate(
	bev.upf 			= rowSums1(across(all_of(paste0("w.", bev)))),
	dairy.upf 		= rowSums1(across(all_of(paste0("w.", dairy)))),
	veg_fruit.upf 	= rowSums1(across(all_of(paste0("w.", veg_fruit)))),
	meat_fish_egg.upf 		= rowSums1(across(all_of(paste0("w.", meat_fish_egg)))),
	soup.upf 		= rowSums1(across(all_of(paste0("w.", soup)))),
	starch.upf 		= rowSums1(across(all_of(paste0("w.", starch)))),
	snack_sugar.upf 		= rowSums1(across(all_of(paste0("w.", snack_sugar)))),
	snack_savory.upf 		= rowSums1(across(all_of(paste0("w.", snack_savory)))),
	diet.total.w 					= rowSums1(across(all_of(paste0("w.", weight_vars)))),
	diet.upf.sum 					= rowSums1(across(ends_with(".upf"))) 
) %>% mutate(diet.upf.prop = diet.upf.sum / diet.total.w * 100	# % of total
) %>% dplyr::select(eid, all_of(paste0(cols, ".upf")), diet.upf.prop) 
saveRDS(dat.upf, paste0(indir, "/Rdata/diet.upf.rds"))


# 🚩 HLCD: Healthy Low-Carbohydrate Diet
score_11grp <- function(x, sex, reverse = FALSE){x <- as.numeric(x); sex <- as.numeric(sex); g <- ave(x, sex, FUN = function(v) if(all(is.na(v))) rep(NA_real_, length(v)) else dplyr::ntile(v, 11) - 1); if(reverse) g <- 10 - g; g}
cols <- 'weight_lowqcarb n.Vegetable_protein n.Added_fat_unsaturated' %>% strsplit(' ') %>% unlist()
dat_field.w <- field.w %>% filter (HLCD > 0) %>% mutate(name = paste0("w.", name))
dat1 <- dat.food
	grp.link <- split(dat_field.w$name, dat_field.w$Category)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE])}	
dat1 <- dat1 %>% mutate(
	energy_kcal 			= energy.kJ * 0.239,
	weight_lowqcarb 		= weight_Refined_grains + weight_Starchy_vegetables + w.p26064 + w.p26095 + w.p26134, 
	weight_lowqcarb 		= (weight_lowqcarb * 4 / energy_kcal) * 100,
	n.Vegetable_protein 	= (w.p26006 * 4 / energy_kcal) * 100,
	n.Added_fat_unsaturated = rowSums1(across(all_of(paste0("w.", c("p26032","p26015","p26016"))))),
	n.Added_fat_unsaturated = (n.Added_fat_unsaturated * 9 / energy_kcal) * 100
) %>% dplyr::select(eid, all_of(cols), sex)
dat.hlcd <- dat1 %>% mutate(
	lowqcarb.hlcd 			   = score_11grp(weight_lowqcarb, sex, reverse = TRUE),
	Vegetable_protein.hlcd	   = score_11grp(n.Vegetable_protein, sex),
	Added_fat_unsaturated.hlcd = score_11grp(n.Added_fat_unsaturated, sex)
) %>% mutate(diet.hlcd.sum = rowSums1(across(ends_with(".hlcd")))
) %>% dplyr::select(eid, diet.hlcd.sum)  # ends_with(".hlcd")
saveRDS(dat.hlcd, paste0(indir, "/Rdata/diet.hlcd.rds"))


# 🚩 HLFD: Healthy Low-Fat Diet
cols <- 'weight_highqcarb n.Vegetable_protein n.Added_fat_saturated' %>% strsplit(' ') %>% unlist()
dat_field.w <- field.w %>% filter(HLFD > 0) %>% mutate(name = paste0("w.", name))
dat1 <- dat.food
	grp.link <- split(dat_field.w$name, dat_field.w$Category)
	for (i in names(grp.link)) {dat1[[i]] <- rowSums1(dat1[, grp.link[[i]], drop = FALSE])}
dat1 <- dat1 %>% mutate(
	energy_kcal 		  = energy.kJ * 0.239,
	weight_highqcarb 	  = weight_Vegetables + weight_Legumes + weight_fruits + rowSums1(across(all_of(paste0("w.", c("p26076","p26077","p26074","p26114"))))) + 0.5 * w.p26105,
	weight_highqcarb 	  = (weight_highqcarb * 4 / energy_kcal) * 100,
	n.Vegetable_protein   = (w.p26006 * 4 / energy_kcal) * 100,
	n.Added_fat_saturated = (w.p26014 * 9 / energy_kcal) * 100
) %>% dplyr::select(eid, all_of(cols), sex)
dat.hlfd <- dat1 %>% mutate(
	highqcarb.hlfd 			 = score_11grp(weight_highqcarb, sex),
	Vegetable_protein.hlfd   = score_11grp(n.Vegetable_protein, sex),
	Added_fat_saturated.hlfd = score_11grp(n.Added_fat_saturated, sex, reverse = TRUE)
) %>% mutate(diet.hlfd.sum = rowSums1(across(c(highqcarb.hlfd, Added_fat_saturated.hlfd, Vegetable_protein.hlfd)))
) %>% dplyr::select(eid, diet.hlfd.sum) #ends_with(".hlfd")
saveRDS(dat.hlfd, paste0(indir, "/Rdata/diet.hlfd.rds"))


# 🚩 rEDII (reversed energy‑adjusted diet inflammatory index)
get_z  <- function(x) {x <- as.numeric(x); s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(NA_real_, length(x)) else (x - mean(x, na.rm = TRUE)) / s}
get_cp <- function(z) {(pnorm(as.numeric(z)) * 2) - 1}
dat_field.w <- field.w %>% filter(rEDII != 0, FiledID != 26002, FiledID != 26142) %>% mutate(name = paste0("w.", name))
vars <- dat_field.w$name
beta <- dat_field.w$rEDII; names(beta) <- vars
dat.redii <- dat.food %>% mutate(
	energy_kcal = energy.kJ * 0.239,
	w.p26141 = w.p26141 + w.p26142,
	across(all_of(vars), ~ .x / energy_kcal * 1000),
	across(all_of(vars), ~ get_cp(get_z(.x)) * beta[cur_column()], .names = "{.col}.redii")
) %>% mutate(
	diet.edii.sum = rowSums1(across(ends_with(".redii"))),
	diet.redii.sum = -diet.edii.sum
) %>% dplyr::select(eid, diet.redii.sum)
saveRDS(dat.redii, paste0(indir, "/Rdata/diet.redii.rds"))


# 🚩 FDS: Flavodiet Score(servings/day) 【PMID: 40456886, 39292460, 38778045】
cols <- 'tea red_wine apples berries grapes oranges grapefruit sweet_peppers onions dark_chocolate' %>% strsplit(' ') %>% unlist()
# dat_field <- field %>% filter (FDS > 0)
dat.fds <- dat.food %>% mutate(
	tea.fds 		   = p100400 + p100420,
	red_wine.fds	   = p26152,
	apples.fds		   = p104450,
	berries.fds 	   = p104470,
	grapes.fds 		   = p104500,
	oranges.fds		   = p104530,
	grapefruit.fds 	   = p104490,
	sweet_peppers.fds  = p104290,
	onions.fds 		   = p104260,
	dark_chocolate.fds = p102290
) %>% mutate(diet.fds.sum = rowSums1(across(all_of(paste0(cols, ".fds"))))) %>% dplyr::select(eid, diet.fds.sum)
saveRDS(dat.fds, paste0(indir, "/Rdata/diet.fds.rds"))
# flavonoid subclass intake(mg/day)


# 🚩 merge all together
lst <- c("medito", "medi24", "dash", "mind", "hpdi", "ahei", "phdi", "modern", "digm", "maha", 
		  "upf", "hlcd", "hlfd", "redii", "fds") 
for (l in lst) assign(l, readRDS(paste0(indir, "/Rdata/diet.", l, ".rds")) %>% mutate(eid = as.character(eid)))
dat0 <- Reduce(function(x, y) merge(x, y, by = "eid", all = TRUE), mget(lst)) %>% mutate(diet.upf.sum = diet.upf.prop) 
sum_cols <- names(dat0)[grepl("\\.sum$", names(dat0))]
lapply(dat0 %>% dplyr::select(matches("\\.sum$")), summary)
lapply(dat0 %>% dplyr::select(matches("\\.(pts|q5)$")), table, useNA = "always")
dat.fudan <- readRDS(paste0(indir, "/Rdata/diet.fudan.rds")) %>% transmute(eid = as.character(eid), energy.kJ, enroll, start_date) # number of dietary assessments
dat0 <- dat0 %>% left_join(dat.fudan, by = "eid")
saveRDS(dat0, paste0(indir, "/Rdata/diet.rds"))


# 🚩 QC
std_10_90 <- function(x){  # 10th→90th
	x <- as.numeric(x); q <- quantile(x, c(0.10, 0.90), na.rm = TRUE, names = FALSE); if(!all(is.finite(q)) || q[2] == q[1]) rep(NA_real_, length(x)) else (x - q[1]) / (q[2] - q[1])
}
score_q5 <- function(x){ # 0/25/50/75/100
	x <- as.numeric(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x); out[ok] <- c(0, 25, 50, 75, 100)[dplyr::ntile(x[ok], 5)]; out
}
score_0_100 <- function(x){
	x <- as.numeric(x); ok <- is.finite(x); if(sum(ok) < 2) rep(NA_real_, length(x)) else {mn <- min(x[ok]); mx <- max(x[ok]); if(mx == mn) rep(NA_real_, length(x)) else {out <- rep(NA_real_, length(x)); out[ok] <- (x[ok] - mn) / (mx - mn) * 100; out}}
}
score_0_100_q5 <- function(x){ # 0-100 scale first, then q5
	x <- as.numeric(x); ok <- is.finite(x); if(sum(ok) < 2) rep(NA_real_, length(x)) else {mn <- min(x[ok]); mx <- max(x[ok]); if(mx == mn) rep(NA_real_, length(x)) else {out <- rep(NA_real_, length(x)); out[ok] <- (x[ok] - mn) / (mx - mn) * 100; as.numeric(as.character(cut(out, breaks = c(-Inf, 20, 40, 60, 80, Inf), labels = c(0, 25, 50, 75, 100), right = FALSE)))}}
}
dat <- dat0 %>% mutate(
	across(all_of(sum_cols), std_10_90, .names = "{sub('\\\\.sum$', '', .col)}.std"),
    across(all_of(sum_cols), score_0_100, .names = "{sub('\\\\.sum$', '', .col)}.s100"),
	across(all_of(sum_cols), score_0_100_q5, .names = "{sub('\\\\.sum$', '', .col)}.pts"),
	across(all_of(sum_cols), score_q5, .names = "{sub('\\\\.sum$', '', .col)}.q5"),
	across(all_of(sum_cols), f3c, .names = "{sub('\\\\.sum$', '', .col)}.3c")
)
lapply(dat %>% dplyr::select(matches("\\.sum$")), summary)
lapply(dat %>% dplyr::select(matches("\\.(pts|q5)$")), table, useNA = "always")
