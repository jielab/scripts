dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
raw_dir <- file.path(dir0, "data", "nhanes", "raw")
maha_outdir <- Sys.getenv("MAHA_OUTDIR", unset = file.path(dir0, "analysis", "maha"))
.this_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.this_file <- if (length(.this_arg)) sub("^--file=", "", .this_arg[[1]]) else file.path(dir0, "scripts", "maha", "f", "MAHA_nhanes.R")
.this_dir <- dirname(normalizePath(.this_file, winslash = "/", mustWork = FALSE))
.helper_dir <- Sys.getenv("MAHA_HELPER_DIR", unset = file.path(dir0, "scripts", "0f"))
.bootstrap_files <- c(file.path(.helper_dir,"0phe.f.R"),file.path(.helper_dir,"assoc.f.R"),file.path(.helper_dir,"pred.f.R"),file.path(.helper_dir,"plot.f.R"),file.path(.this_dir,"comm.f.R"))
.bootstrap_missing <- .bootstrap_files[!file.exists(.bootstrap_files)]
if(length(.bootstrap_missing)){cat("[NHANES INPUT CHECK] Required code/helper file(s) do not exist:\n",paste0("  MISSING: ",.bootstrap_missing,collapse="\n"),"\n",sep="");stop("NHANES bootstrap input check failed.",call.=FALSE)}
pacman::p_load(tidyverse, haven, survey, survival, broom, writexl, patchwork, cowplot)
invisible(lapply(c("0phe.f.R","assoc.f.R","pred.f.R","plot.f.R"),function(f)source(file.path(.helper_dir,f))))
source(file.path(.this_dir,"comm.f.R"))
cohort_prefix <- "nhanes"
maha_auxdir <- Sys.getenv("MAHA_AUXDIR", unset = file.path(maha_outdir, cohort_prefix))
dir.create(maha_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(maha_auxdir, recursive = TRUE, showWarnings = FALSE)
cohort_file <- function(path) file.path(maha_auxdir, basename(path))
pub_file <- function(path) file.path(maha_outdir, basename(path))
.save_plot_unprefixed <- save_plot
save_plot <- function(plot,filename,...) .save_plot_unprefixed(plot,cohort_file(filename),...)
save_pub_plot <- function(plot,filename,...) .save_plot_unprefixed(plot,pub_file(filename),...)
write_xlsx <- function(x,path,...) writexl::write_xlsx(x,cohort_file(path),...)
write_pub_xlsx <- function(x,path,...) writexl::write_xlsx(x,pub_file(path),...)

cycle_info <- tibble::tribble(
	~cycle,      ~suffix, ~prefix, ~period_years,
	"1999-2000", "",     "",      2,
	"2001-2002", "_B",   "",      2,
	"2003-2004", "_C",   "",      2,
	"2005-2006", "_D",   "",      2,
	"2007-2008", "_E",   "",      2,
	"2009-2010", "_F",   "",      2,
	"2011-2012", "_G",   "",      2,
	"2013-2014", "_H",   "",      2,
	"2015-2016", "_I",   "",      2,
	"2017-2018", "_J",   "",      2,
	"2017-2020", "",     "P_",    3.2,
	"2021-2023", "_L",   "",      2
)

main_cycles <- c("1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2020", "2021-2023")
mort_cycles <- c("1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018")

use_mortality <- TRUE

maha_map <- c(
	maha = "MAHA",
	maha_bal = "MAHA-balanced",
	maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy",
	maha_noprotein = "MAHA-no protein"
)

diet.lst <- c(
	medi = "MEDI",
	dash = "DASH",
	mind = "MIND",
	maha = "MAHA",
	maha_bal = "MAHA-balanced",
	maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy",
	maha_noprotein = "MAHA-no protein"
)

dx.lst <- c(
	cad = "Coronary artery disease",
	mi = "Myocardial infarction",
	stroke = "Stroke",
	heart_failure = "Heart failure",
	t2dm = "Type 2 diabetes",
	ckd = "Chronic kidney disease",
	hypertension = "Hypertension",
	copd = "COPD/emphysema/chronic bronchitis",
	asthma = "Asthma",
	depression = "Depression",
	obesity = "Obesity",
	death = "All-cause mortality",
	heart_death = "Heart disease mortality",
	cvd_death = "CVD mortality",
	cancer_death = "Cancer mortality",
	diabetes_death = "Diabetes mortality"
)

Y.inc <- c("cad", "stroke", "heart_failure", "t2dm", "ckd", "death")
Y.mort <- c("death", "heart_death", "cvd_death", "cancer_death", "diabetes_death")
diet.inc <- c("maha", "dash", "mind", "medi")
group_pct_th <- 0.40

zstd <- function(x) { x <- as.numeric(x); s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x))); as.numeric(scale(x)) }
has_cols <- function(d, x) all(x %in% names(d))

pick_col <- function(d, x) {
	x <- x[x %in% names(d)]
	if (length(x) == 0) return(NULL)
	x[1]
}

vnum <- function(d, x) {
	nm <- pick_col(d, x)
	if (is.null(nm)) return(rep(NA_real_, nrow(d)))
	as.numeric(d[[nm]])
}

vchr <- function(d, x) {
	nm <- pick_col(d, x)
	if (is.null(nm)) return(rep(NA_character_, nrow(d)))
	as.character(d[[nm]])
}

yesno <- function(x) {
	x <- suppressWarnings(as.numeric(x))
	ifelse(x == 1, 1L, ifelse(x %in% c(2, 3), 0L, NA_integer_))
}

miss77 <- function(x) {
	x <- suppressWarnings(as.numeric(x))
	x[x %in% c(7, 9, 77, 99, 777, 999)] <- NA_real_
	x
}

rowmean_min <- function(d, vars, prop = 0.60) {
	vars <- vars[vars %in% names(d)]
	if (length(vars) == 0) return(rep(NA_real_, nrow(d)))
	m <- as.matrix(d[, vars, drop = FALSE])
	n_ok <- rowSums(!is.na(m))
	out <- rowMeans(m, na.rm = TRUE)
	out[n_ok < ceiling(length(vars) * prop)] <- NA_real_
	out
}

rowmean2 <- function(...) {
	m <- cbind(...)
	n_ok <- rowSums(!is.na(m))
	out <- rowMeans(m, na.rm = TRUE)
	out[n_ok == 0] <- NA_real_
	out
}

std_10_90 <- function(x) {
	x <- as.numeric(x)
	q <- quantile(x, c(0.1, 0.9), na.rm = TRUE, names = FALSE)
	if (!all(is.finite(q)) || q[1] == q[2]) return(rep(NA_real_, length(x)))
	(x - q[1]) / (q[2] - q[1])
}

score_q5 <- function(x) {
	x <- as.numeric(x)
	out <- rep(NA_real_, length(x))
	ok <- is.finite(x)
	if (sum(ok) < 5 || length(unique(x[ok])) < 2) return(out)
	r <- dplyr::percent_rank(x[ok])
	out[ok] <- as.numeric(as.character(cut(
		r, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
		labels = c(0, 25, 50, 75, 100), right = TRUE
	)))
	out
}

score_0_100 <- function(x) {
	x <- as.numeric(x)
	ok <- is.finite(x)
	if (sum(ok) < 2 || length(unique(x[ok])) < 2) return(rep(NA_real_, length(x)))
	out <- rep(NA_real_, length(x))
	out[ok] <- (x[ok] - min(x[ok])) / diff(range(x[ok])) * 100
	out
}

score_0_100_rank <- function(x) {
	x <- as.numeric(x)
	out <- rep(NA_real_, length(x))
	ok <- is.finite(x)
	if (sum(ok) < 5 || length(unique(x[ok])) < 2) return(out)
	out[ok] <- dplyr::percent_rank(x[ok]) * 100
	out
}

f3c <- function(x) {
	x <- as.numeric(x)
	out <- rep(NA_character_, length(x))
	ok <- is.finite(x)
	if (sum(ok) < 3 || length(unique(x[ok])) < 2) return(out)
	r <- dplyr::percent_rank(x[ok])
	out[ok] <- as.character(cut(r, breaks = c(-Inf, 1/3, 2/3, Inf), labels = c("low", "middle", "high"), right = TRUE))
	out
}

mk_hml <- function(x, p = 0.40) {
	x <- as.numeric(x)
	out <- rep(NA_character_, length(x))
	ok <- is.finite(x)
	if (sum(ok) < 3 || length(unique(x[ok])) < 2) return(out)
	r <- dplyr::percent_rank(x[ok])
	out[ok] <- case_when(r <= p ~ "low", r >= 1 - p ~ "high", TRUE ~ "middle")
	out
}

qscore <- function(x, reverse = FALSE) {
	z <- score_q5(x)
	if (reverse) z <- 100 - z
	z / 10
}

moderate_alcohol_score <- function(g, sex) {
	g <- as.numeric(g)
	sex <- as.numeric(sex)
	out <- rep(NA_real_, length(g))
	out[sex == 2 & g >= 5 & g <= 15] <- 10
	out[sex == 1 & g >= 10 & g <= 25] <- 10
	out[is.na(out) & is.finite(g)] <- 0
	out
}

weighted_mode <- function(x, w = NULL) {
	x <- x[!is.na(x)]
	if (length(x) == 0) return(NA)
	names(sort(table(x), decreasing = TRUE))[1]
}

add_pooled_weights <- function(d, cycles.use, wt_name = "wt") {
	d0 <- d %>% filter(cycle %in% cycles.use)
	if (nrow(d0) == 0) return(d0)
	total_years <- d0 %>% distinct(cycle, period_years) %>% summarise(total = sum(period_years, na.rm = TRUE)) %>% pull(total)
	d0 %>% mutate(
		!!wt_name := wt_diet * period_years / total_years,
		strata = interaction(cycle, strata0, drop = TRUE),
		psu = interaction(cycle, psu0, drop = TRUE)
	)
}

mk_stem <- function(base, cycle, suffix, prefix) {
	if (cycle == "2017-2020") return(paste0("P_", base))
	paste0(base, suffix)
}

stem_candidates <- function(base, cycle, suffix, prefix) {
	if (base == "DR1TOT") {
		if (cycle %in% c("1999-2000", "2001-2002")) return(paste0("DRXTOT", suffix))
		return(mk_stem("DR1TOT", cycle, suffix, prefix))
	}
	if (base == "DR2TOT") {
		if (cycle %in% c("1999-2000", "2001-2002")) return(character(0))
		return(mk_stem("DR2TOT", cycle, suffix, prefix))
	}
	if (base == "DR1IFF") {
		if (cycle %in% c("1999-2000", "2001-2002")) return(c(paste0("DRXIFF", suffix), paste0("DR1IFF", suffix)))
		return(mk_stem("DR1IFF", cycle, suffix, prefix))
	}
	if (base == "DR2IFF") {
		if (cycle %in% c("1999-2000", "2001-2002")) return(character(0))
		return(mk_stem("DR2IFF", cycle, suffix, prefix))
	}
	if (base == "DRXFCD") {
		if (cycle %in% c("1999-2000", "2001-2002")) return(c(paste0("DRXFCD", suffix), paste0("DRXFMT", suffix)))
		return(c(mk_stem("DRXFCD", cycle, suffix, prefix), mk_stem("DRXFMT", cycle, suffix, prefix)))
	}
	if (base == "GHB") return(unique(c(mk_stem("GHB", cycle, suffix, prefix), "LAB10", "L10_B", "L10_C")))
	if (base == "TCHOL") return(unique(c(mk_stem("TCHOL", cycle, suffix, prefix), "LAB13", "L13_B", "L13_C")))
	if (base == "HDL") return(unique(c(mk_stem("HDL", cycle, suffix, prefix), "LAB13", "L13_B", "L13_C")))
	if (base == "TRIGLY") return(unique(c(mk_stem("TRIGLY", cycle, suffix, prefix), "LAB13AM", "L13AM_B", "L13AM_C")))
	if (base == "KIQ") return(unique(c(mk_stem("KIQ_U", cycle, suffix, prefix), mk_stem("KIQ", cycle, suffix, prefix))))
	mk_stem(base, cycle, suffix, prefix)
}

find_xpt <- function(cycle, stems) {
	ff <- list.files(file.path(raw_dir, cycle), pattern = "\\.xpt$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
	if (length(ff) == 0) return(NA_character_)
	stem0 <- tools::file_path_sans_ext(basename(ff))
	ff[match(toupper(stems), toupper(stem0), nomatch = 0)][1] %||% NA_character_
}

read_base <- function(base, cycle, suffix, prefix) {
	stems <- stem_candidates(base, cycle, suffix, prefix)
	if (length(stems) == 0) return(NULL)
	f <- find_xpt(cycle, stems)
	if (is.na(f)) {
		message("missing ", cycle, " : ", base, " [", paste(stems, collapse = ", "), "]")
		return(NULL)
	}
	message("read ", f)
	d <- haven::read_xpt(f) %>% as_tibble()
	names(d) <- toupper(names(d))
	d
}

std_tot <- function(d, day = 1) {
	if (is.null(d) || !"SEQN" %in% names(d)) return(NULL)
	p1 <- if (day == 1) c("DR1", "DRX") else "DR2"
	tibble(
		SEQN = as.numeric(d$SEQN),
		kcal = vnum(d, paste0(p1, "TKCAL")),
		protein_g = vnum(d, paste0(p1, "TPROT")),
		carb_g = vnum(d, paste0(p1, "TCARB")),
		sugar_g = vnum(d, paste0(p1, "TSUGR")),
		fiber_g = vnum(d, paste0(p1, "TFIBE")),
		fat_g = vnum(d, paste0(p1, "TTFAT")),
		sfat_g = vnum(d, paste0(p1, "TSFAT")),
		mufa_g = vnum(d, paste0(p1, "TMFAT")),
		pufa_g = vnum(d, paste0(p1, "TPFAT")),
		chol_mg = vnum(d, paste0(p1, "TCHOL")),
		sodium_mg = vnum(d, paste0(p1, "TSODI")),
		alcohol_g = vnum(d, paste0(p1, "TALCO")),
		wt_diet = vnum(d, if (day == 1) c("WTDRD1PP", "WTDRD1", "WTDRD1_2YR") else c("WTDR2DPP", "WTDR2D"))
	)
}

prep_food_code <- function(fcd) {
	if (is.null(fcd) || !"SEQN" %in% names(fcd) && nrow(fcd) == 0) return(NULL)
	code_col <- names(fcd)[grepl("FDCD|FDCODE|FOODCODE", names(fcd), ignore.case = TRUE)][1]
	desc_col <- names(fcd)[grepl("DESC|DESCRIPTION|FOOD|NAME|DRXFCLD|DRXFCSD", names(fcd), ignore.case = TRUE) & !grepl("FDCD|FDCODE", names(fcd), ignore.case = TRUE)][1]
	if (is.na(code_col) || is.na(desc_col)) return(NULL)
	tibble(food_code = as.numeric(fcd[[code_col]]), food_desc = as.character(fcd[[desc_col]])) %>%
		filter(is.finite(food_code), !is.na(food_desc)) %>%
		distinct(food_code, .keep_all = TRUE)
}

food_group_day <- function(iff, fcd = NULL, day = 1) {
	if (is.null(iff) || !"SEQN" %in% names(iff)) return(NULL)
	code <- vnum(iff, if (day == 1) c("DR1IFDCD", "DRXIFDCD") else c("DR2IFDCD"))
	grams <- vnum(iff, if (day == 1) c("DR1IGRMS", "DRXIGRMS") else c("DR2IGRMS"))
	out <- tibble(SEQN = as.numeric(iff$SEQN), food_code = code, grams = grams)
	desc <- prep_food_code(fcd)
	if (!is.null(desc)) out <- out %>% left_join(desc, by = "food_code")
	if (!"food_desc" %in% names(out)) out$food_desc <- as.character(out$food_code)
	out <- out %>%
		mutate(
			fd = str_to_lower(food_desc),
			fruit = str_detect(fd, "apple|pear|banana|orange|citrus|berry|berries|grape|melon|peach|plum|fruit"),
			berry = str_detect(fd, "berry|berries|strawberry|blueberry|raspberry|blackberry"),
			vegetable = str_detect(fd, "vegetable|broccoli|spinach|lettuce|salad|carrot|tomato|onion|garlic|pepper|greens|cabbage|kale|cauliflower|squash"),
			green_leafy = str_detect(fd, "spinach|lettuce|greens|kale|collard|mustard greens|turnip greens|romaine"),
			allium = str_detect(fd, "onion|garlic|leek|shallot|chive"),
			legumes = str_detect(fd, "bean|lentil|chickpea|garbanzo|soybean|tofu|peas"),
			nuts = str_detect(fd, "nut|almond|walnut|peanut|pecan|cashew|pistachio|seed"),
			dairy = str_detect(fd, "milk|yogurt|cheese|dairy|cream"),
			lowfat_dairy = str_detect(fd, "skim|nonfat|non-fat|low fat|low-fat|reduced fat|reduced-fat|yogurt"),
			whole_grain = str_detect(fd, "whole wheat|whole grain|wholemeal|oat|oatmeal|bran|brown rice|bulgur|barley|quinoa|rye"),
			refined_grain = str_detect(fd, "white bread|white rice|pasta|noodle|cracker|biscuit|roll|bagel|muffin|tortilla|pancake|waffle"),
			fish = str_detect(fd, "fish|salmon|tuna|sardine|cod|trout|seafood|shrimp|crab|lobster"),
			poultry = str_detect(fd, "chicken|turkey"),
			red_processed_meat = str_detect(fd, "beef|pork|lamb|bacon|sausage|ham|hot dog|frankfurter|pepperoni|salami"),
			sweets_pastries = str_detect(fd, "cake|cookie|pie|doughnut|donut|pastry|candy|chocolate|ice cream|dessert|sweet|sugar"),
			ssb = str_detect(fd, "soft drink|soda|cola|fruit drink|sweetened beverage|sports drink|energy drink"),
			fried_fast = str_detect(fd, "fried|french fries|pizza|burger|fast food"),
			coffee_tea = str_detect(fd, "coffee|tea")
		)
	groups <- c("fruit","berry","vegetable","green_leafy","allium","legumes","nuts","dairy","lowfat_dairy",
		"whole_grain","refined_grain","fish","poultry","red_processed_meat","sweets_pastries","ssb","fried_fast","coffee_tea")
	out %>%
		filter(is.finite(SEQN), is.finite(grams)) %>%
		group_by(SEQN) %>%
		summarise(across(all_of(groups), ~ sum(grams * as.numeric(.x), na.rm = TRUE)), n_food_items = n(), .groups = "drop")
}

combine_day <- function(d1, d2, prefix) {
	if (is.null(d1) && is.null(d2)) return(NULL)
	if (is.null(d1)) d1 <- d2 %>% transmute(SEQN)
	if (is.null(d2)) d2 <- d1 %>% transmute(SEQN)
	names(d1)[names(d1) != "SEQN"] <- paste0(names(d1)[names(d1) != "SEQN"], "_1")
	names(d2)[names(d2) != "SEQN"] <- paste0(names(d2)[names(d2) != "SEQN"], "_2")
	d <- full_join(d1, d2, by = "SEQN")
	base <- unique(gsub("_[12]$", "", names(d)[names(d) != "SEQN"]))
	# Nutrient/food amounts are repeated-recall means, but NHANES dietary
	# weights are survey-design variables and must never be averaged.  For
	# cycles with Day 2 data, the two-day analytic sample requires both recalls
	# and uses WTDR2D.  Early one-day-only cycles use WTDRD1.
	mean_fields <- setdiff(base, "wt_diet")
	for (b in mean_fields) d[[b]] <- rowmean2(d[[paste0(b, "_1")]], d[[paste0(b, "_2")]])
	if (identical(prefix, "diet")) {
		w1 <- if ("wt_diet_1" %in% names(d)) as.numeric(d$wt_diet_1) else rep(NA_real_, nrow(d))
		w2 <- if ("wt_diet_2" %in% names(d)) as.numeric(d$wt_diet_2) else rep(NA_real_, nrow(d))
		has1 <- is.finite(w1) & w1 > 0
		has2 <- is.finite(w2) & w2 > 0
		cycle_has_day2_weights <- any(has2)
		d$n_recall_days <- as.integer(has1) + as.integer(has2)
		d$wt_diet <- if (cycle_has_day2_weights) ifelse(has1 & has2, w2, NA_real_) else ifelse(has1, w1, NA_real_)
		d$recall_weight_source <- ifelse(is.finite(d$wt_diet), if (cycle_has_day2_weights) "WTDR2D_two_complete_days" else "WTDRD1_one_day_cycle", NA_character_)
		base <- unique(c(mean_fields, "n_recall_days", "wt_diet", "recall_weight_source"))
	}
	day_cols <- grep("_[12]$", names(d), value = TRUE)
	d %>% dplyr::select(SEQN, all_of(base), all_of(day_cols))
}

read_mortality <- function() {
	mort_dir <- file.path(raw_dir, "mortality_2019_public")
	if (!dir.exists(mort_dir)) return(NULL)
	ff <- list.files(mort_dir, pattern = "NHANES_.*_MORT_2019_PUBLIC\\.dat$", full.names = TRUE)
	if (length(ff) == 0) return(NULL)
	read_one <- function(f) {
		readr::read_fwf(
			f,
			readr::fwf_positions(
				start = c(1, 15, 16, 17, 20, 21, 43, 46),
				end = c(6, 15, 16, 19, 20, 21, 45, 48),
				col_names = c("SEQN", "eligstat", "mortstat", "ucod_leading", "diabetes_mort_flag", "hypertension_mort_flag", "permth_int", "permth_exm")
			),
			col_types = "iiiiiiii",
			progress = FALSE
		) %>% mutate(mort_cycle_file = basename(f))
	}
	map_dfr(ff, read_one) %>% mutate(SEQN = as.numeric(SEQN))
}

read_cycle <- function(cycle, suffix, prefix, period_years) {
	step_header(paste0("Read cycle: ", cycle))
	demo <- read_base("DEMO", cycle, suffix, prefix)
	if (is.null(demo)) return(NULL)

	d1 <- std_tot(read_base("DR1TOT", cycle, suffix, prefix), day = 1)
	d2 <- std_tot(read_base("DR2TOT", cycle, suffix, prefix), day = 2)
	diet_tot <- combine_day(d1, d2, "diet")

	fcd <- read_base("DRXFCD", cycle, suffix, prefix)
	fg1 <- food_group_day(read_base("DR1IFF", cycle, suffix, prefix), fcd, day = 1)
	fg2 <- food_group_day(read_base("DR2IFF", cycle, suffix, prefix), fcd, day = 2)
	food <- combine_day(fg1, fg2, "food")

	bmx <- read_base("BMX", cycle, suffix, prefix)
	bpx <- read_base("BPX", cycle, suffix, prefix)
	diq <- read_base("DIQ", cycle, suffix, prefix)
	mcq <- read_base("MCQ", cycle, suffix, prefix)
	bpq <- read_base("BPQ", cycle, suffix, prefix)
	alq <- read_base("ALQ", cycle, suffix, prefix)
	smq <- read_base("SMQ", cycle, suffix, prefix)
	paq <- read_base("PAQ", cycle, suffix, prefix)
	slq <- read_base("SLQ", cycle, suffix, prefix)
	dpq <- read_base("DPQ", cycle, suffix, prefix)
	kiq <- read_base("KIQ", cycle, suffix, prefix)
	hiq <- read_base("HIQ", cycle, suffix, prefix)
	ocq <- read_base("OCQ", cycle, suffix, prefix)

	ghb <- read_base("GHB", cycle, suffix, prefix)
	glu <- read_base("GLU", cycle, suffix, prefix)
	hdl <- read_base("HDL", cycle, suffix, prefix)
	tchol <- read_base("TCHOL", cycle, suffix, prefix)
	trig <- read_base("TRIGLY", cycle, suffix, prefix)
	biopro <- read_base("BIOPRO", cycle, suffix, prefix)
	albcr <- read_base("ALB_CR", cycle, suffix, prefix)

	base <- tibble(
		SEQN = as.numeric(demo$SEQN),
		cycle = cycle,
		period_years = period_years,
		age = vnum(demo, "RIDAGEYR"),
		sex = vnum(demo, "RIAGENDR"),
		female = ifelse(vnum(demo, "RIAGENDR") == 2, 1L, ifelse(vnum(demo, "RIAGENDR") == 1, 0L, NA_integer_)),
		race = factor(vnum(demo, c("RIDRETH3", "RIDRETH1"))),
		edu = factor(miss77(vnum(demo, c("DMDEDUC2", "DMDEDUC3")))),
		marital = factor(miss77(vnum(demo, c("DMDMARTL", "DMDMARTZ")))),
		pir = vnum(demo, "INDFMPIR"),
		psu0 = vnum(demo, "SDMVPSU"),
		strata0 = vnum(demo, "SDMVSTRA"),
		wt_mec = vnum(demo, c("WTMECPRP", "WTMEC2YR")),
		wt_int = vnum(demo, c("WTINTPRP", "WTINT2YR"))
	)

	add <- list(diet_tot, food)
	add <- add[!vapply(add, is.null, logical(1))]
	for (a in add) base <- left_join(base, a, by = "SEQN")
	base$diet_food_group_source <- "description_keyword_proxy_not_FPED"

	base <- base %>%
		mutate(
			bmi = if (!is.null(bmx)) vnum(bmx, "BMXBMI")[match(SEQN, as.numeric(bmx$SEQN))] else NA_real_,
			weight_kg = if (!is.null(bmx)) vnum(bmx, "BMXWT")[match(SEQN, as.numeric(bmx$SEQN))] else NA_real_,
			waist = if (!is.null(bmx)) vnum(bmx, "BMXWAIST")[match(SEQN, as.numeric(bmx$SEQN))] else NA_real_
		)

	if (!is.null(bpx)) {
		bpdat <- bpx %>%
			transmute(
				SEQN = as.numeric(SEQN),
				sbp = rowmean2(!!!syms(intersect(c("BPXSY1","BPXSY2","BPXSY3","BPXSY4"), names(bpx)))),
				dbp = rowmean2(!!!syms(intersect(c("BPXDI1","BPXDI2","BPXDI3","BPXDI4"), names(bpx))))
			)
		base <- left_join(base, bpdat, by = "SEQN")
	} else {
		base$sbp <- NA_real_; base$dbp <- NA_real_
	}

	add_q <- function(dat, nm, vars) {
		if (is.null(dat)) return(rep(NA_real_, nrow(base)))
		vnum(dat, vars)[match(base$SEQN, as.numeric(dat$SEQN))]
	}

	base <- base %>%
		mutate(
			diabetes_q = yesno(add_q(diq, "diq", "DIQ010")),
			cad_q = pmax(yesno(add_q(mcq, "mcq", "MCQ160C")), yesno(add_q(mcq, "mcq", "MCQ160D")), na.rm = TRUE),
			mi_q = yesno(add_q(mcq, "mcq", "MCQ160E")),
			stroke_q = yesno(add_q(mcq, "mcq", "MCQ160F")),
			hf_q = yesno(add_q(mcq, "mcq", "MCQ160B")),
			asthma_q = yesno(add_q(mcq, "mcq", "MCQ010")),
			copd_q = pmax(yesno(add_q(mcq, "mcq", "MCQ160G")), yesno(add_q(mcq, "mcq", "MCQ160K")), na.rm = TRUE),
			ckd_q = yesno(add_q(kiq, "kiq", c("KIQ022", "KIQ025"))),
			htn_q = yesno(add_q(bpq, "bpq", "BPQ020")),
			smoke_ever = yesno(add_q(smq, "smq", "SMQ020")),
			smoke_now = ifelse(add_q(smq, "smq", "SMQ040") %in% c(1, 2), 1L, ifelse(add_q(smq, "smq", "SMQ040") == 3, 0L, NA_integer_)),
			health_insurance = factor(yesno(add_q(hiq, "hiq", "HIQ011"))),
			employment = factor(miss77(add_q(ocq, "ocq", c("OCQ180", "OCD150")))),
			sleep_h = add_q(slq, "slq", c("SLD012", "SLD010H", "SLD010")),
			hba1c = add_q(ghb, "ghb", c("LBXGH", "LBXGH_L")),
			glu_mgdl = add_q(glu, "glu", c("LBXGLU", "LBDGLUSI")),
			hdl_mgdl = add_q(hdl, "hdl", c("LBDHDD", "LBDHDL")),
			tc_mgdl = add_q(tchol, "tchol", c("LBXTC", "LBXTC_L")),
			tg_mgdl = add_q(trig, "trig", c("LBXTR", "LBXTR_L")),
			creat_mgdl = add_q(biopro, "biopro", c("LBXSCR", "LBDSCR")),
			uacr_mgg = add_q(albcr, "albcr", c("URDACT")),
			# Dietary outcomes use the dietary-recall design weight selected in
			# combine_day(); MEC/interview weights are not valid substitutes.
			wt_diet = as.numeric(wt_diet)
		)

	if (!is.null(dpq)) {
		dpq_vars <- intersect(paste0("DPQ0", c("10","20","30","40","50","60","70","80","90")), names(dpq))
		if (length(dpq_vars) > 0) {
			dp <- dpq %>%
				transmute(SEQN = as.numeric(SEQN), phq9 = rowSums(as.data.frame(lapply(across(all_of(dpq_vars)), miss77)), na.rm = FALSE))
			base <- left_join(base, dp, by = "SEQN")
		} else base$phq9 <- NA_real_
	} else base$phq9 <- NA_real_

	base
}

egfr_2021 <- function(scr, age, female) {
	scr <- as.numeric(scr); age <- as.numeric(age); female <- as.numeric(female)
	k <- ifelse(female == 1, 0.7, 0.9)
	alpha <- ifelse(female == 1, -0.241, -0.302)
	142 * pmin(scr / k, 1)^alpha * pmax(scr / k, 1)^(-1.200) * 0.9938^age * ifelse(female == 1, 1.012, 1)
}

construct_diet_scores <- function(dat.in, score_source = "observed_repeated_recall_mean") {
	d <- dat.in
	for (v in food_vars) if (!v %in% names(d)) d[[v]] <- NA_real_
	for (v in nutrient_vars) if (!v %in% names(d)) d[[v]] <- NA_real_

	d <- d %>%
		mutate(
			protein_g_kg = protein_g / weight_kg,
			healthy_fat_ratio = (mufa_g + pufa_g) / (sfat_g + 0.1),
			upf_proxy = refined_grain + sweets_pastries + ssb + fried_fast + red_processed_meat,
			nuts_legumes = nuts + legumes,
			other_veg = pmax(vegetable - green_leafy, 0),
			alcohol_limit = qscore(alcohol_g, reverse = TRUE),
			alcohol_moderate = moderate_alcohol_score(alcohol_g, sex),

			maha_c_protein = qscore(protein_g_kg),
			maha_c_dairy = qscore(dairy),
			maha_c_veg = qscore(vegetable),
			maha_c_fruit = qscore(fruit),
			maha_c_wholegrain = qscore(whole_grain),
			maha_c_fat = qscore(healthy_fat_ratio),
			maha_c_upf = qscore(upf_proxy, reverse = TRUE),
			maha_c_alcohol = qscore(alcohol_g, reverse = TRUE),
			maha_c_sodium = qscore(sodium_mg, reverse = TRUE),

			dash_c_fruit = qscore(fruit),
			dash_c_veg = qscore(vegetable),
			dash_c_wholegrain = qscore(whole_grain),
			dash_c_lowfatdairy = qscore(lowfat_dairy),
			dash_c_nutslegumes = qscore(nuts_legumes),
			dash_c_sodium = qscore(sodium_mg, reverse = TRUE),
			dash_c_redmeat = qscore(red_processed_meat, reverse = TRUE),
			dash_c_sweets = qscore(sweets_pastries + ssb, reverse = TRUE),

			medi_c_fruit = qscore(fruit),
			medi_c_veg = qscore(vegetable),
			medi_c_legumes = qscore(legumes),
			medi_c_wholegrain = qscore(whole_grain),
			medi_c_nuts = qscore(nuts),
			medi_c_fish = qscore(fish),
			medi_c_fat = qscore(healthy_fat_ratio),
			medi_c_alcohol = alcohol_moderate,
			medi_c_meat = qscore(red_processed_meat, reverse = TRUE),
			medi_c_dairy = qscore(dairy, reverse = TRUE),

			mind_c_green = qscore(green_leafy),
			mind_c_otherveg = qscore(other_veg),
			mind_c_berry = qscore(berry),
			mind_c_nuts = qscore(nuts),
			mind_c_wholegrain = qscore(whole_grain),
			mind_c_fish = qscore(fish),
			mind_c_poultry = qscore(poultry),
			mind_c_beans = qscore(legumes),
			mind_c_fat = qscore(healthy_fat_ratio),
			mind_c_meat = qscore(red_processed_meat, reverse = TRUE),
			mind_c_fried = qscore(fried_fast, reverse = TRUE),
			mind_c_sweets = qscore(sweets_pastries, reverse = TRUE),
			mind_c_cheese = qscore(dairy, reverse = TRUE)
		) %>%
		mutate(
			diet.maha.sum = rowmean_min(., c("maha_c_protein","maha_c_dairy","maha_c_veg","maha_c_fruit","maha_c_wholegrain","maha_c_fat","maha_c_upf","maha_c_alcohol","maha_c_sodium"), 0.60),
			diet.maha_bal.sum = rowmean_min(., c("maha_c_protein","maha_c_dairy","maha_c_veg","maha_c_fruit","maha_c_wholegrain","maha_c_fat","maha_c_upf","maha_c_sodium"), 0.60),
			diet.maha_strict.sum = rowmean_min(., c("maha_c_veg","maha_c_fruit","maha_c_wholegrain","maha_c_fat","maha_c_upf","maha_c_alcohol","maha_c_sodium", "dash_c_redmeat"), 0.60),
			diet.maha_nodairy.sum = rowmean_min(., c("maha_c_protein","maha_c_veg","maha_c_fruit","maha_c_wholegrain","maha_c_fat","maha_c_upf","maha_c_alcohol","maha_c_sodium"), 0.60),
			diet.maha_noprotein.sum = rowmean_min(., c("maha_c_dairy","maha_c_veg","maha_c_fruit","maha_c_wholegrain","maha_c_fat","maha_c_upf","maha_c_alcohol","maha_c_sodium"), 0.60),

			diet.dash.sum = rowmean_min(., c("dash_c_fruit","dash_c_veg","dash_c_wholegrain","dash_c_lowfatdairy","dash_c_nutslegumes","dash_c_sodium","dash_c_redmeat","dash_c_sweets"), 0.60),
			diet.medi.sum = rowmean_min(., c("medi_c_fruit","medi_c_veg","medi_c_legumes","medi_c_wholegrain","medi_c_nuts","medi_c_fish","medi_c_fat","medi_c_alcohol","medi_c_meat","medi_c_dairy"), 0.60),
			diet.mind.sum = rowmean_min(., c("mind_c_green","mind_c_otherveg","mind_c_berry","mind_c_nuts","mind_c_wholegrain","mind_c_fish","mind_c_poultry","mind_c_beans","mind_c_fat","mind_c_meat","mind_c_fried","mind_c_sweets","mind_c_cheese"), 0.60),
			diet_score_source = score_source
		)

	diet.sum.cols0 <- paste0("diet.", names(diet.lst), ".sum")
	diet.sum.cols0 <- diet.sum.cols0[diet.sum.cols0 %in% names(d)]
	for (x in diet.sum.cols0) {
		pre <- sub("\\.sum$", "", x)
		d[[paste0(pre, ".pts")]] <- std_10_90(d[[x]])
		d[[paste0(pre, ".s100")]] <- score_0_100(d[[x]])
		d[[paste0(pre, ".q5")]] <- score_q5(d[[x]])
		d[[paste0(pre, ".3c")]] <- f3c(d[[x]])
		d[[paste0(pre, ".hml")]] <- mk_hml(d[[paste0(pre, ".pts")]], group_pct_th)
	}
	d
}

add_analysis_weight <- function(d, cycles.use) {
	cycles.use <- intersect(cycles.use, unique(d$cycle))
	total_years <- d %>% distinct(cycle, period_years) %>% filter(cycle %in% cycles.use) %>% summarise(total = sum(period_years, na.rm = TRUE)) %>% pull(total)
	if (!is.finite(total_years) || total_years <= 0) stop("No valid cycle years for weight construction.")
	d %>%
		filter(cycle %in% cycles.use) %>%
		mutate(
			wt = wt_diet * period_years / total_years,
			strata = interaction(cycle, strata0, drop = TRUE),
			psu = interaction(cycle, psu0, drop = TRUE)
		)
}

prep_analysis_dataset <- function(d, cycles.use, require_mortality = FALSE) {
	d1 <- add_analysis_weight(d, cycles.use) %>%
		filter(age >= 20) %>%
		filter(is.finite(wt), wt > 0, !is.na(strata), !is.na(psu)) %>%
		filter(if_any(all_of(paste0("diet.", diet.inc, ".pts")), is.finite))
	if (require_mortality) {
		d1 <- d1 %>% filter(mort_eligible == 1, is.finite(death_time_y), death_time_y > 0)
	}
	d1
}

good_covs <- function(d, covs.in, min_nonmiss_prop = 0.55) {
	covs.in <- covs.in[covs.in %in% names(d)]
	keep <- vapply(covs.in, function(v) {
		x <- d[[v]]
		mean(!is.na(x)) >= min_nonmiss_prop && length(unique(x[!is.na(x)])) > 1
	}, logical(1))
	covs.in[keep]
}

make_design <- function(d) survey::svydesign(ids = ~psu, strata = ~strata, weights = ~wt, nest = TRUE, data = d)

rhs_formula <- function(x, covs.use) paste(c(x, covs.use), collapse = " + ")

term_result <- function(fit, term) {
	tt <- broom::tidy(fit, conf.int = TRUE) %>% filter(.data$term == term)
	if (nrow(tt) == 0) return(NULL)
	beta <- as.numeric(tt$estimate[1])
	se <- as.numeric(tt$std.error[1])
	lo <- if ("conf.low" %in% names(tt)) as.numeric(tt$conf.low[1]) else beta - 1.96 * se
	hi <- if ("conf.high" %in% names(tt)) as.numeric(tt$conf.high[1]) else beta + 1.96 * se
	p <- if ("p.value" %in% names(tt) && is.finite(tt$p.value[1])) as.numeric(tt$p.value[1]) else 2 * pnorm(-abs(beta / se))
	tibble(beta = beta, se = se, estimate = exp(beta), conf.low = exp(lo), conf.high = exp(hi), p.value = p)
}

fit_one <- function(dat.in, Y, X, covs.in = covs, family = c("auto", "logistic", "cox")) {
	family <- match.arg(family)
	if (!Y %in% names(dat.in) || !X %in% names(dat.in)) return(NULL)
	is_surv <- family == "cox" || (family == "auto" && Y %in% Y.mort.plot && "death_time_y" %in% names(dat.in))
	covs.use <- good_covs(dat.in, covs.in)
	if (Y == "heart_death") covs.use <- setdiff(covs.use, c("cad", "heart_failure", "stroke"))
	if (Y == "cvd_death") covs.use <- setdiff(covs.use, c("cad", "heart_failure", "stroke"))
	if (Y == "diabetes_death") covs.use <- setdiff(covs.use, "t2dm")
	need <- c(Y, X, "wt", "psu", "strata", covs.use, if (is_surv) "death_time_y")
	d0 <- dat.in %>% dplyr::select(all_of(unique(need))) %>% drop_na()
	if (is_surv) d0 <- d0 %>% filter(is.finite(death_time_y), death_time_y > 0)
	if (nrow(d0) < 200 || sum(d0[[Y]] == 1, na.rm = TRUE) < 10 || sum(d0[[Y]] == 0, na.rm = TRUE) < 10) return(NULL)
	d0[[X]] <- as.numeric(scale(d0[[X]]))
	des <- make_design(d0)
	fm <- if (is_surv) as.formula(paste0("Surv(death_time_y, ", Y, ") ~ ", rhs_formula(X, covs.use))) else as.formula(paste0(Y, " ~ ", rhs_formula(X, covs.use)))
	fit <- tryCatch(if (is_surv) survey::svycoxph(fm, design = des) else survey::svyglm(fm, design = des, family = quasibinomial()), error = function(e) NULL)
	if (is.null(fit)) return(NULL)
	tt <- term_result(fit, X)
	if (is.null(tt)) return(NULL)
	tt %>% transmute(Outcome = Y, Exposure = X, beta, se, estimate, conf.low, conf.high, p.value, N_total = nrow(d0), N_event = sum(d0[[Y]] == 1, na.rm = TRUE), model = ifelse(is_surv, "survey-weighted Cox", "survey-weighted logistic"))
}

fit_disc <- function(dat.in, Y, trad_nm, trad_lab, maha_nm = "maha", maha_lab = "MAHA", covs.in = covs, family = c("auto", "logistic", "cox")) {
	family <- match.arg(family)
	trad_var <- paste0("diet.", trad_nm, ".hml")
	maha_var <- paste0("diet.", maha_nm, ".hml")
	if (!all(c(Y, trad_var, maha_var) %in% names(dat.in))) return(NULL)
	is_surv <- family == "cox" || (family == "auto" && Y %in% Y.mort.plot && "death_time_y" %in% names(dat.in))
	covs.use <- good_covs(dat.in, covs.in)
	if (Y == "heart_death") covs.use <- setdiff(covs.use, c("cad", "heart_failure", "stroke"))
	if (Y == "cvd_death") covs.use <- setdiff(covs.use, c("cad", "heart_failure", "stroke"))
	if (Y == "diabetes_death") covs.use <- setdiff(covs.use, "t2dm")
	need <- c(Y, trad_var, maha_var, "wt", "psu", "strata", covs.use, if (is_surv) "death_time_y")
	d0 <- dat.in %>%
		dplyr::select(all_of(unique(need))) %>%
		rename(trad = all_of(trad_var), maha = all_of(maha_var)) %>%
		filter(trad %in% c("low", "high"), maha %in% c("low", "high")) %>%
		drop_na()
	if (is_surv) d0 <- d0 %>% filter(is.finite(death_time_y), death_time_y > 0)
	g.ref <- paste0(trad_lab, " high + ", maha_lab, " low")
	g.alt <- paste0(trad_lab, " low + ", maha_lab, " high")
	d0 <- d0 %>% mutate(discord = factor(case_when(trad == "high" & maha == "low" ~ g.ref, trad == "low" & maha == "high" ~ g.alt, TRUE ~ NA_character_), levels = c(g.ref, g.alt))) %>% filter(!is.na(discord))
	if (nrow(d0) < 100 || sum(d0[[Y]] == 1, na.rm = TRUE) < 10) return(NULL)
	des <- make_design(d0)
	rhs_disc <- paste(c("discord", covs.use), collapse = " + ")
	fm <- if (is_surv) as.formula(paste0("Surv(death_time_y, ", Y, ") ~ ", rhs_disc)) else as.formula(paste0(Y, " ~ ", rhs_disc))
	fit <- tryCatch(if (is_surv) survey::svycoxph(fm, design = des) else survey::svyglm(fm, design = des, family = quasibinomial()), error = function(e) NULL)
	if (is.null(fit)) return(NULL)
	tt <- broom::tidy(fit, conf.int = TRUE) %>% filter(grepl("^discord", term))
	if (nrow(tt) == 0) return(NULL)
	beta <- as.numeric(tt$estimate[1]); se <- as.numeric(tt$std.error[1])
	lo <- if ("conf.low" %in% names(tt)) as.numeric(tt$conf.low[1]) else beta - 1.96 * se
	hi <- if ("conf.high" %in% names(tt)) as.numeric(tt$conf.high[1]) else beta + 1.96 * se
	p <- if ("p.value" %in% names(tt) && is.finite(tt$p.value[1])) as.numeric(tt$p.value[1]) else 2 * pnorm(-abs(beta / se))
	tibble(MAHA = maha_lab, Pattern = trad_lab, Outcome = Y, contrast = paste0(g.alt, " vs ", g.ref), estimate = exp(beta), conf.low = exp(lo), conf.high = exp(hi), p.value = p, N_total = nrow(d0), N_event = sum(d0[[Y]] == 1, na.rm = TRUE), model = ifelse(is_surv, "survey-weighted Cox", "survey-weighted logistic"))
}

fmt_p <- function(p) ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

run_fig2_data <- function(dat.in, Ys, covs.in, family = c("logistic", "cox")) {
	family <- match.arg(family)
	Xs <- paste0("diet.", diet.inc, ".pts")
	empty_assoc <- tibble(Outcome = character(), Exposure = character(), beta = numeric(), se = numeric(), estimate = numeric(), conf.low = numeric(), conf.high = numeric(), p.value = numeric(), N_total = integer(), N_event = integer(), model = character(), Outcome_label = factor(), Diet = factor())
	empty_disc <- tibble(MAHA = character(), Pattern = factor(), Outcome = character(), contrast = character(), estimate = numeric(), conf.low = numeric(), conf.high = numeric(), p.value = numeric(), N_total = integer(), N_event = integer(), model = character(), Outcome_label = factor())
	assoc <- map_dfr(Ys, function(y) map_dfr(Xs, function(x) fit_one(dat.in, y, x, covs.in, family = family)))
	if (nrow(assoc) > 0) {
		assoc <- assoc %>% mutate(Outcome_label = factor(unname(dx.lst[Outcome]), levels = rev(unname(dx.lst[Ys]))), Diet = factor(recode(gsub("^diet\\.|\\.pts$", "", Exposure), !!!diet.lst), levels = c("MAHA", "DASH", "MIND", "MEDI")))
	} else assoc <- empty_assoc
	disc <- map_dfr(c(dash = "DASH", medi = "MEDI", mind = "MIND"), function(lab) {
		nm <- names(which(c(dash = "DASH", medi = "MEDI", mind = "MIND") == lab))
		map_dfr(Ys, function(y) fit_disc(dat.in, y, nm, lab, maha_nm = "maha", maha_lab = "MAHA", covs.in = covs.in, family = family))
	})
	if (nrow(disc) > 0) {
		disc <- disc %>% mutate(Outcome_label = factor(unname(dx.lst[Outcome]), levels = rev(unname(dx.lst[Ys]))), Pattern = factor(Pattern, levels = c("DASH", "MEDI", "MIND")))
	} else disc <- empty_disc
	list(assoc = assoc, disc = disc)
}

plot_fig2 <- function(res, Ys, title_left, title_right, xlab_left, xlab_right = NULL, file_png = NULL, file_xlsx = NULL,
	xlim_left = c(0.45, 1.35), xlim_right = c(0.25, 4.0), width = 14, height = 7.5,
	save_outputs = TRUE) {

	outcome_levels <- rev(unname(dx.lst[Ys]))
	y_axis_title <- if (length(outcome_levels) == 1) outcome_levels else NULL
	y_axis_labels <- if (length(outcome_levels) == 1) rep("", length(outcome_levels)) else outcome_levels
	diet_cols <- c(MAHA = "#F26D60", DASH = "#16B9C0", MIND = "#B97AF7", MEDI = "#7CAE00")

	clip_for_plot <- function(d, xlim) {
		d %>%
			mutate(
				conf.low.plot = pmax(conf.low, xlim[1], na.rm = TRUE),
				conf.high.plot = pmin(conf.high, xlim[2], na.rm = TRUE),
				estimate.plot = pmin(pmax(estimate, xlim[1]), xlim[2]),
				clipped = conf.low < xlim[1] | conf.high > xlim[2] | estimate < xlim[1] | estimate > xlim[2]
			)
	}

	pA.dat <- res$assoc %>% filter(Outcome %in% Ys, !is.na(Outcome_label), !is.na(Diet))
	if (nrow(pA.dat) > 0) {
		shift <- c(MAHA = -0.27, DASH = -0.09, MIND = 0.09, MEDI = 0.27)
		pA.dat <- pA.dat %>%
			mutate(y = as.numeric(Outcome_label) + unname(shift[as.character(Diet)])) %>%
			clip_for_plot(xlim_left)

		pA <- ggplot(pA.dat, aes(color = Diet)) +
			geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
			geom_segment(aes(x = conf.low.plot, xend = conf.high.plot, y = y, yend = y), linewidth = 1) +
			geom_point(aes(x = estimate.plot, y = y, shape = clipped), size = 3) +
			scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
			scale_color_manual(values = diet_cols, name = NULL) +
			scale_y_continuous(breaks = seq_along(outcome_levels), labels = y_axis_labels, expand = expansion(add = c(0.5, 0.5))) +
			scale_x_continuous(limits = xlim_left, breaks = pretty(xlim_left, n = 5), expand = expansion(mult = c(0.01, 0.03))) +
			labs(title = title_left, x = xlab_left, y = y_axis_title) +
			theme_classic(base_size = 15) +
			theme(
				plot.title = element_text(face = "bold", hjust = 0.5),
				axis.text.y = element_text(face = "bold"),
				axis.title.y = element_text(face = "bold", margin = margin(r = 4)),
				axis.title.x = element_text(face = "bold"),
				legend.position = "bottom"
			)
	} else {
		pA <- ggplot() + theme_void() + labs(title = paste0(title_left, "\nNo estimable models"))
	}

	pB.dat <- res$disc %>% filter(Outcome %in% Ys, !is.na(Outcome_label), !is.na(Pattern))
	if (nrow(pB.dat) > 0) {
		shift.b <- c(DASH = -0.20, MEDI = 0, MIND = 0.20)
		pB.dat <- pB.dat %>%
			mutate(y = as.numeric(Outcome_label) + unname(shift.b[as.character(Pattern)])) %>%
			clip_for_plot(xlim_right)

		right_xlab <- xlab_right %||% ifelse(any(grepl("Cox", pB.dat$model, ignore.case = TRUE)), "HR", "OR")
		pB <- ggplot(pB.dat, aes(color = Pattern)) +
			geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
			geom_segment(aes(x = conf.low.plot, xend = conf.high.plot, y = y, yend = y), linewidth = 1) +
			geom_point(aes(x = estimate.plot, y = y, shape = clipped), size = 3) +
			scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
			scale_color_manual(values = diet_cols[c("DASH", "MEDI", "MIND")], name = NULL) +
			scale_y_continuous(breaks = seq_along(outcome_levels), labels = rep("", length(outcome_levels)), expand = expansion(add = c(0.5, 0.5))) +
			scale_x_continuous(limits = xlim_right, breaks = pretty(xlim_right, n = 5), expand = expansion(mult = c(0.01, 0.03))) +
			labs(title = title_right, x = right_xlab, y = NULL) +
			theme_classic(base_size = 15) +
			theme(
				plot.title = element_text(face = "bold", hjust = 0.5),
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.title.x = element_text(face = "bold"),
				legend.position = "bottom"
			)
	} else {
		pB <- ggplot() + theme_void() + labs(title = paste0(title_right, "\nNo estimable discordant models"))
	}

	Fig <- pA | pB
	if (save_outputs && !is.null(file_png) && !is.null(file_xlsx)) {
		save_plot(Fig, file_png, width, height)
		write_xlsx(list(associations = res$assoc, discordant = res$disc), file_xlsx)
	}
	invisible(list(fig = Fig, A = pA, B = pB))
}

make_typical_grid_value <- function(x, w = NULL) {
	x2 <- x[!is.na(x)]
	if (length(x2) == 0) return(NA)
	if (is.numeric(x)) {
		if (!is.null(w)) return(weighted.mean(x, w, na.rm = TRUE))
		return(mean(x2, na.rm = TRUE))
	}
	if (is.factor(x)) return(factor(names(sort(table(x), decreasing = TRUE))[1], levels = levels(x)))
	names(sort(table(x), decreasing = TRUE))[1]
}

weighted_cox_joint <- function(dat.in, outcome = "death", x1 = "diet.dash.3c", x2 = "diet.maha.3c", covs.in = covs.mort, t0 = 10) {
	covs.use <- good_covs(dat.in, covs.in)
	if (outcome == "death") covs.use <- setdiff(covs.use, c("death"))
	need <- c(outcome, "death_time_y", x1, x2, "wt", "psu", covs.use)
	d0 <- dat.in %>% dplyr::select(all_of(unique(need))) %>% drop_na() %>% filter(is.finite(death_time_y), death_time_y > 0)
	if (nrow(d0) < 500 || sum(d0[[outcome]] == 1, na.rm = TRUE) < 30) return(NULL)
	d0[[x1]] <- factor(d0[[x1]], levels = c("low", "middle", "high"))
	d0[[x2]] <- factor(d0[[x2]], levels = c("low", "middle", "high"))
	fm <- as.formula(paste0("Surv(death_time_y, ", outcome, ") ~ ", x1, " * ", x2, " + ", paste(covs.use, collapse = " + ")))
	fit <- tryCatch(coxph(fm, data = d0, weights = wt, robust = TRUE, cluster = psu), error = function(e) NULL)
	if (is.null(fit)) return(NULL)
	nd <- expand.grid(a = factor(c("low", "middle", "high"), levels = c("low", "middle", "high")), b = factor(c("low", "middle", "high"), levels = c("low", "middle", "high")))
	names(nd) <- c(x1, x2)
	for (v in covs.use) nd[[v]] <- make_typical_grid_value(d0[[v]], d0$wt)
	tt <- seq(0, t0, by = 0.1)
	curves0 <- bind_rows(lapply(seq_len(nrow(nd)), function(i) {
		s <- summary(survfit(fit, newdata = nd[i, , drop = FALSE]), times = tt, extend = TRUE)
		tibble(row = i, time = s$time, risk = 1 - s$surv, lower = pmax(0, 1 - s$upper), upper = pmin(1, 1 - s$lower))
	}))
	curves <- bind_cols(curves0, nd[curves0$row, c(x1, x2)] %>% as_tibble()) %>% dplyr::select(-row)
	risk10 <- curves %>% group_by(across(all_of(c(x1, x2)))) %>% filter(abs(time - t0) == min(abs(time - t0))) %>% slice(1) %>% ungroup()
	risk10 <- risk10 %>% left_join(d0 %>% count(across(all_of(c(x1, x2))), name = "n"), by = c(x1, x2))
	list(fit = fit, curves = curves, risk10 = risk10, data = d0)
}

make_mortality_bar_panels <- function(dat.in, outcome = "death", t0 = 10) {
	res1 <- weighted_cox_joint(dat.in, outcome, "diet.dash.3c", "diet.maha.3c", covs.mort, t0)
	res2 <- weighted_cox_joint(dat.in, outcome, "diet.maha.3c", "diet.dash.3c", covs.mort, t0)
	if (is.null(res1) || is.null(res2)) stop("Fig3 mortality bars skipped: no estimable model.")
	prep_risk_plot_dat <- function(d, x_var) {
		d %>%
			mutate(
				mean = risk * 100,
				lower_plot = lower * 100,
				upper_plot = upper * 100,
				lbl = sprintf("%.1f%%", mean)
			) %>%
			group_by(.data[[x_var]]) %>%
			mutate(bg_mean = weighted.mean(mean, n, na.rm = TRUE)) %>%
			ungroup()
	}
	make_bg <- function(d, x_var) {
		d %>%
			group_by(.data[[x_var]]) %>%
			summarise(mean = weighted.mean(mean, n, na.rm = TRUE), .groups = "drop") %>%
			mutate(bg_lbl = sprintf("%.1f%%", mean), y_bg = max(d$upper_plot, mean, na.rm = TRUE) * 0.01)
	}
	plot_risk_bars <- function(d, bg, x_var, fill_var, ttl, xlab, fill_lab) {
		ggplot(d, aes(x = .data[[x_var]], y = mean, fill = .data[[fill_var]])) +
			geom_col(data = bg, aes(x = .data[[x_var]], y = mean), inherit.aes = FALSE, fill = "grey85", width = 0.92) +
			geom_col(position = position_dodge(0.8), width = 0.7) +
			geom_errorbar(aes(ymin = lower_plot, ymax = upper_plot), width = 0.2, position = position_dodge(0.8)) +
			geom_text(aes(label = lbl, y = upper_plot), vjust = -0.5, position = position_dodge(0.8), size = 3, fontface = "bold") +
			geom_text(data = bg, aes(x = .data[[x_var]], y = y_bg, label = bg_lbl), inherit.aes = FALSE, fontface = "bold", size = 3, vjust = 0) +
			scale_y_continuous(labels = \(x) sprintf("%g%%", x), expand = expansion(mult = c(0, 0.2))) +
			scale_fill_brewer(palette = "Set2") +
			labs(title = ttl, x = xlab, y = paste0(t0, "-year Risk"), fill = fill_lab) +
			theme_classic(base_size = 13) +
			theme(
				plot.title = element_text(face = "bold", hjust = 0.5),
				axis.title = element_text(face = "bold"),
				axis.title.y = element_text(face = "bold", margin = margin(r = -6)),
				axis.text = element_text(face = "bold"),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				legend.position = "top",
				legend.direction = "horizontal",
				legend.justification = "center",
				legend.box.just = "center",
				legend.text = element_text(face = "bold", color = "grey30", margin = margin(l = 4, r = 10)),
				legend.title = element_text(face = "bold", color = "grey30", margin = margin(r = 10)),
				legend.key.width = grid::unit(0.44, "cm"),
				legend.spacing.x = grid::unit(0.14, "cm")
			)
	}
	r1 <- prep_risk_plot_dat(res1$risk10, "diet.dash.3c")
	r2 <- prep_risk_plot_dat(res2$risk10, "diet.maha.3c")
	bg1 <- make_bg(r1, "diet.dash.3c")
	bg2 <- make_bg(r2, "diet.maha.3c")
	p1 <- plot_risk_bars(r1, bg1, "diet.dash.3c", "diet.maha.3c", "b. DASH gradient", "DASH", "MAHA")
	p2 <- plot_risk_bars(r2, bg2, "diet.maha.3c", "diet.dash.3c", "c. MAHA gradient", "MAHA", "DASH")
	list(C = p1, D = p2, risk10_DASH_within_MAHA = r1, risk10_MAHA_within_DASH = r2, curves_DASH_within_MAHA = res1$curves, curves_MAHA_within_DASH = res2$curves)
}

make_shift <- function(levels.in, max_shift = 0.42) {
	levels.in <- as.character(levels.in)
	if (length(levels.in) == 0) return(setNames(numeric(0), character(0)))
	if (length(levels.in) == 1) return(setNames(0, levels.in))
	setNames(seq(-max_shift, max_shift, length.out = length(levels.in)), levels.in)
}

run_fig3_assoc_data <- function(dat.in, Ys, covs.in, diet_names, diet_label_map, family = c("logistic", "cox"), min_nonmiss = 200) {
	family <- match.arg(family)
	diet_names <- intersect(diet_names, names(diet_label_map))
	pts_all <- paste0("diet.", diet_names, ".pts")
	names(pts_all) <- diet_names
	ok <- vapply(pts_all, function(x) {
		x %in% names(dat.in) &&
			sum(is.finite(dat.in[[x]]), na.rm = TRUE) >= min_nonmiss &&
			length(unique(dat.in[[x]][is.finite(dat.in[[x]])])) > 1
	}, logical(1))
	diet_use <- names(pts_all)[ok]
	pts_use <- unname(pts_all[ok])
	diag <- tibble(
		diet = diet_names,
		label = unname(diet_label_map[diet_names]),
		variable = unname(pts_all[diet_names]),
		available = unname(ok[diet_names]),
		n_nonmissing = vapply(unname(pts_all[diet_names]), function(x) {
			if (x %in% names(dat.in)) as.integer(sum(is.finite(dat.in[[x]]), na.rm = TRUE)) else 0L
		}, integer(1))
	)
	empty_assoc <- tibble(
		Outcome = character(), Exposure = character(), beta = numeric(), se = numeric(), estimate = numeric(),
		conf.low = numeric(), conf.high = numeric(), p.value = numeric(), N_total = integer(),
		N_event = integer(), model = character(), Outcome_label = factor(), Diet = factor()
	)
	if (length(pts_use) == 0) return(list(assoc = empty_assoc, diagnostics = diag))
	assoc <- map_dfr(Ys, function(y) map_dfr(pts_use, function(x) fit_one(dat.in, y, x, covs.in, family = family)))
	if (nrow(assoc) > 0) {
		lab_use <- unname(diet_label_map[diet_use])
		assoc <- assoc %>%
			mutate(
				Outcome_label = factor(unname(dx.lst[Outcome]), levels = rev(unname(dx.lst[Ys]))),
				Diet = factor(unname(diet_label_map[gsub("^diet\\.|\\.pts$", "", Exposure)]), levels = lab_use)
			)
	} else assoc <- empty_assoc
	list(assoc = assoc, diagnostics = diag)
}

run_fig3_discord_data <- function(dat.in, Ys, covs.in, compare_names, compare_label_map,
	ref_nm = "maha", ref_lab = "MAHA", family = c("logistic", "cox"), min_n_pair = 50) {
	family <- match.arg(family)
	compare_names <- intersect(compare_names, names(compare_label_map))
	ref_var <- paste0("diet.", ref_nm, ".hml")
	cmp_vars <- paste0("diet.", compare_names, ".hml")
	names(cmp_vars) <- compare_names
	pair_counts <- function(cmp_var) {
		if (!all(c(ref_var, cmp_var) %in% names(dat.in))) return(c(n_high_low = 0L, n_low_high = 0L, n_pair = 0L))
		cmp <- dat.in[[cmp_var]]
		ref <- dat.in[[ref_var]]
		c(
			n_high_low = sum(cmp == "high" & ref == "low", na.rm = TRUE),
			n_low_high = sum(cmp == "low" & ref == "high", na.rm = TRUE),
			n_pair = sum(cmp %in% c("low", "high") & ref %in% c("low", "high"), na.rm = TRUE)
		)
	}
	cnt <- t(vapply(unname(cmp_vars), pair_counts, numeric(3)))
	ok <- (ref_var %in% names(dat.in)) & (unname(cmp_vars) %in% names(dat.in)) & cnt[, "n_high_low"] >= min_n_pair & cnt[, "n_low_high"] >= min_n_pair
	diag <- tibble(
		contrast_diet = compare_names,
		label = unname(compare_label_map[compare_names]),
		compare_variable = unname(cmp_vars[compare_names]),
		ref_variable = ref_var,
		available = unname(ok),
		n_high_compare_low_MAHA = as.integer(cnt[, "n_high_low"]),
		n_low_compare_high_MAHA = as.integer(cnt[, "n_low_high"]),
		n_pair = as.integer(cnt[, "n_pair"])
	)
	empty_disc <- tibble(
		MAHA = character(), Pattern = factor(), Outcome = character(), contrast = character(),
		estimate = numeric(), conf.low = numeric(), conf.high = numeric(), p.value = numeric(),
		N_total = integer(), N_event = integer(), model = character(), Outcome_label = factor()
	)
	if (!any(ok)) return(list(disc = empty_disc, diagnostics = diag))
	cmp_use <- compare_names[ok]
	disc <- map_dfr(cmp_use, function(nm) {
		lab <- unname(compare_label_map[nm])
		map_dfr(Ys, function(y) fit_disc(dat.in, y, nm, lab, maha_nm = ref_nm, maha_lab = ref_lab, covs.in = covs.in, family = family))
	})
	if (nrow(disc) > 0) {
		level_use <- unname(compare_label_map[cmp_use])
		disc <- disc %>%
			mutate(
				Outcome_label = factor(unname(dx.lst[Outcome]), levels = rev(unname(dx.lst[Ys]))),
				Pattern = factor(Pattern, levels = level_use)
			)
	} else disc <- empty_disc
	list(disc = disc, diagnostics = diag)
}

plot_fig3_line_panel <- function(dat.plot, Ys, item_col, item_levels, item_cols, title, xlab,
	xlim = c(0.80, 1.10), show_y = TRUE, show_legend = TRUE, max_shift = 0.42, y_expand = 0.18,
	show_hr_text = FALSE, show_missing_levels = FALSE) {
	outcome_levels <- rev(unname(dx.lst[Ys]))
	y_axis_title <- if (show_y && length(outcome_levels) == 1) outcome_levels else NULL
	y_axis_labels <- if (show_y) {
		if (length(outcome_levels) == 1) rep("", length(outcome_levels)) else outcome_levels
	} else rep("", length(outcome_levels))
	item_levels <- as.character(item_levels)
	item_cols <- item_cols[item_levels]
	item_legend_labels <- sub("^MAHA-", "", item_levels)
	item_shapes <- ifelse(grepl("^MAHA-", item_levels), 21, 16)
	names(item_shapes) <- item_levels
	item_fills <- ifelse(grepl("^MAHA-", item_levels), "white", unname(item_cols))
	names(item_fills) <- item_levels
	clip_for_plot <- function(d) {
		d %>% mutate(
			conf.low.plot = pmax(conf.low, xlim[1], na.rm = TRUE),
			conf.high.plot = pmin(conf.high, xlim[2], na.rm = TRUE),
			estimate.plot = pmin(pmax(estimate, xlim[1]), xlim[2]),
			clipped = conf.low < xlim[1] | conf.high > xlim[2] | estimate < xlim[1] | estimate > xlim[2]
		)
	}
	pdat <- dat.plot %>%
		filter(Outcome %in% Ys, !is.na(Outcome_label), !is.na(.data[[item_col]])) %>%
		mutate(Item = as.character(.data[[item_col]]))
	present_levels <- item_levels[item_levels %in% unique(pdat$Item)]
	plot_levels <- if (show_missing_levels) item_levels else present_levels
	missing_levels <- setdiff(plot_levels, present_levels)
	if (nrow(pdat) == 0 && !show_missing_levels) return(ggplot() + theme_void() + labs(title = paste0(title, "\nNo estimable models")))
	pdat <- pdat %>%
		mutate(Item = factor(Item, levels = plot_levels), missing_estimate = FALSE)
	if (length(missing_levels) > 0) {
		pdat <- bind_rows(
			pdat,
			tibble(
				Outcome = Ys[1],
				Outcome_label = factor(unname(dx.lst[Ys[1]]), levels = outcome_levels),
				Item = factor(missing_levels, levels = plot_levels),
				estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
				p.value = NA_real_, N_total = NA_integer_, N_event = NA_integer_,
				model = "not estimable", missing_estimate = TRUE
			)
		)
	}
	pdat <- pdat %>%
		mutate(
			y = as.numeric(Outcome_label) - unname(make_shift(plot_levels, max_shift)[as.character(Item)]),
			hr_label = ifelse(is.finite(estimate), sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high), "not estimable")
		) %>%
		clip_for_plot()
	x_text <- if (show_hr_text) 1.07 else xlim[2] + diff(xlim) * 0.05
	x_missing <- xlim[1] + diff(xlim) * 0.03
	xlim_plot <- if (show_hr_text) c(xlim[1], xlim[2] + diff(xlim) * 0.12) else xlim
	p <- ggplot(pdat, aes(color = Item)) +
		geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
		geom_segment(data = pdat %>% filter(!missing_estimate), aes(x = conf.low.plot, xend = conf.high.plot, y = y, yend = y), linewidth = 1) +
		geom_point(data = pdat %>% filter(!missing_estimate), aes(x = estimate.plot, y = y, shape = Item, fill = Item), size = 3.2, stroke = 1.1) +
		geom_point(data = pdat %>% filter(missing_estimate), aes(x = x_missing, y = y, shape = Item, fill = Item), size = 3.2, stroke = 1.1) +
		geom_text(data = pdat %>% filter(missing_estimate), aes(x = x_missing + diff(xlim) * 0.04, y = y, label = "not estimable"), hjust = 0, size = 3, fontface = "bold", show.legend = FALSE) +
		geom_text(data = if (show_hr_text) pdat %>% filter(!missing_estimate) else pdat[0, ], aes(x = x_text, y = y, label = hr_label, color = Item), hjust = 0, size = 3.5, show.legend = FALSE) +
		scale_color_manual(values = item_cols, limits = item_levels, breaks = item_levels, labels = item_legend_labels, drop = FALSE, name = NULL) +
		scale_fill_manual(values = item_fills, limits = item_levels, breaks = item_levels, labels = item_legend_labels, drop = FALSE, name = NULL) +
		scale_shape_manual(values = item_shapes, limits = item_levels, breaks = item_levels, labels = item_legend_labels, drop = FALSE, name = NULL) +
		scale_y_continuous(breaks = seq_along(outcome_levels), labels = y_axis_labels, expand = expansion(add = c(y_expand, y_expand))) +
		scale_x_continuous(breaks = pretty(xlim, n = 5), expand = expansion(mult = c(0.01, 0.03))) +
		coord_cartesian(xlim = xlim_plot, clip = "off") +
		labs(title = title, x = xlab, y = y_axis_title) +
		theme_classic(base_size = 15) +
		theme(
			plot.title = element_text(face = "bold", hjust = 0.5),
			axis.text.y = if (show_y) element_text(face = "bold") else element_blank(),
			axis.ticks.y = if (show_y) element_line() else element_blank(),
			axis.title.y = element_text(face = "bold", margin = margin(r = 4)),
			axis.title.x = element_text(face = "bold"),
			legend.position = if (show_legend) "bottom" else "none",
			legend.text = element_text(face = "bold"),
			legend.margin = margin(0, 0, 0, 0),
			plot.margin = if (show_hr_text) margin(5.5, 18, 5.5, 5.5) else margin(5.5, 8, 5.5, 5.5)
		) +
		guides(
			color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(shape = unname(item_shapes), fill = unname(item_fills), stroke = 1.1, size = 4)),
			shape = "none",
			fill = "none"
		)
	p
}

ap_from_assoc <- function(assoc_tbl, source_label, margins = c(1.05, 1.075, 1.10), priors = c(1.05, 1.10, 1.20)) {
	if (is.null(assoc_tbl) || nrow(assoc_tbl) == 0) return(tibble())
	num_col <- function(d, nm) {
		if (!nm %in% names(d)) return(rep(NA_real_, nrow(d)))
		suppressWarnings(as.numeric(as.character(d[[nm]])))
	}
	chr_col <- function(d, nm, default = NA_character_) {
		if (!nm %in% names(d)) return(rep(default, nrow(d)))
		as.character(d[[nm]])
	}
	hr0  <- num_col(assoc_tbl, "estimate")
	lo0  <- num_col(assoc_tbl, "conf.low")
	hi0  <- num_col(assoc_tbl, "conf.high")
	b0   <- num_col(assoc_tbl, "beta")
	se0  <- num_col(assoc_tbl, "se")
	p0   <- num_col(assoc_tbl, "p.value")
	ntot <- num_col(assoc_tbl, "N_total")
	nevt <- num_col(assoc_tbl, "N_event")
	beta2 <- ifelse(is.finite(hr0) & hr0 > 0, log(hr0), b0)
	se2 <- ifelse(
		is.finite(lo0) & lo0 > 0 & is.finite(hi0) & hi0 > 0 & hi0 > lo0,
		(log(hi0) - log(lo0)) / (2 * 1.96),
		se0
	)
	base <- tibble(
		Source = source_label,
		Outcome = chr_col(assoc_tbl, "Outcome"),
		Diet = chr_col(assoc_tbl, "Diet"),
		beta = beta2,
		se = se2,
		HR = exp(beta2),
		LCI95 = exp(beta2 - 1.96 * se2),
		UCI95 = exp(beta2 + 1.96 * se2),
		P = p0,
		N_total = ntot,
		N_event = nevt
	) %>%
		filter(is.finite(beta), is.finite(se), se > 0, is.finite(HR), HR > 0)
	if (nrow(base) == 0) return(tibble())
	purrr::map_dfr(seq_len(nrow(base)), function(i) {
		d <- base[i, ]
		tost <- purrr::map_dfr(margins, function(m) ap_tost_one(d$beta, d$se, m))
		pow <- purrr::map_dfr(margins, function(m) ap_power_one(d$se, m)) %>%
			dplyr::select(margin_hr, MDE_logHR, MDE_HR_lower, MDE_HR_upper, powered_for_margin)
		bf <- tibble(
			prior_95_hr = priors,
			BF01 = vapply(priors, function(p) ap_bf01_one(d$beta, d$se, p), numeric(1)),
			BF10 = 1 / BF01,
			BF_interpretation = case_when(
				BF01 >= 3 ~ "evidence_for_point_null",
				BF01 <= 1/3 ~ "evidence_against_point_null",
				TRUE ~ "inconclusive"
			)
		)
		tost_pow <- tost %>% left_join(pow, by = "margin_hr")
		tidyr::crossing(tost_pow, bf) %>%
			mutate(
				Source = d$Source,
				Outcome = d$Outcome,
				Diet = d$Diet,
				beta = d$beta,
				se = d$se,
				HR = d$HR,
				LCI95 = d$LCI95,
				UCI95 = d$UCI95,
				P = d$P,
				N_total = d$N_total,
				N_event = d$N_event
			) %>%
			dplyr::select(Source, Outcome, Diet, N_total, N_event, beta, se, HR, LCI95, UCI95, P, everything())
	})
}

ap_fit_head2head_nhanes <- function(dat.in, Y = "death", trad_nm, trad_lab, covs.in = covs.mort, maha_nm = "maha", maha_lab = "MAHA") {
	maha_var <- paste0("diet.", maha_nm, ".pts")
	trad_var <- paste0("diet.", trad_nm, ".pts")
	if (!all(c(Y, "death_time_y", "wt", "psu", "strata", maha_var, trad_var) %in% names(dat.in))) return(tibble())
	covs.use <- good_covs(dat.in, covs.in)
	need <- unique(c(Y, "death_time_y", "wt", "psu", "strata", maha_var, trad_var, covs.use))
	d0 <- dat.in %>%
		dplyr::select(all_of(need)) %>%
		drop_na() %>%
		filter(is.finite(death_time_y), death_time_y > 0) %>%
		mutate(
			maha = as.numeric(scale(.data[[maha_var]])),
			trad = as.numeric(scale(.data[[trad_var]]))
		)
	if (nrow(d0) < 200 || sum(d0[[Y]] == 1, na.rm = TRUE) < 20) return(tibble())
	des <- make_design(d0)
	fb <- as.formula(paste0("Surv(death_time_y, ", Y, ") ~ maha + trad + ", paste(covs.use, collapse = " + ")))
	fitb <- tryCatch(survey::svycoxph(fb, design = des), error = function(e) NULL)
	if (is.null(fitb)) return(tibble())
	co <- coef(fitb); vv <- vcov(fitb)
	if (!all(c("maha", "trad") %in% names(co))) return(tibble())
	diff <- unname(co["maha"] - co["trad"])
	se_diff <- sqrt(vv["maha", "maha"] + vv["trad", "trad"] - 2 * vv["maha", "trad"])
	p_maha_mutual <- 2 * pnorm(-abs(unname(co["maha"]) / sqrt(vv["maha", "maha"])))
	p_trad_mutual <- 2 * pnorm(-abs(unname(co["trad"]) / sqrt(vv["trad", "trad"])))
	tibble(
		Outcome = unname(dx.lst[Y]),
		Comparison = paste0(maha_lab, " vs ", trad_lab),
		N_total = nrow(d0),
		N_event = sum(d0[[Y]] == 1, na.rm = TRUE),
		HR_maha_mutual = exp(unname(co["maha"])),
		HR_trad_mutual = exp(unname(co["trad"])),
		P_maha_given_trad = p_maha_mutual,
		P_trad_given_maha = p_trad_mutual,
		logHR_difference_maha_minus_trad = diff,
		HR_ratio_maha_vs_trad = exp(diff),
		SE_difference = se_diff,
		P_equal_coefficients = 2 * pnorm(-abs(diff / se_diff))
	)
}

ap_clean_equiv <- function(x, dataset_label = "NHANES", margin = ap_margin_main, prior = ap_prior_main) {
	if (is.null(x) || nrow(x) == 0) return(tibble())
	x %>%
		filter(abs(margin_hr - margin) < 1e-8, abs(prior_95_hr - prior) < 1e-8) %>%
		mutate(
			Dataset = dataset_label,
			Outcome = ifelse(Outcome %in% names(dx.lst), unname(dx.lst[Outcome]), as.character(Outcome)),
			Diet = as.character(Diet),
			HR_95CI = ap_ci_lab(HR, LCI95, UCI95),
			HR_90CI = ap_ci_lab(HR, CI90_low_HR, CI90_high_HR),
			TOST_P_label = ap_p_lab(TOST_p),
			BF01_label = ap_bf_lab(BF01),
			MDE80_HR_label = sprintf("%.2f-%.2f", MDE_HR_lower, MDE_HR_upper),
			Equivalence_5pct = ifelse(equivalent, "Equivalent within HR 0.95-1.05", "Not equivalent/inconclusive"),
			Powered_5pct = ifelse(powered_for_margin, "Powered for HR 1.05", "Not powered for HR 1.05"),
			BF_summary = case_when(
				BF01 >= 3 ~ "BF favors near-zero",
				BF01 <= 1/3 ~ "BF favors non-zero",
				TRUE ~ "BF inconclusive"
			)
		)
}

ap_clean_head2head <- function(x, dataset_label = "NHANES") {
	if (is.null(x) || nrow(x) == 0) return(tibble())
	x %>%
		mutate(
			Dataset = dataset_label,
			Outcome = as.character(Outcome),
			Comparison = as.character(Comparison),
			Comp2 = factor(gsub("^MAHA vs ", "", Comparison), levels = c("DASH", "MIND", "MEDI")),
			LCI_ratio = exp(logHR_difference_maha_minus_trad - 1.96 * SE_difference),
			UCI_ratio = exp(logHR_difference_maha_minus_trad + 1.96 * SE_difference),
			Ratio_95CI = ap_ci_lab(HR_ratio_maha_vs_trad, LCI_ratio, UCI_ratio),
			P_equal_label = ap_p_lab(P_equal_coefficients)
		)
}

ap_plot_maha_equiv <- function(d, title = "d. MAHA equivalence", margin = ap_margin_main) {
	d <- d %>% filter(Diet %in% ap_maha_levels) %>% mutate(Diet = factor(Diet, levels = rev(ap_maha_levels)))
	if (nrow(d) == 0) return(ggplot() + theme_void() + labs(title = title))
	xlim <- c(0.90, 1.20)
	x_txt <- 1.105
	ggplot(d, aes(y = Diet, color = Diet)) +
		annotate("rect", xmin = 1 / margin, xmax = margin, ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "grey70") +
		geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.7) +
		geom_vline(xintercept = c(1 / margin, margin), linetype = 3, color = "grey55", linewidth = 0.55) +
		geom_segment(aes(x = CI90_low_HR, xend = CI90_high_HR, yend = Diet), linewidth = 1) +
		geom_point(aes(x = HR, shape = Equivalence_5pct), size = 3, stroke = 1.0, fill = "white") +
		geom_text(aes(x = x_txt, label = paste0(HR_90CI, "; TOST P=", TOST_P_label)), hjust = 0, size = 3.0, fontface = "bold", show.legend = FALSE) +
		scale_color_manual(values = ap_maha_cols, guide = "none") +
		scale_shape_manual(values = c("Equivalent within HR 0.95-1.05" = 16, "Not equivalent/inconclusive" = 1), drop = FALSE) +
		scale_x_continuous(limits = xlim, breaks = c(0.9, 1.0, 1.1, 1.2), expand = expansion(mult = c(0, 0))) +
		coord_cartesian(xlim = c(xlim[1], 1.235), clip = "off") +
		labs(title = title, x = "MAHA hazard ratio with 90% CI", y = NULL, shape = NULL) +
		ap_theme(13) +
		theme(
			legend.position = "bottom",
			plot.title = element_text(face = "bold", hjust = 0.5),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		)
}

ap_plot_bf <- function(d, title = "e. Bayes evidence for near-zero") {
	d <- d %>% filter(Diet %in% ap_maha_levels) %>% mutate(Diet = factor(Diet, levels = rev(ap_maha_levels)), log10BF01 = log10(BF01))
	if (nrow(d) == 0) return(ggplot() + theme_void() + labs(title = title))
	lo_thr <- log10(1 / 3)
	hi_thr <- log10(3)
	d <- d %>% mutate(
		label_x = dplyr::case_when(
			Diet == "MAHA-balanced" ~ pmax(log10BF01 - 0.02, -0.48),
			TRUE ~ pmin(log10BF01 + 0.02, 0.47)
		),
		label_hjust = dplyr::if_else(Diet == "MAHA-balanced", 1, 0)
	)
	ggplot(d, aes(y = Diet, color = Diet)) +
		annotate("rect", xmin = lo_thr, xmax = hi_thr, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.10) +
		geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.7) +
		geom_vline(xintercept = c(lo_thr, hi_thr), linetype = 3, color = "grey55", linewidth = 0.55) +
		geom_segment(aes(x = 0, xend = log10BF01, yend = Diet), linewidth = 1) +
		geom_point(aes(x = log10BF01, shape = BF_summary), size = 3, fill = "white", stroke = 1.0) +
		geom_text(aes(x = label_x, label = paste0("BF01=", BF01_label), hjust = label_hjust), size = 3.0, fontface = "bold", show.legend = FALSE) +
		scale_color_manual(values = ap_maha_cols, guide = "none") +
		scale_shape_manual(values = c("BF favors near-zero" = 17, "BF inconclusive" = 1, "BF favors non-zero" = 16), drop = FALSE) +
		scale_x_continuous(limits = c(-0.50, 0.50), breaks = c(-0.5, -0.25, 0, 0.25, 0.5), expand = expansion(mult = c(0, 0))) +
		coord_cartesian(xlim = c(-0.50, 0.50), clip = "off") +
		labs(title = title, x = expression(log[10](BF[01])), y = NULL, shape = NULL) +
		ap_theme(13) +
		theme(
			legend.position = "bottom",
			plot.title = element_text(face = "bold", hjust = 0.5),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		)
}

ap_plot_power <- function(d, title = "f. MAHA detectable effect", margin = ap_margin_main) {
	d <- d %>% filter(Diet %in% ap_maha_levels) %>% mutate(Diet = factor(Diet, levels = rev(ap_maha_levels)))
	if (nrow(d) == 0) return(ggplot() + theme_void() + labs(title = title))
	xmin <- max(1.045, min(c(margin, d$MDE_HR_upper), na.rm = TRUE) - 0.005)
	xmax <- min(1.105, max(c(margin, d$MDE_HR_upper), na.rm = TRUE) + 0.010)
	ggplot(d, aes(y = Diet, color = Diet)) +
		geom_vline(xintercept = margin, linetype = 3, color = "grey55", linewidth = 0.65) +
		geom_segment(aes(x = margin, xend = MDE_HR_upper, yend = Diet), linewidth = 1.2, alpha = 0.9) +
		geom_point(aes(x = MDE_HR_upper, shape = Powered_5pct), size = 3) +
		geom_text(aes(x = pmin(MDE_HR_upper + 0.0013, xmax - 0.001), label = MDE80_HR_label), hjust = 0, size = 3.2, fontface = "bold", show.legend = FALSE) +
		scale_color_manual(values = ap_maha_cols, guide = "none") +
		scale_shape_manual(values = c("Powered for HR 1.05" = 16, "Not powered for HR 1.05" = 17), drop = FALSE) +
		scale_x_continuous(limits = c(xmin, xmax), breaks = seq(1.05, 1.10, by = 0.01), expand = expansion(mult = c(0, 0))) +
		coord_cartesian(xlim = c(xmin, xmax), clip = "off") +
		labs(title = title, x = "Minimum detectable HR at 80% power (threshold = 1.05)", y = NULL, shape = NULL) +
		ap_theme(13) +
		theme(
			legend.position = "bottom",
			plot.title = element_text(face = "bold", hjust = 0.5),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank()
		)
}

ap_plot_head2head <- function(d, title = "g. Direct attenuation versus established scores") {
	if (nrow(d) == 0) return(ggplot() + theme_void() + labs(title = title))
	d <- d %>% mutate(Comparison = factor(Comparison, levels = rev(c("MAHA vs DASH", "MAHA vs MIND", "MAHA vs MEDI"))))
	comp_cols <- c(DASH = fig3.shared_cols[["DASH"]], MIND = fig3.shared_cols[["MIND"]], MEDI = fig3.shared_cols[["MEDI"]])
	xmax <- 1.70
	ggplot(d, aes(y = Comparison, color = Comp2)) +
		annotate("rect", xmin = 1.00, xmax = xmax, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.10) +
		geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.75) +
		geom_segment(aes(x = LCI_ratio, xend = UCI_ratio, yend = Comparison), linewidth = 1.3) +
		geom_point(aes(x = HR_ratio_maha_vs_trad, size = -log10(pmax(P_equal_coefficients, 1e-12))), shape = 16) +
		geom_text(aes(x = xmax - 0.015, label = paste0(Ratio_95CI, "; Pdiff=", P_equal_label)), hjust = 1, size = 3.2, fontface = "bold", show.legend = FALSE) +
		scale_color_manual(values = comp_cols, breaks = c("DASH", "MIND", "MEDI")) +
		scale_size_continuous(range = c(3, 5), guide = "none") +
		scale_x_continuous(limits = c(1.00, xmax), breaks = seq(1.0, xmax, by = 0.1), expand = expansion(mult = c(0, 0))) +
		coord_cartesian(xlim = c(1.00, xmax), clip = "off") +
		labs(title = title, x = "HR ratio: MAHA coefficient / established-score coefficient", y = NULL, color = NULL) +
		ap_theme(13) +
		theme(
			legend.position = "bottom",
			plot.title = element_text(face = "bold", hjust = 0.5),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			plot.margin = margin(5.5, 16, 5.5, 5.5)
		)
}

nhanes_component_keys <- c("protein", "dairy", "veg", "fruit", "wholegrain", "fat", "upf", "alcohol", "sodium")
nhanes_component_cols <- paste0("maha_c_", nhanes_component_keys)
nhanes_component_labels <- c(
  protein="Protein", dairy="Dairy", veg="Vegetables", fruit="Fruit", wholegrain="Whole grains",
  fat="Healthy fats", upf="Low UPF/refined foods", alcohol="Low alcohol", sodium="Low sodium"
)
nhanes_jobs <- suppressWarnings(as.integer(Sys.getenv("MAHA_JOBS", unset=ifelse(.Platform$OS.type=="windows","1","4"))))
if (!is.finite(nhanes_jobs) || nhanes_jobs < 1) nhanes_jobs <- 1L
nhanes_null_max <- suppressWarnings(as.integer(Sys.getenv("MAHA_NULL_MAX", unset="0")))
if (!is.finite(nhanes_null_max) || nhanes_null_max < 0) nhanes_null_max <- 0L

nhanes_make_sign_grid <- function(p, max_null=0L) {
  g <- expand.grid(rep(list(c(-1,1)), p), KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
  names(g) <- nhanes_component_keys[seq_len(p)]
  obs <- which(rowSums(g==1)==p); nul <- setdiff(seq_len(nrow(g)), obs)
  if (max_null>0L && length(nul)>max_null) nul <- nul[unique(round(seq(1,length(nul),length.out=max_null)))]
  g[c(obs,nul),,drop=FALSE]
}
nhanes_parallel_lapply <- function(X,FUN) {
  if (.Platform$OS.type!="windows" && nhanes_jobs>1L) parallel::mclapply(X,FUN,mc.cores=nhanes_jobs,mc.preschedule=TRUE) else lapply(X,FUN)
}

nhanes_validate_maha <- function(dat.in, covs.in) {
  step_header("Reviewer validation before NHANES main analysis: MAHA is not an arbitrary score")
  need <- c("SEQN","death","death_time_y","wt","psu","strata",nhanes_component_cols,
            "diet.dash.sum","diet.mind.sum","diet.medi.sum")
  miss <- setdiff(need,names(dat.in)); if(length(miss)) stop("NHANES validation missing columns: ",paste(miss,collapse=", "),call.=FALSE)
  cov.use <- good_covs(dat.in,covs.in)
  d0 <- dat.in %>% dplyr::select(all_of(unique(c(need,cov.use)))) %>%
    filter(is.finite(death_time_y),death_time_y>0,is.finite(wt),wt>0) %>% drop_na()
  if(nrow(d0)<1000 || sum(d0$death==1)<100) stop("NHANES validation sample unexpectedly small: N=",nrow(d0)," deaths=",sum(d0$death==1),call.=FALSE)
  message("[NHANES VALIDATE] common complete-case N=",nrow(d0),"; deaths=",sum(d0$death==1),"; covariates=",paste(cov.use,collapse=","))
  X <- as.matrix(d0[,nhanes_component_cols,drop=FALSE]); storage.mode(X)<-"double"
  primary <- rowMeans(X)
  anchor <- rowMeans(cbind(zstd(d0$diet.dash.sum),zstd(d0$diet.mind.sum),zstd(d0$diet.medi.sum)))
  patterns <- nhanes_make_sign_grid(ncol(X),nhanes_null_max)
  fit_pattern <- function(i) {
    sg <- as.numeric(patterns[i,]); raw_score <- as.numeric(X %*% sg)/ncol(X)
    dd <- d0; dd$score <- zstd(raw_score); des <- make_design(dd)
    fm <- as.formula(paste0("Surv(death_time_y, death) ~ score",if(length(cov.use)) paste0(" + ",paste(cov.use,collapse=" + ")) else ""))
    fit <- tryCatch(survey::svycoxph(fm,design=des),error=function(e)NULL)
    beta<-se<-p<-NA_real_
    if(!is.null(fit)) {tt<-broom::tidy(fit)%>%filter(term=="score");if(nrow(tt)){beta<-tt$estimate[1];se<-tt$std.error[1];p<-tt$p.value[1]}}
    tibble(pattern_id=i,pattern=paste(ifelse(sg>0,"+","-"),collapse=""),is_prespecified_MAHA=all(sg==1),N=nrow(dd),deaths=sum(dd$death==1),
           beta=beta,se=se,HR=exp(beta),CI_low=exp(beta-1.96*se),CI_high=exp(beta+1.96*se),P=p,protective_Z=-beta/se,
           rho_established_anchor=suppressWarnings(cor(raw_score,anchor,method="spearman",use="complete.obs")))
  }
  message("[NHANES VALIDATE] fitting ",nrow(patterns)," matched direction patterns; workers=",nhanes_jobs)
  pattern_res <- bind_rows(nhanes_parallel_lapply(seq_len(nrow(patterns)),fit_pattern))
  obs <- pattern_res%>%filter(is_prespecified_MAHA); nul<-pattern_res%>%filter(!is_prespecified_MAHA)
  summary_tbl <- tibble(Cohort="NHANES",N=obs$N,Deaths=obs$deaths,MAHA_HR_per_SD=obs$HR,MAHA_CI_low=obs$CI_low,MAHA_CI_high=obs$CI_high,
    Direction_null_n=nrow(nul),MAHA_protective_Z=obs$protective_Z,
    Mortality_percentile_vs_direction_null=mean(nul$protective_Z<=obs$protective_Z,na.rm=TRUE),
    Mortality_empirical_one_sided_P=(1+sum(nul$protective_Z>=obs$protective_Z,na.rm=TRUE))/(1+sum(is.finite(nul$protective_Z))),
    MAHA_anchor_Spearman=obs$rho_established_anchor,
    Anchor_percentile_vs_direction_null=mean(nul$rho_established_anchor<=obs$rho_established_anchor,na.rm=TRUE),
    Anchor_empirical_one_sided_P=(1+sum(nul$rho_established_anchor>=obs$rho_established_anchor,na.rm=TRUE))/(1+sum(is.finite(nul$rho_established_anchor))))
  loo_specs <- c("Primary MAHA",paste0("Without ",unname(nhanes_component_labels))); loo_drop<-c(NA_character_,nhanes_component_keys)
  loo <- purrr::map2_dfr(loo_specs,loo_drop,function(lab,drop_key){
    keep<-if(is.na(drop_key))nhanes_component_keys else setdiff(nhanes_component_keys,drop_key); raw_score<-rowMeans(as.matrix(d0[,paste0("maha_c_",keep),drop=FALSE]))
    dd<-d0;dd$score<-zstd(raw_score);des<-make_design(dd);fm<-as.formula(paste0("Surv(death_time_y, death) ~ score",if(length(cov.use))paste0(" + ",paste(cov.use,collapse=" + "))else""))
    fit<-survey::svycoxph(fm,design=des);tt<-broom::tidy(fit)%>%filter(term=="score")
    tibble(label=lab,omitted=ifelse(is.na(drop_key),"None",drop_key),N=nrow(dd),deaths=sum(dd$death==1),beta=tt$estimate[1],se=tt$std.error[1],HR=exp(tt$estimate[1]),
      CI_low=exp(tt$estimate[1]-1.96*tt$std.error[1]),CI_high=exp(tt$estimate[1]+1.96*tt$std.error[1]),P=tt$p.value[1],rho_primary=suppressWarnings(cor(raw_score,primary,method="spearman")))
  })
  p1<-ggplot(nul,aes(protective_Z))+geom_histogram(bins=35,boundary=0)+geom_vline(xintercept=obs$protective_Z,linewidth=.9)+
    labs(title="a. Mortality criterion validity",subtitle="Matched scores use identical 9 components; only directions differ",x="Protective Z statistic (-beta/SE)",y="Matched-null scores")+theme_classic(base_size=11)
  p2<-ggplot(nul,aes(rho_established_anchor))+geom_histogram(bins=35,boundary=0)+geom_vline(xintercept=obs$rho_established_anchor,linewidth=.9)+
    labs(title="b. Convergent construct validity",subtitle="Anchor = mean z(DASH, MIND, Mediterranean diet)",x="Spearman rho with established-diet anchor",y="Matched-null scores")+theme_classic(base_size=11)
  lp<-loo%>%mutate(label=factor(label,levels=rev(label)))
  p3<-ggplot(lp,aes(HR,label))+geom_vline(xintercept=1,linetype=2)+geom_segment(aes(x=CI_low,xend=CI_high,yend=label),linewidth=.6)+geom_point(size=2)+
    labs(title="c. Leave-one-component-out mortality",x="Survey-weighted HR per 1 SD",y=NULL)+theme_classic(base_size=11)
  p4<-ggplot(lp,aes(rho_primary,label))+geom_segment(aes(x=0,xend=rho_primary,yend=label),linewidth=.45)+geom_point(size=2)+coord_cartesian(xlim=c(0,1))+
    labs(title="d. Rank stability after omitting one domain",x="Spearman rho with primary MAHA",y=NULL)+theme_classic(base_size=11)
  fig<-(p1|p2)/(p3|p4)+plot_annotation(title="NHANES: prespecified MAHA versus matched arbitrary scores")
  save_plot(fig,"FigS1.validate.png",width=14,height=10.5,dpi=500,bg="white")
  write_xlsx(list(summary=summary_tbl,direction_null=pattern_res,leave_one_component_out=loo),"FigS1.validate.out.xlsx")
  print(summary_tbl); invisible(list(summary=summary_tbl,pattern=pattern_res,loo=loo))
}

nhanes_construct_profile <- function(dat.in) {
  step_header("Reviewer validation before NHANES main analysis: MAHA construct profile")
  need<-c(nhanes_component_cols,"diet.maha.sum","diet.dash.sum","diet.mind.sum","diet.medi.sum")
  miss<-setdiff(need,names(dat.in));if(length(miss))stop("NHANES construct profile missing: ",paste(miss,collapse=", "),call.=FALSE)
  d<-dat.in%>%dplyr::select(all_of(need))%>%drop_na(); X<-as.data.frame(d[,nhanes_component_cols,drop=FALSE]);names(X)<-nhanes_component_keys
  X$maha<-d$diet.maha.sum;X$decile<-dplyr::ntile(X$maha,10)
  profile<-X%>%group_by(decile)%>%summarise(across(all_of(nhanes_component_keys),mean),.groups="drop")%>%pivot_longer(all_of(nhanes_component_keys),names_to="component",values_to="mean_score")%>%mutate(Component=unname(nhanes_component_labels[component]))
  metrics<-purrr::map_dfr(nhanes_component_keys,function(k){others<-setdiff(nhanes_component_keys,k);rest<-rowMeans(as.matrix(X[,others,drop=FALSE]));
    tibble(component=k,Component=nhanes_component_labels[[k]],item_rest_rho=suppressWarnings(cor(X[[k]],rest,method="spearman",use="complete.obs")),rho_primary_MAHA=suppressWarnings(cor(X[[k]],X$maha,method="spearman",use="complete.obs")),
      mean_decile1=mean(X[[k]][X$decile==1]),mean_decile10=mean(X[[k]][X$decile==10]),D10_minus_D1=mean(X[[k]][X$decile==10])-mean(X[[k]][X$decile==1]))})
  conv<-tibble(Comparator=c("DASH","MIND","Mediterranean diet","Established-diet anchor"),Spearman_rho=c(
    cor(d$diet.maha.sum,d$diet.dash.sum,method="spearman"),cor(d$diet.maha.sum,d$diet.mind.sum,method="spearman"),cor(d$diet.maha.sum,d$diet.medi.sum,method="spearman"),
    cor(zstd(d$diet.maha.sum),rowMeans(cbind(zstd(d$diet.dash.sum),zstd(d$diet.mind.sum),zstd(d$diet.medi.sum))),method="spearman")))
  p1<-ggplot(profile,aes(decile,mean_score,group=Component,color=Component))+geom_line(linewidth=.7)+geom_point(size=1.2)+scale_x_continuous(breaks=1:10)+
    labs(title="a. All nine domains across MAHA deciles",x="Primary MAHA decile",y="Mean component score (0-10)",color=NULL)+theme_classic(base_size=11)
  mp<-metrics%>%mutate(Component=factor(Component,levels=rev(Component)))
  p2<-ggplot(mp,aes(D10_minus_D1,Component))+geom_vline(xintercept=0,linetype=2)+geom_segment(aes(x=0,xend=D10_minus_D1,yend=Component))+geom_point(size=2)+labs(title="b. Component separation",x="Mean score difference: decile 10 - decile 1",y=NULL)+theme_classic(base_size=11)
  p3<-ggplot(mp,aes(item_rest_rho,Component))+geom_vline(xintercept=0,linetype=2)+geom_point(size=2)+labs(title="c. Corrected item-rest correlation",subtitle="Descriptive only: MAHA is a formative index",x="Spearman rho",y=NULL)+theme_classic(base_size=11)
  p4<-ggplot(conv,aes(Spearman_rho,reorder(Comparator,Spearman_rho)))+geom_segment(aes(x=0,xend=Spearman_rho,yend=reorder(Comparator,Spearman_rho)))+geom_point(size=2.3)+coord_cartesian(xlim=c(0,1))+labs(title="d. Convergent validity",x="Spearman rho",y=NULL)+theme_classic(base_size=11)
  fig<-(p1|p2)/(p3|p4)+plot_annotation(title="NHANES: MAHA construct profile before outcome benchmarking")
  save_plot(fig,"FigS2.construct_profile.png",width=14,height=10.5,dpi=500,bg="white")
  write_xlsx(list(decile_profile=profile,component_metrics=metrics,convergent_validity=conv),"FigS2.construct_profile.out.xlsx")
  invisible(list(profile=profile,metrics=metrics,convergent=conv))
}

nhanes_check_inputs <- function() {
  cat("\n[NHANES INPUT CHECK] Required raw inputs:\n")
  if(!dir.exists(raw_dir)){cat("  MISSING: NHANES raw directory does not exist: ",raw_dir,"\n",sep="");stop("NHANES raw directory missing.",call.=FALSE)}
  req_cycles<-unique(c(main_cycles,if(use_mortality)mort_cycles else character()))
  req_bases<-c("DEMO","BMX","DR1TOT","DR1IFF","DRXFCD")
  missing<-character()
  for(cyc in req_cycles){info<-cycle_info%>%filter(cycle==cyc);if(!nrow(info)){missing<-c(missing,paste0(cyc,": cycle metadata"));next}
    for(base in req_bases){st<-stem_candidates(base,cyc,info$suffix[1],info$prefix[1]);f<-find_xpt(cyc,st)
      if(is.na(f)){lab<-paste0(cyc,"/",base," [",paste(st,collapse=" | "),"]");missing<-c(missing,lab);cat("  MISSING: ",lab," 文件不存在\n",sep="")}
      else cat("  OK: ",cyc,"/",base," -> ",f,"\n",sep="")
    }}
  if(use_mortality){md<-file.path(raw_dir,"mortality_2019_public");mf<-if(dir.exists(md))list.files(md,pattern="NHANES_.*_MORT_2019_PUBLIC\\.dat$",full.names=TRUE)else character()
    if(!length(mf)){missing<-c(missing,paste0("mortality files under ",md));cat("  MISSING: linked mortality .dat 文件不存在: ",md,"\n",sep="")} else cat("  OK: linked mortality files = ",length(mf),"\n",sep="")}
  if(length(missing)){cat("[NHANES INPUT CHECK] STOP: ",length(missing)," required input(s) missing.\n",sep="");stop("NHANES required input file(s) missing. See log above.",call.=FALSE)}
  cat("[NHANES INPUT CHECK] PASS\n\n")
  invisible(TRUE)
}

setwd2(maha_outdir)

options(width = 200, warn = 1, survey.lonely.psu = "adjust")

log_file <- cohort_file("maha.log")
if (file.exists(log_file)) invisible(file.remove(log_file))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
on.exit({ while (sink.number() > 0) sink(); close(log_con) }, add = TRUE)

nhanes_check_inputs()

maha_configure_steps(
	c("data", "fig3", "fig4", "figS1"),
	dependencies = list(
		fig3 = "data",
		fig4 = c("data", "fig3"),
		figS1 = "data"
	)
)

maha_run_step("data", {
dat <- pmap_dfr(cycle_info, read_cycle)

if (nrow(dat) == 0) stop("No NHANES data were read. Check raw_dir: ", raw_dir)

if (use_mortality) {
	mort <- read_mortality()
	if (!is.null(mort)) dat <- left_join(dat, mort, by = "SEQN")
}

dat <- dat %>%
	mutate(
		nonhdl_mgdl = tc_mgdl - hdl_mgdl,
		egfr = egfr_2021(creat_mgdl, age, female),
		t2dm = as.integer(diabetes_q == 1 | hba1c >= 6.5 | glu_mgdl >= 126),
		cad = as.integer(cad_q == 1 | mi_q == 1),
		mi = as.integer(mi_q == 1),
		stroke = as.integer(stroke_q == 1),
		heart_failure = as.integer(hf_q == 1),
		ckd = as.integer(ckd_q == 1 | egfr < 60 | uacr_mgg >= 30),
		hypertension = as.integer(htn_q == 1 | sbp >= 130 | dbp >= 80),
		copd = as.integer(copd_q == 1),
		asthma = as.integer(asthma_q == 1),
		depression = as.integer(phq9 >= 10),
		obesity = as.integer(bmi >= 30),
		mort_eligible = if ("eligstat" %in% names(.)) as.integer(eligstat == 1) else NA_integer_,
		death = if ("mortstat" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1), NA_integer_) else NA_integer_,
		death_time_y = if ("permth_exm" %in% names(.)) as.numeric(permth_exm) / 12 else NA_real_,
		heart_death = if ("ucod_leading" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1 & ucod_leading == 1), NA_integer_) else NA_integer_,
		cvd_death = if ("ucod_leading" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1 & ucod_leading %in% c(1, 5)), NA_integer_) else NA_integer_,
		cancer_death = if ("ucod_leading" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1 & ucod_leading == 2), NA_integer_) else NA_integer_,
		diabetes_death = if ("ucod_leading" %in% names(.)) ifelse(eligstat == 1, as.integer(mortstat == 1 & ucod_leading == 7), NA_integer_) else NA_integer_
	)

for (v in c("cad_q", "copd_q")) dat[[v]][is.infinite(dat[[v]])] <- NA_real_

food_vars <- c("fruit","berry","vegetable","green_leafy","allium","legumes","nuts","dairy","lowfat_dairy",
	"whole_grain","refined_grain","fish","poultry","red_processed_meat","sweets_pastries","ssb","fried_fast","coffee_tea")
nutrient_vars <- c("kcal", "protein_g", "carb_g", "sugar_g", "fiber_g", "fat_g", "sfat_g", "mufa_g", "pufa_g", "chol_mg", "sodium_mg", "alcohol_g")
for (v in food_vars) if (!v %in% names(dat)) dat[[v]] <- NA_real_
for (v in nutrient_vars) if (!v %in% names(dat)) dat[[v]] <- NA_real_

step_header("Diet score construction: observed repeated-recall mean")
dat_observed <- construct_diet_scores(dat, score_source = "observed_repeated_recall_mean")

dat1_obs <- prep_analysis_dataset(dat_observed, main_cycles, require_mortality = FALSE)
dat_mort_obs <- prep_analysis_dataset(dat_observed, mort_cycles, require_mortality = TRUE)

dat1 <- dat1_obs
dat_mort1 <- dat_mort_obs

diet.sum.cols <- paste0("diet.", names(diet.lst), ".sum")
diet.sum.cols <- diet.sum.cols[diet.sum.cols %in% names(dat1_obs)]
diet.inc.pts <- paste0("diet.", diet.inc, ".pts")
diet.inc.s100 <- paste0("diet.", diet.inc, ".s100")
diet.inc.q5 <- paste0("diet.", diet.inc, ".q5")
diet.inc.3c <- paste0("diet.", diet.inc, ".3c")

message("N after QC, observed main analysis: ", nrow(dat1_obs))
message("N after QC, mortality analysis: ", nrow(dat_mort_obs))
message("Main cycles retained: ", paste(sort(unique(dat1_obs$cycle)), collapse = ", "))
message("Mortality cycles retained: ", paste(sort(unique(dat_mort_obs$cycle)), collapse = ", "))

covs.base <- c("age", "female", "race", "edu", "pir", "cycle", "kcal")
covs.add <- c("bmi", "smoke_ever", "smoke_now")
covs <- c(covs.base, covs.add)

covs.mort <- c("age", "female", "race", "edu", "pir", "marital", "health_insurance", "cycle", "kcal", "bmi", "smoke_ever", "smoke_now", "hypertension", "cad", "stroke", "heart_failure", "t2dm", "ckd")
Y.cross <- c("cad", "stroke", "heart_failure", "t2dm", "ckd")
Y.mort.plot <- c("death", "heart_death", "cvd_death", "cancer_death", "diabetes_death")

nhanes_validation <- nhanes_validate_maha(dat_mort_obs, covs.mort)
nhanes_construct <- nhanes_construct_profile(dat1_obs)
})

maha_run_step("fig3", {
step_header("Figure 3: NHANES all-cause mortality")

Y.mort.main <- "death"
fig3.all_map <- c(
	maha = "MAHA", dash = "DASH", mind = "MIND", medi = "MEDI",
	maha_bal = "MAHA-balanced", maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy", maha_noprotein = "MAHA-no protein"
)
fig3.discord_map <- c(
	dash = "DASH", medi = "MEDI", mind = "MIND",
	maha_bal = "MAHA-balanced", maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy", maha_noprotein = "MAHA-no protein"
)
fig3.shared_cols <- c(
	"MAHA" = "#F26D60", "DASH" = "#16B9C0", "MIND" = "#B97AF7", "MEDI" = "#7CAE00",
	"MAHA-balanced" = "#E76F51", "MAHA-strict" = "#C77DFF",
	"MAHA-no dairy" = "#2A9D8F", "MAHA-no protein" = "#577590"
)
fig3.legend_levels <- unname(fig3.all_map)
fig3.discord_levels <- unname(fig3.discord_map)

Fig3.all <- run_fig3_assoc_data(dat_mort_obs, Y.mort.main, covs.mort, names(fig3.all_map), fig3.all_map, family = "cox")
Fig3.disc <- run_fig3_discord_data(dat_mort_obs, Y.mort.main, covs.mort, names(fig3.discord_map), fig3.discord_map, ref_nm = "maha", ref_lab = "MAHA", family = "cox", min_n_pair = 20)

Fig3.ab <- list(
	A = plot_fig3_line_panel(
		Fig3.all$assoc, Y.mort.main, "Diet", fig3.legend_levels, fig3.shared_cols,
		"a. Dietary patterns", "Hazard Ratio (per SD)", xlim = c(0.80, 1.10),
		show_y = TRUE, show_legend = TRUE, max_shift = 0.58, y_expand = 0.08, show_hr_text = TRUE
	)
)

Fig3.cd <- make_mortality_bar_panels(dat_mort_obs, outcome = "death", t0 = 10)
legend_ab <- cowplot::get_legend(Fig3.ab$A + theme(legend.position = "bottom"))
Fig3.top <- cowplot::plot_grid(
	Fig3.ab$A + theme(legend.position = "none"),
	legend_ab,
	ncol = 1, rel_heights = c(1, 0.16)
)
Fig3.right <- cowplot::plot_grid(Fig3.cd$C, Fig3.cd$D, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr")
Fig3 <- cowplot::plot_grid(
	Fig3.top, NULL, Fig3.right,
	nrow = 1, rel_widths = c(0.69, 0.08, 0.60),
	align = "h", axis = "tb"
)

print(Fig3.all$diagnostics)
print(Fig3.all$assoc)
print(Fig3.disc$diagnostics)
print(Fig3.disc$disc)
})

maha_run_step("fig4", {
step_header("Figure 4: NHANES mortality validation and precision checks")
pacman::p_load(tidyverse, patchwork, writexl, scales, cowplot)

ap_margin_main <- 1.05
ap_prior_main <- 1.10

ap_equiv_main <- ap_from_assoc(Fig3.all$assoc, "NHANES mortality")
ap_head2head <- purrr::map2_dfr(c("dash", "mind", "medi"), c("DASH", "MIND", "MEDI"), function(nm, lab) {
	ap_fit_head2head_nhanes(dat_mort_obs, Y = "death", trad_nm = nm, trad_lab = lab, covs.in = covs.mort, maha_nm = "maha", maha_lab = "MAHA")
})

ap_equiv_clean <- ap_clean_equiv(ap_equiv_main, "NHANES")
ap_head_clean <- ap_clean_head2head(ap_head2head, "NHANES")
ap_maha_levels <- c("MAHA", "MAHA-balanced", "MAHA-strict", "MAHA-no dairy", "MAHA-no protein")
ap_maha_cols <- fig3.shared_cols[ap_maha_levels]

p4a <- Fig3.ab$A +
	labs(title = "a. Dietary patterns", y = "\n\n\nAll-cause mortality") +
	guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE)) +
	theme(
		legend.position = "top",
		legend.box = "vertical",
		legend.margin = margin(0, 0, 2, 0),
		legend.box.margin = margin(0, 0, 2, 0),
		plot.title = element_text(margin = margin(b = 6)),
		axis.title.y = element_text(face = "bold", angle = 90, vjust = 1.8, margin = margin(r = -28)),
		plot.margin = margin(18, 18, 5.5, -6)
	)
p4b <- Fig3.cd$C +
	labs(title = "b. DASH gradient") +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.title.y = element_text(face = "bold", margin = margin(r = -5))
	)
p4c <- Fig3.cd$D +
	labs(title = "c. MAHA gradient") +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.title.y = element_text(face = "bold", margin = margin(r = -4)),
		plot.margin = margin(5.5, 5.5, 5.5, 0)
	)
p4d <- ap_plot_maha_equiv(ap_equiv_clean, "d. MAHA equivalence")
p4e <- ap_plot_bf(ap_equiv_clean, "e. Bayes evidence for near-zero")
p4f <- ap_plot_power(ap_equiv_clean, "f. MAHA detectable effect")
p4g <- ap_plot_head2head(ap_head_clean, "g. Direct attenuation versus established scores")

Fig4_top <- p4a + patchwork::plot_spacer() + p4b + patchwork::plot_spacer() + p4c + plot_layout(widths = c(1.08, 0.07, 1.00, 0.07, 1.00))
Fig4_bottom <- (p4d + p4e) / patchwork::plot_spacer() / (p4f + p4g) + plot_layout(heights = c(1.00, 0.08, 1.00))
Fig4_nhanes <- Fig4_top / patchwork::plot_spacer() / Fig4_bottom + plot_layout(heights = c(0.86, 0.06, 1.44))

save_plot(Fig4_nhanes, "Fig4.mortality_validation.png", width = 17.4, height = 13.0, dpi = 320)
write_xlsx(
	list(
		Fig4A_assoc = Fig3.all$assoc,
		Fig4A_diag = Fig3.all$diagnostics,
		Fig4_disc = Fig3.disc$disc,
		Fig4_disc_diag = Fig3.disc$diagnostics,
		Fig4B_risk_DASH = Fig3.cd$risk10_DASH_within_MAHA,
		Fig4C_risk_MAHA = Fig3.cd$risk10_MAHA_within_DASH,
		Fig4D_equiv = ap_equiv_clean %>% filter(Diet %in% ap_maha_levels) %>% dplyr::select(Dataset, Outcome, Diet, N_total, N_event, HR_90CI, TOST_p, TOST_P_label, Equivalence_5pct, margin_hr, CI90_low_HR, CI90_high_HR),
		Fig4E_BF = ap_equiv_clean %>% filter(Diet %in% ap_maha_levels) %>% dplyr::select(Dataset, Outcome, Diet, N_total, N_event, BF01, BF01_label, BF_summary, prior_95_hr),
		Fig4F_power = ap_equiv_clean %>% filter(Diet %in% ap_maha_levels) %>% dplyr::select(Dataset, Outcome, Diet, N_total, N_event, MDE80_HR_label, Powered_5pct, MDE_HR_lower, MDE_HR_upper, margin_hr),
		Fig4G_head2head = ap_head_clean,
		Fig4_raw_equiv = ap_equiv_main,
		Fig4_raw_h2h = ap_head2head,
		curves_DASH_MAHA = Fig3.cd$curves_DASH_within_MAHA,
		curves_MAHA_DASH = Fig3.cd$curves_MAHA_within_DASH
	),
	"Fig4.mortality_validation.out.xlsx"
)

print(ap_equiv_clean %>% filter(Diet %in% ap_maha_levels))
print(ap_head_clean)
})

maha_run_step("figS1", {
step_header("Figure S3: NHANES diet score distributions and correlation")
pacman::p_load(tidyverse, patchwork, cowplot)
dat1 <- dat1_obs

figS1_score_map <- c(
	maha = "MAHA", dash = "DASH", mind = "MIND", medi = "MEDI",
	maha_bal = "MAHA-balanced", maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy", maha_noprotein = "MAHA-no protein"
)
figS1_cols <- c(
	"MAHA" = "#F26D60", "DASH" = "#16B9C0", "MIND" = "#B97AF7", "MEDI" = "#7CAE00",
	"MAHA-balanced" = "#E76F51", "MAHA-strict" = "#C77DFF",
	"MAHA-no dairy" = "#2A9D8F", "MAHA-no protein" = "#577590"
)
figS1_dist_axis_text_size <- 8.2
figS1_top_height <- 0.76 * 1.30
figS1_gap_height <- 0.10
figS1_bottom_height <- 1.24
figS1_plot_height <- 9.6 * (figS1_top_height + figS1_gap_height + figS1_bottom_height) / (0.76 + 1.24)

for (nm in names(figS1_score_map)) {
	sumv <- paste0("diet.", nm, ".sum"); s100v <- paste0("diet.", nm, ".s100")
	if (!s100v %in% names(dat1) && sumv %in% names(dat1)) dat1[[s100v]] <- score_0_100(dat1[[sumv]])
}

figS1_vars <- paste0("diet.", names(figS1_score_map), ".s100")
missing_figS1 <- setdiff(figS1_vars, names(dat1))
if (length(missing_figS1)) stop("Missing FigS1 score variables: ", paste(missing_figS1, collapse = ", "))
figS1_var_labs <- setNames(unname(figS1_score_map), figS1_vars)

plot_dat <- dat1 %>%
	dplyr::select(all_of(figS1_vars)) %>%
	pivot_longer(everything(), names_to = "var", values_to = "score") %>%
	mutate(score = as.numeric(score), Diet_score = factor(figS1_var_labs[var], levels = rev(unname(figS1_score_map)))) %>%
	filter(is.finite(score), !is.na(Diet_score))

plot_sum <- plot_dat %>%
	group_by(Diet_score) %>%
	summarise(
		N = n(), Mean = mean(score), SD = sd(score), Median = median(score),
		P25 = quantile(score, 0.25) %>% as.numeric(), P75 = quantile(score, 0.75) %>% as.numeric(),
		Min = min(score), Max = max(score), .groups = "drop"
	) %>%
	mutate(across(c(Mean, SD, Median, P25, P75, Min, Max), ~ round(.x, 1))) %>%
	as.data.frame()
density_dat <- plot_dat %>% mutate(Diet_score = as.character(Diet_score)) %>%
	group_by(Diet_score) %>% group_modify(~ {
		dn <- density(.x$score, from = 0, to = 100, n = 256, na.rm = TRUE)
		tibble(score = dn$x, density = dn$y)
	}) %>% ungroup()

p_dist <- ggplot(plot_dat, aes(x = score, y = Diet_score, fill = Diet_score)) +
	geom_violin(width = 0.82, scale = "width", trim = TRUE, alpha = 0.72, color = NA) +
	geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.92, color = "grey25") +
	scale_fill_manual(values = figS1_cols, guide = "none") +
	scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0.01, 0.01))) +
	labs(x = "Harmonized dietary score (0-100)", y = NULL, title = "b. Score distribution (NHANES)") +
	theme_classic(base_size = 10) +
	theme(
		legend.position = "none",
		plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
		axis.text.y = element_text(face = "bold", size = figS1_dist_axis_text_size, color = "black"),
		axis.text.x = element_text(face = "bold", size = figS1_dist_axis_text_size),
		axis.title.x = element_text(face = "bold", size = figS1_dist_axis_text_size),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		plot.margin = margin(5.5, 12, 5.5, 5.5)
	)

score_mat <- dat1 %>% dplyr::select(all_of(figS1_vars)) %>% mutate(across(everything(), as.numeric))
cor_mat <- cor(score_mat, use = "pairwise.complete.obs", method = "spearman")
dimnames(cor_mat) <- rep(list(unname(figS1_score_map)), 2)
cor_dat <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE) %>%
	mutate(Var1 = factor(Var1, levels = rev(unname(figS1_score_map))), Var2 = factor(Var2, levels = unname(figS1_score_map)), lab = sprintf("%.2f", Freq))

p_cor <- ggplot(cor_dat, aes(Var2, Var1, fill = Freq)) +
	geom_tile(color = "white", linewidth = 0.35) +
	geom_text(aes(label = lab), size = 2.55, fontface = "bold") +
	scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0.5, limits = c(0, 1), breaks = c(0, 0.5, 1), oob = scales::squish, name = NULL) +
	coord_fixed() +
	labs(x = NULL, y = NULL, title = "d. Spearman correlation (NHANES)") +
	theme_minimal(base_size = 10) +
	theme(
		plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
		axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, face = "bold", size = 8.2),
		axis.text.y = element_text(face = "bold", size = figS1_dist_axis_text_size, color = "black"),
		panel.grid = element_blank(),
		legend.position = "right",
		legend.key.height = unit(0.55, "cm"),
		plot.margin = margin(5.5, 5.5, 5.5, 5.5)
	)

FigS1_right_nhanes <- p_dist / plot_spacer() / p_cor +
	plot_layout(heights = c(figS1_top_height, figS1_gap_height, figS1_bottom_height))
save_plot(FigS1_right_nhanes, "FigS3.score_distribution_concordance.png", width = 6.9, height = figS1_plot_height, dpi = 500, bg = "white")

write_xlsx(
	list(
		score_distribution = plot_sum,
		density_curve = density_dat,
		spearman_matrix = as.data.frame(cor_mat, check.names = FALSE)
	),
	"FigS3.score_distribution_concordance.out.xlsx"
)
print(plot_sum)
print(round(cor_mat, 2))
})
