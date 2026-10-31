pacman::p_load(tidyverse, data.table, survival)

dir0 <- Sys.getenv("DIR0", unset = ifelse(Sys.info()[["sysname"]] == "Windows", "D:/", "/mnt/d"))
indir <- Sys.getenv("PHEDIR", unset = paste0(dir0, "/data/ukb"))
outdir <- Sys.getenv("UKB_OUT", unset = paste0(dir0, "/analysis/ukb/phe"))
script_dir <- Sys.getenv("SCRIPT_DIR", unset = file.path(dir0, "scripts"))
gen_dir <- Sys.getenv("GEN_DIR", unset = if (Sys.info()[["sysname"]] == "Windows") "F:/gen/ukb/37" else "/mnt/f/gen/ukb/37")
pgs_dir <- Sys.getenv("PGS_DIR", unset = if (Sys.info()[["sysname"]] == "Windows") "F:/gen/ukb/pgs" else "/mnt/f/gen/ukb/pgs")
map_tile_cache <- normalizePath(
	file.path(dir0, "data", "ukb", "map", "maha", "rosm.cache"),
	winslash = "/", mustWork = FALSE
)
invisible(lapply(c("phenotype.R", "association.R", "plotting.R", "prediction.R"), \(f) source(file.path(script_dir, "0f", f))))
source(file.path(script_dir, "ukb", "f", "date_utils.R"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd2(outdir)

step_order <- c(
	"data_qc", "vip_phe", "ses", "move", "death", "srd", "icd", "fod_ref", "gp_prep", "gp",
	"biom", "hla_pca", "pgs", "merge",
	"qc", "drug", "phe4gwas"
)

start_step <- Sys.getenv("START_STEP", unset = "data_qc")
end_step <- Sys.getenv("END_STEP", unset = step_order[length(step_order)])
run_steps <- Sys.getenv("RUN_STEPS", unset = "")

show_steps <- function() cat(paste0(seq_along(step_order), ". ", step_order, collapse = "\n"), "\n")
if (start_step %in% c("help", "--help", "-h")) { show_steps(); quit(save = "no") }
if (start_step %in% c("all", "ALL", "")) start_step <- step_order[1]
if (end_step %in% c("all", "ALL", "")) end_step <- step_order[length(step_order)]
if (!start_step %in% step_order && run_steps == "") stop("Unknown START_STEP: ", start_step, "\nAvailable steps: ", paste(step_order, collapse = ", "))
if (!end_step %in% step_order && run_steps == "") stop("Unknown END_STEP: ", end_step, "\nAvailable steps: ", paste(step_order, collapse = ", "))

run_vec <- trimws(unlist(strsplit(run_steps, ",|;|\\s+")))
run_vec <- run_vec[nzchar(run_vec)]
bad_run_vec <- setdiff(run_vec, step_order)
if (length(bad_run_vec)) stop("Unknown RUN_STEPS entry: ", paste(bad_run_vec, collapse = ", "), "\nAvailable steps: ", paste(step_order, collapse = ", "))
run_step <- function(step) {
	if (length(run_vec)) return(step %in% run_vec)
	i <- match(step, step_order); i0 <- match(start_step, step_order); i1 <- match(end_step, step_order)
	!is.na(i) && !is.na(i0) && !is.na(i1) && i >= i0 && i <= i1
}
run_if <- function(step, expr) {
	if (run_step(step)) {
		cat("\n========== RUN STEP:", step, "==========\n")
		eval(substitute(expr), envir = parent.frame())
		cat("========== DONE STEP:", step, "==========\n")
	} else {
		cat("Skip step:", step, "\n")
	}
}

asDate2 <- function(x) {
	if (inherits(x, "Date")) return(x)
	if (inherits(x, "IDate")) return(as.Date(x))
	if (is.character(x)) return(as.Date(x))
	y <- suppressWarnings(as.numeric(unclass(x)))
	y[!is.finite(y)] <- NA_real_
	as.Date(y, origin = "1970-01-01")
}
pminDate2 <- function(...) {
	lst <- list(...)
	if (length(lst) == 1 && is.data.frame(lst[[1]])) lst <- as.list(lst[[1]])
	lst <- Filter(Negate(is.null), lst)
	if (!length(lst)) return(as.Date(NA))
	lst <- lapply(lst, asDate2)
	z <- do.call(pmin, c(lst, list(na.rm = TRUE)))
	z <- asDate2(z)
	z[!is.finite(as.numeric(z))] <- as.Date(NA)
	z
}
read_rds <- function(name) readRDS(paste0(indir, "/Rdata/", name, ".rds"))
get_existing_date_cols <- function(dat, trait, prefixes = c("", "fod_", "date_")) {
	trait2 <- gsub("^cvd_", "", trait)
	cand <- unique(unlist(lapply(prefixes, function(p) c(paste0(p, trait), paste0(p, trait2)))))
	intersect(cand, names(dat))
}
get_all_trait_date_cols <- function(dat, traits, prefixes) {
	unique(unlist(lapply(traits, function(tr) get_existing_date_cols(dat, tr, prefixes))))
}
trait_prevalent <- function(dat, trait, date_attend_col = "date_attend", prefixes = c("fod_icd10_", "fod_icd9_", "fod_srd_", "fod_gp_", "fod_ref_", "date_", "")) {
	cols <- get_existing_date_cols(dat, trait, prefixes)
	if (!length(cols)) return(rep(FALSE, nrow(dat)))
	d <- pminDate2(dat[, cols, drop = FALSE])
	!is.na(d) & !is.na(dat[[date_attend_col]]) & d <= asDate2(dat[[date_attend_col]])
}
get_date_cols <- function(dat, cols) {
	cols <- intersect(cols, names(dat))
	if (!length(cols)) return(data.frame(eid = as.character(dat$eid)))
	dat %>% dplyr::select(eid, all_of(cols)) %>% mutate(eid = as.character(eid))
}
coalesce_date_from_join <- function(dat, base) {
	cand <- intersect(c(base, paste0(base, ".phe"), paste0(base, ".death"), paste0(base, ".x"), paste0(base, ".y")), names(dat))
	if (!length(cand)) return(rep(as.Date(NA), nrow(dat)))
	do.call(coalesce, lapply(cand, function(cc) asDate2(dat[[cc]])))
}


phe_num <- function(x) {
	z <- suppressWarnings(as.numeric(as.character(x)))
	z[!is.finite(z)] <- NA_real_
	z
}

phe_date <- function(x) {
	if (inherits(x, "Date")) return(x)
	if (inherits(x, "IDate")) return(as.Date(x))
	if (inherits(x, "POSIXt")) return(as.Date(x))
	if (is.character(x)) return(suppressWarnings(as.Date(x)))
	z <- suppressWarnings(as.numeric(unclass(x)))
	z[!is.finite(z)] <- NA_real_
	as.Date(z, origin = "1970-01-01")
}

# Compute baseline UACR from numeric urine albumin and creatinine values.
derive_ukb_uacr <- function(dat,
			alb_col = "bb_UALB",
			cr_col = "bb_UCRE") {
	stopifnot("eid" %in% names(dat))
	n <- nrow(dat)
	get_col <- function(nm, default = NA) {
		if (nm %in% names(dat)) dat[[nm]] else rep(default, n)
	}

	alb <- phe_num(get_col(alb_col))
	cr <- phe_num(get_col(cr_col))
	assessed <- is.finite(alb) & alb >= 0 & is.finite(cr) & cr > 0
	uacr <- rep(NA_real_, n)
	uacr[assessed] <- alb[assessed] * 8840 / cr[assessed]

	category <- dplyr::case_when(
		is.finite(uacr) & uacr < 30 ~ "A1",
		is.finite(uacr) & uacr < 300 ~ "A2",
		is.finite(uacr) ~ "A3",
		TRUE ~ NA_character_
	)

	data.frame(
		eid = as.character(dat$eid),
		bb_UACR = uacr,
		bb_UACR_exact = uacr,
		bb_UACR_lo = uacr,
		bb_UACR_hi = uacr,
		bb_UACR_category = category,
		bb_UACR_assessed = assessed,
		bb_UALB_below_lod = FALSE,
		bb_UCRE_below_lod = FALSE,
		bb_UACR_interval_uncertain = FALSE,
		bb_UACR_prevent_sens = uacr,
		stringsAsFactors = FALSE
	)
}


run_if("data_qc", {


# 🚩 data QC
fn <- fread(paste0(indir, "/common/vip.lst"), sep = "\t", select = 1:3, fill = Inf, header = FALSE) %>% as.data.frame()
fn$V1[duplicated(fn$V1)]; fn$V2[duplicated(fn$V2)]
phe0 <- fread(paste0(indir, "/rap/vip.tab.gz"), sep = "\t", header = TRUE) %>% as.data.frame()
names(phe0) <- names(phe0) %>% gsub("_i0$|_i0_a0$", "", .)
setdiff(fn$V1, names(phe0)); setdiff(names(phe0), fn$V1)
names(phe0) <- sapply(names(phe0), replace_code, sep = "_", mapping_df = fn) %>% as.character()

phe0 <- phe0 %>% mutate(across(grep("date", names(phe0), value = TRUE), as.Date), edu = ifelse(edu ==  - 7, 7, edu))
#special negatvie values
phe0 <- neg2NA(phe0, c( - 1, - 2, - 3, - 7, - 10, - 121, - 818)); table(phe0$social) # 
#special dates
date_4na <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))
date_cols <- grep("date", names(phe0), ignore.case = TRUE, value = TRUE)
purrr::walk(date_cols, \(cc) print(summary(phe0[[cc]])))
phe0 <- phe0 %>% mutate(across(all_of(date_cols), \(x) { x[x %in% date_4na] <- NA; x }))
#combo variables
phe0 %>% dplyr::select(where( ~ any(stringr::str_detect(as.character(.x), stringr::fixed("|")), na.rm = TRUE))) %>% names()
saveRDS(phe0, paste0(indir, "/Rdata/phe0.rds"))
})


run_if("vip_phe", {


# 🚩 VIP Phe
if (!exists("phe0")) phe0 <- read_rds("phe0")
phe <- phe0 %>% mutate( # age_reception [p21003], not age_recruit [p21022]
	center = paste0("c", center),
	abo = recode_factor(blood_type, "AA" = "A", "BB" = "B", "OO" = "O", "AO" = "A", "BO" = "B"),
	ethnic.c = ifelse(ethnicity == 2002, 12, ifelse(ethnicity == 2003, 13, ifelse(ethnicity %in% c(2004, 3004), 6, ifelse(grepl("^1", ethnicity), 1, ifelse(grepl("^2", ethnicity), 2, ifelse(grepl("^3", ethnicity), 3, ifelse(grepl("^4", ethnicity), 4, ethnicity))))))),
	ethnic.c = factor(ethnic.c, levels = c(1:6, 12:13), labels = c("White", "Mixed", "Asian", "Black", "Chinese", "Other", "White-Black", "White-Asian")),
	tdi.i = scale(inormal( - tdi)),
	age.c = cut(age, breaks = seq(35, 75, 5)), age2 = age * age,
	birth_date = as.Date(paste0(birth_year, "-", birth_month, "-15")),
	sex.b = ifelse(sex == 0, "female", ifelse(sex == 1, "male", NA)),
	sex.f = factor(sex, c(0, 1), c("Female", "Male")),
	sex_first.c = factor(ifelse(sex_first < 12, NA, ifelse(sex_first < 16, "early", ifelse(sex_first < 22, "average", "late"))), levels = c("early", "average", "late")),
	sex_par.c = factor(ifelse(sex_par <= 1, "low", ifelse(sex_par <= 10, "average", "high")), levels = c("low", "average", "high")),
	whr = waist / hip,
	leg = height - height_sit, lhr = leg / height, lhr = residuals(lm(lhr ~ height, na.action = na.exclude)),
	bald.c = factor(bald, levels = 1:4, labels = c("ctrl", "bald2", "bald3", "bald4")),
	bald12 = ifelse(bald == 1, 0, ifelse(bald == 2, 1, NA)), bald13 = ifelse(bald == 1, 0, ifelse(bald == 3, 1, NA)), bald14 = ifelse(bald == 1, 0, ifelse(bald == 4, 1, NA)),
	# SES 
	edu.b = factor(ifelse(is.na(edu), NA, ifelse(grepl("1", edu), 1, 0))),
		edu.years = unname(c("1" = 20, "2" = 13, "3" = 10, "4" = 10, "5" = 19, "6" = 15, "7" = 7)[as.character(split_multiple(edu, fun = "max"))]),
		edu.sco = if_else(is.na(edu.years), NA_integer_, as.integer(findInterval(edu.years, vec = c(0, 3, 6, 9, 12, 15, 17, 19, 21, 23)))),
		edu.c = factor(case_when(edu == 1 ~ "high", edu %in% c(2, 5, 6) ~ "middle", edu %in% c(3, 4, 7) ~ "low", TRUE ~ NA_character_), levels = c("low", "middle", "high"), ordered = TRUE),
		job.b = ifelse(is.na(job), NA, ifelse(grepl("[1267]", job), 1, 0)),
		job.c = factor(case_when(job == 1 ~ "employed", job %in% c(2, 7, 6) ~ "inactive", job %in% c(3, 4, 5, 8) ~ "unemployed", TRUE ~ NA_character_), levels = c("unemployed", "inactive", "employed"), ordered = TRUE),
		inc.c = factor(case_when(inc %in% c(1, 2) ~ "low", inc == 3 ~ "middle", inc %in% c(4, 5) ~ "high", TRUE ~ NA_character_), levels = c("low", "middle", "high"), ordered = TRUE),
	# happy  
	mental_bad = (depress_freq - 1) + (disinterest_freq - 1) + tense_freq,
	mental.pts = case_match(mental_bad, 0 ~ 100, c(1, 2) ~ 75, c(3, 4) ~ 50, c(5, 6) ~ 25, 7 ~ 0, .default = NA),
	happy = ifelse(happy_general %in% 1:2, 1, ifelse(happy_general %in% 3:6, 0, NA)),
	haphealth = ifelse(happy_health %in% 1:2, 1, ifelse(happy_health %in% 3:6, 0, NA)),
	haplife = ifelse(happy_life %in% c(4, 5), 1, ifelse(happy_life %in% 1:3, 0, NA)),
	# social 
	social.n = ifelse(is.na(social), NA_real_, lengths(strsplit(social, "\\|"))),
	social.n = case_when(is.na(social.n) ~ NA_real_, social.n == 0 ~ 0, social.n == 1 ~ 70, social.n == 2 ~ 85, social.n >= 3 ~ 100, TRUE ~ NA_real_),
	ffvisit.n = case_when(ffvisit %in% 1:3 ~ 100, ffvisit == 4 ~ 80, ffvisit == 5 ~ 60, ffvisit == 6 ~ 20, ffvisit == 7 ~ 0, TRUE ~ NA_real_),
	confide.n = case_when(confide == 5 ~ 100, confide == 4 ~ 90, confide == 3 ~ 80, confide == 2 ~ 60, confide == 1 ~ 40, confide == 0 ~ 0, TRUE ~ NA_real_),
	lone.n = case_when(lone == 0 ~ 100, lone == 1 ~ 0, TRUE ~ NA_real_),
	social.raw = rowMeans(cbind(social.n, ffvisit.n, confide.n, lone.n), na.rm = TRUE),
	social.pts = cut(social.raw, breaks = unique(quantile(social.raw, probs = c(0, .2, .4, .6, .8, 1), na.rm = TRUE)), labels = c(0, 25, 50, 75, 100)[seq_len(length(unique(quantile(social.raw, probs = c(0, .2, .4, .6, .8, 1), na.rm = TRUE))) - 1)], include.lowest = TRUE) %>% as.character() %>% as.numeric(),

)
saveRDS(phe, paste0(indir, "/Rdata/phe.rds"))
})


run_if("ses", {


# 🚩 SES (PMID: 39715877)
if (!exists("phe")) phe <- read_rds("phe")
# 22601: Job coding history; 22617 uses the first four of 22601; 22620: job shift
library(ukbjobs); ls("package:ukbjobs")
job <- phe %>% dplyr::select(eid, sex, job_soc_a0) %>%
	add_occ_status(soc2000_var = "job_soc_a0", sex_var = "sex") %>%
	add_traits(soc2000_var = "job_soc_a0") %>%
	add_isco88(soc2000_var = "job_soc_a0") %>%
	dplyr::select(!matches("sex|^job_soc_a0$"))
saveRDS(job, paste0(indir, "/Rdata/ses.rds"))
})


run_if("move", {


# 🚩 move
if (!exists("phe")) phe <- read_rds("phe")
pacman::p_load(dplyr, tidyr, sf, RANN, ggplot2, patchwork, rnaturalearth, ggspatial, scales)
# phe <- readRDS(paste0(indir, "/Rdata/phe.rds"))
dat <- phe %>% dplyr::select(eid, tdi, center, birth_east, birth_north, home_east, home_north) %>%
	mutate(center0 = sub("\\.0$", "", trimws(as.character(center))),
		center = case_when(
			center0 %in% names(center.lst) ~ center0,
			paste0("c", center0) %in% names(center.lst) ~ paste0("c", center0),
			tolower(center0) %in% tolower(unname(center.lst)) ~ names(center.lst)[match(tolower(center0), tolower(unname(center.lst)))],
			TRUE ~ NA_character_
		)
	) %>% dplyr::select( - center0) %>%
	drop_na(home_east, home_north, birth_east, birth_north, center)
gb <- st_as_sf(dat, coords = c("birth_east", "birth_north"), crs = 27700, remove = FALSE)
gh <- st_as_sf(dat, coords = c("home_east", "home_north"), crs = 27700, remove = FALSE)
cb <- st_coordinates(st_transform(st_geometry(gb), 4326))
ch <- st_coordinates(st_transform(st_geometry(gh), 4326))

# tdi at birth 
ref <- phe %>% dplyr::select(home_east, home_north, tdi) %>%
	drop_na(home_east, home_north, tdi) %>% group_by(home_east, home_north) %>% summarise(tdi.ref = median(tdi), .groups = "drop")
nn <- RANN::nn2(data = as.matrix(ref[, c("home_east", "home_north")]), query = as.matrix(dat[, c("birth_east", "birth_north")]), k = 1)
dat1 <- dat %>% mutate(
		lon.birth = cb[, 1], lat.birth = cb[, 2],
		lon.home = ch[, 1], lat.home = ch[, 2],
		move = as.numeric(st_distance(st_geometry(gb), st_geometry(gh), by_element = TRUE)) / 1000,
		tdi.birth = ref$tdi.ref[nn$nn.idx[, 1]],
		tdi.birth.nn = nn$nn.dists[, 1]
	)

# draw maps 
p1 <- plot_map(dat1, separate = FALSE, transparency = "tdi", center.lst = center.lst, basemap = TRUE, basemap_type = "osm", tile_zoom = 6, tile_cachedir = map_tile_cache)
	ggsave("Fig1_UK_all_centers_tdi.png", p1, width = 10, height = 12, dpi = 300, bg = "white")
p2 <- plot_map(dat1, separate = FALSE, transparency = "move", center.lst = center.lst, basemap = TRUE, basemap_type = "osm", tile_zoom = 6, tile_cachedir = map_tile_cache)
	ggsave("Fig2_UK_all_centers_move.png", p2, width = 10, height = 12, dpi = 300, bg = "white")
p3 <- plot_map(dat1, centers = c("Croydon", "Hounslow", "Bristol", "Edinburgh"), separate = TRUE, width = 50, nrow = 2, transparency = "tdi", center.lst = center.lst, basemap = TRUE, basemap_type = "osm", tile_zoom = 11, tile_cachedir = map_tile_cache)
	ggsave("Fig3_4centers_tdi.png", p3, width = 10, height = 10, dpi = 300, bg = "white")

dat1 <- dat1 %>% dplyr::select(eid, lon.birth, lat.birth, lon.home, lat.home, move, tdi.birth, tdi.birth.nn)
saveRDS(dat1, paste0(indir, "/Rdata/move.rds"))
})


run_if("death", {


# 🚩 Death
fn <- fread(paste0(indir, "/common/death.lst"), sep = "\t", header = FALSE, select = 1:2) %>% as.data.frame()
dat0 <- fread(paste0(indir, "/rap/death.tab.gz"), sep = "\t", header = TRUE) %>% as.data.frame()
death_dates <- ukb_death_dates(dat0)
concat <- function(x) {
	x <- as.character(x); x <- x[!is.na(x) & x != ""]
	if (!length(x)) NA_character_ else paste(x, collapse = " ")
}
death_concat <- function(field) {
	cc <- grep(paste0("^", field, "(_|$)"), names(dat0), value = TRUE)
	if (!length(cc)) return(rep(NA_character_, nrow(dat0)))
	apply(dat0[, cc, drop = FALSE], 1, concat)
}
death_first <- function(field) {
	cc <- grep(paste0("^", field, "(_|$)"), names(dat0), value = TRUE)
	if (!length(cc)) return(rep(NA, nrow(dat0)))
	do.call(dplyr::coalesce, as.list(dat0[, cc, drop = FALSE]))
}

dat <- dat0 %>% mutate(
	p40000 = death_dates$date_death,
	p40001 = death_concat("p40001"),
	p40002 = death_concat("p40002"),
	p40007 = death_first("p40007"),
	p40010 = as.character(death_first("p40010"))
) %>% dplyr::select(eid, any_of(fn$V1)) %>% rename_with( ~ fn$V2[match(.x, fn$V1)], any_of(fn$V1)) %>% dplyr::select( - death_cause) %>% mutate(
	# Keep true dates even after the analysis end; t2e() applies censoring.
	date_cvd_death = death_dates$date_cvd_death
)
saveRDS(dat, paste0(indir, "/Rdata/death.rds"))
})


run_if("srd", {


# 🚩 SRD: p20008 is extracted as srd.year.tab.gz; keep its instance/array names.
# The icd step below pairs p20002 with these slots by eid + instance + array.
# p87 remains a clearly labelled approximate fallback when p20008 was not exported.
srd_interp <- file.path(indir, "rap", "srd.year.tab.gz")
use_interp <- file.exists(srd_interp)
srd_time_file <- if (use_interp) srd_interp else file.path(indir, "rap", "srd.time.tab.gz")
if (!file.exists(srd_time_file)) stop("Missing rap/srd.year.tab.gz (p20008) and legacy rap/srd.time.tab.gz (p87). Run extract_multi first.")
if (!use_interp) warning("rap/srd.year.tab.gz is unavailable; using legacy p87 age/year approximation. Extract p20008 to use UKB interpolated diagnosis years.", call. = FALSE)
date0 <- fread(srd_time_file, sep = "\t", header = TRUE) %>% as.data.frame()
expected_field <- if (use_interp) "p20008" else "p87"
if (!"eid" %in% names(date0) || !length(setdiff(names(date0), "eid")) ||
	any(!grepl(paste0("^", expected_field, "(?:_i[0-9]+)?(?:_a[0-9]+)?$"), setdiff(names(date0), "eid"), perl = TRUE)))
	stop("Unexpected columns in ", srd_time_file, "; preserve eid and the original ", expected_field, " instance/array column names.")
date0$eid <- as.character(date0$eid)
birth <- readRDS(paste0(indir, "/Rdata/phe.rds")) %>% dplyr::select(eid, birth_year)
birth$eid <- as.character(birth$eid)
stopifnot(!anyNA(date0$eid), !anyDuplicated(date0$eid), !anyNA(birth$eid), !anyDuplicated(birth$eid))
birth_year <- birth$birth_year[match(date0$eid, birth$eid)]
# Keep the SRD year limit tied to this dataset's declared follow-up end.
# The day-level event/censor decision remains in t2e(); do not move July-1
# approximations to the administrative cutoff.
year_cap <- as.integer(format(date_follow_end, "%Y"))
if (length(year_cap) != 1L || is.na(year_cap) || year_cap < 1900) stop("Invalid date_follow_end for SRD")
for (cc in setdiff(names(date0), "eid")) {
	date0[[cc]] <- ukb_srd_date(date0[[cc]], birth_year, interpolated = use_interp, year_cap = year_cap)
}
saveRDS(list(field = if (use_interp) "p20008" else "p87", input = srd_time_file,
	method = if (use_interp) "interpolated decimal year; floor fractional actual calendar-year days" else "calendar year or birth year plus reported age; July 1 approximation",
	precision = "approximate", year_cap = year_cap), file.path(indir, "Rdata", "srd.date.meta.rds"))
message("SRD date method: ", if (use_interp) "p20008 interpolated year" else "p87 approximate midyear fallback")
fwrite(date0, paste0(indir, "/rap/srd.date.tab.gz"), sep = "\t", na = "NA", quote = FALSE)
})


run_if("icd", {


# 🚩 ICD10, ICD9, OPCS4, SRD
# Scan each ICD source once and build all trait columns.
icd_sources <- trimws(unlist(strsplit(Sys.getenv("ICD_SOURCES", unset = "icd10,icd9,opcs4,srd"), ",|;|\\s+")))
icd_sources <- icd_sources[nzchar(icd_sources)]
save_aod_text <- toupper(Sys.getenv("SAVE_AOD_TEXT", unset = "TRUE")) %in% c("TRUE", "T", "1", "YES", "Y")
save_aod_list <- toupper(Sys.getenv("SAVE_AOD_LIST", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")
cat("ICD sources:", paste(icd_sources, collapse = ","), "\n")
cat("FAST_ICD=TRUE; SAVE_AOD_TEXT=", save_aod_text, "; SAVE_AOD_LIST=", save_aod_list, "\n", sep = "")
if (!exists("get_icd_fast", mode = "function")) stop("get_icd_fast() is missing. Update scripts/0f/phenotype.R first.", call. = FALSE)

for (src in icd_sources) {
	date0 <- fread(paste0(indir, "/rap/", src, ".date.tab.gz"), header = TRUE) %>% as.data.frame()
	date0$eid <- as.character(date0$eid)
	# p20002's canonical extract is srd.tab.gz. Never prefer an older srd.code
	# export when both exist; its rows/instances may belong to another extraction.
	code_file <- file.path(indir, "rap", if (src == "srd") "srd.tab.gz" else paste0(src, ".code.tab.gz"))
	if (!file.exists(code_file)) stop("Missing ", code_file, "; run extract_multi first.")
	code_raw <- fread(code_file, sep = "\t", header = TRUE, colClasses = "character") %>% as.data.frame()
	code0 <- ukb_expand_codes(code_raw, date0, src)
	date_fields <- setdiff(names(date0), "eid")
	names(date0) <- c("eid", paste0("d", seq_along(date_fields)))
	attr(date0, "ukb_date_fields") <- date_fields
	attr(code0, "ukb_code_fields") <- setdiff(names(code_raw), "eid")
	if (src == "srd" && file.exists(file.path(indir, "Rdata", "srd.date.meta.rds")))
		attr(date0, "date_method") <- readRDS(file.path(indir, "Rdata", "srd.date.meta.rds"))
	saveRDS(date0, paste0(indir, "/Rdata/", src, ".date0.rds"))
	saveRDS(code0, paste0(indir, "/Rdata/", src, ".code0.rds"))
	stopifnot(identical(code0$eid, date0$eid))
	fn <- fread(paste0(indir, "/common/", src, ".lst"), sep = "\t", select = 1:2, fill = TRUE, header = FALSE) %>% as.data.frame()
	dat <- get_icd_fast(src, code0, date0, fn, save_aod_list = save_aod_list, save_aod_text = save_aod_text)
	saveRDS(dat, paste0(indir, "/Rdata/fod.", src, ".rds"))
}

# Phecodes 
# Optional only. RUN_STEPS=icd,merge does not need PheWAS.rds.
# Default is FALSE because the PheWAS package is often unavailable on newer R versions.
make_phewas <- toupper(Sys.getenv("MAKE_PHEWAS", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")
cat("MAKE_PHEWAS=", make_phewas, "\n", sep = "")

if (isTRUE(make_phewas)) {
	icd10_code_rds <- paste0(indir, "/Rdata/icd10.code0.rds")
	icd10_date_rds <- paste0(indir, "/Rdata/icd10.date0.rds")
	if (!file.exists(icd10_code_rds) || !file.exists(icd10_date_rds)) {
		stop("MAKE_PHEWAS=TRUE but icd10.code0.rds/icd10.date0.rds is missing; PheWAS.rds was not rebuilt.", call. = FALSE)
	} else if (!requireNamespace("PheWAS", quietly = TRUE)) {
		stop(
			"MAKE_PHEWAS=TRUE but package PheWAS is not installed; PheWAS.rds was not rebuilt. ",
			"This does not affect fod/lod/cnt/aodtxt or all.rds.",
			call. = FALSE
		)
	} else {
        # The UK hospital export is ICD10, not ICD10CM. A full, documented
        # linkage-eligible roster is required so missing export/coverage is not
        # silently represented as disease-free controls.
        source(file.path(script_dir, "ukb", "f", "phewas.R"))
        phewas_roster <- Sys.getenv("PHEWAS_ROSTER", "")
        if (!nzchar(phewas_roster) || !file.exists(phewas_roster))
            stop("MAKE_PHEWAS=TRUE requires PHEWAS_ROSTER: documented id, phewas_observable, sex roster; see PHENOTYPE_REBUILD.md. Existing PheWAS.rds was not changed.")
        phewas_coverage <- Sys.getenv("PHEWAS_COVERAGE_DEFINITION", "")
        if (!nzchar(phewas_coverage)) phewas_coverage <- NULL
        phewas_out <- Sys.getenv("PHEWAS_REBUILD_OUTDIR", file.path(indir, "Rdata"))
        phewas_replace <- toupper(Sys.getenv("PHEWAS_REPLACE", "FALSE")) %in% c("TRUE", "1", "YES")
        manifest <- ukb_rebuild_phewas(icd10_code_rds, icd10_date_rds,
            phewas_roster, phewas_out, coverage_definition=phewas_coverage,
            overwrite=phewas_replace)
        cat("Wrote validated UK ICD10 PheWAS panel: ", file.path(phewas_out,"PheWAS.rds"),
            " | N=",manifest$n_roster,"; phenotypes=",manifest$n_phenotypes,"\n",sep="")
	}
} else {
	message("MAKE_PHEWAS=FALSE; skipped PheWAS.rds. This is not needed for icd/merge/all.rds.")
}
})


run_if("fod_ref", {


# 🚩 FOD coming with UKB
date0 <- fread(paste0(indir, "/rap/fod.date.tab.gz"), header = TRUE) %>% as.data.frame()
date0$eid <- as.character(date0$eid)
if (anyNA(date0$eid) || anyDuplicated(date0$eid)) stop("Invalid eid in fod.date.tab.gz")
# Keep the date fields exactly as exported, including coding-819 sentinels.
# Companion source fields can be separate or included in a combined export.
field_id <- function(x) suppressWarnings(as.integer(sub("^p([0-9]+).*$", "\\1", x)))
source_cols <- names(date0)[grepl("^p13[0-9]{4}(_|$)", names(date0)) & field_id(names(date0)) %% 2L == 1L]
source0 <- date0[, c("eid", source_cols), drop = FALSE]
date0 <- date0[, setdiff(names(date0), source_cols), drop = FALSE]
source_file <- file.path(indir, "rap", "fod.source.tab.gz")
if (file.exists(source_file)) {
	source_extra <- fread(source_file, sep = "\t", header = TRUE) %>% as.data.frame()
	source_extra$eid <- as.character(source_extra$eid)
	if (anyNA(source_extra$eid) || anyDuplicated(source_extra$eid)) stop("Invalid eid in fod.source.tab.gz")
	if (!setequal(source0$eid, source_extra$eid)) stop("FOD date/source exports have different participant sets")
	for (cc in setdiff(names(source_extra), "eid"))
		source0[[cc]] <- source_extra[[cc]][match(source0$eid, source_extra$eid)]
}
saveRDS(date0, paste0(indir, "/Rdata/fod.date0.rds"))
saveRDS(source0, paste0(indir, "/Rdata/fod.ref.source0.rds"))
if (ncol(source0) == 1L) message("No UKB source fields exported; date-check will mark source provenance unavailable. Re-run extract_fod with companion source fields.")
fn <- fread(paste0(indir, "/common/fod.lst"), sep = "\t", header = FALSE, select = 1:2) %>% as.data.frame()
if (anyNA(fn$V1) || anyNA(fn$V2) || any(!nzchar(trimws(fn$V2))))
	stop("fod.lst requires two tab-separated fields: field pattern and trait name")
dat <- data.frame(eid = date0$eid)
for (i in seq_len(nrow(fn))) {
		print(i); dat <- get_fod(dat.fod = date0, dat.this = dat, code.this = fn$V1[i], label.this = fn$V2[i])
}
saveRDS(dat, paste0(indir, "/Rdata/fod.ref.rds"))
})


run_if("gp_prep", {


# 🚩 GP prep: generate GP phecode lists before GP step
# This step only creates mapping lists. It does not depend on external phe.8gp.f.R.

pacman::p_load(tidyverse, data.table, stringr)

mapdir <- file.path(indir, "common", "map")
lstdir <- file.path(outdir, "gp_prep")
dir.create(lstdir, showWarnings = FALSE, recursive = TRUE)

r2 <- fread(file.path(mapdir, "read2_to_phecode.csv"))[, .(code = read2_code, phecode)]
r3 <- fread(file.path(mapdir, "ctv3_to_phecode.csv"))[, .(code = ctv3_code, phecode)]
r2 <- r2[nzchar(code) & nzchar(phecode)]
r3 <- r3[nzchar(code) & nzchar(phecode)]

write.table(
	r2[, .(codes = paste(sort(unique(code)), collapse = ",")), by = phecode],
	file.path(lstdir, "gp2.phecode.lst"),
	sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)
write.table(
	r3[, .(codes = paste(sort(unique(code)), collapse = ",")), by = phecode],
	file.path(lstdir, "gp3.phecode.lst"),
	sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
)

})


run_if("gp", {
# GP disease first-occurrence dates
pacman::p_load(data.table, fst)

gp_file <- paste0(indir, "/rap/raw/gp_clinical.tsv.gz")
# The old gp_clinical.fst omitted value1/value2/value3 and cannot support ABI,
# CAC, LVEF, BNP, or troponin extraction. Use a new full cache.
fst_file <- paste0(indir, "/rap/raw/gp_clinical.full.fst")
out_dir <- paste0(indir, "/Rdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rebuild_gp_fst <- toupper(Sys.getenv("GP_REBUILD_FULL_FST", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")
if (!file.exists(fst_file) || rebuild_gp_fst) {
	keep_gp <- c("eid", "data_provider", "event_dt", "read_2", "read_3",
						"value1", "value2", "value3")
	dt <- fread(gp_file, select = keep_gp, showProgress = TRUE)
	dt[, event_dt := as.integer(ukb_date_parse(event_dt))]
	write_fst(dt, fst_file, compress = 100)
	rm(dt); invisible(gc())
} else message("Skip full GP convert; fst exists: ", fst_file)

make_fod_from_fst <- function(src, fst_file) {
	code_col <- if (src == "gp2") "read_2" else "read_3"
	fn <- fread(paste0(indir, "/common/", src, ".lst"), sep = "\t",
						header = FALSE, fill = TRUE)[, .(codes_csv = V1, trait = V2)]
	long <- fn[!is.na(codes_csv) & !is.na(trait), .(
		code = trimws(unlist(strsplit(codes_csv, ",", fixed = TRUE)))
	), by = trait]
	long <- unique(long[!is.na(code) & nzchar(code)])
	setkey(long, code)
	dt <- as.data.table(read_fst(fst_file, columns = c("eid", "event_dt", code_col)))
	dt[, event_dt := ukb_date_valid(event_dt)]
	setnames(dt, code_col, "code")
	dt <- dt[!is.na(code) & nzchar(code)]
	x <- long[dt, on = "code", nomatch = 0L, allow.cartesian = TRUE]
	rm(dt); invisible(gc())
	x[, eid := as.character(eid)]
	# Missing and special dates remain coded disease evidence. AOD is the set
	# of distinct available recorded dates, not a count of independent episodes.
	ans <- x[, {
		d <- sort(unique(event_dt[!is.na(event_dt)]))
		list(fod = if (length(d)) min(d) else as.Date(NA),
			 lod = if (length(d)) max(d) else as.Date(NA),
			 cnt = .N, cnt_undated = sum(is.na(event_dt)),
			 aodtxt = if (length(d)) paste(format(d, "%Y-%m-%d"), collapse = "|") else NA_character_)
	}, by = .(eid, trait)]
	wide <- data.table(eid = unique(x$eid))
	traits0 <- unique(long$trait)
	for (metric in c("fod", "lod", "cnt", "cnt_undated", "aodtxt")) {
		prefix <- paste0(metric, "_", src, "_")
		if (nrow(ans)) {
			w <- dcast(ans, eid ~ trait, value.var = metric)
			old <- setdiff(names(w), "eid"); setnames(w, old, paste0(prefix, old))
			wide <- merge(wide, w, by = "eid", all.x = TRUE, sort = FALSE)
		}
		for (tr in traits0) {
			cc <- paste0(prefix, tr)
			if (!cc %in% names(wide)) wide[[cc]] <- rep(if (metric %in% c("fod", "lod")) as.Date(NA) else if (grepl("^cnt", metric)) 0L else NA_character_, nrow(wide))
			if (grepl("^cnt", metric)) wide[is.na(get(cc)), (cc) := 0L]
		}
	}
	for (cc in grep("^(fod|lod)_", names(wide), value = TRUE)) wide[, (cc) := ukb_date_valid(get(cc))]
	wide
}

fod_gp2 <- make_fod_from_fst("gp2", fst_file)
fod_gp3 <- make_fod_from_fst("gp3", fst_file)
saveRDS(fod_gp2, paste0(out_dir, "/fod.gp2.rds"))
saveRDS(fod_gp3, paste0(out_dir, "/fod.gp3.rds"))

tmp <- merge_check(fod_gp2 = as.data.frame(fod_gp2), fod_gp3 = as.data.frame(fod_gp3), by = "eid")
fod_gp <- as.data.table(merge(tmp$fod_gp2, tmp$fod_gp3, by = "eid", all = TRUE))
traits <- union(
	sub("^fod_gp2_", "", grep("^fod_gp2_", names(fod_gp), value = TRUE)),
	sub("^fod_gp3_", "", grep("^fod_gp3_", names(fod_gp), value = TRUE))
)
for (t in traits) {
	get_gp <- function(prefix, default) {
		cc <- paste0(prefix, t)
		if (cc %in% names(fod_gp)) fod_gp[[cc]] else rep(default, nrow(fod_gp))
	}
	a <- ukb_date_valid(get_gp("fod_gp2_", as.Date(NA)))
	b <- ukb_date_valid(get_gp("fod_gp3_", as.Date(NA)))
	z <- pmin(as.numeric(a), as.numeric(b), na.rm = TRUE)
	z[!is.finite(z)] <- NA_real_
	fod_gp[[paste0("fod_gp_", t)]] <- as.Date(z, origin = "1970-01-01")
	a <- ukb_date_valid(get_gp("lod_gp2_", as.Date(NA)))
	b <- ukb_date_valid(get_gp("lod_gp3_", as.Date(NA)))
	z <- pmax(as.numeric(a), as.numeric(b), na.rm = TRUE); z[!is.finite(z)] <- NA_real_
	fod_gp[[paste0("lod_gp_", t)]] <- as.Date(z, origin = "1970-01-01")
	for (metric in c("cnt", "cnt_undated")) {
		a <- get_gp(paste0(metric, "_gp2_"), 0L); b <- get_gp(paste0(metric, "_gp3_"), 0L)
		a[is.na(a)] <- 0L; b[is.na(b)] <- 0L
		fod_gp[[paste0(metric, "_gp_", t)]] <- a + b
	}
	a <- get_gp("aodtxt_gp2_", NA_character_); b <- get_gp("aodtxt_gp3_", NA_character_)
	fod_gp[[paste0("aodtxt_gp_", t)]] <- vapply(seq_len(nrow(fod_gp)), function(i) {
		v <- c(a[i], b[i]); v <- v[!is.na(v) & nzchar(v)]
		if (!length(v)) return(NA_character_)
		paste(sort(unique(unlist(strsplit(v, "|", fixed = TRUE)))), collapse = "|")
	}, character(1))
}
attr(fod_gp, "count_unit") <- "GP counts are Read2 plus CTV3 code matches; a raw record containing both can count twice. AOD dates are deduplicated."
for (cc in grep("^fod_", names(fod_gp), value = TRUE)) fod_gp[, (cc) := phe_date(get(cc))]
saveRDS(fod_gp, paste0(out_dir, "/fod.gp.rds"))

})


run_if("biom", {


# 🚩 BBC, Prot, Met, Img
fn <- fread(paste0(indir, "/common/bbc.lst"), sep = "\t",
					select = 1:2, header = FALSE, fill = TRUE) %>% as.data.frame()
fn$V1[duplicated(fn$V1)]; fn$V2[duplicated(fn$V2)]
bbc_raw <- fread(paste0(indir, "/rap/bbc.tab.gz"), sep = "\t",
					header = TRUE) %>% as.data.frame()
names(bbc_raw) <- gsub("_i0$|_i0_a0$", "", names(bbc_raw))
setDT(bbc_raw)
setnames(bbc_raw, old = fn$V1, new = fn$V2, skip_absent = TRUE)
bbc_raw <- as.data.frame(bbc_raw)
bbc_raw$eid <- as.character(bbc_raw$eid)
for (cc in intersect(c("date_bb_UALB", "date_bb_UCRE"), names(bbc_raw))) {
	bbc_raw[[cc]] <- phe_date(bbc_raw[[cc]])
}

# clean_biom() coerces non-numeric results such as "<6.7" to NA. Dates are
# kept outside numeric biomarker QC and merged back unchanged.
bbc_date_cols <- grep("^date_", names(bbc_raw), value = TRUE)
bbc_dates <- bbc_raw[, c("eid", bbc_date_cols), drop = FALSE]
uacr <- derive_ukb_uacr(bbc_raw)
bbc_main <- bbc_raw[, setdiff(names(bbc_raw), bbc_date_cols), drop = FALSE]
bbc_main <- clean_biom(bbc_main, miss_col = 0.2, nonnegative = TRUE)
bbc <- Reduce(function(x, y) merge(x, y, by = "eid", all = TRUE, sort = FALSE),
					list(bbc_main, bbc_dates, uacr))
saveRDS(bbc, paste0(indir, "/Rdata/bbc.rds"))
cat("UACR assessed from numeric values:",
		sum(bbc$bb_UACR_assessed %in% TRUE, na.rm = TRUE), "\n")

prot <- fread(paste0(indir, "/rap/raw/prot.tab.gz"), sep = "\t",
					header = TRUE) %>% as.data.frame()
names(prot)[ - 1] <- toupper(names(prot)[ - 1])
# Olink NPX is a signed log2 scale, not an absolute concentration.
prot <- clean_biom(prot, miss_col = 0.2, nonnegative = FALSE)
saveRDS(prot, paste0(indir, "/Rdata/prot.rds"))

fn <- fread(paste0(indir, "/common/met.lst"), sep = "\t",
					select = 1:2, header = FALSE, fill = TRUE) %>% as.data.frame()
fn$V1[duplicated(fn$V1)]; fn$V2[duplicated(fn$V2)]
met <- fread(paste0(indir, "/rap/met.tab.gz"), sep = "\t",
					header = TRUE) %>% as.data.frame()
names(met) <- gsub("_i0$", "", names(met))
setDT(met)
col2calc <- fn[grepl("/", fn$V1), ]
col2rename <- fn[!grepl("/", fn$V1), ]
for (i in seq_len(nrow(col2calc))) {
	met[, (col2calc$V2[i]) := eval(parse(text = col2calc$V1[i]))]
}
setnames(met, old = col2rename$V1, new = col2rename$V2, skip_absent = TRUE)
met <- clean_biom(as.data.frame(met), nonnegative = TRUE)
saveRDS(met, paste0(indir, "/Rdata/met.rds"))

# Imaging uses a reusable raw reader; keep the established cleaned img.rds too.


# 🚩 img.f
# Preserve numeric imaging measurements without imputation or missingness filters.
# Metadata is kept separately from measured brain phenotypes.
prepare_ukb_imaging <- function(indir) {
	file <- file.path(indir, 'rap/img.tab.gz')
	x <- data.table::fread(file, showProgress = FALSE)
	x[, eid := as.character(eid)];stopifnot(!anyDuplicated(x$eid))
	extra <- file.path(indir, 'rap/img.extra.tab.gz')
	if(file.exists(extra)){
		e <- data.table::fread(extra, showProgress = FALSE);e[, eid := as.character(eid)]
		missing <- setdiff(names(e), names(x))
		if(length(missing))x <- merge(x, e[, c('eid', missing), with = FALSE], by = 'eid', all = TRUE, sort = FALSE)
	}
	map <- data.table::fread(file.path(indir, 'common/img.lst'), header = FALSE, fill = TRUE)[, 1:2]
	data.table::setnames(map, c('field', 'label'))
	cols <- grep('_i2$', names(x), value = TRUE)
	nums <- cols[vapply(x[, ..cols], is.numeric, logical(1))]
	nums <- setdiff(nums, c('p54_i2'))
	img <- x[, c('eid', nums), with = FALSE]
	data.table::setnames(img, nums, paste0('img_', sub('_i2$', '', nums)))
	img <- img[rowSums(!is.na(img[, - 1])) > 0]
	dir.create(file.path(indir, 'Rdata'), recursive = TRUE, showWarnings = FALSE)
	saveRDS(as.data.frame(img), file.path(indir, 'Rdata/img.raw.rds'), compress = FALSE)
	visit_cols <- intersect(c('p53_i2', 'p53_i3', 'p54_i2', 'p54_i3'), names(x))
	visits <- x[, c('eid', visit_cols), with = FALSE]
	for(cc in intersect(c('p53_i2', 'p53_i3'), names(visits)))visits[, (cc) := as.Date(get(cc))]
	saveRDS(as.data.frame(visits), file.path(indir, 'Rdata/img.visits.rds'))
	map <- map[field %in% sub('^img_', '', names(img))];map[, variable := paste0('img_', field)]
	data.table::fwrite(map, file.path(indir, 'Rdata/img.dictionary.csv'))
	message('Unimputed imaging: ', nrow(img), ' participants, ', ncol(img) - 1, ' variables')
	as.data.frame(img)
}


img <- prepare_ukb_imaging(indir)
img <- clean_biom(img)
saveRDS(img, paste0(indir, "/Rdata/img.rds"))

})


run_if("hla_pca", {


# 🚩 HLA & PCA
hla <- fread(paste0(indir, "/rap/hla.tab.gz"), sep = "\t", header = TRUE) %>% as.data.frame() # 362 commas
	alleles <- fread(paste0(indir, "/common/ukb_hla_v2.txt"), header = TRUE) %>% as.data.frame() # 1 row, 362 columns
	hla <- hla %>% separate(p22182, into = names(alleles), remove = TRUE, convert = TRUE, sep = ",") %>% drop_na()
	names(hla)[ - 1] = paste0("hla_", names(hla)[ - 1])
	pc <- prcomp(hla[, - 1], center = TRUE, scale. = TRUE)
	pc <- pc$x[, 1:4]; colnames(pc) <- paste0("hla_PC", 1:4)
	hla <- cbind(hla, pc)
	saveRDS(hla, paste0(indir, "/Rdata/hla.rds"))
pca <- fread(paste0(indir, "/rap/pca.tab.gz"), header = TRUE) %>% as.data.frame()
	colnames(pca) <- c("eid", paste0("PC", 1:40))
	saveRDS(pca, paste0(indir, "/Rdata/pca.rds"))
})


run_if("pgs", {


# 🚩 GEN & PGS
imp <- fread(file.path(gen_dir, "imp/snp.raw.gz"), header = TRUE) %>% as.data.frame() %>% rename(eid = IID)
typ <- fread(file.path(gen_dir, "typ/snp.raw.gz"), header = TRUE) %>% as.data.frame() %>% rename(eid = IID)
apoe <- fread(file.path(gen_dir, "hap/apoe.hap"), header = TRUE) %>% as.data.frame() %>% rename(eid = IID)
tmp <- merge_check(imp = imp, typ = typ, by = "eid")
gen0 <- merge(merge(tmp$imp, tmp$typ, by = "eid", all = TRUE), apoe, by = "eid", all = TRUE)
gen <- gen0 %>% mutate(
	se = ifelse(fut2.rs601338_G %in% 0:1, 0, ifelse(fut2.rs601338_G == 2, 1, NA)), # [PMID: 30345375]
	apoe = factor(apoe, levels = c("e3e3", "e2e2", "e2e3", "e2e4", "e3e4", "e4e4")),
	apoe.e4 = stringr::str_count(apoe, "4"),
	sp1 = ifelse(sp1.Z.rs28929474_C == 0, "ZZ", ifelse (sp1.S.rs17580_T == 0, "SS",
		ifelse( (sp1.Z.rs28929474_C == 1 & sp1.S.rs17580_T == 1), "SZ", ifelse( (sp1.Z.rs28929474_C == 1 & sp1.S.rs17580_T == 2), "MZ",
		ifelse( (sp1.Z.rs28929474_C == 2 & sp1.S.rs17580_T == 1), "MS", ifelse( (sp1.Z.rs28929474_C == 2 & sp1.S.rs17580_T == 2), "MM",
		NA)))))),
	sp1.M = stringr::str_count(sp1, "M"), sp1.S = stringr::str_count(sp1, "S"), sp1.Z = stringr::str_count(sp1, "Z")
)
saveRDS(gen, paste0(indir, "/Rdata/gen.rds"))

pgs <- fread(file.path(pgs_dir, "all.pgs.gz"), header = TRUE) %>% as.data.frame()
names(pgs) <- sub("\\.score_sum$", ".pgs", names(pgs))
saveRDS(pgs, paste0(indir, "/Rdata/pgs.rds"))

#rel <- fread(file.path(gen_dir, "typ/ukb_rel.dat"), header = TRUE) %>% as.data.frame()
#rel.id <- ukbtools::ukb_gen_samples_to_remove(data = rel, ukb_with_data = dat$eid, cutoff = 0.0882) 
#rel.res <- rel_seek(rel = rel, eid = unique(c(rel$ID1, rel$ID2)), cutoff = 0.0882)
})


run_if("merge", {


# 🚩 merge (diet.rds must exist)
pacman::p_load(matrixStats)
lst <- c("phe", "bbc", "diet", "move", "ses", "death", "fod.icd10", "fod.icd9", "fod.opcs4", "fod.srd", "fod.gp", "fod.ref", "pca", "gen", "pgs")
rds_files <- setNames(file.path(indir, "Rdata", paste0(lst, ".rds")), lst)
missing_rds <- names(rds_files)[!file.exists(rds_files)]
if (length(missing_rds)) {
	stop("Missing required RDS file(s) for merge: ", paste(missing_rds, collapse = ", "),
		". Run the corresponding earlier step(s) before RUN_STEPS=merge.", call. = FALSE)
}
for (l in lst) {
	x <- readRDS(rds_files[[l]])
	if (!"eid" %in% names(x)) stop("Missing eid column in ", rds_files[[l]], call. = FALSE)
	assign(l, x %>% mutate(eid = as.character(eid)))
}
dat_list <- merge_check(mget(lst), by = "eid")
# Build the participant-level outer join column-wise to control memory.
dup_eid <- names(dat_list)[vapply(dat_list, function(x) anyDuplicated(x$eid) > 0L, logical(1))]
if (length(dup_eid)) {
	stop("Duplicate eid values in merge input(s): ", paste(dup_eid, collapse = ", "), call. = FALSE)
}
all_eid <- unique(unlist(lapply(dat_list, function(x) as.character(x$eid)), use.names = FALSE))
all_eid <- all_eid[!is.na(all_eid) & !grepl("^-", all_eid)]
merged_cols <- list(eid = all_eid)
for (nm in names(dat_list)) {
	x <- dat_list[[nm]]
	idx <- match(all_eid, x$eid)
	for (cc in setdiff(names(x), "eid")) merged_cols[[cc]] <- x[[cc]][idx]
}
dat0 <- as.data.frame(merged_cols, check.names = FALSE, optional = TRUE)
rm(list = lst)
rm(merged_cols, dat_list, x, idx, all_eid)
gc()
aodtxt_cols0 <- grep("^aodtxt_", names(dat0), value = TRUE)
cat("Merged aodtxt columns before derived variables: ", length(aodtxt_cols0), "\n", sep = "")
if (length(aodtxt_cols0)) cat("First aodtxt columns: ", paste(head(aodtxt_cols0, 5), collapse = ", "), "\n", sep = "")
	idate_cols <- vapply(dat0, inherits, logical(1), "IDate"); idate_cols
	bad_idate_cols <- idate_cols & !vapply(dat0, function(x) is.integer(unclass(x)), logical(1)); bad_idate_cols
	# data.table::IDate objects with double storage are rejected by vctrs/dplyr.
	# Normalize every IDate to base Date before derived variables and saveRDS.
	for (cc in names(dat0)[idate_cols]) {
		dat0[[cc]] <- as.Date(as.numeric(unclass(dat0[[cc]])), origin = "1970-01-01")
	}

dat0 <- dat0 %>% mutate(
	eid = as.character(eid),
	fod_icd10_death = date_death,
	drug.lipid = ifelse(matrixStats::rowSums2(as.matrix(across(starts_with("drug.big3"), ~ grepl("1", as.character(.))))) > 0, 1, 0),
	drug.dm = ifelse(matrixStats::rowSums2(as.matrix(across(starts_with("drug.big3"), ~ grepl("3", as.character(.))))) > 0, 1, 0),
	drug.htn = ifelse(matrixStats::rowSums2(as.matrix(across(starts_with("drug.big3"), ~ grepl("2", as.character(.))))) > 0, 1, 0),
	date_cad = do.call(pminDate2, across(any_of(c("fod_icd10_cvd_cad", "fod_icd9_cvd_cad", "fod_srd_cvd_cad")))),
	date_mi = ukb_mi_date(dat0),
	date_stroke_i = do.call(pminDate2, across(matches("^fod_icd.*_cvd_stroke_i$"))),
	date_stroke_h = do.call(pminDate2, across(matches("^fod_icd.*_cvd_stroke_(i|s)h$"))),
	date_stroke = pminDate2(date_stroke_i, date_stroke_h),
	date_mace3 = pminDate2(date_mi, date_stroke, date_cvd_death)
) %>% mutate( # ⓸⓸ = LE8
	# 
	diet.sum = diet.medito.sum, diet.pts = diet.medito.pts,
	#  
	across(c(days_pa_mod, dura_pa_mod, days_pa_vig, dura_pa_vig), ~ ifelse(.x < 0, NA, .x)),
	pa_mod_w = days_pa_mod * pmin(dura_pa_mod, 180),
	pa_vig_w = days_pa_vig * pmin(dura_pa_vig, 180),
	pa_mins = rowSums1(data.frame(pa_mod_w, 2 * pa_vig_w), per = 0.5),
	pa.pts = cut(pa_mins, breaks = c(0, 1, 30, 60, 90, 120, 150, Inf), labels = c(0, 20, 40, 60, 80, 90, 100), right = FALSE) %>% as.character() %>% as.numeric(),
	# 
	smoke.c = factor(smoke_status, levels = 0:2, labels = c("never", "previous", "current")),
	smoke_quit_years = ifelse(smoke_quit_age > 0 & age >= smoke_quit_age, age - smoke_quit_age, NA_real_),
	smoke_base = case_when(smoke_status == 0 ~ 100, smoke_status == 2 ~ 0, smoke_status == 1 & smoke_quit_years >= 5 ~ 75, smoke_status == 1 & smoke_quit_years >= 1 ~ 50, smoke_status == 1 ~ 25, TRUE ~ NA_real_),
	smoke_shs_home = if_else(smoke_2nd_home > 0, 1L, 0L),
	smoke.pts = pmax(smoke_base - if_else(smoke_shs_home == 1L & smoke_base > 0, 20, 0), 0),
	# 
	sleep.pts = cut(sleep_duration, breaks = c(0, 4, 5, 6, 7, 9, 10, 24), labels = c(0, 20, 40, 70, 100, 90, 40), right = FALSE) %>% as.character() %>% as.numeric(),
	# 
	bmi.c = cut(bmi, breaks = c(10, 18.5, 25, 30, 100), labels = c("lean", "healthy", "overweight", "obese")),
	bmi.c = factor(bmi.c, levels = c("healthy", "lean", "overweight", "obese")),
	bmi.pts = cut(bmi, breaks = c(0, 25, 30, 35, 40, Inf), labels = c(100, 70, 30, 15, 0), right = FALSE) %>% as.character() %>% as.numeric(),
	# 
	nonhdl = pmax((bb_TC - bb_HDL) * 38.67, 0), # Units of measurement are mmol/L, convert to mg/dL
	nonhdl.c = cut(nonhdl, breaks = c(0, 130, 160, 190, 220, Inf), labels = c(100, 60, 40, 20, 0), right = FALSE) %>% as.character() %>% as.numeric(),
	nonhdl.pts = ifelse(drug.lipid == 1 & nonhdl.c > 0, pmax(nonhdl.c - 20, 0), nonhdl.c),
	# 
	hba1c_ngsp = bb_HBA1C * 0.0915 + 2.15, # NGSP%
	dm.yes = if_else(!is.na(fod_srd_t2dm) | !is.na(fod_icd10_t2dm) | !is.na(fod_ref_t2dm) | drug.dm == 1 | (!is.na(hba1c_ngsp) & hba1c_ngsp >= 6.5), 1, 0),
	hba1c.pts = case_when(is.na(hba1c_ngsp) ~ NA_real_, dm.yes == 0 & hba1c_ngsp < 5.7 ~ 100, dm.yes == 0 & hba1c_ngsp < 6.5 ~ 60, hba1c_ngsp < 7 ~ 40, hba1c_ngsp < 8 ~ 30, hba1c_ngsp < 9 ~ 20, hba1c_ngsp < 10 ~ 10, TRUE ~ 0),
	#  sbp_man[93], dbp_man[94], dbp_auto[4079], sbp_auto[4080]
	sbp = rowMeans(across(matches("^sbp_.*_i0.*")), na.rm = TRUE),
	dbp = rowMeans(across(matches("^dbp_.*_i0.*")), na.rm = TRUE),
	htn.yes = if_else(!is.na(fod_srd_cvd_htn) | !is.na(fod_icd10_cvd_htn) | !is.na(fod_ref_cvd_htn) | drug.htn == 1, 1L, 0L),
	bp_base = case_when(sbp >= 160 | dbp >= 100 ~ 0, (sbp >= 140 & sbp < 160) | (dbp >= 90 & dbp < 100) ~ 25, (sbp >= 130 & sbp < 140) | (dbp >= 80 & dbp < 90) ~ 50, (sbp >= 120 & sbp < 130) & (dbp > 0 & dbp < 80) ~ 75, (sbp > 0 & sbp < 120) & (dbp > 0 & dbp < 80) ~ 100, TRUE ~ NA_real_),
	bp.pts = if_else(coalesce(drug.htn, 0) == 1 & bp_base > 0, pmax(bp_base - 20, 0), bp_base),
	# bp1.pts = case_when(is.na(sbp) | is.na(dbp) ~ NA_real_, TRUE ~ pmax(pmin(100 - 2 * pmax(sbp - 115, 0) - 2 * pmax(dbp - 75, 0), 100), 0)),
	# bp2.pts = case_when(is.na(htn.yes) & is.na(drug.htn) ~ NA_real_, TRUE ~ pmax(100 - 30 * coalesce(htn.yes, 0) - 10 * coalesce(drug.htn, 0), 0))	
	# 
	alcohol.c = factor(alcohol_status, levels = 0:2, labels = c("never", "previous", "current")),
	alcohol_days = case_match(alcohol_freq, 1 ~ 7, 2 ~ 3.5, 3 ~ 1.5, 4 ~ 0.5, 5 ~ 0.1, 6 ~ 0, .default = NA),
	alcohol_units = case_match(alcohol_amount, 1 ~ 1.5, 2 ~ 3.5, 3 ~ 5.5, 4 ~ 8, 5 ~ 12, .default = NA),
	# alcohol.pts = case_when(alcohol_status %in% c(0, 1) ~ 100, sex == 0 ~ as.numeric(as.character(cut(alcohol_days * alcohol_units * 8 / 7, c(0, 7, 14, 21, 28, Inf), c(100, 80, 50, 20, 0), right = FALSE))), sex == 1 ~ as.numeric(as.character(cut(alcohol_days * alcohol_units * 8 / 7, c(0, 14, 28, 42, 56, Inf), c(100, 80, 50, 20, 0), right = FALSE))), TRUE ~ NA_real_),
	alcohol.pts = case_when(alcohol_status %in% c(0, 1) ~ 100, is.na(alcohol_status) ~ NA_real_, sex == 0 ~ as.numeric(as.character(cut(coalesce(alcohol_days * alcohol_units * 8 / 7, 14.1), c(0, 7, 14, 21, 28, Inf), c(100, 80, 50, 20, 0), right = FALSE))), sex == 1 ~ as.numeric(as.character(cut(coalesce(alcohol_days * alcohol_units * 8 / 7, 28.1), c(0, 14, 28, 42, 56, Inf), c(100, 80, 50, 20, 0), right = FALSE))), TRUE ~ NA_real_)
) %>% mutate(
	le8.sco = rowMeans1(across(all_of(vars.le8)), 0.7), # 
	le8r.sco =  - le8.sco, le8r.3c = f3c(le8r.sco),
	cvd_cad.pgs.3c = f3c(cvd_cad.pgs), # 
	cvh.3c = cut(le8.sco, breaks = c(0, 50, 80, Inf), right = FALSE, labels = c("unhealthy", "average", "healthy")),
	cvh.3c = factor(cvh.3c, levels = c("healthy", "average", "unhealthy"))
)

# Save the large merged baseline table once. Downstream steps should read this as dat_base/dat_all,
# and should avoid reassigning dat0 after this point.
saveRDS(dat0, paste0(indir, "/Rdata/all.rds"))
dup_cols <- grep("\\.x$|\\.y$", names(dat0), value = TRUE)
if (length(dup_cols)) warning("Columns with .x/.y after merge: ", paste(dup_cols, collapse = ", "), call. = FALSE)
aodtxt_cols <- grep("^aodtxt_", names(dat0), value = TRUE)
cat("Merged aodtxt columns saved in all.rds: ", length(aodtxt_cols), "\n", sep = "")
if (length(aodtxt_cols)) cat("Example aodtxt columns: ", paste(head(aodtxt_cols, 10), collapse = ", "), "\n", sep = "")
dat0 %>% dplyr::select(matches("^(date_|fod_|lod_|cnt_|aodtxt_).*")) %>%
	sapply(function(x) format(sum(!is.na(x)), big.mark = ",")) %>%
	as.matrix() %>% print(quote = FALSE)
})


run_if("drug", {


# 🚩 drug⭕
prot <- fread(paste0('/mnt/f/gwas/prot/map_3k_v1.tsv'), sep = '\t', select = c(1, 3, 5), header = TRUE, fill = TRUE) %>% as.data.frame()
	drug <- fread(paste0(dir0, '/files/prot-drug.txt'), sep = '\t', header = TRUE, fill = TRUE) %>% as.data.frame()
	tmp <- merge_check(prot = prot, drug = drug, by.x = 'UniProt', by.y = 'UNIPROT_ACCESSION')
	drug <- merge(tmp$prot, tmp$drug, by.x = 'UniProt', by.y = 'UNIPROT_ACCESSION') %>%
	mutate(drug = ifelse(APPROVED_DRUG_TARG_CONF != '', 'approved', ifelse(CLINICAL_DRUG_TARG_CONF != '', 'clinical', ifelse(DRUG_POTENTIAL_TARGET != '', 'potential', NA)))) %>%
	dplyr::select(ID, olink_target_fullname, drug) %>% arrange('drug') %>% setNames(c("X", "olink_target", "drug")) %>% distinct(X, .keep_all = TRUE)
	saveRDS(drug, paste0(indir, "/Rdata/drug.rds"))
})


run_if("phe4gwas", {


# 🚩 GWAS data
dat <- read_rds("all")
bbc_gwas <- read_rds("bbc")
stopifnot(all(c("eid", "bb_LDL") %in% names(bbc_gwas)),
					!anyNA(dat$eid), !anyDuplicated(dat$eid),
					!anyNA(bbc_gwas$eid), !anyDuplicated(bbc_gwas$eid))
# Match by ID without changing the full cohort's order or duplicating rows.
dat$ldl <- bbc_gwas$bb_LDL[match(as.character(dat$eid), as.character(bbc_gwas$eid))]
message("input LDL: ", file.path(indir, "Rdata", "bbc.rds"),
				"; bb_LDL -> ldl; nonmissing=", sum(!is.na(dat$ldl)))
rm(bbc_gwas)
for (Y in c("cvd_cad", "t2dm", "stroke", "cvd_stroke_i", "stroke_o")) {
	# Diabetes evidence must refer to diabetes, not cardiovascular disease.
	domain <- if (Y == "t2dm") "cnt_icd10_t2dm" else "cvd"
	stopifnot(paste0("fod_icd10_", Y) %in% names(dat))
	dat <- t2e(dat, domain, paste0("fod_icd10_", Y), "birth_date", "date_attend", "date_lost", "date_death", date_follow_end, Y, "year")
}
dat1 <- dat %>% filter(ethnic.c == "White") %>% rename(IID = eid) %>% mutate(FID = IID)
# ADuLT must use the original dates, before select() drops them. Existing
# .t2e/.Yt2e are retained; ADuLT includes dated prevalent cases as well.
# Enable with ADULT_ENABLE=TRUE; the default new phenotype is cvd_cad.adu.
if (adu_bool(Sys.getenv("ADULT_ENABLE", "FALSE")) &&
		adu_bool(Sys.getenv("ADULT_CIP_FROM_UKB", "FALSE"))) {
	adu_cohort_cip(dat1, "cvd_cad", date_follow_end,
		Sys.getenv("ADULT_CIP_FILE", file.path(indir, "common", "adu_cip.tsv")),
		audit_dir = file.path(outdir, "adu"),
		input_source = paste(file.path(indir, "Rdata", "all.rds"), "White GWAS cohort"))
}
dat1 <- ukb_add_adult(dat1, indir = indir, outdir = outdir,
					end_date = date_follow_end, default_traits = "cvd_cad")
dat1 <- dat1 %>% mutate(center = factor(center)) %>%
	dplyr::select(FID, IID, ethnic.c, center, tdi, edu.sco, age, sex, bmi, height, ldl,
					drug.lipid, drug.dm, drug.htn, bald,
					matches("^bald1|^cvd_cad|^t2dm|^stroke_|^hap|t2e$|^PC[1-4]$"), ends_with(".adu"), - starts_with("happy_"))

# mm <- model.matrix(~ center - 1, data = dat1)
# colnames(mm) <- paste0("CR", seq_len(ncol(mm)))
# dat1 <- cbind(dat1, as.data.table(mm))
dir.create(file.path(indir, "common"), recursive = TRUE, showWarnings = FALSE)
# Stage the completed table beside the destination; errors cannot truncate the
# existing phenotype file. Rename only after ADuLT and the entire write succeed.
phe_destination <- file.path(indir, "common", "ukb.phe")
phe_staged <- tempfile("ukb.phe.", tmpdir = dirname(phe_destination))
write.table(dat1, phe_staged, na = "NA", append = FALSE, quote = FALSE, row.names = FALSE)
all_destination <- file.path(indir, "Rdata", "all.rds")
all_staged <- tempfile("all.rds.", tmpdir = dirname(all_destination))
# Save the full cohort, not the White-only GWAS export.
saveRDS(dat, all_staged)
if (!file.rename(all_staged, all_destination)) {
	unlink(c(all_staged, phe_staged))
	stop("Could not replace ", all_destination)
}
if (!file.rename(phe_staged, phe_destination)) {
	unlink(phe_staged)
	stop("Could not replace ", phe_destination)
}
message("output full phenotype: ", all_destination, "; N=", nrow(dat))
message("output GWAS phenotype: ", phe_destination, "; N=", nrow(dat1),
				"; height=", sum(!is.na(dat1$height)), "; ldl=", sum(!is.na(dat1$ldl)),
				"; t2dm events=", sum(dat1$t2dm.Yt2e == 1, na.rm = TRUE))
})
