# EMS120 FINAL: TEST ML/CV/data front-end + CURRENT main/supplementary figure suite.
ems120_prime_conda_ssl <- function() {
	sys <- Sys.info()[["sysname"]]
	if (identical(sys, "Windows")) return(invisible(FALSE))
	py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
	candidates <- c(
		if (nzchar(py)) normalizePath(file.path(dirname(py), ".."), winslash = "/", mustWork = FALSE) else "",
		file.path(Sys.getenv("HOME"), "anaconda3", "envs", "ai"),
		file.path(Sys.getenv("HOME"), "miniconda3", "envs", "ai"),
		"/home/huangj/anaconda3/envs/ai",
		"/home/jiehuang001/anaconda3/envs/ai"
	)
	candidates <- unique(candidates[nzchar(candidates)])
	env_prefix <- candidates[dir.exists(candidates)][1]
	if (is.na(env_prefix) || !nzchar(env_prefix)) return(invisible(FALSE))
	lib_dir <- file.path(env_prefix, "lib")
	if (!dir.exists(lib_dir)) return(invisible(FALSE))
	old_ld <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
	if (!grepl(paste0("(^|:)", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", lib_dir), "(:|$)"), old_ld)) {
		Sys.setenv(LD_LIBRARY_PATH = paste(c(lib_dir, old_ld[nzchar(old_ld)]), collapse = ":"))
	}
	for (lib in file.path(lib_dir, c("libcrypto.so.3", "libssl.so.3"))) {
		if (file.exists(lib)) try(dyn.load(lib, local = FALSE, now = TRUE), silent = TRUE)
	}
	invisible(TRUE)
}
ems120_prime_conda_ssl()

ems120_cli_args <- commandArgs(trailingOnly = TRUE)
ems120_arg_value <- function(flag, default = "") {
	i <- match(flag, ems120_cli_args)
	if (is.na(i) || i >= length(ems120_cli_args)) return(default)
	ems120_cli_args[[i + 1L]]
}
ems120_step_rank <- c(fig1 = 1L, fig2 = 2L, fig3 = 3L, fig4 = 4L, fig5 = 5L, fig6 = 6L, figs = 7L)
ems120_start_step_requested <- tolower(trimws(ems120_arg_value("--start_step", Sys.getenv("START_STEP", unset = "all"))))
ems120_end_step <- tolower(trimws(ems120_arg_value("--end_step", Sys.getenv("END_STEP", unset = Sys.getenv("STOP_AFTER_STEP", unset = "")))))
if (!ems120_start_step_requested %in% c("all", names(ems120_step_rank))) {
	stop("START_STEP must be one of: all, fig1, fig2, fig3, fig4, fig5, fig6, figS.", call. = FALSE)
}
if (nzchar(ems120_end_step) && !ems120_end_step %in% names(ems120_step_rank)) {
	stop("END_STEP/STOP_AFTER_STEP must be one of: fig1, fig2, fig3, fig4, fig5, fig6, figS.", call. = FALSE)
}
# Supplementary figures depend on objects created in the main-figure blocks.
# Treat --start-step figs as a cached-data redraw from fig1 onward so those
# dependencies are available; partial starts at fig2-fig6 skip supplementary output.
ems120_start_step <- if (identical(ems120_start_step_requested, "figs")) "fig1" else ems120_start_step_requested
ems120_should_run <- function(step) {
	step <- tolower(trimws(step))
	if (identical(ems120_start_step, "all")) return(TRUE)
	unname(ems120_step_rank[[step]]) >= unname(ems120_step_rank[[ems120_start_step]])
}
ems120_maybe_exit_after <- function(step) {
	step <- tolower(trimws(step))
	if (nzchar(ems120_end_step) && identical(step, ems120_end_step)) {
		cat(sprintf("END_STEP=%s reached; stopping EMS120 pipeline.\n", ems120_end_step))
		quit(save = "no", status = 0, runLast = FALSE)
	}
	invisible(FALSE)
}

pacman::p_load(readxl, writexl, data.table, tidyverse, scales, RColorBrewer, reticulate, lubridate, patchwork, zoo, broom, forcats, circlize, jpeg, png, grid)

dir0 <- if (Sys.info()[["sysname"]] == "Windows") "D:" else "/mnt/d"
dir.analysis <- Sys.getenv("EMS120_ANALYSIS_DIR", unset = file.path(dir0, "analysis", "ems120"))
script.dir <- Sys.getenv("EMS120_SCRIPT_DIR", unset = file.path(dir0, "scripts", "ems120"))

dir.cache <- file.path(dir.analysis, "dat"); 
dir.create(dir.cache, showWarnings = FALSE, recursive = TRUE)
file.dat.phone.rds <- file.path(dir.cache, "dat.ml_phone.rds")
file.dat.geo.rds <- file.path(dir.cache, "dat.ml_geo.rds")
file.dat.rds <- file.path(dir.cache, "dat.list.rds")

invisible(lapply(c("0phe.f.R", "plot.f.R"), function(f) { ff <- file.path(dir0, "scripts", "0f", f); if (file.exists(ff)) source(ff) }))
source(file.path(script.dir, "f", "ems120.f.R"))

ml_python_ready <- FALSE
ensure_ml_python <- function() {
	if (isTRUE(ml_python_ready)) return(invisible(TRUE))
	source(file.path(dir0, "scripts", "0f", "0conf_ML.R"))
	py_file <- file.path(script.dir, "f", "ems120.py")
	if (!file.exists(py_file)) stop("Missing Python helper file: ", py_file, call. = FALSE)
	# source_python() is called inside this function, so its default envir would be
	# the local function frame.  Use .GlobalEnv explicitly; otherwise ml_phone(),
	# ml_geo(), and ml_dx() disappear after ensure_ml_python() returns.
	reticulate::source_python(py_file, envir = .GlobalEnv)
	required_py_funs <- c("ml_phone", "ml_geo", "release_geo_model", "ml_dx", "run_dx_stratified_cv", "refresh_dx_cv_derived_outputs")
	missing_py_funs <- required_py_funs[!vapply(required_py_funs, exists, logical(1), envir = .GlobalEnv, mode = "function")]
	if (length(missing_py_funs)) {
		stop(sprintf("Python helper functions not loaded from %s: %s", py_file, paste(missing_py_funs, collapse = ", ")), call. = FALSE)
	}
	ml_python_ready <<- TRUE
	invisible(TRUE)
}

setwd2(dir.analysis)

.ems120_cleanup <- new.env(parent = emptyenv())
.ems120_cleanup$temp_rds <- character()
.ems120_register_temp_rds <- function(path) {
	path <- normalizePath(path, winslash = "/", mustWork = FALSE)
	.ems120_cleanup$temp_rds <- unique(c(.ems120_cleanup$temp_rds, path))
	path
}
.ems120_remove_temp_rds <- function(paths = .ems120_cleanup$temp_rds) {
	paths <- unique(paths[!is.na(paths) & nzchar(paths)])
	if (!length(paths)) return(invisible(paths))
	existing <- file.exists(paths)
	if (any(existing)) unlink(paths[existing], force = TRUE)
	invisible(paths)
}
reg.finalizer(.ems120_cleanup, function(e) {
	paths <- unique(e$temp_rds[!is.na(e$temp_rds) & nzchar(e$temp_rds)])
	if (length(paths)) {
		existing <- file.exists(paths)
		if (any(existing)) unlink(paths[existing], force = TRUE)
	}
}, onexit = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Constants and phenotype mapping
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.dat <- file.path(dir0, "data", "ems120")
years <- 2013:2024
yrs_trend <- 2017:2024
phone.grp.use <- c("low", "high")
phone_cols <- c(low = "grey55", high = "#D55E00")

phone.score.use <- "sco0" #📍 Options: "sco0", "sco1", "sco2", "sco3".
phone.luck.use <- "phone.luck.sco0" #📍 Options: "phone.luck.quarter", "phone.luck.sco0", "phone.luck.sco1", "phone.luck.sco2", "phone.luck.sco3".
geo.type.use <- "geo.type1" #📍 Options: "地址类型", "geo.type1"
geo.map.use <- "keplergl_3d" #📍 Options: "deckgl_satellite", "deckgl_hex_satellite", "keplergl_3d", "none".
geo.map.transparent <- 0  #📍
house_price.use <-  "住宅区" # 📍 Options: "住宅区", "工作区", NA
house_price_min <- 10^4.25 # 📍 Options
house_price_scale <- "log10" #📍 Options: "log" or "log10".
house_price_bar <- "sqrt" #📍 Fig2c/Fig2d bar height uses sqrt-capped housing price, matching the working map code.
dx.type.use <- "dx.type1" # 📍 Options: "dx.type0" or "dx.type1".
dxs.vip <- c("Violence", "Fall", "Traffic", "Intoxication", "CVD", "Respiratory", "Psychiatric", "Death") # 📍
call_volume <- "adjusted" # raw
DX_PIPELINE_VERSION <- "macbert_dx_v3_fieldtag_cv5_kw4_confidence_20260809"
if (grepl("^phone\\.luck\\.sco[0-3]$", phone.luck.use)) phone.score.use <- sub("^phone\\.luck\\.", "", phone.luck.use)
geo.type.use <- match.arg(geo.type.use, c("地址类型", "geo.type1"))
geo.map.use <- match.arg(geo.map.use, c("deckgl_satellite", "deckgl_hex_satellite", "keplergl_3d", "none"))
geo.map.transparent <- suppressWarnings(as.numeric(geo.map.transparent)[1])
if (!is.finite(geo.map.transparent)) stop("geo.map.transparent must be numeric between 0 and 1.", call. = FALSE)
geo.map.transparent <- pmin(pmax(geo.map.transparent, 0), 1)
# Fig2 map rendering is configured by geo.map.use; these flags are internal.
draw_2D_map <- !identical(geo.map.use, "none")  # clean coordinate map used in submitted Fig2
draw_satellite_map <- FALSE  # avoid web-tile/Chrome failures; HTML satellite map is not required for submission
draw_3D_map <- FALSE  # interactive 3D map is not part of the submission figure set
if (length(house_price.use) != 1L) stop("house_price.use must be one of: 住宅区, 工作区, or NA.", call. = FALSE)
if (is.character(house_price.use) && identical(toupper(trimws(house_price.use)), "NA")) house_price.use <- NA_character_
if (!is.na(house_price.use)) house_price.use <- match.arg(house_price.use, c("住宅区", "工作区"))
dx.type.use <- match.arg(dx.type.use, c("dx.type0", "dx.type1"))
call_volume <- match.arg(call_volume, c("adjusted", "raw"))
house_price_bar <- match.arg(house_price_bar, c("sqrt", "log", "log10"))
fmt_p_compact <- function(p, sig_threshold = 0.05, nonsig_digits = 2) {
	p <- suppressWarnings(as.numeric(p))
	vapply(p, function(x) {
		if (!is.finite(x)) return("NA")
		if (x <= .Machine$double.xmin) return("<2.23E-308")
		if (x < sig_threshold) return(formatC(x, format = "E", digits = 2))
		formatC(x, format = "f", digits = nonsig_digits)
	}, character(1))
}
fmt_p_x_math <- function(p, digits = 2) {
	p <- suppressWarnings(as.numeric(p))
	vapply(p, function(x) {
		if (!is.finite(x)) return("NA")
		if (x <= .Machine$double.xmin) return("'<2.23'~'x'~10^{-308}")
		if (x > 0 && x < 0.001) {
			e <- floor(log10(x)); m <- x / 10^e
			return(sprintf(paste0("%.", digits, "f~'x'~10^{%d}"), m, as.integer(e)))
		}
		sprintf("%.2f", x)
	}, character(1))
}

cache_complete <- function(dat_list, required_cols, year_set = years) {
	all(vapply(as.character(year_set), function(y) {
		if (!y %in% names(dat_list)) return(FALSE)
		d <- dat_list[[y]]
		is.null(d) || all(required_cols %in% names(d))
	}, logical(1)))
}
format_count_table <- function(x, levels = NULL) {
	if (is.null(x)) return("missing")
	if (!length(x)) return("empty")
	if (!is.null(levels)) x <- factor(as.character(x), levels = levels)
	tab <- table(x, useNA = "ifany")
	nm <- names(tab)
	nm[is.na(nm) | nm == "" | nm == "<NA>"] <- "NA"
	paste(sprintf("%s %s", nm, as.integer(tab)), collapse = ", ")
}
format_quintile_summary <- function(x) {
	x <- suppressWarnings(as.numeric(x))
	ok <- is.finite(x)
	na_n <- sum(!ok)
	if (!any(ok)) return(sprintf("quintile unavailable; NA %d", na_n))
	qs <- stats::quantile(x[ok], probs = seq(0, 1, 0.2), na.rm = TRUE, type = 7)
	q_txt <- paste(sprintf("%s=%s", names(qs), formatC(qs, format = "fg", digits = 4)), collapse = ", ")
	bin <- paste0("Q", dplyr::ntile(x[ok], 5))
	bin_txt <- format_count_table(factor(bin, levels = paste0("Q", 1:5)))
	paste0("quintile cutpoints ", q_txt, "; bins ", bin_txt, "; NA ", na_n)
}
vars.basic.ems <- c("电话", "地址", "地址类型", "开始受理时刻", "派车时间", "去程时间", "现场时间", "返程时间", "急救时间", "疾病类型", "接车地址经度", "接车地址纬度")
vars.basic.alias <- c("^病人电话号码|^联系电话.1|^联系电话", "^接车地址$|^接车地点$|^现场地址$", "地址类型", "^开始受理时刻|^开始时刻|^摘机时刻|^收到指令时刻", "^派车时间|^受理调度时间", "^去程时间|^去程在途时间", "^现场时间|^现场救援时间|^现场治疗时间|^现场急救时间", "^返程时间|^返程在途时间", "^急救时间|^急救反应时间", "疾病类型", "接车地址经度", "接车地址纬度")
vars.phone <- setNames(c(rep("病人电话号码", 7), rep("联系电话.1", 3), rep("联系电话", 2)), as.character(years))

vars.dxs <- c("性别", "年龄", "呼救原因", "病种判断", "病因", "伤病程度", "症状", "主诉", "病史", "初步诊断", "补充诊断")
vars.dxs.alias <- c("^性别$|病人性别|患者性别", "^年龄$|病人年龄|患者年龄", "^呼救原因|^呼叫原因", "病种判断", "^病因$|辅助诊断", "伤病程度|病情分级", "^症状$|患者症状", "^主诉$|病情\\\\(主诉\\\\)", "^病史$|现病史", "^初步诊断$", "初步诊断2|补充诊断")
vars <- c(vars.basic.ems, vars.dxs); vars.alias <- c(vars.basic.alias, vars.dxs.alias); vars.geo <- c("ID", "地址", "地址类型"); vars.time <- c("开始受理时刻", "驶向现场时刻", "到达现场时刻", "病人上车时刻", "到达医院时刻")
ems120_text_vars <- unique(c(
	vars.basic.ems[c(1, 2, 3, 10)], vars.dxs[-2],
	"geo.type", "geo.type1",
	"dx.type0", "dx.type0.reason", "dx.type1", "dx.type1.reason", "dx_raw",
	"phone.luck", "phone.luck.quarter",
	"phone.sco.reason", "phone.sco0.reason", "phone.sco1.reason", "phone.sco2.reason", "phone.sco3.reason"
))
normalize_ems120_year_types <- function(d) {
	if (is.null(d) || !nrow(d)) return(d)
	d %>% mutate(across(any_of(ems120_text_vars), as.character))
}

dxs.type.list <- list(
	Violence = "创伤-暴力事件",
	Fall = "创伤-跌倒",
	Traffic = "创伤-交通事故",
	Trauma = c("创伤-其他原因", "创伤-高处坠落"),
	Respiratory = "呼吸系统疾病",
	Psychiatric = "精神病",
	Intoxication = "理化中毒",
	Urinary = "泌尿系统疾病",
	Endocrine = "内分泌系统疾病",
	Death = "其他-死亡",
	Neurological = c("神经系统疾病-脑卒中", "神经系统疾病-其他疾病"),
	Digestive = "消化系统疾病",
	CVD = c("心血管系统疾病-其他疾病", "心血管系统疾病-胸痛", "其他-胸闷"),
	"Ob/Gyn" = "妇产科",
	Pediatric = "儿科",
	Other = c("其他-其他症状", "其他-昏迷")
)
dxs.all0 <- names(dxs.type.list)
dxs.all <- setdiff(dxs.all0, "Other")
dxs.all.color <- setNames(c(
	"#E64B35", "#F39B7F", "#D62728", "#B07AA1",
	"#2ECED0", "#2F55D4", "#D9A51E",
	"#4E79A7", "#59A14F", "#D62DB5",
	"#8CD17D", "#00A087", "#28C85A",
	"#FF9DA7", "#76B7B2", "#9D9D9D"
), dxs.all0)
dxs.vip <- intersect(dxs.vip, dxs.all)
dxs.vip.color <- dxs.all.color[dxs.vip]
map_grp <- stack(dxs.type.list) %>% setNames(c("dx_raw", "dx_grp")) %>% mutate(dx_raw = trimws(as.character(dx_raw)), dx_grp = trimws(as.character(dx_grp)))
apply_dx_type_use <- function(d) {
	if (is.null(d) || !nrow(d)) return(d)
	if (!dx.type.use %in% names(d)) stop(sprintf("Selected dx.type.use column is missing: %s", dx.type.use), call. = FALSE)
	d %>%
		mutate(dx_raw = trimws(as.character(.data[[dx.type.use]]))) %>%
		dplyr::select(-any_of("dx_grp")) %>%
		left_join(map_grp, by = "dx_raw") %>%
		mutate(dx_grp = factor(as.character(dx_grp), levels = dxs.all0))
}
write_dat_summary_log <- function(dat_list, cache_file, file = "dat.summary.log", year_set = years) {
	lines <- c(
		"EMS120 processed data summary",
		sprintf("source cache: %s", normalizePath(cache_file, winslash = "/", mustWork = FALSE)),
		sprintf("phone.score.use = %s", phone.score.use),
		sprintf("phone.luck.use = %s", phone.luck.use),
		sprintf("dx.type.use = %s", dx.type.use),
		""
	)
	score_vars <- c("phone.sco0", "phone.sco1", "phone.sco2", "phone.sco3")
	for (year in year_set) {
		key <- as.character(year)
		d <- dat_list[[key]]
		lines <- c(lines, sprintf("Year %s:", key))
		if (is.null(d) || !nrow(d)) {
			lines <- c(lines, "  no rows", "")
			next
		}
		missing_fig2_vars <- setdiff(c("房价指数", "1km渔网建筑密度", "建筑高度", "到最近湖泊距离", "到最近水系距离", "接车地址经度", "接车地址纬度"), names(d))
		if (length(missing_fig2_vars)) {
			for (v_miss in missing_fig2_vars) lines <- c(lines, sprintf("❌ 缺少 %s 变量", v_miss))
		} else {
			lines <- c(lines, "✅ Fig2 房价/坐标变量存在")
			for (v_ok in c("房价指数", "接车地址经度", "接车地址纬度")) {
				x_ok <- suppressWarnings(as.numeric(d[[v_ok]]))
				lines <- c(lines, sprintf("%s: finite %s / %s", v_ok, sum(is.finite(x_ok), na.rm = TRUE), length(x_ok)))
			}
		}
		lines <- c(lines, sprintf("geo.type1: %s", format_count_table(if ("geo.type1" %in% names(d)) d$geo.type1 else NULL)))
		for (v in score_vars) {
			lines <- c(lines, sprintf("%s: %s", v, if (v %in% names(d)) format_quintile_summary(d[[v]]) else "missing"))
		}
		lines <- c(
			lines,
			sprintf("phone.luck (based on %s): %s", phone.luck.use, format_count_table(if ("phone.luck" %in% names(d)) d$phone.luck else NULL, levels = c("low", "middle", "high"))),
			sprintf("dx.type0: %s", format_count_table(if ("dx.type0" %in% names(d)) d$dx.type0 else NULL)),
			sprintf("dx.type1: %s", format_count_table(if ("dx.type1" %in% names(d)) d$dx.type1 else NULL)),
			""
		)
	}
	writeLines(lines, file, useBytes = TRUE)
	invisible(lines)
}


# English display labels used only for plots; raw Chinese values are kept in the data.
vars.dx.eng <- c(
	setNames(dxs.all0, dxs.all0),
	setNames(rep(names(dxs.type.list), lengths(dxs.type.list)), unlist(dxs.type.list, use.names = FALSE))
)
vars.geo.eng <- c("Housing price" = "房价指数", "Building density" = "1km渔网建筑密度", "Building height" = "建筑高度", "Distance to nearest lake" = "到最近湖泊距离", "Distance to nearest water body" = "到最近水系距离")
to_eng <- function(x, dict) { x0 <- as.character(x); hit <- match(x0, unname(dict)); ifelse(!is.na(hit), names(dict)[hit], x0) }
dx_to_eng <- function(x) to_eng(x, vars.dx.eng)
dx_raw_to_type_label <- function(x) {
	x0 <- trimws(as.character(x))
	hit <- match(x0, map_grp$dx_raw)
	out <- x0
	ok <- !is.na(hit)
	out[ok] <- dx_to_eng(map_grp$dx_grp[hit[ok]])
	out
}
geo_to_eng <- function(x) to_eng(x, vars.geo.eng)

# Required variables that must be present in dat.list.rds for Fig2.
# These should be created by the Xia geo merge step and carried into dat.list.rds.
fig2_required_geo_vars <- c("房价指数", "1km渔网建筑密度", "建筑高度", "到最近湖泊距离", "到最近水系距离")
fig2_required_coord_vars <- c("接车地址经度", "接车地址纬度")
fig2_required_vars <- c(fig2_required_geo_vars, fig2_required_coord_vars)
fig2_num_vars <- fig2_required_geo_vars
fig2_coord_vars <- fig2_required_coord_vars

check_required_vars_in_dat_list <- function(dat_list, required_vars = fig2_required_vars, year_set = years, stop_on_missing = FALSE, context = "dat.list.rds") {
	miss_tbl <- purrr::imap_dfr(dat_list[as.character(year_set)], function(d, y) {
		if (is.null(d) || !nrow(d)) {
			return(tibble(year = as.integer(y), variable = required_vars, status = "no rows"))
		}
		missing_vars <- setdiff(required_vars, names(d))
		if (!length(missing_vars)) return(tibble())
		tibble(year = as.integer(y), variable = missing_vars, status = "missing")
	})
	if (nrow(miss_tbl) && isTRUE(stop_on_missing)) {
		msg <- miss_tbl %>%
			group_by(year) %>%
			summarise(missing = paste(variable, collapse = ", "), .groups = "drop") %>%
			mutate(line = sprintf("Year %s: %s", year, missing)) %>%
			pull(line) %>%
			paste(collapse = "\n")
		stop(sprintf("%s is missing required Fig2 variables. Rebuild dat.list.rds after Xia geo merge.\n%s", context, msg), call. = FALSE)
	}
	invisible(miss_tbl)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Data harmonization, ml_phone, ml_geo, and ml_dx
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (file.exists("dat0.summary.log") && file.info("dat0.summary.log")$size > 0) {
	cat("dat0.summary.log exists; skipping phone-variable audit.\n")
} else {
	phone_var_patterns <- strsplit(vars.basic.alias[1], "\\|")[[1]]
	con <- file("dat0.summary.log", open = "wt", encoding = "UTF-8"); sink(con, split = TRUE)
	tryCatch({ phone_audit <- lapply(years, function(year) audit_phone_vars(file.path(dir.dat, "clean", paste0(year, ".xlsx")), year, phone_var_patterns)) }, finally = { sink(); close(con) })
}
cat("vars.phone default:\n"); for (yy in names(vars.phone)) cat(sprintf("  %s = %s\n", yy, vars.phone[[yy]])); cat(sprintf("phone.score.use = %s; phone.luck.use = %s; dx.type.use = %s\n", phone.score.use, phone.luck.use, dx.type.use))
cat(c("vars.phone default:", sprintf("  %s = %s", names(vars.phone), unname(vars.phone)), sprintf("phone.score.use = %s", phone.score.use), sprintf("phone.luck.use = %s", phone.luck.use), sprintf("dx.type.use = %s", dx.type.use)), file = "dat0.summary.log", sep = "\n", append = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 1) Raw yearly data harmonization.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file.dat0.rds <- file.path(dir.cache, "dat0.list.rds")
if (file.exists(file.dat0.rds)) {
	dat0.list <- readRDS(file.dat0.rds)
} else {
	dat0.list <- list()
	for (year in years) {
		cat(sprintf("========== %d ==========", year), "\n"); key <- as.character(year)
		dat <- suppressMessages(readxl::read_excel(file.path(dir.dat, "clean", paste0(year, ".xlsx")))); names(dat) <- trimws(names(dat))
		col_names <- names(dat); new_names <- col_names
		for (i in seq_along(vars)) {
			all_m <- grep(vars.alias[i], col_names, value = TRUE); if (!length(all_m)) next
			best_m <- all_m[1]
			if (length(all_m) > 1) for (p in strsplit(vars.alias[i], "\\|")[[1]]) { m <- grep(p, all_m, value = TRUE); if (length(m)) { best_m <- if (vars[i] == "电话" && !is.null(vars.phone[[key]]) && vars.phone[[key]] %in% all_m) vars.phone[[key]] else m[1]; break } }
			new_names[match(best_m, col_names)] <- vars[i]
		}
		names(dat) <- new_names
		dat <- dat %>% mutate(across(any_of(vars.time), ~ lubridate::ymd_hms(trimws(.x))),
			!!vars.basic.ems[5] := if (vars.basic.ems[5] %in% names(.)) as.numeric(.data[[vars.basic.ems[5]]]) else if (all(c(vars.time[2], vars.time[1]) %in% names(.))) as.numeric(.data[[vars.time[2]]] - .data[[vars.time[1]]], units = "secs") else NA_real_,
			!!vars.basic.ems[6] := if (vars.basic.ems[6] %in% names(.)) as.numeric(.data[[vars.basic.ems[6]]]) else if (all(c(vars.time[3], vars.time[2]) %in% names(.))) as.numeric(.data[[vars.time[3]]] - .data[[vars.time[2]]], units = "secs") else NA_real_,
			!!vars.basic.ems[7] := if (vars.basic.ems[7] %in% names(.)) as.numeric(.data[[vars.basic.ems[7]]]) else if (all(c(vars.time[4], vars.time[3]) %in% names(.))) as.numeric(.data[[vars.time[4]]] - .data[[vars.time[3]]], units = "secs") else NA_real_,
			!!vars.basic.ems[8] := if (vars.basic.ems[8] %in% names(.)) as.numeric(.data[[vars.basic.ems[8]]]) else if (all(c(vars.time[5], vars.time[4]) %in% names(.))) as.numeric(.data[[vars.time[5]]] - .data[[vars.time[4]]], units = "secs") else NA_real_,
			!!vars.basic.ems[9] := if (vars.basic.ems[9] %in% names(.)) as.numeric(.data[[vars.basic.ems[9]]]) else if (all(vars.basic.ems[5:8] %in% names(.))) .data[[vars.basic.ems[5]]] + .data[[vars.basic.ems[6]]] + .data[[vars.basic.ems[7]]] + .data[[vars.basic.ems[8]]] else NA_real_)
		missing_vars <- setdiff(vars, names(dat)); if (length(missing_vars)) { cat(sprintf("-> Warning: %d missing vars: %s\n", year, paste(missing_vars, collapse = ", "))); for (v in missing_vars) dat[[v]] <- NA }
		dat0.list[[key]] <- dat %>% dplyr::select(all_of(vars)) %>% mutate(!!vars.basic.ems[1] := clean_phone_value(.data[[vars.basic.ems[1]]]), phone = ifelse(is.na(.data[[vars.basic.ems[1]]]), NA_character_, substring(.data[[vars.basic.ems[1]]], 4, 11)), !!vars.dxs[2] := suppressWarnings(as.numeric(gsub("岁", "", .data[[vars.dxs[2]]])))) %>% mutate(ID = row_number(), .before = 1)
	}
	saveRDS(dat0.list, file.dat0.rds)
}
dat0.list <- lapply(dat0.list, normalize_ems120_year_types)
file.train2019 <- file.path(dir.cache, "2019.train_dx.xlsx")
if (!file.exists(file.train2019) && "2019" %in% names(dat0.list)) {
	train2019 <- dat0.list[["2019"]][
		1:min(10000, nrow(dat0.list[["2019"]])),
		intersect(
			c("ID", "电话", "phone", vars.dxs, "地址", "地址类型", "疾病类型"),
			names(dat0.list[["2019"]])
		),
		drop = FALSE
	]

	train2019_phone <- score_phone_simple(train2019$phone)
	train2019$phone.sco <- train2019_phone$phone.sco
	train2019$phone.sco.reason <- train2019_phone$phone.sco.reason

	writexl::write_xlsx(train2019, file.train2019)

	cat(
		"Created 2019.train_dx.xlsx from the first",
		nrow(train2019),
		"expert-annotated rows of the 2019 EMS file.\n"
	)
}

train_check <- suppressMessages(readxl::read_excel(file.train2019))
cat(sprintf(
	"Training workbook: %s rows; %s raw disease classes.\n",
	nrow(train_check),
	dplyr::n_distinct(train_check$疾病类型)
))
rm(train_check)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 2) phone: score phone number luckiness; keep only RDS cache, not dat.ml_phone/[year].xlsx.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file.dat.phone.rds <- file.path(dir.cache, "dat.ml_phone.rds")
phone.required <- c("phone.sco0", "phone.sco1", "phone.sco2", "phone.sco3", "phone.sco0.reason", "phone.sco1.reason", "phone.sco2.reason", "phone.sco3.reason", "phone.luck.sco0", "phone.luck.sco1", "phone.luck.sco2", "phone.luck.sco3")
phone.ready <- FALSE
if (file.exists(file.dat.phone.rds)) {
	dat.phone.list <- readRDS(file.dat.phone.rds)
	dat.phone.list <- lapply(dat.phone.list, apply_phone_luck_use)
	phone.ready <- TRUE
	cat("dat.ml_phone.rds exists; skipping ml_phone:", file.dat.phone.rds, "\n")
}
if (!phone.ready) {
	ensure_ml_python()
	dat.phone.list <- list()
	for (year in years) {
		key <- as.character(year)
		dat <- dat0.list[[key]] %>% filter(!is.na(电话)); if (!nrow(dat)) { dat.phone.list[[key]] <- NULL; next }
		cat("Processing ml_phone year:", year, "\n"); res_phone <- ml_phone(dat$phone); res0 <- t(vapply(dat$phone, score_phone_sco0_one, character(2)))
		dat.phone.list[[key]] <- dat %>% mutate(phone.sco0 = as.numeric(res0[, "score"]), phone.sco1 = as.numeric(unlist(res_phone$sco1)), phone.sco2 = as.numeric(unlist(res_phone$sco2)), phone.sco3 = as.numeric(unlist(res_phone$sco3)), phone.sco0.reason = as.character(res0[, "reason"]), phone.sco1.reason = as.character(unlist(res_phone$reason1)), phone.sco2.reason = as.character(unlist(res_phone$reason2)), phone.sco3.reason = as.character(unlist(res_phone$reason3))) %>% apply_phone_luck_use()
		print_phone_score_bins(dat.phone.list[[key]], key, phone.required[1:4])
	}
	saveRDS(dat.phone.list, file.dat.phone.rds)
}
dat.phone.list <- lapply(dat.phone.list, normalize_ems120_year_types)
phone_score_vars <- c("phone.sco0", "phone.sco1", "phone.sco2", "phone.sco3")
phone_score_summary <- purrr::imap_dfr(dat.phone.list, function(d, y) { if (is.null(d) || !nrow(d)) return(tibble()); bind_rows(lapply(phone_score_vars, function(v) as_tibble(table(score = d[[v]], useNA = "always")) %>% transmute(Year = as.integer(y), score_name = v, score = as.character(score), n = as.integer(n)))) })
phone_score_bin_summary <- purrr::imap_dfr(dat.phone.list, function(d, y) { if (is.null(d) || !nrow(d)) return(tibble()); bind_rows(lapply(phone_score_vars, function(v) tibble(Year = as.integer(y), score_name = v, score_bin = score_bin_label(d[[v]])) %>% count(Year, score_name, score_bin, name = "n"))) })
phone_quarter_summary <- bind_rows(lapply(phone_score_vars, function(v) { x <- unlist(lapply(dat.phone.list, function(d) if (is.null(d) || !nrow(d)) numeric(0) else suppressWarnings(as.numeric(d[[v]]))), use.names = FALSE); q <- quasi_quarter(x, .25, .5, .25); tibble(score_name = v, low_cutoff = ifelse(any(q == "low", na.rm = TRUE), max(x[q == "low"], na.rm = TRUE), NA_real_), low_n = sum(q == "low", na.rm = TRUE), high_cutoff = ifelse(any(q == "high", na.rm = TRUE), min(x[q == "high"], na.rm = TRUE), NA_real_), high_n = sum(q == "high", na.rm = TRUE)) }))
writexl::write_xlsx(list(score_distribution = phone_score_summary, score_distribution_integer_bins = phone_score_bin_summary, quasi_quarter_cutoffs = phone_quarter_summary), "phone.sco.summary.xlsx")

# dat.phone.list supersedes dat0.list for all downstream stages.  Keeping both
# multi-year objects can exceed the WSL memory limit while Python trains a model.
rm(dat0.list)
invisible(gc())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 3) geo: classify address type after phone scoring, then merge yearly [year].geo.xlsx files.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file.dat.geo.rds <- file.path(dir.cache, "dat.ml_geo.rds")
geo.ready <- FALSE
if (file.exists(file.dat.geo.rds)) {
	dat.geo.list <- readRDS(file.dat.geo.rds)
	geo.ready <- TRUE
	cat("dat.ml_geo.rds exists; skipping ml_geo and Xia geo merge:", file.dat.geo.rds, "\n")
}
if (!geo.ready) {
	ensure_ml_python()
	dat.geo.list <- list()
	for (year in years) {
		key <- as.character(year)
		dat <- dat.phone.list[[key]]; if (is.null(dat) || !nrow(dat)) { dat.geo.list[[key]] <- NULL; next }
		cat("Processing ml_geo year:", year, "\n"); res_geo <- ml_geo(dat %>% dplyr::select(any_of(vars.geo[-1])), batch_size = 1024)
		dat.geo.list[[key]] <- dat %>% mutate(geo.type1 = as.character(unlist(res_geo)))
	}
	dir.geo.xia <- file.path(dir.dat, "szu_xia")
	file.geo.log <- "geo_xia_merge.log"
	geo_xia_results <- purrr::imap(dat.geo.list[as.character(years)], function(d, y) merge_xia_geo_year(d, y, dir.geo.xia, vars.basic.ems = vars.basic.ems))
	dat.geo.list <- purrr::map(geo_xia_results, "data") %>% lapply(apply_phone_luck_use)
	geo_xia_merge_log <- purrr::map_dfr(geo_xia_results, "log")
	saveRDS(dat.geo.list, file.dat.geo.rds)
	data.table::fwrite(geo_xia_merge_log, file.geo.log, sep = "\t", na = "")
} else {
	file.geo.log <- "geo_xia_merge.log"
	geo_xia_merge_log <- if (file.exists(file.geo.log)) data.table::fread(file.geo.log, sep = "\t", data.table = FALSE) %>% tibble::as_tibble() else tibble::tibble()
}
dat.geo.list <- lapply(dat.geo.list, normalize_ems120_year_types)

if (isTRUE(ml_python_ready)) release_geo_model()

# dat.geo.list contains the phone columns needed downstream; release the earlier
# full copy before CV/model training and 1.8-million-record disease inference.
rm(dat.phone.list)
invisible(gc())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 4) dx CV: stratified 5-fold validation immediately before ml_dx.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reuse an existing valid OOF run whenever possible.  Reviewer-ready bootstrap
# confidence intervals and confidence-calibration summaries can be regenerated
# from the saved OOF predictions without retraining the five MacBERT folds.
run_dx_cv <- !tolower(Sys.getenv("EMS120_RUN_DX_CV", unset = "1")) %in% c("0", "false", "no", "off")
file.dx.cv.rds <- file.path(dir.cache, "dat.ml_dx_cv.rds")
dir_dx_cv <- file.path(dir.analysis, "raw", "hfl_cv5")
# Backward compatibility: some earlier runs wrote cv_*.csv directly under raw/.
# Copy those files once into the versioned hfl_cv5 directory so the current
# Python and R helpers share one location without forcing an expensive rerun.
dir_dx_cv_legacy <- file.path(dir.analysis, "raw")
dx_cv_legacy_names <- c("cv_summary.csv", "cv_fold_metrics.csv", "cv_per_class_raw.csv", "cv_per_class_grouped.csv", "cv_confusion_raw.csv", "cv_confusion_grouped.csv", "cv_oof_predictions.csv", "cv_manifest.json")
dx_cv_nested_min <- file.path(dir_dx_cv, dx_cv_legacy_names[c(1, 2, 4, 5, 6, 7)])
dx_cv_legacy_min <- file.path(dir_dx_cv_legacy, dx_cv_legacy_names[c(1, 2, 4, 5, 6, 7)])
if (!all(file.exists(dx_cv_nested_min)) && all(file.exists(dx_cv_legacy_min))) {
	dir.create(dir_dx_cv, recursive = TRUE, showWarnings = FALSE)
	legacy_present <- file.path(dir_dx_cv_legacy, dx_cv_legacy_names[file.exists(file.path(dir_dx_cv_legacy, dx_cv_legacy_names))])
	invisible(file.copy(legacy_present, dir_dx_cv, overwrite = FALSE))
	cat("Copied legacy raw/cv_* outputs into raw/hfl_cv5 for reuse.\n")
}
dx_cv_base_files <- file.path(dir_dx_cv, c(
	"cv_summary.csv", "cv_fold_metrics.csv", "cv_per_class_raw.csv", "cv_per_class_grouped.csv",
	"cv_confusion_raw.csv", "cv_confusion_grouped.csv", "cv_oof_predictions.csv"
))
dx_cv_derived_files <- file.path(dir_dx_cv, c("cv_bootstrap_metrics.csv", "cv_confidence_calibration.csv"))

if (run_dx_cv) {
	if (all(file.exists(dx_cv_base_files))) {
		# Existing OOF predictions are the expensive part.  Refresh only derived
		# uncertainty/calibration outputs if this is the first run of the v3 pipeline.
		if (!all(file.exists(dx_cv_derived_files))) {
			ensure_ml_python()
			cat("Existing MacBERT OOF predictions found; refreshing bootstrap/calibration outputs without retraining folds.\n")
			refresh_dx_cv_derived_outputs()
		}
		dx_cv_files <- c(
			summary = "cv_summary.csv", fold_metrics = "cv_fold_metrics.csv",
			per_class_raw = "cv_per_class_raw.csv", per_class_grouped = "cv_per_class_grouped.csv",
			confusion_raw = "cv_confusion_raw.csv", confusion_grouped = "cv_confusion_grouped.csv",
			oof_predictions = "cv_oof_predictions.csv", bootstrap_metrics = "cv_bootstrap_metrics.csv",
			confidence_calibration = "cv_confidence_calibration.csv"
		)
		dx_cv_tables <- lapply(file.path(dir_dx_cv, dx_cv_files), function(f) {
			if (file.exists(f)) data.table::fread(f, data.table = FALSE) else NULL
		})
		names(dx_cv_tables) <- names(dx_cv_files)
		dx_cv_info <- list(status = "cached_oof", summary = file.path(dir_dx_cv, "cv_summary.csv"), dir = dir_dx_cv, tables = dx_cv_tables, generated_at = Sys.time())
		saveRDS(dx_cv_info, file.dx.cv.rds)
		cat("Reused stratified 5-fold OOF validation:", dir_dx_cv, "\n")
	} else {
		ensure_ml_python()
		cat("Complete OOF validation outputs are absent; running stratified 5-fold MacBERT validation.\n")
		dx_cv_run <- run_dx_stratified_cv(force = TRUE)
		dir_dx_cv_run <- as.character(dx_cv_run$dir)
		dx_cv_files <- c(
			summary = "cv_summary.csv", fold_metrics = "cv_fold_metrics.csv",
			per_class_raw = "cv_per_class_raw.csv", per_class_grouped = "cv_per_class_grouped.csv",
			confusion_raw = "cv_confusion_raw.csv", confusion_grouped = "cv_confusion_grouped.csv",
			oof_predictions = "cv_oof_predictions.csv", bootstrap_metrics = "cv_bootstrap_metrics.csv",
			confidence_calibration = "cv_confidence_calibration.csv"
		)
		dx_cv_tables <- lapply(file.path(dir_dx_cv_run, dx_cv_files), function(f) {
			if (file.exists(f)) data.table::fread(f, data.table = FALSE) else NULL
		})
		names(dx_cv_tables) <- names(dx_cv_files)
		dx_cv_info <- c(dx_cv_run, list(tables = dx_cv_tables, generated_at = Sys.time()))
		saveRDS(dx_cv_info, file.dx.cv.rds)
	}
	cat("MacBERT CV status:", as.character(dx_cv_info$status), "; summary:", as.character(dx_cv_info$summary), "; directory:", as.character(dx_cv_info$dir), "\n")
} else {
	cat("EMS120_RUN_DX_CV=0; skipping CV execution. Existing raw/hfl_cv5 files remain available for figure redraws.\n")
}


# 🚩 5) dx: classify disease after phone, geo, and dx CV; keep only RDS cache.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file.dat.rds <- file.path(dir.cache, "dat.list.rds")
final.required <- c(
	"geo.type1", "dx.type0", "dx.type0.reason", "dx.type1", "dx.type1.reason",
	"dx.type1.confidence", "dx.pipeline.version", "phone.repeated_gt5",
	"phone.sco", "phone.sco.reason", phone.required, "phone.luck", "phone.luck.quarter"
)
final.ready <- FALSE
if (file.exists(file.dat.rds)) {
	dat1.list <- readRDS(file.dat.rds)
	dat1.list <- lapply(dat1.list, apply_phone_luck_use)
	cache_cols_ok <- cache_complete(dat1.list, final.required)
	cache_version_ok <- all(vapply(dat1.list, function(d) {
		is.null(d) || !nrow(d) || ("dx.pipeline.version" %in% names(d) && all(as.character(d$dx.pipeline.version) == DX_PIPELINE_VERSION, na.rm = TRUE))
	}, logical(1)))
	if (isTRUE(cache_cols_ok) && isTRUE(cache_version_ok)) {
		final.ready <- TRUE
		cat(basename(file.dat.rds), "is current; skipping ml_dx:", file.dat.rds, "\n")
	} else {
		cat(basename(file.dat.rds), "is stale/incomplete; rebuilding ml_dx while reusing phone/geo caches.\n")
		rm(dat1.list)
		invisible(gc())
	}
}
if (!final.ready) {
	ensure_ml_python()
	dat1.list <- list()
	for (year in years) {
		key <- as.character(year)
		file.year.xlsx <- file.path(dir.cache, paste0(key, ".xlsx"))
		dat <- dat.geo.list[[key]]; if (is.null(dat) || !nrow(dat)) { dat1.list[[key]] <- NULL; next }
		# Count repeated telephone values in the complete harmonized year BEFORE the
		# disease-analysis eligibility filter.  This makes the >5/year sensitivity
		# definition correspond to the source-year archive rather than only to the
		# final analytic subset.
		dat <- dat %>%
			add_count(.data[[vars.basic.ems[1]]], name = "n_tel_year") %>%
			mutate(phone.repeated_gt5 = !is.na(.data[[vars.basic.ems[1]]]) & n_tel_year > 5) %>%
			filter(!is.na(电话), !is.na(性别), !if_all(c(主诉, 病史, 初步诊断, 补充诊断), is.na))
		cat("Processing ml_dx year:", year, "\n"); res_dx <- ml_dx(dat %>% dplyr::select(any_of(vars.dxs)), data_name = key, batch_size = 1024)
		dat1.list[[key]] <- dat %>%
			mutate(
				dx.type0 = trimws(as.character(unlist(res_dx$kw))),
				dx.type0.reason = as.character(unlist(res_dx$kw_reason)),
				dx.type1 = trimws(as.character(unlist(res_dx$ml))),
				dx.type1.reason = as.character(unlist(res_dx$ml_reason)),
				dx.type1.confidence = suppressWarnings(as.numeric(unlist(res_dx$ml_confidence))),
				dx.pipeline.version = DX_PIPELINE_VERSION
			) %>%
			mutate(
				!!vars.basic.ems[1] := ifelse(phone.repeated_gt5, NA_character_, .data[[vars.basic.ems[1]]]),
				.time0 = .data[[vars.time[1]]], .date0 = as.Date(.time0), hour = lubridate::hour(.time0)
			) %>%
			dplyr::select(-n_tel_year, -.time0) %>%
			rename(!!"日期" := .date0) %>%
			apply_phone_luck_use() %>% apply_dx_type_use()
		writexl::write_xlsx(dat1.list[[key]], file.year.xlsx)
	}
	saveRDS(dat1.list, file.dat.rds)
}
dat1.list <- lapply(dat1.list, function(d) apply_dx_type_use(apply_phone_luck_use(normalize_ems120_year_types(d))))
write_dat_summary_log(dat1.list, file.dat.rds)
# Final data files: only dat/[year].xlsx plus RDS caches are written.
invisible(purrr::iwalk(dat1.list[as.character(years)], function(d, y) if (!is.null(d) && nrow(d)) writexl::write_xlsx(d, file.path(dir.cache, paste0(y, ".xlsx")))))
cat("Phenotype counts by year:\n"); print(lapply(dat1.list, function(dat) if (is.null(dat)) NA else table(dat$dx_grp, useNA = "ifany")))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 TabS1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (ems120_should_run("fig1")) {
the_table1 <- purrr::imap_dfr(dat1.list, function(d, y) {
	if (is.null(d) || !nrow(d)) return(tibble())
	tibble(
		Year = as.integer(y),
		`EMS calls, N` = format(nrow(d), big.mark = ","),
		`Unique caller numbers, N` = format(n_distinct(d[[vars.basic.ems[1]]], na.rm = TRUE), big.mark = ","),
		`Unique / total, %` = sprintf("%.2f%%", 100 * n_distinct(d[[vars.basic.ems[1]]], na.rm = TRUE) / nrow(d)),
		`Age, mean (SD)` = sprintf("%.1f (%.2f)", mean(d[[vars.dxs[2]]], na.rm = TRUE), sd(d[[vars.dxs[2]]], na.rm = TRUE)),
		`Female, n (%)` = sprintf("%s (%.2f%%)", format(sum(d[[vars.dxs[1]]] == "女", na.rm = TRUE), big.mark = ","), 100 * mean(d[[vars.dxs[1]]] == "女", na.rm = TRUE)),
		`Low-luck, n (%)` = sprintf("%s (%.2f%%)", format(sum(d$phone.luck == "low", na.rm = TRUE), big.mark = ","), 100 * mean(d$phone.luck == "low", na.rm = TRUE)),
		`High-luck, n (%)` = sprintf("%s (%.2f%%)", format(sum(d$phone.luck == "high", na.rm = TRUE), big.mark = ","), 100 * mean(d$phone.luck == "high", na.rm = TRUE))
	)
})
writexl::write_xlsx(the_table1, "TabS1.xlsx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig1. Overall EMS call spectrum and hourly distribution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!exists("fig_hourly_all")) {
	fig_hourly_all <- bind_rows(lapply(years, function(y) {
		d0 <- dat1.list[[as.character(y)]]; if (is.null(d0) || !nrow(d0)) return(tibble())
		d0 %>% filter(dx_grp %in% dxs.all, phone.luck %in% c("low", "middle", "high"), is.finite(hour), between(as.integer(hour), 0, 23)) %>% transmute(year = as.integer(y), dx_grp = factor(as.character(dx_grp), levels = dxs.all), phone.luck = factor(as.character(phone.luck), levels = c("low", "middle", "high")), hour = as.integer(hour)) %>% count(year, dx_grp, phone.luck, hour, name = "call_n") %>% complete(year = as.integer(y), dx_grp = dxs.all, phone.luck = c("low", "middle", "high"), hour = 0:23, fill = list(call_n = 0)) %>% group_by(year, dx_grp, phone.luck) %>% mutate(total_n = sum(call_n, na.rm = TRUE), pct = ifelse(total_n > 0, call_n / total_n, NA_real_), pct_smooth = (call_n + 0.5) / (total_n + 24 * 0.5)) %>% ungroup()
	}))
}

datS1.week <- bind_rows(lapply(yrs_trend, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(!is.na(dx_grp), !is.na(.data[["日期"]]), dx_grp %in% dxs.all) %>%
		transmute(year = y, date = as.Date(.data[["日期"]]), dx = factor(as.character(dx_grp), levels = dxs.all)) %>%
		mutate(week_start = floor_date(date, "week", week_start = 1)) %>%
		group_by(year, week_start, dx) %>%
		summarise(call_count = n(), days = n_distinct(date), .groups = "drop") %>%
		filter(days == 7)
}))

datS1.week <- datS1.week %>%
	group_by(year, week_start) %>%
	mutate(total_calls_week = sum(call_count), pct_week = call_count / total_calls_week) %>%
	ungroup()
fig1B_y <- if (call_volume == "raw") "call_count" else "pct_week"
fig1B_y_scale <- if (call_volume == "raw") {
	scale_y_continuous(labels = scales::comma, breaks = scales::breaks_pretty(n = 3))
} else {
	scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = scales::breaks_pretty(n = 3))
}
fig1B_title <- if (call_volume == "raw") "Raw weekly call volume" else "Call-volume-adjusted weekly proportion"

fig1B_dat <- datS1.week %>% filter(dx %in% dxs.all)

pFig1 <- ggplot(fig1B_dat, aes(week_start, .data[[fig1B_y]], color = dx, group = dx)) +
	geom_line(linewidth = .75) +
	facet_wrap(~dx, scales = "free_y", ncol = 3, labeller = as_labeller(setNames(dx_to_eng(dxs.all), dxs.all))) +
	scale_color_manual(values = dxs.all.color[dxs.all], guide = "none") +
	fig1B_y_scale +
	scale_x_date(breaks = as.Date(paste0(yrs_trend, "-01-01")), labels = yrs_trend) +
	labs(title = fig1B_title, x = NULL, y = NULL) +
	fig_theme(base_size = 8.5) +
	theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 7), axis.text.y = element_text(size = 7), strip.text = element_text(size = 8.2, face = "bold"), plot.title = element_text(size = 10.5, face = "bold"), axis.title = element_text(size = 8.5, face = "bold"))

Fig1 <- pFig1
save_plot(Fig1, "Fig1.png", width = 13.26, height = 12.2, dpi = 600)
unlink(c("Fig1c.png", "Fig1c.out.xlsx", "Fig1ac.png", "Fig1ac.out.xlsx"), force = TRUE)

# Deep-dive tables for the post-2023 increase in Death-coded EMS calls.
# These tables do not alter Fig1 layout; they only add reviewer-ready diagnostics to Fig1.out.xlsx.
fig1_dx_base <- bind_rows(lapply(yrs_trend, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(!is.na(dx_grp), !is.na(.data[["日期"]]), dx_grp %in% dxs.all) %>%
		transmute(
			year = as.integer(y),
			date = as.Date(.data[["日期"]]),
			dx_grp = factor(as.character(dx_grp), levels = dxs.all),
			sex = as.character(.data[[vars.dxs[1]]]),
			age = suppressWarnings(as.numeric(.data[[vars.dxs[2]]])),
			phone.luck = factor(as.character(phone.luck), levels = c("low", "middle", "high")),
			geo.type1 = as.character(geo.type1)
		)
}))
fig1_death_daily <- fig1_dx_base %>%
	mutate(is_death = dx_grp == "Death") %>%
	group_by(year, date) %>%
	summarise(death_calls = sum(is_death, na.rm = TRUE), major_calls = n(), death_pct = death_calls / major_calls, .groups = "drop") %>%
	arrange(date) %>%
	group_by(year) %>%
	mutate(death_calls_7d = roll7(death_calls), death_pct_7d = roll7(death_pct)) %>%
	ungroup()
fig1_death_monthly <- fig1_dx_base %>%
	mutate(month_start = floor_date(date, "month"), is_death = dx_grp == "Death") %>%
	group_by(year, month_start) %>%
	summarise(death_calls = sum(is_death, na.rm = TRUE), major_calls = n(), death_pct = death_calls / major_calls, mean_age_death = mean(age[is_death], na.rm = TRUE), female_pct_death = mean(sex[is_death] == "女", na.rm = TRUE), high_luck_pct_death = mean(phone.luck[is_death] == "high", na.rm = TRUE), .groups = "drop") %>%
	arrange(month_start)
fig1_death_yearly <- fig1_dx_base %>%
	mutate(is_death = dx_grp == "Death") %>%
	group_by(year) %>%
	summarise(death_calls = sum(is_death, na.rm = TRUE), major_calls = n(), death_pct = death_calls / major_calls, mean_age_death = mean(age[is_death], na.rm = TRUE), female_pct_death = mean(sex[is_death] == "女", na.rm = TRUE), high_luck_pct_death = mean(phone.luck[is_death] == "high", na.rm = TRUE), .groups = "drop") %>%
	arrange(year) %>%
	mutate(death_pct_yoy_ratio = death_pct / lag(death_pct), death_calls_yoy_ratio = death_calls / lag(death_calls), excess_death_calls_vs_previous_year = death_calls - lag(death_calls))
fig1_compare_death_years <- function(tbl, target_year, ref_year) {
	a <- tbl$death_calls[tbl$year == target_year][1]; A <- tbl$major_calls[tbl$year == target_year][1]
	b <- tbl$death_calls[tbl$year == ref_year][1]; B <- tbl$major_calls[tbl$year == ref_year][1]
	if (!all(is.finite(c(a, A, b, B))) || min(a, b, A - a, B - b, na.rm = TRUE) < 1) return(tibble(target_year = target_year, ref_year = ref_year, RR = NA_real_, RR_lo = NA_real_, RR_hi = NA_real_, p = NA_real_))
	RR <- (a / A) / (b / B); se <- sqrt(1 / a - 1 / A + 1 / b - 1 / B)
	tibble(target_year = target_year, ref_year = ref_year, death_calls_target = a, major_calls_target = A, death_pct_target = a / A, death_calls_ref = b, major_calls_ref = B, death_pct_ref = b / B, RR = RR, RR_lo = exp(log(RR) - 1.96 * se), RR_hi = exp(log(RR) + 1.96 * se), p = suppressWarnings(chisq.test(matrix(c(a, A - a, b, B - b), nrow = 2))$p.value))
}
fig1_death_year_contrasts <- bind_rows(
	fig1_compare_death_years(fig1_death_yearly, 2020, 2019),
	fig1_compare_death_years(fig1_death_yearly, 2021, 2019),
	fig1_compare_death_years(fig1_death_yearly, 2022, 2019),
	fig1_compare_death_years(fig1_death_yearly, 2023, 2022),
	fig1_compare_death_years(fig1_death_yearly, 2024, 2022),
	fig1_compare_death_years(fig1_death_yearly, 2024, 2023)
) %>% mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj), label = ifelse(is.finite(RR), sprintf("%.2f (%.2f–%.2f)%s", RR, RR_lo, RR_hi, sig05), NA_character_))
fig1_death_peak_weeks <- fig1_death_daily %>%
	mutate(week_start = floor_date(date, "week", week_start = 1)) %>%
	group_by(week_start) %>%
	summarise(death_calls = sum(death_calls), major_calls = sum(major_calls), death_pct = death_calls / major_calls, .groups = "drop") %>%
	arrange(desc(death_pct), desc(death_calls)) %>%
	slice_head(n = 30)
fig1_death_context_by_period <- fig1_dx_base %>%
	mutate(period = case_when(year <= 2019 ~ "2017-2019 baseline", year == 2020 ~ "2020 first COVID wave", year %in% 2021:2022 ~ "2021-2022 pre-reopening", year == 2023 ~ "2023 post-reopening", year == 2024 ~ "2024 sustained post-reopening", TRUE ~ as.character(year)), is_death = dx_grp == "Death") %>%
	group_by(period, is_death) %>%
	summarise(n = n(), mean_age = mean(age, na.rm = TRUE), female_pct = mean(sex == "女", na.rm = TRUE), high_luck_pct = mean(phone.luck == "high", na.rm = TRUE), low_luck_pct = mean(phone.luck == "low", na.rm = TRUE), residential_pct = mean(geo.type1 == "住宅区", na.rm = TRUE), .groups = "drop")

writexl::write_xlsx(list(
	weekly_data_plotted = fig1B_dat,
	configuration = tibble(call_volume = call_volume, panel_y = fig1B_y, n_phenotypes = dplyr::n_distinct(fig1B_dat$dx), layout = "5 rows x 3 columns", death_deep_dive_note = "Additional Death sheets quantify the post-2023 increase and compare it with 2019/2022 references."),
	weekly_summary = datS1.week %>% group_by(dx) %>% summarise(weeks = n(), mean_weekly_calls = round(mean(call_count), 1), sd_weekly_calls = round(sd(call_count), 1), mean_weekly_pct = round(mean(pct_week), 4), sd_weekly_pct = round(sd(pct_week), 4), .groups = "drop"),
	yearly_summary = datS1.week %>% group_by(year, dx) %>% summarise(year_total_calls = sum(call_count), mean_weekly_calls = round(mean(call_count), 1), mean_weekly_pct = round(mean(pct_week), 4), .groups = "drop"),
	death_daily_2017_2024 = fig1_death_daily,
	death_monthly_2017_2024 = fig1_death_monthly,
	death_yearly_deep_dive = fig1_death_yearly,
	death_year_contrasts = fig1_death_year_contrasts,
	death_peak_weeks = fig1_death_peak_weeks,
	death_context_by_period = fig1_death_context_by_period
), "Fig1.out.xlsx")
ems120_maybe_exit_after("fig1")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig2. NLP / transformer phenotyping pipeline and validation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main computational figure.  This is deliberately built from the existing
# stratified 5-fold OOF outputs; no second model is trained here.
if (ems120_should_run("fig2")) {
fig2nlp_dir <- file.path(dir.analysis, "raw", "hfl_cv5")
fig2nlp_files <- c(
	summary = "cv_summary.csv",
	per_class = "cv_per_class_grouped.csv",
	confusion = "cv_confusion_grouped.csv",
	oof = "cv_oof_predictions.csv",
	bootstrap = "cv_bootstrap_metrics.csv",
	calibration = "cv_confidence_calibration.csv",
	folds = "cv_fold_metrics.csv"
)
fig2nlp_paths <- setNames(file.path(fig2nlp_dir, unname(fig2nlp_files)), names(fig2nlp_files))
if (!all(file.exists(fig2nlp_paths[c("summary", "per_class", "confusion", "oof")]))) {
	stop("Fig2 requires the MacBERT stratified 5-fold OOF outputs under raw/hfl_cv5. Run ./ems120.sh without --skip-cv.", call. = FALSE)
}

fig2nlp_summary <- read.csv(fig2nlp_paths[["summary"]], check.names = FALSE, stringsAsFactors = FALSE)
fig2nlp_perclass <- read.csv(fig2nlp_paths[["per_class"]], check.names = FALSE, stringsAsFactors = FALSE)
fig2nlp_cm0 <- read.csv(fig2nlp_paths[["confusion"]], check.names = FALSE, stringsAsFactors = FALSE)
fig2nlp_oof <- read.csv(fig2nlp_paths[["oof"]], check.names = FALSE, stringsAsFactors = FALSE)
fig2nlp_boot <- if (file.exists(fig2nlp_paths[["bootstrap"]])) read.csv(fig2nlp_paths[["bootstrap"]], check.names = FALSE, stringsAsFactors = FALSE) else tibble()
fig2nlp_cal <- if (file.exists(fig2nlp_paths[["calibration"]])) read.csv(fig2nlp_paths[["calibration"]], check.names = FALSE, stringsAsFactors = FALSE) else tibble()
fig2nlp_folds <- if (file.exists(fig2nlp_paths[["folds"]])) read.csv(fig2nlp_paths[["folds"]], check.names = FALSE, stringsAsFactors = FALSE) else tibble()

fig2nlp_n <- fig2nlp_summary %>% filter(level == "grouped", system == "MacBERT") %>% pull(n) %>% dplyr::first(default = nrow(fig2nlp_oof))
fig2nlp_n <- suppressWarnings(as.integer(fig2nlp_n))
if (!is.finite(fig2nlp_n)) fig2nlp_n <- nrow(fig2nlp_oof)

# a. Workflow schematic: source fields -> expert labels -> OOF branches -> deployment.
fig2nlp_nodes <- tibble(
	id = c("raw", "gold", "kw", "bert", "oof", "deploy"),
	x = c(0.10, 0.31, 0.53, 0.53, 0.75, 0.93),
	y = c(0.50, 0.50, 0.70, 0.30, 0.50, 0.50),
	w = c(0.16, 0.16, 0.17, 0.17, 0.17, 0.13),
	h = c(0.29, 0.29, 0.25, 0.25, 0.30, 0.30),
	label = c(
		"Raw EMS records\n11 clinical/text fields",
		sprintf("Expert-labelled 2019 set\nN = %s", scales::comma(fig2nlp_n)),
		"Keyword baseline\ntransparent rules + log-odds",
		"MacBERT\nChinese transformer fine-tuning",
		"Stratified 5-fold OOF\neach record tested once",
		"Confidence-aware\n16-group reconstruction"
	)
) %>% mutate(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2)
fig2nlp_edges <- tibble(
	x = c(.18, .39, .39, .615, .615, .835),
	y = c(.50, .50, .50, .70, .30, .50),
	xend = c(.23, .445, .445, .665, .665, .865),
	yend = c(.50, .70, .30, .53, .47, .50)
)
p2NLP_A <- ggplot() +
	geom_segment(data = fig2nlp_edges, aes(x = x, y = y, xend = xend, yend = yend),
		arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed"), linewidth = .55, color = "grey45") +
	geom_rect(data = fig2nlp_nodes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
		fill = "white", color = "grey25", linewidth = .55) +
	geom_text(data = fig2nlp_nodes, aes(x = x, y = y, label = label), size = 3.0, fontface = "bold", lineheight = .95) +
	annotate("text", x = .53, y = .94, label = "Two OOF comparators", size = 2.7, fontface = "bold", color = "grey35") +
	coord_cartesian(xlim = c(.01, 1), ylim = c(.08, .98), clip = "off") +
	labs(title = "a. EMS text-to-phenotype transformer workflow") +
	theme_void(base_size = 9) +
	theme(plot.title = element_text(face = "bold", size = 10.5, hjust = 0), plot.margin = margin(6, 5, 5, 5))

# b. Expert-labelled OOF cohort composition.
fig2nlp_comp <- fig2nlp_perclass %>%
	filter(level == "grouped", system == "MacBERT") %>%
	transmute(class = as.character(class), support = as.integer(support)) %>%
	filter(is.finite(support), support > 0) %>%
	arrange(support) %>%
	mutate(class_f = factor(class, levels = class), pct = support / sum(support))
p2NLP_B <- ggplot(fig2nlp_comp, aes(support, class_f)) +
	geom_col(width = .68, fill = "#4C78A8") +
	geom_text(aes(label = scales::comma(support)), hjust = -.12, size = 2.65, fontface = "bold") +
	scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0, .18))) +
	labs(title = "b. Expert-labelled OOF cohort composition", x = "Expert-labelled records", y = NULL) +
	fig_theme(base_size = 8.4) +
	theme(plot.title = element_text(size = 10.5, hjust = 0), axis.text.y = element_text(size = 7.4), panel.grid.major.y = element_blank())

# c. Overall OOF performance with bootstrap 95% CIs.
fig2nlp_metric_labels <- c(accuracy = "Accuracy", macro_f1 = "Macro-F1", weighted_f1 = "Weighted-F1")
if (nrow(fig2nlp_boot)) {
	fig2nlp_overall <- fig2nlp_boot %>%
		filter(level == "grouped", system %in% c("Keyword", "MacBERT"), metric %in% names(fig2nlp_metric_labels)) %>%
		mutate(metric_label = factor(fig2nlp_metric_labels[metric], levels = rev(unname(fig2nlp_metric_labels))),
			system = factor(system, levels = c("Keyword", "MacBERT")))
} else {
	fig2nlp_overall <- fig2nlp_summary %>%
		filter(level == "grouped", system %in% c("Keyword", "MacBERT")) %>%
		select(system, accuracy, macro_f1, weighted_f1) %>%
		pivot_longer(c(accuracy, macro_f1, weighted_f1), names_to = "metric", values_to = "estimate") %>%
		mutate(lo = estimate, hi = estimate, metric_label = factor(fig2nlp_metric_labels[metric], levels = rev(unname(fig2nlp_metric_labels))),
			system = factor(system, levels = c("Keyword", "MacBERT")))
}
fig2nlp_gain_overall <- fig2nlp_overall %>%
	select(metric, system, estimate) %>%
	pivot_wider(names_from = system, values_from = estimate) %>%
	mutate(delta = MacBERT - Keyword)
p2NLP_C <- ggplot(fig2nlp_overall, aes(estimate, metric_label, color = system, shape = system)) +
	geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .12, position = position_dodge(width = .34), linewidth = .6) +
	geom_point(position = position_dodge(width = .34), size = 2.8) +
	geom_text(aes(label = sprintf("%.3f", estimate)), position = position_dodge(width = .34), hjust = -.22, size = 2.55, show.legend = FALSE) +
	scale_color_manual(values = c(Keyword = "grey50", MacBERT = "#2166AC")) +
	scale_shape_manual(values = c(Keyword = 16, MacBERT = 17)) +
	scale_x_continuous(limits = c(.62, .94), breaks = seq(.65, .90, .05), expand = expansion(mult = c(0, .02))) +
	labs(title = "c. MacBERT OOF performance", subtitle = "Points show complete OOF estimates; bars are bootstrap 95% CIs", x = "Performance", y = NULL, color = NULL, shape = NULL) +
	fig_theme(base_size = 8.5) +
	theme(plot.title = element_text(size = 10.5, hjust = 0), plot.subtitle = element_text(size = 7.3), legend.position = "top", panel.grid.major.y = element_blank())

# d. Per-phenotype F1 gain: dumbbell rather than a MacBERT-only lollipop.
fig2nlp_gain <- fig2nlp_perclass %>%
	filter(level == "grouped", system %in% c("Keyword", "MacBERT")) %>%
	select(system, class, f1, support) %>%
	pivot_wider(names_from = system, values_from = f1) %>%
	filter(is.finite(Keyword), is.finite(MacBERT)) %>%
	mutate(delta = MacBERT - Keyword) %>%
	arrange(delta) %>%
	mutate(class_f = factor(class, levels = class))
p2NLP_D <- ggplot(fig2nlp_gain, aes(y = class_f)) +
	geom_segment(aes(x = Keyword, xend = MacBERT, yend = class_f), color = "grey70", linewidth = 1.0) +
	geom_point(aes(x = Keyword), color = "grey45", size = 2.0) +
	geom_point(aes(x = MacBERT), color = "#2166AC", size = 2.35) +
	geom_text(aes(x = pmax(Keyword, MacBERT) + .025, label = sprintf("%+.2f", delta)), hjust = 0, size = 2.25, fontface = "bold") +
	scale_x_continuous(limits = c(0, 1.16), breaks = seq(0, 1, .2)) +
	labs(title = "d. Per-phenotype F1 gain over keyword baseline", x = "F1  (grey: keyword; blue: MacBERT)", y = NULL) +
	fig_theme(base_size = 8.0) +
	theme(plot.title = element_text(size = 10.5, hjust = 0), axis.text.y = element_text(size = 7.2), panel.grid.major.y = element_blank())

# e. Summarize the confusion structure as the largest off-diagonal OOF errors.
names(fig2nlp_cm0)[1] <- "true_group"
fig2nlp_cm_long <- fig2nlp_cm0 %>%
	pivot_longer(-true_group, names_to = "pred_group", values_to = "n") %>%
	group_by(true_group) %>%
	mutate(row_total = sum(n), row_pct = ifelse(row_total > 0, n / row_total, NA_real_)) %>%
	ungroup()
fig2nlp_errors <- fig2nlp_cm_long %>%
	filter(true_group != pred_group, n > 0, is.finite(row_pct)) %>%
	arrange(desc(n), desc(row_pct)) %>%
	slice_head(n = 10) %>%
	mutate(error = paste0(true_group, "  →  ", pred_group), error_f = factor(error, levels = rev(error)))
p2NLP_E <- ggplot(fig2nlp_errors, aes(row_pct, error_f)) +
	geom_col(width = .65, fill = "#7B9EBD") +
	geom_text(aes(label = sprintf("%.1f%%  (n=%s)", 100 * row_pct, scales::comma(n))), hjust = -.08, size = 2.45, fontface = "bold") +
	scale_x_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, .34))) +
	labs(title = "e. Largest OOF misclassification pathways", x = "% of expert-labelled row", y = NULL) +
	fig_theme(base_size = 8.0) +
	theme(plot.title = element_text(size = 10.5, hjust = 0), axis.text.y = element_text(size = 7.3), panel.grid.major.y = element_blank())

# f. Confidence threshold utility: coverage and retained-set performance.
fig2nlp_macro_f1 <- function(true, pred, lev) {
	z <- vapply(lev, function(cl) {
		tp <- sum(true == cl & pred == cl, na.rm = TRUE)
		fp <- sum(true != cl & pred == cl, na.rm = TRUE)
		fn <- sum(true == cl & pred != cl, na.rm = TRUE)
		pr <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
		rc <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
		if (!is.finite(pr) || !is.finite(rc) || (pr + rc) == 0) return(0)
		2 * pr * rc / (pr + rc)
	}, numeric(1))
	mean(z, na.rm = TRUE)
}
fig2nlp_oof2 <- fig2nlp_oof %>%
	transmute(
		true_group = as.character(true_group),
		pred_group = as.character(macbert_pred_group),
		confidence = suppressWarnings(as.numeric(macbert_confidence))
	) %>% filter(is.finite(confidence), !is.na(true_group), !is.na(pred_group))
fig2nlp_levels <- unique(c(dxs.all0, fig2nlp_oof2$true_group))
fig2nlp_levels <- fig2nlp_levels[fig2nlp_levels %in% unique(fig2nlp_oof2$true_group)]
fig2nlp_thresholds <- sort(unique(c(seq(.50, .98, by = .02), .80)))
fig2nlp_utility <- purrr::map_dfr(fig2nlp_thresholds, function(th) {
	d <- fig2nlp_oof2 %>% filter(confidence >= th)
	if (!nrow(d)) return(tibble(threshold = th, coverage = 0, accuracy = NA_real_, macro_f1 = NA_real_, n = 0L))
	tibble(
		threshold = th,
		coverage = nrow(d) / nrow(fig2nlp_oof2),
		accuracy = mean(d$true_group == d$pred_group),
		macro_f1 = fig2nlp_macro_f1(d$true_group, d$pred_group, fig2nlp_levels),
		n = nrow(d)
	)
})
fig2nlp_u80 <- fig2nlp_utility %>% filter(abs(threshold - .80) < 1e-8) %>% slice(1)
fig2nlp_u80_sub <- if (nrow(fig2nlp_u80)) sprintf(">=0.80 retains %.1f%%; accuracy %.1f%%; macro-F1 %.3f", 100*fig2nlp_u80$coverage, 100*fig2nlp_u80$accuracy, fig2nlp_u80$macro_f1) else NULL
fig2nlp_utility_long <- fig2nlp_utility %>%
	select(threshold, Coverage = coverage, Accuracy = accuracy, `Macro-F1` = macro_f1) %>%
	pivot_longer(-threshold, names_to = "metric", values_to = "value") %>%
	mutate(metric = factor(metric, levels = c("Accuracy", "Macro-F1", "Coverage")))
p2NLP_F <- ggplot(fig2nlp_utility_long, aes(threshold, value, color = metric, linetype = metric)) +
	geom_vline(xintercept = .80, linetype = "dotted", color = "grey45", linewidth = .7) +
	geom_line(linewidth = .9) +
	geom_point(data = fig2nlp_utility_long %>% filter(abs(threshold - .80) < 1e-8), size = 2.2) +
	scale_color_manual(values = c(Accuracy = "#2166AC", `Macro-F1` = "#008B8B", Coverage = "grey45")) +
	scale_linetype_manual(values = c(Accuracy = "solid", `Macro-F1` = "solid", Coverage = "dashed")) +
	scale_x_continuous(breaks = c(.5, .6, .7, .8, .9, .98)) +
	scale_y_continuous(limits = c(0, 1.02), labels = scales::percent_format(accuracy = 1)) +
	labs(title = "f. Confidence-aware deployment utility", subtitle = fig2nlp_u80_sub, x = "Minimum MacBERT confidence", y = "Retained fraction / performance", color = NULL, linetype = NULL) +
	fig_theme(base_size = 8.3) +
	theme(plot.title = element_text(size = 10.5, hjust = 0), plot.subtitle = element_text(size = 7.4), legend.position = "top", legend.text = element_text(size = 7.2))



# Publication-oriented main Fig2.  The analytical objects above are unchanged;
# only the visual presentation is rebuilt here.  The figure now focuses on the
# key publishable messages rather than on a workflow cartoon: (1) MacBERT beats
# the keyword baseline, (2) performance gains are broad but heterogeneous across
# phenotypes, (3) residual errors cluster in a small number of pathways, and
# (4) confidence filtering is practically useful.

fig2_pub_theme <- function(base_size = 9) {
	theme_classic(base_size = base_size) +
		theme(
			panel.grid = element_blank(),
			axis.title = element_text(face = "bold", color = "grey15"),
			axis.text = element_text(color = "grey20"),
			plot.title = element_text(face = "bold", size = base_size + 1.4, hjust = 0),
			plot.subtitle = element_text(size = base_size - .6, color = "grey35", hjust = 0),
			legend.title = element_text(face = "bold"),
			plot.margin = margin(5, 7, 5, 5)
		)
}

# Match the visual language of Fig5 panel b across all six Fig2 panels.
fig2_panel_theme <- theme(
	plot.title = element_text(size = 10, face = "bold", hjust = .5),
	plot.subtitle = element_blank(),
	axis.text = element_text(size = 8.4, face = "bold", color = "grey15"),
	axis.title = element_text(size = 8.5, face = "bold", color = "grey15"),
	axis.title.x = element_text(margin = margin(t = 8)),
	axis.title.y = element_text(margin = margin(r = 8)),
	# Keep legends inside the panel viewport.  External top legends change the
	# title-grob height and make titles in the same row sit at different levels.
	legend.position = "inside",
	legend.position.inside = c(.5, .985),
	legend.justification = c(.5, 1),
	legend.direction = "horizontal",
	legend.text = element_text(size = 7.2, face = "bold"),
	legend.title = element_text(size = 7.2, face = "bold"),
	legend.key.width = grid::unit(16, "pt"),
	legend.spacing.x = grid::unit(8, "pt"),
	legend.margin = margin(0, 0, 0, 0),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	plot.margin = margin(6, 7, 6, 7)
)

# a. Overall grouped performance gain as a compact dumbbell panel.
fig2nlp_overall_wide <- fig2nlp_overall %>%
	select(system, metric_label, estimate, lo, hi) %>%
	mutate(metric_label = factor(metric_label, levels = rev(unname(fig2nlp_metric_labels)))) %>%
	pivot_wider(names_from = system, values_from = c(estimate, lo, hi), names_sep = "__") %>%
	mutate(gain = estimate__MacBERT - estimate__Keyword)

p2NLP_mainA <- ggplot(fig2nlp_overall_wide) +
	geom_segment(aes(y = estimate__Keyword, yend = estimate__MacBERT, x = metric_label, xend = metric_label), linewidth = 1.0, color = "grey76") +
	geom_errorbar(aes(ymin = lo__Keyword, ymax = hi__Keyword, x = metric_label, color = "Keyword"), width = .08, linewidth = .55) +
	geom_errorbar(aes(ymin = lo__MacBERT, ymax = hi__MacBERT, x = metric_label, color = "MacBERT"), width = .08, linewidth = .55) +
	geom_point(aes(y = estimate__Keyword, x = metric_label, color = "Keyword"), size = 2.4) +
	geom_point(aes(y = estimate__MacBERT, x = metric_label, color = "MacBERT"), size = 2.8) +
	geom_text(aes(y = estimate__Keyword - .014, x = metric_label, label = sprintf("%.3f", estimate__Keyword)), vjust = 1, size = 2.65, color = "grey35", fontface = "bold") +
	geom_text(aes(y = estimate__MacBERT + .014, x = metric_label, label = sprintf("%.3f", estimate__MacBERT)), vjust = 0, size = 2.65, color = "#D55E00", fontface = "bold") +
	scale_color_manual(values = c(Keyword = "grey48", MacBERT = "#D55E00"), breaks = c("Keyword", "MacBERT"), name = NULL) +
	scale_y_continuous(limits = c(.50, 1.00), breaks = seq(.50, 1.00, .10), expand = expansion(mult = c(.01, .01))) +
	labs(title = "a. Grouped held-out performance", x = NULL, y = "Performance") +
	fig2_pub_theme(base_size = 8.6) + fig2_panel_theme

# b. Per-phenotype MacBERT precision-recall map with bubble size = support and
# colour = F1 gain over keyword.
fig2nlp_pr <- fig2nlp_perclass %>%
	filter(level == "grouped", system %in% c("Keyword", "MacBERT")) %>%
	select(system, class, precision, recall, f1, support) %>%
	pivot_wider(names_from = system, values_from = c(precision, recall, f1, support), names_sep = "__") %>%
	mutate(
		delta_f1 = f1__MacBERT - f1__Keyword,
		class = as.character(class),
		lab = ifelse(class == "Other", "Other", gsub("/", "/\n", class)),
		label_x = pmin(.965, pmax(.535, recall__MacBERT + rep(c(-.055, .055), length.out = n()))),
		label_y = pmin(.965, pmax(.535, precision__MacBERT + rep(c(.040, -.040, .055, -.055), length.out = n())))
	)

p2NLP_mainB <- ggplot(fig2nlp_pr, aes(recall__MacBERT, precision__MacBERT)) +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey72", linewidth = .55) +
	geom_vline(xintercept = .80, linetype = "dotted", color = "grey82", linewidth = .45) +
	geom_hline(yintercept = .80, linetype = "dotted", color = "grey82", linewidth = .45) +
	geom_segment(aes(xend = label_x, yend = label_y), linewidth = .32, color = "grey48") +
	geom_point(aes(size = support__MacBERT, fill = delta_f1), shape = 21, color = "grey25", stroke = .35, alpha = .92) +
	geom_text(aes(x = label_x, y = label_y, label = lab), size = 2.65, fontface = "bold", lineheight = .88, check_overlap = TRUE) +
	scale_fill_gradient2(low = "#F1E6FF", mid = "#B89AE6", high = "#6A3D9A", midpoint = 0.15, name = "Δ F1 vs keyword") +
	scale_size_continuous(range = c(2.3, 7.0), breaks = c(500, 1000, 2000), labels = scales::comma, name = "Expert-labelled records") +
	# Fill the same column width as the panels above and below.  coord_fixed()
	# forced this panel into a narrow square inside a wide grid cell, which made
	# the second row appear horizontally misaligned.
	coord_cartesian(xlim = c(.50, 1.00), ylim = c(.50, 1.00), expand = FALSE) +
	labs(title = "c. Per-phenotype precision-recall profile", x = "MacBERT recall", y = "MacBERT precision") +
	fig2_pub_theme(base_size = 7.9) +
	fig2_panel_theme +
	guides(size = guide_legend(nrow = 1, order = 1), fill = "none") +
	theme(legend.spacing.x = grid::unit(13, "pt"), legend.key.width = grid::unit(19, "pt"))

# Keep the sample-size legend across the top, but give the colour key its own
# compact vertical inset in the lower-right corner so the two guides cannot crowd.
fig2_c_colorbar_vals <- range(fig2nlp_pr$delta_f1, na.rm = TRUE)
fig2_c_colorbar_dat <- tibble(
	y = 3:1,
	level = factor(c("High", "Mid", "Low"), levels = c("Low", "Mid", "High")),
	label = sprintf("%.2f", c(fig2_c_colorbar_vals[2], mean(fig2_c_colorbar_vals), fig2_c_colorbar_vals[1]))
)
fig2_c_colorbar <- ggplot(fig2_c_colorbar_dat, aes(1, y)) +
	geom_tile(aes(fill = level), width = .48, height = .72) +
	geom_text(aes(x = 1.38, label = label), hjust = 0, size = 2.3, fontface = "bold", color = "grey20") +
	scale_fill_manual(values = c(Low = "#F1E6FF", Mid = "#B89AE6", High = "#6A3D9A"), guide = "none") +
	scale_x_continuous(breaks = NULL, labels = NULL) +
	scale_y_continuous(breaks = NULL, labels = NULL) +
	coord_cartesian(xlim = c(.68, 1.95), ylim = c(.5, 3.5), expand = FALSE, clip = "off") +
	labs(title = "Δ F1", x = NULL, y = NULL) +
	theme_void(base_size = 7) +
	theme(plot.title = element_text(size = 6.6, face = "bold", hjust = .5), plot.background = element_rect(fill = "white", color = "grey70", linewidth = .25), plot.margin = margin(2, 3, 2, 2))
p2NLP_mainB <- p2NLP_mainB + patchwork::inset_element(fig2_c_colorbar, left = .84, bottom = .15, right = .965, top = .45, align_to = "panel")

# c. Full row-normalised confusion matrix.
fig2nlp_cm_plot <- fig2nlp_cm_long %>%
	mutate(
		true_group = factor(true_group, levels = rev(dxs.all0)),
		pred_group = factor(pred_group, levels = dxs.all0),
		is_diag = as.character(true_group) == as.character(pred_group),
		cell_lab = ifelse(is.finite(row_pct) & (is_diag | row_pct >= .05), sprintf("%.0f", 100 * row_pct), "")
	)

p2NLP_mainC <- ggplot(fig2nlp_cm_plot, aes(pred_group, true_group, fill = row_pct)) +
	geom_tile(color = "white", linewidth = .20) +
	geom_tile(data = fig2nlp_cm_plot %>% filter(is_diag), fill = NA, color = "grey18", linewidth = .42) +
	geom_text(aes(label = cell_lab, fontface = ifelse(is_diag, "bold", "plain")), size = 1.85, show.legend = FALSE) +
	scale_fill_gradient(low = "grey98", high = "#2166AC", limits = c(0, 1), labels = scales::percent_format(accuracy = 1), name = "Row %", na.value = "grey96") +
	coord_fixed() +
	labs(title = "c. OOF confusion matrix", subtitle = "Errors are sparse and structured; diagonal concentration is strongest for traffic, intoxication, and Ob/Gyn", x = "Predicted phenotype", y = "Expert label") +
	fig2_pub_theme(base_size = 7.6) +
	theme(
		axis.line = element_blank(), axis.ticks = element_blank(),
		axis.text.x = element_text(angle = 52, hjust = 1, vjust = 1, size = 6.2),
		axis.text.y = element_text(size = 6.4),
		legend.position = "right", legend.key.height = grid::unit(.45, "cm")
	)

# d. Largest misclassification pathways.
fig2nlp_toperr <- fig2nlp_errors %>%
	arrange(desc(row_pct), desc(n)) %>%
	slice_head(n = 12) %>%
	mutate(error_f = factor(error, levels = rev(error)))

p2NLP_mainD <- ggplot(fig2nlp_toperr, aes(row_pct, error_f)) +
	geom_col(aes(fill = row_pct), width = .66) +
	geom_text(aes(label = sprintf("%.1f%%  (n=%s)", 100 * row_pct, scales::comma(n))), hjust = -0.10, size = 2.75, fontface = "bold", color = "grey18") +
	scale_fill_gradient(low = "#F6C1B5", high = "#B2182B", guide = "none") +
	scale_x_continuous(limits = c(0, max(.08, max(fig2nlp_toperr$row_pct, na.rm = TRUE) * 1.18)), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, .02))) +
	labs(title = "d. Largest off-diagonal pathways", x = "Row-normalised error share", y = NULL) +
	fig2_pub_theme(base_size = 7.7) +
	fig2_panel_theme + theme(axis.text.y = element_text(size = 7.2, face = "bold"))

# e. Confidence-aware deployment utility.
fig2nlp_utility_plot <- fig2nlp_utility %>%
	select(threshold, Coverage = coverage, Accuracy = accuracy, `Macro-F1` = macro_f1) %>%
	pivot_longer(-threshold, names_to = "metric", values_to = "value") %>%
	mutate(metric = factor(metric, levels = c("Coverage", "Accuracy", "Macro-F1")))
fig2nlp_u80_label <- if (nrow(fig2nlp_u80)) {
	sprintf("Threshold: 0.80\nRetained: %.1f%%\nAccuracy: %.1f%%\nMacro-F1: %.3f",
		100 * fig2nlp_u80$coverage, 100 * fig2nlp_u80$accuracy, fig2nlp_u80$macro_f1)
} else "0.80 threshold summary unavailable"

p2NLP_mainE <- ggplot(fig2nlp_utility_plot, aes(threshold, value, color = metric, linetype = metric)) +
	geom_vline(xintercept = .80, linetype = "dotted", color = "grey42", linewidth = .6) +
	geom_line(linewidth = .90) +
	geom_point(data = fig2nlp_utility_plot %>% filter(abs(threshold - .80) < 1e-8), size = 2.2) +
	annotate("label", x = .515, y = .575, label = fig2nlp_u80_label, hjust = 0, vjust = 0, size = 2.75, fontface = "bold", lineheight = 1.25, linewidth = .28, label.padding = grid::unit(.34, "lines"), fill = "white", color = "grey20") +
	scale_color_manual(values = c(Coverage = "#7A7A7A", Accuracy = "#2A9D8F", `Macro-F1` = "#8AB17D")) +
	scale_linetype_manual(values = c(Coverage = "dashed", Accuracy = "solid", `Macro-F1` = "solid")) +
	scale_x_continuous(breaks = c(.5, .6, .7, .8, .9, .98)) +
	scale_y_continuous(limits = c(.50, 1.00), labels = scales::percent_format(accuracy = 1)) +
	labs(title = "e. Confidence-threshold utility", x = "Minimum MacBERT confidence", y = "Retained fraction / performance", color = NULL, linetype = NULL) +
	fig2_pub_theme(base_size = 7.9) +
	fig2_panel_theme

# Cross-validation diagnostic panels are built here so the main figure can be
# assembled before the supplementary staging block is reached.
fig2_fold_cols <- intersect(names(fig2nlp_metric_labels), names(fig2nlp_folds))
fig2_fold_long <- if (nrow(fig2nlp_folds) && length(fig2_fold_cols) == length(fig2nlp_metric_labels)) {
	fig2nlp_folds %>%
		filter(level == "grouped", system %in% c("Keyword", "MacBERT")) %>%
		select(any_of(c("fold", "system", fig2_fold_cols))) %>%
		pivot_longer(cols = all_of(fig2_fold_cols), names_to = "metric", values_to = "value") %>%
		mutate(metric = factor(fig2nlp_metric_labels[metric], levels = unname(fig2nlp_metric_labels)), system = factor(system, levels = c("Keyword", "MacBERT")))
} else tibble()
p2NLP_fold_stability <- if (nrow(fig2_fold_long)) {
	ggplot(fig2_fold_long, aes(metric, value, color = system)) +
		geom_boxplot(aes(group = interaction(metric, system)), position = position_dodge(width = .62), width = .46, outlier.shape = NA, alpha = .12) +
		geom_point(position = position_jitterdodge(jitter.width = .08, dodge.width = .62), size = 1.8, alpha = .85) +
		scale_color_manual(values = c(Keyword = "grey50", MacBERT = "#7B2CBF")) +
		scale_y_continuous(limits = c(.50, 1.00), breaks = seq(.5, 1, .1), expand = expansion(mult = c(.02, .04))) +
		labs(title = "b. Fold-to-fold stability", x = NULL, y = "Performance", color = NULL) +
		fig2_pub_theme(base_size = 8.3) + fig2_panel_theme +
		theme(legend.spacing.x = grid::unit(16, "pt"), legend.key.width = grid::unit(24, "pt"))
} else empty_panel("b. Fold-to-fold stability", "Fold metrics unavailable")
p2NLP_calibration <- if (nrow(fig2nlp_oof)) {
	fig2nlp_cal_multi <- fig2nlp_oof %>%
		transmute(confidence = suppressWarnings(as.numeric(macbert_confidence)), correct = as.character(macbert_pred_group) == as.character(true_group)) %>%
		filter(is.finite(confidence), !is.na(correct)) %>%
		arrange(confidence)
	fig2nlp_cal_multi <- purrr::map_dfr(c(500L, 1000L, 2000L), function(bin_n) {
		fig2nlp_cal_multi %>%
			mutate(confidence_bin = ceiling(row_number() / bin_n)) %>%
			group_by(confidence_bin) %>%
			summarise(n = n(), mean_confidence = mean(confidence), observed_accuracy = mean(correct), .groups = "drop") %>%
			filter(n >= max(100L, floor(bin_n * .50))) %>%
			mutate(bin_size = factor(paste(bin_n, "records/bin"), levels = c("500 records/bin", "1000 records/bin", "2000 records/bin")))
	})
	ggplot(fig2nlp_cal_multi, aes(mean_confidence, observed_accuracy, color = bin_size, group = bin_size)) +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = .55) +
		geom_line(linewidth = .90) +
		geom_point(size = 1.7, alpha = .92) +
		scale_color_manual(values = c(`500 records/bin` = "#9ECAE1", `1000 records/bin` = "#3182BD", `2000 records/bin` = "#08519C"), name = NULL) +
		coord_cartesian(xlim = c(.50, 1.00), ylim = c(.50, 1.00), expand = FALSE) +
		labs(title = "f. Confidence calibration", x = "Mean confidence", y = "Observed accuracy") +
		fig2_pub_theme(base_size = 7.8) + fig2_panel_theme
} else empty_panel("f. Confidence calibration", "Calibration output unavailable")

p2NLP_mainB_fig2 <- p2NLP_mainB
p2NLP_mainE_fig2 <- p2NLP_mainE + labs(title = "d. Confidence-threshold utility")
p2NLP_mainD_fig2 <- p2NLP_mainD + labs(title = "e. Largest off-diagonal pathways")
Fig2_NLP <- align_panel_rows(
	rows = list(
		list(p2NLP_mainA, p2NLP_fold_stability),
		list(p2NLP_mainB_fig2, p2NLP_mainE_fig2),
		list(p2NLP_mainD_fig2, p2NLP_calibration)
	),
	rel_widths = c(1, 1),
	rel_heights = c(1, 1, 1),
	row_gap = .07,
	side_pad = .018,
	axis_text_size = 7.9,
	axis_title_size = 7.9
)

save_plot(Fig2_NLP, "Fig2.png", width = 13.2, height = 13.2, dpi = 600, bg = "white")

writexl::write_xlsx(list(
	panel_a_overall_performance = fig2nlp_overall_wide,
	panel_a_bootstrap_CI = fig2nlp_overall,
	panel_b_fold_stability = fig2_fold_long,
	panel_c_precision_recall = fig2nlp_pr,
	panel_d_threshold_utility = fig2nlp_utility,
	panel_e_error_pathways = fig2nlp_toperr,
	panel_f_calibration = fig2nlp_cal_multi,
	configuration = tibble(
		figure = "Fig2. MacBERT disease phenotyping validation",
		layout = "3 rows x 2 equal-width columns",
		panel_order = "a overall performance; b fold stability; c precision-recall; d threshold utility; e error pathways; f confidence calibration",
		CV = "stratified 5-fold out-of-fold",
		OOF_n = fig2nlp_n,
		confidence_threshold_highlight = .80,
		note = "Main Fig2 is a six-panel 3-by-2 layout with equal column widths and includes fold stability and confidence calibration."
	)
), "Fig2.out.xlsx")
ems120_maybe_exit_after("fig2")


}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig3. Housing-price and environmental validation for phone/geo/dx
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (ems120_should_run("fig3")) {
geo_env_vars <- fig2_required_geo_vars[-1]
# Fig2c/Fig2d coordinates are the original pickup-address longitude/latitude.
fig2_coord_vars <- fig2_required_coord_vars
fig2_num_vars <- fig2_required_geo_vars
fig2_house_var <- "房价指数"
fig2_geo_type_var <- geo.type.use
fig2_house_price_filter <- if (is.na(house_price.use)) NA_character_ else as.character(house_price.use)
fig2_house_price_filter_label <- if (is.na(fig2_house_price_filter)) "ALL" else fig2_house_price_filter

# Strict check: dat.list.rds should already contain Xia geo payload columns.
# Do not fall back to dat/[year].xlsx, because missing columns indicate a broken cache.
check_required_vars_in_dat_list(dat1.list, required_vars = c(fig2_geo_type_var, fig2_num_vars, fig2_coord_vars), stop_on_missing = TRUE, context = basename(file.dat.rds))
fig2_keep_vars <- unique(c("year", "ID", "phone.sco", "phone.luck", "dx_grp", "fig2_geo_type", "地址类型", "geo.type", "geo.type1", fig2_num_vars, fig2_coord_vars))
fig2_text_vars <- setdiff(fig2_keep_vars, c("year", "ID", "phone.sco", fig2_num_vars, fig2_coord_vars))

fig2_base <- bind_rows(lapply(years, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		mutate(
			across(any_of(fig2_text_vars), as.character),
			year = as.integer(y),
			fig2_geo_type = as.character(.data[[fig2_geo_type_var]])
		) %>%
		{ if (is.na(fig2_house_price_filter)) . else dplyr::filter(., .data$fig2_geo_type == fig2_house_price_filter) } %>%
		dplyr::select(any_of(fig2_keep_vars)) %>%
		mutate(across(all_of(c(fig2_num_vars, fig2_coord_vars)), ~ suppressWarnings(as.numeric(.x))))
})) %>%
	mutate(
		phone.sco = suppressWarnings(as.numeric(phone.sco)),
		phone.luck = factor(as.character(phone.luck), levels = c("low", "middle", "high")),
		dx_grp = factor(as.character(dx_grp), levels = dxs.all)
	)

fig2_n_before_house_price_min <- nrow(fig2_base)
if (fig2_house_var %in% names(fig2_base) && length(house_price_min) == 1L && is.finite(house_price_min)) {
	fig2_base <- fig2_base %>%
		filter(is.finite(.data[[fig2_house_var]]), .data[[fig2_house_var]] >= house_price_min)
}
fig2_n_after_house_price_min <- nrow(fig2_base)
if (!nrow(fig2_base)) {
	stop(sprintf("Fig2 has no rows after applying house_price_min = %s to %s. Lower house_price_min or check 房价指数.", house_price_min, fig2_house_var), call. = FALSE)
}

fig2_num_vars <- intersect(fig2_num_vars, names(fig2_base))
fig2_house_plot_var <- "house_price_capped"
fig2_house_price_cap <- c(low = 10000, high = 300000)
fig2_density_var <- "1km渔网建筑密度"
fig2_base <- fig2_base %>%
	mutate(across(any_of(c(fig2_num_vars, fig2_coord_vars)), ~ suppressWarnings(as.numeric(.x)))) %>%
	mutate(!!fig2_house_plot_var := if (fig2_house_var %in% names(.)) clamp(.data[[fig2_house_var]], fig2_house_price_cap["low"], fig2_house_price_cap["high"]) else NA_real_)
saveRDS(fig2_base, .ems120_register_temp_rds("Fig2.rds"))

top_luck_th <- c(0.20, 0.10, 0.05, 0.01)
top_luck_th <- unique(top_luck_th[is.finite(top_luck_th) & top_luck_th > 0 & top_luck_th < 1])
fig2_top_luck_cuts <- tibble(
	threshold = top_luck_th,
	threshold_label = paste0("Top ", scales::percent(top_luck_th, accuracy = 1)),
	cutoff = as.numeric(stats::quantile(fig2_base$phone.sco[is.finite(fig2_base$phone.sco)], probs = 1 - top_luck_th, na.rm = TRUE, type = 7))
)
fig2_top10_cut <- fig2_top_luck_cuts$cutoff[which.min(abs(fig2_top_luck_cuts$threshold - 0.10))][1]

fit_fig2_top_trend <- function(raw) {
	if (!nrow(raw) || !all(c("phone_top", "x") %in% names(raw))) return(tibble(beta = NA_real_, OR_per_x = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, p_adj = NA_real_, sig05 = "", n = nrow(raw)))
	d <- raw %>% filter(is.finite(x), !is.na(phone_top)) %>% mutate(phone_top = as.integer(phone_top))
	if (nrow(d) < 100 || length(unique(d$phone_top)) < 2) return(tibble(beta = NA_real_, OR_per_x = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, p_adj = NA_real_, sig05 = "", n = nrow(d)))
	form <- if ("year" %in% names(d) && n_distinct(d$year, na.rm = TRUE) > 1) phone_top ~ x + factor(year) else phone_top ~ x
	fit <- tryCatch(glm(form, family = binomial(), data = d), error = function(e) NULL)
	if (is.null(fit)) return(tibble(beta = NA_real_, OR_per_x = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, p_adj = NA_real_, sig05 = "", n = nrow(d)))
	td <- broom::tidy(fit) %>% filter(term == "x")
	if (!nrow(td)) return(tibble(beta = NA_real_, OR_per_x = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, p_adj = NA_real_, sig05 = "", n = nrow(d)))
	tibble(beta = td$estimate[1], OR_per_x = exp(td$estimate[1]), lo = exp(td$estimate[1] - 1.96 * td$std.error[1]), hi = exp(td$estimate[1] + 1.96 * td$std.error[1]), p = td$p.value[1], p_adj = td$p.value[1], sig05 = sig_star05(td$p.value[1]), n = nrow(d))
}

make_fig2_top10_panel <- function(base, value_col, value_label, panel_letter, transform = c("log", "log10", "log10p1", "identity"), trim_low_q = NULL, bins = 12) {
	transform <- match.arg(transform)
	if (!value_col %in% names(base)) {
		return(list(plot = empty_panel(sprintf("%s. Top phone luck vs %s", panel_letter, tolower(value_label)), msg = paste("Missing variable:", value_label)), raw = tibble(), bin = tibble(), hist = tibble(), trend = tibble(), trim_cut = NA_real_))
	}
	raw <- base %>%
		filter(is.finite(.data[[value_col]]), is.finite(phone.sco)) %>%
		mutate(value_raw = as.numeric(.data[[value_col]]))
	if (transform == "log") raw <- raw %>% filter(value_raw > 0) %>% mutate(x = log(value_raw))
	if (transform == "log10") raw <- raw %>% filter(value_raw > 0) %>% mutate(x = log10(value_raw))
	if (transform == "log10p1") raw <- raw %>% filter(value_raw >= 0) %>% mutate(x = log10(value_raw + 1))
	if (transform == "identity") raw <- raw %>% mutate(x = value_raw)
	raw <- raw %>% filter(is.finite(x))
	trim_cut <- NA_real_
	if (!is.null(trim_low_q) && is.finite(trim_low_q) && trim_low_q > 0 && nrow(raw)) {
		trim_cut <- as.numeric(stats::quantile(raw$value_raw, trim_low_q, na.rm = TRUE, type = 7))
		raw <- raw %>% filter(value_raw >= trim_cut)
	}
	if (nrow(raw) < 100 || !nrow(fig2_top_luck_cuts) || !is.finite(fig2_top10_cut)) {
		return(list(plot = empty_panel(sprintf("%s. Top phone luck vs %s", panel_letter, tolower(value_label))), raw = raw, bin = tibble(), hist = tibble(), trend = tibble(), trim_cut = trim_cut))
	}
	h4 <- hist(raw$x, breaks = bins, plot = FALSE)
	hist_df <- tibble(
		bin_id = seq_along(h4$counts),
		xmin = head(h4$breaks, -1),
		xmax = tail(h4$breaks, -1),
		x = (xmin + xmax) / 2,
		count = as.integer(h4$counts)
	)
	raw <- raw %>%
		mutate(bin_id = pmin(pmax(findInterval(x, h4$breaks, rightmost.closed = TRUE), 1L), length(h4$counts)))
	bin_df <- purrr::map_dfr(seq_len(nrow(fig2_top_luck_cuts)), function(i) {
		cut_i <- fig2_top_luck_cuts$cutoff[i]
		lab_i <- fig2_top_luck_cuts$threshold_label[i]
		raw %>%
			mutate(phone_top = phone.sco >= cut_i) %>%
			group_by(bin_id) %>%
			summarise(n = n(), top_pct = 100 * mean(phone_top, na.rm = TRUE), top_se = 100 * sqrt(pmax(mean(phone_top, na.rm = TRUE) * (1 - mean(phone_top, na.rm = TRUE)) / n(), 0)), .groups = "drop") %>%
			left_join(hist_df, by = "bin_id") %>%
			mutate(threshold_label = lab_i, threshold = fig2_top_luck_cuts$threshold[i], cutoff = cut_i)
	}) %>%
		arrange(threshold, bin_id)
	sec_axis_ylim <- c(0, max(5, ceiling(max(bin_df$top_pct + 1.96 * bin_df$top_se, na.rm = TRUE) / 5) * 5))
	scale_factor <- max(hist_df$count, na.rm = TRUE) / max(sec_axis_ylim[2], 1)
	if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1
	bin_df <- bin_df %>%
		mutate(
			threshold_label = factor(threshold_label, levels = fig2_top_luck_cuts$threshold_label),
			top_pct_plot = top_pct * scale_factor,
			top_lo_plot = pmax(0, top_pct - 1.96 * top_se) * scale_factor,
			top_hi_plot = pmax(0, top_pct + 1.96 * top_se) * scale_factor
		) %>%
		group_by(threshold_label) %>%
		mutate(is_endpoint = row_number() %in% c(1L, n()), endpoint_label = ifelse(is_endpoint, sprintf("%.2f", top_pct), NA_character_), endpoint_y = top_pct_plot + 0.045 * max(hist_df$count, na.rm = TRUE)) %>%
		ungroup()
	raw_top10 <- raw %>% mutate(phone_top10 = phone.sco >= fig2_top10_cut)
	trend <- purrr::map_dfr(seq_len(nrow(fig2_top_luck_cuts)), function(i) {
		cut_i <- fig2_top_luck_cuts$cutoff[i]
		lab_i <- fig2_top_luck_cuts$threshold_label[i]
		fit_fig2_top_trend(raw %>% mutate(phone_top = phone.sco >= cut_i)) %>%
			mutate(threshold = fig2_top_luck_cuts$threshold[i], threshold_label = lab_i, cutoff = cut_i, .before = 1)
	})
	fig2A_cols <- setNames(c("#009E73", "#0072B2", "#D55E00", "#CC79A7")[seq_along(fig2_top_luck_cuts$threshold_label)], fig2_top_luck_cuts$threshold_label)
	x_lab <- switch(transform, log = paste0("log(", value_label, ")"), log10 = paste0("log10(", value_label, ")"), log10p1 = paste0("log10(", value_label, " + 1)"), identity = value_label)
	p <- ggplot() +
		geom_rect(data = hist_df %>% filter(count > 0), aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count), fill = "grey86", color = "grey50", linewidth = .25) +
		geom_errorbar(data = bin_df, aes(x = x, ymin = top_lo_plot, ymax = top_hi_plot, color = threshold_label), width = 0, linewidth = .28, alpha = .72) +
		geom_line(data = bin_df, aes(x = x, y = top_pct_plot, color = threshold_label, group = threshold_label), linewidth = .62) +
		geom_point(data = bin_df, aes(x = x, y = top_pct_plot, color = threshold_label), size = 1.0) +
		geom_text(data = dplyr::filter(bin_df, .data$is_endpoint), aes(x = x, y = endpoint_y, label = endpoint_label), color = "blue", size = 2.35, fontface = "bold", vjust = 0) +
		scale_color_manual(values = fig2A_cols, name = NULL, drop = FALSE) +
		guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
		scale_y_continuous(name = NULL, sec.axis = sec_axis(~ . / scale_factor, name = "% top phone-luck score", breaks = seq(sec_axis_ylim[1], sec_axis_ylim[2], 5))) +
		coord_cartesian(ylim = c(0, sec_axis_ylim[2] * scale_factor * 1.08), clip = "off") +
		labs(title = sprintf("%s. Top phone luck vs %s", panel_letter, tolower(value_label)), x = x_lab) +
		fig_theme(base_size = 8.2) +
		theme(axis.title.y.right = element_text(color = "blue"), axis.text.y.right = element_text(color = "blue"), legend.position = "right", legend.direction = "vertical", legend.margin = margin(0, 0, 0, 4), legend.box.margin = margin(0, 0, 0, 0), legend.key.height = grid::unit(.62, "cm"), legend.key.width = grid::unit(.42, "cm"), legend.spacing.y = grid::unit(.26, "cm"), legend.text = element_text(size = 7.0, lineheight = 1.25), plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5, margin = margin(t = 0, b = 1)), plot.margin = margin(0, 2, 1, 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	list(plot = p, raw = raw_top10, bin = bin_df, hist = hist_df, trend = trend, trim_cut = trim_cut)
}
make_fig2_dx_quintile_panel <- function(base, value_col, value_label, panel_letter, trim_low_q = NULL) {
	if (!value_col %in% names(base)) {
		return(list(plot = empty_panel(sprintf("%s. Disease mix by %s quintile", panel_letter, tolower(value_label)), msg = paste("Missing variable:", value_label)), heat = tibble(), trim_cut = NA_real_))
	}
	dat <- base %>%
		filter(dx_grp %in% dxs.all, is.finite(.data[[value_col]])) %>%
		mutate(value_raw = as.numeric(.data[[value_col]]))
	trim_cut <- NA_real_
	if (!is.null(trim_low_q) && is.finite(trim_low_q) && trim_low_q > 0 && nrow(dat)) {
		trim_cut <- as.numeric(stats::quantile(dat$value_raw, trim_low_q, na.rm = TRUE, type = 7))
		dat <- dat %>% filter(value_raw >= trim_cut)
	}
	if (nrow(dat) < 100) {
		return(list(plot = empty_panel(sprintf("%s. Disease mix by %s quintile", panel_letter, tolower(value_label))), heat = tibble(), trim_cut = trim_cut))
	}
	dat <- dat %>%
		mutate(env_q = ntile(value_raw, 5), env_q = factor(paste0("Q", env_q), levels = paste0("Q", 1:5)))
	total_all <- nrow(dat)
	heat <- dat %>%
		count(env_q, dx_grp, name = "n") %>%
		complete(env_q = paste0("Q", 1:5), dx_grp = dxs.all, fill = list(n = 0)) %>%
		group_by(env_q) %>% mutate(q_total = sum(n), pct_q = n / q_total) %>% ungroup() %>%
		left_join(dat %>% count(dx_grp, name = "n_all") %>% mutate(pct_all = n_all / sum(n_all)), by = "dx_grp") %>%
		mutate(
			enrichment = pct_q / pct_all,
			n_rest = pmax(n_all - n, 0),
			total_rest = pmax(total_all - q_total, 0),
			RR_q_vs_rest = div1(n / q_total, n_rest / total_rest),
			se_logRR = sqrt(1 / pmax(n, 1) - 1 / pmax(q_total, 1) + 1 / pmax(n_rest, 1) - 1 / pmax(total_rest, 1)),
			RR_lo = exp(log(RR_q_vs_rest) - 1.96 * se_logRR),
			RR_hi = exp(log(RR_q_vs_rest) + 1.96 * se_logRR),
			p = purrr::pmap_dbl(list(n, q_total, n_rest, total_rest), function(a, A, b, B) suppressWarnings(chisq.test(matrix(c(a, A - a, b, B - b), nrow = 2))$p.value)),
			p_adj = p.adjust(p, "BH"),
			sig05 = sig_star05(p_adj),
			label = sprintf("%.2f%s", enrichment, sig05),
			dx_grp = factor(dx_grp, levels = rev(dxs.all)), variable = value_col, variable_label = value_label
		)
	p <- ggplot(heat, aes(env_q, dx_grp, fill = log(enrichment))) +
		geom_tile(color = "white", linewidth = .25) +
		geom_text(aes(label = label), size = 2.5) +
		scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = NULL, guide = guide_colorbar(direction = "vertical", barheight = grid::unit(2.3, "cm"), barwidth = grid::unit(.25, "cm"))) +
		labs(title = sprintf("%s. Disease mix by %s quintile", panel_letter, tolower(value_label)), x = paste0(value_label, " quintile"), y = NULL) +
		fig_theme(base_size = 8.2) +
		theme(plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5, margin = margin(t = 0, b = 6)), plot.margin = margin(0, 1, 1, 2), legend.position = "right", legend.direction = "vertical", legend.key.height = grid::unit(.36, "cm"), legend.key.width = grid::unit(.42, "cm"), legend.margin = margin(t = 0, r = 0, b = 1, l = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	list(plot = p, heat = heat, trim_cut = trim_cut)
}

fig2A_obj <- make_fig2_top10_panel(fig2_base, fig2_house_plot_var, "Housing price", "a", transform = house_price_scale)
p2A <- fig2A_obj$plot; fig2A_raw <- fig2A_obj$raw; fig2A_bin <- fig2A_obj$bin; fig2A_hist <- fig2A_obj$hist; fig2A_trend <- fig2A_obj$trend

fig2B_obj <- make_fig2_dx_quintile_panel(fig2_base, fig2_house_plot_var, "Housing price", "b")
p2B <- fig2B_obj$plot; fig2B_heat <- fig2B_obj$heat

# Phone-luck group difference across housing and environmental variables is kept for the output workbook.
fig2B_dat <- purrr::map_dfr(fig2_num_vars, function(v) {
	d <- fig2_base %>% filter(phone.luck %in% c("low", "high"), is.finite(.data[[v]]))
	if (nrow(d) < 100) return(tibble(variable = v, smd = NA_real_, low_mean = NA_real_, high_mean = NA_real_, p = NA_real_, n = nrow(d)))
	tibble(
		variable = v,
		smd = calc_smd(d[[v]], d$phone.luck),
		low_mean = mean(d[[v]][d$phone.luck == "low"], na.rm = TRUE),
		high_mean = mean(d[[v]][d$phone.luck == "high"], na.rm = TRUE),
		p = suppressWarnings(tryCatch(wilcox.test(d[[v]] ~ d$phone.luck)$p.value, error = function(e) NA_real_)),
		n = nrow(d)
	)
})
if (!"variable" %in% names(fig2B_dat)) {
	fig2B_dat <- tibble(variable = character(), smd = numeric(), low_mean = numeric(), high_mean = numeric(), p = numeric(), n = integer())
}
fig2B_dat <- fig2B_dat %>%
	mutate(variable_eng = geo_to_eng(as.character(.data$variable)), variable_eng = factor(variable_eng, levels = rev(geo_to_eng(fig2_num_vars))), p_adj = p.adjust(p, "BH"), sig = sig_star(p_adj))

fig2_map_tile_root <- file.path(dir.analysis, "map_tile")
fig2_house_bar_label <- switch(house_price_bar, sqrt = "sqrt housing price", log = "log housing price", log10 = "log10 housing price")

# c. Shenzhen map: static satellite mosaic with 3D housing-price columns.
fig2_map_dat <- tibble()
fig2c_dat <- tibble()
fig2_tile_status <- if (isTRUE(draw_2D_map)) "requested: static satellite mosaic" else "skipped: geo.map.use is none"
fig2c_html_file <- file.path(dir.analysis, "Fig2c.html")
fig2c_png_file <- file.path(dir.analysis, "Fig2c.png")
fig2c_html_width <- 4200L
fig2c_html_height <- 2200L
fig2c_render_mode <- NA_character_
fig2c_status <- fig2_tile_status
fig2c_map_alpha <- 1 - geo.map.transparent
fig2c_bar_height_scale <- 1.75
fig2c_bar_width_scale <- 2.40
p2C <- NULL
if (isTRUE(draw_2D_map)) {
	fig2_map_raw <- fig2_base %>%
		filter(is.finite(.data[[fig2_coord_vars[1]]]), is.finite(.data[[fig2_coord_vars[2]]]), is.finite(.data[[fig2_house_plot_var]]), is.finite(phone.sco))
	if (nrow(fig2_map_raw) >= 100 && is.finite(fig2_top10_cut)) {
		fig2_map_dat <- fig2_map_raw %>%
			mutate(phone_top10 = phone.sco >= fig2_top10_cut, lon = .data[[fig2_coord_vars[1]]], lat = .data[[fig2_coord_vars[2]]], lon_bin = round(lon, 3), lat_bin = round(lat, 3)) %>%
			filter(between(lon, 113.6, 114.8), between(lat, 22.35, 22.95)) %>%
			group_by(lon = lon_bin, lat = lat_bin) %>%
			summarise(n = n(), housing_price = median(.data[[fig2_house_plot_var]], na.rm = TRUE), top10_pct = mean(phone_top10, na.rm = TRUE), .groups = "drop") %>%
			mutate(
				housing_price_sqrt = if_else(housing_price > 0, sqrt(housing_price), NA_real_),
				housing_price_log = if_else(housing_price > 0, log(housing_price), NA_real_),
				housing_price_log10 = if_else(housing_price > 0, log10(housing_price), NA_real_),
				housing_price_bar = if (identical(house_price_bar, "log")) {
					housing_price_log
				} else if (identical(house_price_bar, "log10")) {
					housing_price_log10
				} else {
					housing_price_sqrt
				}
			) %>%
			filter(n >= 10, is.finite(housing_price), is.finite(housing_price_bar), is.finite(top10_pct)) %>%
			slice_max(n, n = 1500, with_ties = FALSE)
		if (nrow(fig2_map_dat) >= 20) {
			fig2c_dat <- fig2_map_dat %>%
				transmute(
					lon,
					lat,
					n,
					housing_price = round(housing_price, 2),
					housing_price_sqrt = round(housing_price_sqrt, 4),
					housing_price_log = round(housing_price_log, 4),
					housing_price_log10 = round(housing_price_log10, 4),
					housing_price_bar = round(housing_price_bar, 4),
					height_metric = fig2_house_bar_label,
					top10_pct = round(top10_pct, 5),
					top10_pct_pct = round(100 * top10_pct, 2),
					height = round(pmin(housing_price, fig2_house_price_cap["high"]) / 1000, 2)
				)
			fig2c_status <- "prepared: static satellite mosaic"
		} else {
			fig2_map_dat <- tibble()
			fig2c_status <- "skipped: insufficient coordinates"
		}
	} else {
		fig2c_status <- "skipped: insufficient coordinates"
	}
}

# d. Optional Fig2d HTML map export, controlled by geo.map.use.
fig2d_data_file <- NA_character_
fig2d_py_file <- NA_character_
fig2d_html_file <- NA_character_
fig2d_html_width <- 1600L
fig2d_html_height <- 900L
fig2d_status <- if (isTRUE(draw_3D_map)) "skipped: insufficient Fig2 map data" else "skipped: geo.map.use is none"
fig2d_python <- NA_character_
fig2d_dat <- tibble()
fig2d_map_dat <- tibble()
fig2d_config_file <- NA_character_
fig2d_render_mode <- if (isTRUE(draw_3D_map)) geo.map.use else NA_character_
py_quote <- function(x) paste0("r'", gsub("'", "\\\\'", gsub("\\\\", "/", x)), "'")
check_fig2d_python <- function() {
	pick_python <- function(paths) {
		paths <- paths[nzchar(paths) & file.exists(paths)]
		if (length(paths)) normalizePath(paths[1], winslash = "/", mustWork = TRUE) else NA_character_
	}
	if (!reticulate::py_available(initialize = FALSE)) {
		py_ai <- switch(
			Sys.info()[["sysname"]],
			Windows = pick_python(c(
				Sys.getenv("RETICULATE_PYTHON", ""),
				file.path(Sys.getenv("USERPROFILE"), "anaconda3", "envs", "ai", "python.exe"),
				file.path(Sys.getenv("USERPROFILE"), "miniconda3", "envs", "ai", "python.exe"),
				"C:/ProgramData/anaconda3/envs/ai/python.exe",
				"C:/ProgramData/miniconda3/envs/ai/python.exe",
				"D:/anaconda3/envs/ai/python.exe",
				"D:/miniconda3/envs/ai/python.exe"
			)),
			Linux = pick_python(c(
				Sys.getenv("RETICULATE_PYTHON", ""),
				file.path(Sys.getenv("HOME"), "anaconda3", "envs", "ai", "bin", "python"),
				file.path(Sys.getenv("HOME"), "miniconda3", "envs", "ai", "bin", "python"),
				"/home/huangj/anaconda3/envs/ai/bin/python",
				"/home/jiehuang001/anaconda3/envs/ai/bin/python"
			)),
			Darwin = pick_python(c(
				Sys.getenv("RETICULATE_PYTHON", ""),
				file.path(Sys.getenv("HOME"), "anaconda3", "envs", "ai", "bin", "python"),
				file.path(Sys.getenv("HOME"), "miniconda3", "envs", "ai", "bin", "python")
			)),
			NA_character_
		)
		if (is.na(py_ai) || !nzchar(py_ai)) {
			stop("geo.map.use requires Fig2d rendering, but conda env 'ai' Python could not be found. Activate conda env 'ai' or set RETICULATE_PYTHON.", call. = FALSE)
		}
		reticulate::use_python(py_ai, required = TRUE)
	}
	if (!reticulate::py_available(initialize = TRUE)) {
		stop("geo.map.use requires Fig2d rendering, but reticulate could not initialize Python. Activate the conda environment or set RETICULATE_PYTHON.", call. = FALSE)
	}
	py <- tryCatch(reticulate::py_config()$python, error = function(e) NA_character_)
	miss <- c()
	if (!reticulate::py_module_available("pandas")) miss <- c(miss, "pandas")
	if (!reticulate::py_module_available("keplergl")) miss <- c(miss, "keplergl")
	if (length(miss)) {
		stop(paste0(
			"geo.map.use requires Fig2d rendering with Python packages ", paste(miss, collapse = ", "),
			" in the active reticulate environment. Python: ", py,
			"\nInstall if needed, for example: conda activate ai && pip install keplergl"
		), call. = FALSE)
	}
	normalizePath(py, winslash = "/", mustWork = FALSE)
}

find_fig2_kepler_python <- function() {
	py_launcher_python <- function(version = "-3.12") {
		py_launcher <- unname(Sys.which("py"))
		if (is.na(py_launcher) || !nzchar(py_launcher)) return(character())
		out <- suppressWarnings(system2(py_launcher, c(version, "-c", shQuote("import sys; print(sys.executable)")), stdout = TRUE, stderr = FALSE))
		status <- attr(out, "status")
		if (is.null(status)) status <- 0L
		if (identical(status, 0L)) out[nzchar(out)] else character()
	}
	candidates <- c(
		Sys.getenv("FIG2D_PYTHON", ""),
		Sys.getenv("RETICULATE_PYTHON", ""),
		py_launcher_python("-3.12"),
		file.path(Sys.getenv("LOCALAPPDATA", unset = ""), "Programs", "Python", "Python312", "python.exe"),
		file.path(Sys.getenv("LOCALAPPDATA", unset = ""), "Programs", "Python", "Python314", "python.exe"),
		file.path(Sys.getenv("USERPROFILE", unset = ""), "anaconda3", "envs", "ai", "python.exe"),
		file.path(Sys.getenv("USERPROFILE", unset = ""), "miniconda3", "envs", "ai", "python.exe"),
		"C:/ProgramData/anaconda3/envs/ai/python.exe",
		"C:/ProgramData/miniconda3/envs/ai/python.exe",
		"D:/anaconda3/envs/ai/python.exe",
		"D:/miniconda3/envs/ai/python.exe",
		file.path(Sys.getenv("HOME", unset = ""), "anaconda3", "envs", "ai", "bin", "python"),
		file.path(Sys.getenv("HOME", unset = ""), "miniconda3", "envs", "ai", "bin", "python"),
		"/home/huangj/anaconda3/envs/ai/bin/python",
		"/home/jiehuang001/anaconda3/envs/ai/bin/python",
		unname(Sys.which(c("python3.12", "python3", "python")))
	)
	candidates <- unique(normalizePath(candidates[nzchar(candidates) & file.exists(candidates)], winslash = "/", mustWork = FALSE))
	failures <- character()
	for (py in candidates) {
		out <- suppressWarnings(system2(py, c("-c", shQuote("import pandas, keplergl")), stdout = TRUE, stderr = TRUE))
		status <- attr(out, "status")
		if (is.null(status)) status <- 0L
		if (identical(status, 0L)) return(py)
		failures <- c(failures, sprintf("%s: %s", py, paste(out, collapse = " | ")))
	}
	stop(paste0(
		"Kepler.gl map rendering requires Python packages pandas and keplergl. ",
		"Set FIG2D_PYTHON or RETICULATE_PYTHON to a Python executable with those packages installed.",
		if (length(failures)) paste0("\nChecked:\n", paste(failures, collapse = "\n")) else ""
	), call. = FALSE)
}

write_fig2d_keplergl_files <- function(dat, data_file, py_file, config_file, html_file, width = 1600L, height = 900L) {
	if (!requireNamespace("jsonlite", quietly = TRUE)) pacman::p_load(jsonlite)
	data_id <- "house"
	empty_object <- structure(list(), names = character(0))
	color_range <- list(
		name = "Fig2d teal yellow rose",
		type = "sequential",
		category = "Custom",
		colors = c("#009392", "#39B185", "#9CCB86", "#E9E29C", "#EEB479", "#E88471", "#CF597E"),
		reversed = FALSE
	)
	vis_config <- list(
		opacity = 0.78,
		worldUnitSize = 0.35,
		resolution = 8,
		colorRange = color_range,
		coverage = 0.82,
		sizeRange = c(0, 500),
		percentile = c(0, 100),
		elevationPercentile = c(0, 100),
		elevationScale = 38,
		colorAggregation = "average",
		sizeAggregation = "average",
		enable3d = TRUE
	)
	vis_config[["hi-precision"]] <- FALSE
	fig2d_config <- list(
		version = "v1",
		config = list(
			visState = list(
				filters = list(),
				layers = list(list(
					id = "fig2d_housing_hex",
					type = "hexagon",
					config = list(
						dataId = data_id,
						label = "Housing price columns",
						color = c(18, 147, 154),
						columns = list(lat = "lat", lng = "lon"),
						isVisible = TRUE,
						visConfig = vis_config,
						hidden = FALSE,
						textLabel = list(list(field = NA, color = c(255, 255, 255), size = 18, offset = c(0, 0), anchor = "start", alignment = "center"))
					),
					visualChannels = list(
						colorField = list(name = "top10_pct_pct", type = "real"),
						colorScale = "quantile",
						sizeField = list(name = "housing_price_bar", type = "real"),
						sizeScale = "linear"
					)
				)),
				interactionConfig = list(
					tooltip = list(fieldsToShow = setNames(list(c("housing_price", "height_metric", "top10_pct_pct", "n")), data_id), enabled = TRUE),
					brush = list(size = 0.5, enabled = FALSE),
					geocoder = list(enabled = FALSE),
					coordinate = list(enabled = TRUE)
				),
				layerBlending = "normal",
				splitMaps = list()
			),
			mapState = list(
				bearing = -28,
				dragRotate = TRUE,
				latitude = mean(dat$lat, na.rm = TRUE),
				longitude = mean(dat$lon, na.rm = TRUE),
				pitch = 56,
				zoom = 10.05,
				isSplit = FALSE
			),
			mapStyle = list(
				styleType = "carto_positron",
				topLayerGroups = empty_object,
				visibleLayerGroups = list(label = TRUE, road = TRUE, border = TRUE, building = FALSE, water = TRUE, land = TRUE, `3d building` = FALSE),
				mapStyles = list(carto_positron = list(
					id = "carto_positron",
					label = "CARTO Positron",
					url = "https://basemaps.cartocdn.com/gl/positron-gl-style/style.json",
					icon = "https://a.basemaps.cartocdn.com/light_all/0/0/0.png"
				))
			)
		)
	)
	dir.create(dirname(config_file), showWarnings = FALSE, recursive = TRUE)
	writeLines(jsonlite::toJSON(fig2d_config, auto_unbox = TRUE, pretty = TRUE, digits = 8, na = "null"), config_file, useBytes = TRUE)
	py_lines <- c(
		"from pathlib import Path",
		"import json",
		"import pandas as pd",
		"from keplergl import KeplerGl",
		sprintf("DATA_FILE = Path(%s)", py_quote(data_file)),
		sprintf("CONFIG_FILE = Path(%s)", py_quote(config_file)),
		sprintf("HTML_FILE = Path(%s)", py_quote(html_file)),
		sprintf("DATA_ID = %s", py_quote(data_id)),
		sprintf("MAP_HEIGHT = %d", as.integer(height)),
		"df = pd.read_csv(DATA_FILE, sep='\\t')",
		"df = df.dropna(subset=['lon', 'lat', 'housing_price_bar', 'top10_pct_pct'])",
		"with CONFIG_FILE.open('r', encoding='utf-8') as f:",
		"    config = json.load(f)",
		"m = KeplerGl(height=MAP_HEIGHT, config=config, show_docs=False)",
		"m.add_data(data=df.copy(), name=DATA_ID)",
		"m.save_to_html(file_name=str(HTML_FILE), read_only=False, config=config)",
		"fit_css = '''<style id=\"fig2d-fit-viewport\">",
		"html, body { width: 100%; height: 100%; margin: 0; overflow: hidden; }",
		"#app-content { position: fixed !important; inset: 0 !important; width: 100vw !important; height: 100vh !important; }",
		"#app-content > div, #app-content .keplergl-widget-container, #kepler-gl__map, #kepler-gl__map .maps, #kepler-gl__map [id^=\"view-\"], #kepler-gl__map .maplibregl-map, #default-deckgl-overlay-wrapper, #default-deckgl-overlay { width: 100vw !important; height: 100vh !important; max-width: none !important; max-height: none !important; }",
		"#kepler-gl__map canvas { width: 100vw !important; height: 100vh !important; }",
		"</style>'''",
		"html = HTML_FILE.read_text(encoding='utf-8')",
		"if 'fig2d-fit-viewport' not in html:",
		"    html = html.replace('</head>', fit_css + '</head>', 1)",
		"    HTML_FILE.write_text(html, encoding='utf-8')",
		"if (not HTML_FILE.exists()) or HTML_FILE.stat().st_size < 1024:",
		"    raise RuntimeError(f'KeplerGl did not create a usable HTML file: {HTML_FILE}')",
		"print(f'Map saved to {HTML_FILE} ({HTML_FILE.stat().st_size} bytes)')"
	)
	writeLines(py_lines, py_file, useBytes = TRUE)
	invisible(list(config_file = config_file, py_file = py_file, html_file = html_file, data_id = data_id))
}

write_fig2_satellite_canvas <- function(dat, file, width = 1600L, height = 900L,
	tile_cache_dir = file.path(dirname(file), "map_tile", "3D", "Fig2d_satellite_tiles"),
	mosaic_stem = sub("\\.[^.]*$", "", basename(file)),
	html_title = "Fig2d satellite 3D housing price",
	note_text = "Fig2d: 3D housing price on Esri satellite imagery",
	show_legend = TRUE,
	show_note = TRUE,
	tile_pad = 1L,
	map_alpha = 1, bar_height_scale = 1.05, bar_width_scale = 1.70) {
	if (!requireNamespace("jsonlite", quietly = TRUE)) pacman::p_load(jsonlite)
	if (!requireNamespace("base64enc", quietly = TRUE)) pacman::p_load(base64enc)
	if (!requireNamespace("jpeg", quietly = TRUE)) pacman::p_load(jpeg)
	html_escape <- function(x) {
		x <- gsub("&", "&amp;", x, fixed = TRUE)
		x <- gsub("<", "&lt;", x, fixed = TRUE)
		x <- gsub(">", "&gt;", x, fixed = TRUE)
		x <- gsub('"', "&quot;", x, fixed = TRUE)
		x
	}
	dat_json <- jsonlite::toJSON(
		dat %>% transmute(lon, lat, housing_price, height, top10_pct_pct, n),
		dataframe = "rows",
		auto_unbox = TRUE,
		digits = 7,
		na = "null"
	)
	map_alpha <- suppressWarnings(as.numeric(map_alpha)[1])
	if (!is.finite(map_alpha)) map_alpha <- 1
	map_alpha <- pmin(pmax(map_alpha, 0), 1)
	bar_height_scale <- suppressWarnings(as.numeric(bar_height_scale)[1])
	if (!is.finite(bar_height_scale) || bar_height_scale <= 0) bar_height_scale <- 1.05
	bar_width_scale <- suppressWarnings(as.numeric(bar_width_scale)[1])
	if (!is.finite(bar_width_scale) || bar_width_scale <= 0) bar_width_scale <- 1.70
	z <- 10L
	tile_size <- 256L
	world <- tile_size * 2^z
	mx_r <- function(lon) (lon + 180) / 360 * world
	my_r <- function(lat) {
		s <- sin(lat * pi / 180)
		(0.5 - log((1 + s) / (1 - s)) / (4 * pi)) * world
	}
	xs <- mx_r(dat$lon)
	ys <- my_r(dat$lat)
	tile_pad <- max(0L, as.integer(tile_pad)[1])
	tx_range <- seq.int(max(0L, floor(min(xs, na.rm = TRUE) / tile_size) - tile_pad), min(2^z - 1L, floor(max(xs, na.rm = TRUE) / tile_size) + tile_pad))
	ty_range <- seq.int(max(0L, floor(min(ys, na.rm = TRUE) / tile_size) - tile_pad), min(2^z - 1L, floor(max(ys, na.rm = TRUE) / tile_size) + tile_pad))
	dir.create(tile_cache_dir, showWarnings = FALSE, recursive = TRUE)
	tile_ok <- function(path) file.exists(path) && !is.na(file.info(path)$size) && file.info(path)$size >= 1024
	ps_quote <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
	ps_local_path <- function(x) {
		x <- gsub("\\\\", "/", x)
		m <- regexec("^/mnt/([A-Za-z])/(.*)$", x)
		hit <- regmatches(x, m)[[1]]
		if (length(hit) == 3L) return(paste0(toupper(hit[2]), ":\\", gsub("/", "\\\\", hit[3])))
		x
	}
	download_tile <- function(url, dest) {
		ok <- FALSE
		ps_bins <- unname(Sys.which(c("powershell.exe", "powershell", "pwsh.exe", "pwsh")))
		ps_bin <- ps_bins[nzchar(ps_bins)][1]
		if (!is.na(ps_bin) && nzchar(ps_bin)) {
			cmd <- sprintf(
				"$ProgressPreference='SilentlyContinue'; try { Invoke-WebRequest -Uri %s -OutFile %s -TimeoutSec 30; exit 0 } catch { exit 1 }",
				ps_quote(url), ps_quote(ps_local_path(dest))
			)
			cmd_encoded <- base64enc::base64encode(iconv(cmd, from = "", to = "UTF-16LE", toRaw = TRUE)[[1]], linewidth = 0)
			ok <- identical(suppressWarnings(system2(ps_bin, c("-NoProfile", "-ExecutionPolicy", "Bypass", "-EncodedCommand", cmd_encoded), stdout = FALSE, stderr = FALSE)), 0L)
		}
		if (!ok) ok <- !inherits(suppressWarnings(try(utils::download.file(url, dest, mode = "wb", quiet = TRUE), silent = TRUE)), "try-error")
		tile_ok(dest)
	}
	tile_files <- matrix(NA_character_, nrow = length(ty_range), ncol = length(tx_range), dimnames = list(as.character(ty_range), as.character(tx_range)))
	for (tx in tx_range) {
		for (ty in ty_range) {
			tile_file <- file.path(tile_cache_dir, sprintf("z%d_x%d_y%d.jpg", z, tx, ty))
			if (!tile_ok(tile_file)) {
				tile_url <- sprintf("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/%d/%d/%d", z, ty, tx)
				download_tile(tile_url, tile_file)
			}
			if (tile_ok(tile_file)) tile_files[as.character(ty), as.character(tx)] <- tile_file
		}
	}
	mosaic <- array(48 / 255, dim = c(length(ty_range) * tile_size, length(tx_range) * tile_size, 3))
	mosaic_tiles <- 0L
	for (tx_i in seq_along(tx_range)) {
		for (ty_i in seq_along(ty_range)) {
			tile_file <- tile_files[ty_i, tx_i]
			if (is.na(tile_file)) next
			img <- try(jpeg::readJPEG(tile_file), silent = TRUE)
			if (inherits(img, "try-error")) next
			if (length(dim(img)) == 2L) img <- array(rep(img, 3), dim = c(dim(img), 3))
			if (dim(img)[3] > 3L) img <- img[, , 1:3, drop = FALSE]
			y0 <- (ty_i - 1L) * tile_size + 1L
			x0 <- (tx_i - 1L) * tile_size + 1L
			mosaic[y0:(y0 + tile_size - 1L), x0:(x0 + tile_size - 1L), ] <- img[seq_len(tile_size), seq_len(tile_size), 1:3]
			mosaic_tiles <- mosaic_tiles + 1L
		}
	}
	mosaic_file <- file.path(tile_cache_dir, sprintf("%s_satellite_mosaic_z%d.jpg", mosaic_stem, z))
	mosaic_src <- ""
	if (mosaic_tiles > 0L) {
		jpeg::writeJPEG(mosaic, target = mosaic_file, quality = 0.9)
		mosaic_src <- paste0("data:image/jpeg;base64,", base64enc::base64encode(mosaic_file))
	}
	satellite_json <- jsonlite::toJSON(list(
		has = mosaic_tiles > 0L,
		tx0 = as.integer(min(tx_range)),
		ty0 = as.integer(min(ty_range)),
		nx = as.integer(length(tx_range)),
		ny = as.integer(length(ty_range))
	), auto_unbox = TRUE)
	html <- c(
		"<!doctype html>",
		"<html>",
		"<head>",
		"  <meta charset=\"utf-8\"/>",
		"  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>",
		sprintf("  <title>%s</title>", html_escape(html_title)),
		"  <style>",
		"    html, body { margin: 0; width: 100%; height: 100%; overflow: hidden; background: #101418; font-family: Arial, sans-serif; }",
		sprintf("    #wrap { position: relative; width: %dpx; height: %dpx; max-width: 100vw; max-height: 100vh; overflow: hidden; background: linear-gradient(135deg, #20313a, #53646b); }", as.integer(width), as.integer(height)),
		"    #sat { position: absolute; display: none; max-width: none; pointer-events: none; user-select: none; }",
		"    canvas { position: absolute; left: 0; top: 0; display: block; width: 100%; height: 100%; }",
		"    .note, .legend, .tip { position: absolute; color: white; background: rgba(0,0,0,.50); padding: 10px 14px; font-size: 24px; border-radius: 4px; pointer-events: none; }",
		"    .note { left: 18px; bottom: 16px; font-size: 20px; }",
		"    .legend { right: 18px; top: 18px; line-height: 1.35; min-width: 430px; }",
		"    .legend-title { font-weight: bold; margin-bottom: 8px; }",
		"    .legend-bar { height: 18px; border: 2px solid rgba(255,255,255,.65); background: linear-gradient(90deg, #009392, #39B185, #9CCB86, #E9E29C, #EEB479, #E88471, #CF597E); }",
		"    .legend-ticks { display: flex; justify-content: space-between; font-size: 20px; margin-top: 4px; }",
		"    .legend-height { margin-top: 8px; }",
		"    .tip { display: none; white-space: pre; }",
		"  </style>",
		"</head>",
		"<body>",
		"  <div id=\"wrap\">",
		sprintf("    <img id=\"sat\" src=\"%s\" alt=\"\"/>", mosaic_src),
		"    <canvas id=\"map\"></canvas>",
		if (isTRUE(show_legend)) "    <div class=\"legend\"><div class=\"legend-title\">% top-decile phone luck</div><div class=\"legend-bar\"></div><div class=\"legend-ticks\"><span>0%</span><span>10%</span><span>20%</span><span>30%</span></div><div class=\"legend-height\">Height: housing price / 1000</div></div>" else character(),
		if (isTRUE(show_note)) sprintf("    <div class=\"note\">%s</div>", html_escape(note_text)) else character(),
		"    <div id=\"tip\" class=\"tip\"></div>",
		"  </div>",
		"  <script>",
		paste0("    const DATA = ", dat_json, ";"),
		paste0("    const SATELLITE = ", satellite_json, ";"),
		sprintf("    const MAP_ALPHA = %.6f;", map_alpha),
		sprintf("    const BAR_HEIGHT_SCALE = %.6f;", bar_height_scale),
		sprintf("    const BAR_WIDTH_SCALE = %.6f;", bar_width_scale),
		"    const canvas = document.getElementById('map');",
		"    const wrap = document.getElementById('wrap');",
		"    const sat = document.getElementById('sat');",
		"    const tip = document.getElementById('tip');",
		"    const ctx = canvas.getContext('2d');",
		"    let placed = [];",
		"    function resizeCanvas() {",
		"      const dpr = window.devicePixelRatio || 1;",
		"      const rect = wrap.getBoundingClientRect();",
		"      canvas.width = Math.max(1, Math.round(rect.width * dpr));",
		"      canvas.height = Math.max(1, Math.round(rect.height * dpr));",
		"      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);",
		"    }",
		"    function colorFor(v, alpha = 220) {",
		"      if (v == null || Number.isNaN(v)) return [110, 110, 130, 190];",
		"      const stops = [[0, [0,147,146]], [5, [57,177,133]], [10, [156,203,134]], [15, [233,226,156]], [20, [238,180,121]], [25, [232,132,113]], [30, [207,89,126]]];",
		"      for (let i = 1; i < stops.length; i++) {",
		"        if (v <= stops[i][0]) {",
		"          const a = stops[i - 1], b = stops[i], t = (v - a[0]) / Math.max(b[0] - a[0], 1e-9);",
		"          return [0, 1, 2].map(k => Math.round(a[1][k] + t * (b[1][k] - a[1][k]))).concat(alpha);",
		"        }",
		"      }",
		"      return [207, 89, 126, alpha];",
		"    }",
		"    function rgba(c) { return `rgba(${c[0]},${c[1]},${c[2]},${(c[3] || 255) / 255})`; }",
		"    const Z = 10, TILE = 256, WORLD = TILE * Math.pow(2, Z);",
		"    function mx(lon) { return (lon + 180) / 360 * WORLD; }",
		"    function my(lat) { const s = Math.sin(lat * Math.PI / 180); return (0.5 - Math.log((1 + s) / (1 - s)) / (4 * Math.PI)) * WORLD; }",
		"    function layout() {",
		"      const rect = wrap.getBoundingClientRect();",
		"      const pts = DATA.map(d => ({d, wx: mx(d.lon), wy: my(d.lat)}));",
		"      const xs = pts.map(p => p.wx), ys = pts.map(p => p.wy);",
		"      const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);",
		"      const margin = 46;",
		"      const scale = Math.min((rect.width - 2 * margin) / (maxX - minX), (rect.height - 2 * margin) / (maxY - minY)) * 0.96;",
		"      const ox = (rect.width - (maxX - minX) * scale) / 2;",
		"      const oy = (rect.height - (maxY - minY) * scale) / 2;",
		"      pts.forEach(p => { p.x = ox + (p.wx - minX) * scale; p.y = oy + (p.wy - minY) * scale; });",
		"      return {rect, pts, minX, maxX, minY, maxY, scale, ox, oy};",
		"    }",
		"    function drawBase(L) {",
		"      ctx.clearRect(0, 0, L.rect.width, L.rect.height);",
		"      if (SATELLITE.has) {",
		"        sat.style.display = 'block';",
		"        sat.style.opacity = MAP_ALPHA;",
		"        sat.style.left = `${L.ox + (SATELLITE.tx0 * TILE - L.minX) * L.scale}px`;",
		"        sat.style.top = `${L.oy + (SATELLITE.ty0 * TILE - L.minY) * L.scale}px`;",
		"        sat.style.width = `${SATELLITE.nx * TILE * L.scale}px`;",
		"        sat.style.height = `${SATELLITE.ny * TILE * L.scale}px`;",
		"      } else {",
		"        sat.style.display = 'none';",
		"        const g = ctx.createLinearGradient(0, 0, L.rect.width, L.rect.height);",
		"        g.addColorStop(0, '#20313a'); g.addColorStop(1, '#53646b');",
		"        ctx.fillStyle = g; ctx.fillRect(0, 0, L.rect.width, L.rect.height);",
		"      }",
		"    }",
		"    function drawColumns(L) {",
		"      placed = [];",
		"      const pts = L.pts.slice().sort((a, b) => a.y - b.y);",
		"      ctx.save();",
		"      ctx.lineCap = 'round';",
		"      pts.forEach(p => {",
		"        const d = p.d;",
		"        const h = Math.max(8, Math.min(240, d.height * 0.52 * BAR_HEIGHT_SCALE));",
		"        const c = colorFor(d.top10_pct_pct, 238);",
		"        const sx = p.x + h * 0.10, sy = p.y - h;",
		"        ctx.strokeStyle = 'rgba(20,25,35,.38)'; ctx.lineWidth = 4.3 * BAR_WIDTH_SCALE; ctx.beginPath(); ctx.moveTo(p.x + 2, p.y + 2); ctx.lineTo(sx + 2, sy + 2); ctx.stroke();",
		"        ctx.strokeStyle = rgba(c); ctx.lineWidth = 2.6 * BAR_WIDTH_SCALE; ctx.beginPath(); ctx.moveTo(p.x, p.y); ctx.lineTo(sx, sy); ctx.stroke();",
		"        ctx.fillStyle = rgba(colorFor(d.top10_pct_pct, 255)); ctx.beginPath(); ctx.ellipse(sx, sy, 2.7 * BAR_WIDTH_SCALE, 1.6 * BAR_WIDTH_SCALE, 0, 0, Math.PI * 2); ctx.fill();",
		"        placed.push({x: sx, y: sy, baseX: p.x, baseY: p.y, d});",
		"      });",
		"      ctx.restore();",
		"    }",
		"    function render() { resizeCanvas(); const L = layout(); drawBase(L); drawColumns(L); }",
		"    canvas.addEventListener('mousemove', ev => {",
		"      const r = canvas.getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top;",
		"      let best = null, bestDist = 999;",
		"      for (const p of placed) { const dist = Math.hypot(p.x - x, p.y - y); if (dist < bestDist) { best = p; bestDist = dist; } }",
		"      if (best && bestDist < 10) {",
		"        tip.style.display = 'block'; tip.style.left = `${x + 12}px`; tip.style.top = `${y + 12}px`;",
		"        tip.textContent = `Housing price: ${Math.round(best.d.housing_price).toLocaleString()}\\nHeight: ${best.d.height}\\nTop-decile phone luck: ${best.d.top10_pct_pct}%\\nN: ${best.d.n}`;",
		"      } else { tip.style.display = 'none'; }",
		"    });",
		"    canvas.addEventListener('mouseleave', () => { tip.style.display = 'none'; });",
		"    window.addEventListener('resize', render);",
		"    render();",
		"  </script>",
		"</body>",
		"</html>"
	)
	writeLines(html, file, useBytes = TRUE)
	invisible(file)
}
write_fig2d_canvas_fallback <- function(dat, file, width = 1600L, height = 900L,
	tile_cache_dir = file.path(dirname(file), "map_tile", "3D", "Fig2d_satellite_tiles")) {
	write_fig2_satellite_canvas(
		dat = dat,
		file = file,
		width = width,
		height = height,
		tile_cache_dir = tile_cache_dir,
		mosaic_stem = "Fig2d",
		html_title = "Fig2d satellite 3D housing price",
		note_text = "Fig2d: 3D housing price on Esri satellite imagery"
	)
}
write_fig2d_deckgl_html <- function(dat, file, width = 1600L, height = 900L,
	tile_cache_dir = NULL, mosaic_stem = "Fig2c", tile_pad = 1L) {
	if (!requireNamespace("jsonlite", quietly = TRUE)) pacman::p_load(jsonlite)
	if (!requireNamespace("base64enc", quietly = TRUE)) pacman::p_load(base64enc)
	html_escape <- function(x) {
		x <- gsub("&", "&amp;", x, fixed = TRUE)
		x <- gsub("<", "&lt;", x, fixed = TRUE)
		x <- gsub(">", "&gt;", x, fixed = TRUE)
		x <- gsub('"', "&quot;", x, fixed = TRUE)
		x
	}
	tile_x_to_lon <- function(x, z) x / 2^z * 360 - 180
	tile_y_to_lat <- function(y, z) atan(sinh(pi * (1 - 2 * y / 2^z))) * 180 / pi
	satellite_src <- ""
	satellite_bounds <- c(NA_real_, NA_real_, NA_real_, NA_real_)
	if (!is.null(tile_cache_dir) && dir.exists(tile_cache_dir)) {
		z <- 10L
		tile_size <- 256L
		world <- tile_size * 2^z
		mx_r <- function(lon) (lon + 180) / 360 * world
		my_r <- function(lat) {
			s <- sin(lat * pi / 180)
			(0.5 - log((1 + s) / (1 - s)) / (4 * pi)) * world
		}
		xs <- mx_r(dat$lon)
		ys <- my_r(dat$lat)
		tile_pad <- max(0L, as.integer(tile_pad)[1])
		tx_range <- seq.int(max(0L, floor(min(xs, na.rm = TRUE) / tile_size) - tile_pad), min(2^z - 1L, floor(max(xs, na.rm = TRUE) / tile_size) + tile_pad))
		ty_range <- seq.int(max(0L, floor(min(ys, na.rm = TRUE) / tile_size) - tile_pad), min(2^z - 1L, floor(max(ys, na.rm = TRUE) / tile_size) + tile_pad))
		mosaic_file <- file.path(tile_cache_dir, sprintf("%s_satellite_mosaic_z%d.jpg", mosaic_stem, z))
		if (file.exists(mosaic_file) && !is.na(file.info(mosaic_file)$size) && file.info(mosaic_file)$size > 1024) {
			satellite_src <- paste0("data:image/jpeg;base64,", base64enc::base64encode(mosaic_file))
			west <- tile_x_to_lon(min(tx_range), z)
			east <- tile_x_to_lon(max(tx_range) + 1L, z)
			north <- tile_y_to_lat(min(ty_range), z)
			south <- tile_y_to_lat(max(ty_range) + 1L, z)
			satellite_bounds <- c(west, south, east, north)
		}
	}
	satellite_json <- jsonlite::toJSON(list(
		has = nzchar(satellite_src),
		src = satellite_src,
		bounds = unname(satellite_bounds)
	), auto_unbox = TRUE, digits = 10, na = "null")
	dat_json <- jsonlite::toJSON(
		dat %>% transmute(lon, lat, housing_price, height, height_metric, top10_pct_pct, n),
		dataframe = "rows",
		auto_unbox = TRUE,
		digits = 7,
		na = "null"
	)
	center_lon <- mean(dat$lon, na.rm = TRUE)
	center_lat <- mean(dat$lat, na.rm = TRUE)
	html <- c(
		"<!doctype html>",
		"<html>",
		"<head>",
		"  <meta charset=\"utf-8\"/>",
		"  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>",
		"  <title>Fig2d interactive 3D satellite map</title>",
		"  <script src=\"https://unpkg.com/deck.gl@8.9.36/dist.min.js\"></script>",
		"  <style>",
		"    html, body { margin: 0; width: 100%; height: 100%; overflow: hidden; background: #07090b; font-family: Arial, sans-serif; }",
		"    #map { position: fixed; inset: 0; width: 100vw; height: 100vh; }",
		"    #loading { position: fixed; left: 16px; bottom: 14px; color: white; background: rgba(0,0,0,.55); padding: 8px 11px; font-size: 13px; border-radius: 3px; pointer-events: none; }",
		"    .deck-tooltip { font-family: Arial, sans-serif; font-size: 12px; line-height: 1.35; }",
		"  </style>",
		"</head>",
		"<body>",
		"  <div id=\"map\"></div>",
		"  <div id=\"loading\">Fig2d: interactive 3D housing price on Esri satellite imagery</div>",
		"  <script>",
		paste0("    const DATA = ", dat_json, ";"),
		paste0("    const SATELLITE = ", satellite_json, ";"),
		sprintf("    const CENTER_LON = %.8f;", center_lon),
		sprintf("    const CENTER_LAT = %.8f;", center_lat),
		"    const TILE_URL = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}';",
		"    const colorStops = [[0, [0,147,146]], [5, [57,177,133]], [10, [156,203,134]], [15, [233,226,156]], [20, [238,180,121]], [25, [232,132,113]], [30, [207,89,126]]];",
		"    function colorFor(v, alpha = 205) {",
		"      if (v == null || Number.isNaN(v)) return [160, 165, 145, alpha];",
		"      for (let i = 1; i < colorStops.length; i++) {",
		"        if (v <= colorStops[i][0]) {",
		"          const a = colorStops[i - 1], b = colorStops[i];",
		"          const t = (v - a[0]) / Math.max(b[0] - a[0], 1e-9);",
		"          return [0, 1, 2].map(k => Math.round(a[1][k] + t * (b[1][k] - a[1][k]))).concat(alpha);",
		"        }",
		"      }",
		"      return [207, 89, 126, alpha];",
		"    }",
		"    function embeddedSatelliteLayer() {",
		"      if (!SATELLITE.has) return null;",
		"      return new deck.BitmapLayer({",
		"        id: 'embedded-esri-satellite-mosaic',",
		"        image: SATELLITE.src,",
		"        bounds: SATELLITE.bounds",
		"      });",
		"    }",
		"    function satelliteLayer() {",
		"      return new deck.TileLayer({",
		"        id: 'esri-world-imagery',",
		"        data: TILE_URL,",
		"        minZoom: 0,",
		"        maxZoom: 19,",
		"        tileSize: 256,",
		"        refinementStrategy: 'best-available',",
		"        getTileData: ({index}) => {",
		"          const {x, y, z} = index;",
		"          return TILE_URL.replace('{z}', z).replace('{y}', y).replace('{x}', x);",
		"        },",
		"        renderSubLayers: props => {",
		"          if (!props.data) return null;",
		"          const bbox = props.tile.bbox;",
		"          const bounds = [bbox.west, bbox.south, bbox.east, bbox.north];",
		"          return new deck.BitmapLayer(props, {",
		"            id: `${props.id}-bitmap`,",
		"            image: props.data,",
		"            bounds",
		"          });",
		"        }",
		"      });",
		"    }",
		"    function layers() {",
		"      return [embeddedSatelliteLayer(), satelliteLayer(),",
		"        new deck.ColumnLayer({",
		"          id: 'housing-price-columns',",
		"          data: DATA,",
		"          pickable: true,",
		"          diskResolution: 12,",
		"          radius: 115,",
		"          extruded: true,",
		"          elevationScale: 55,",
		"          opacity: 0.82,",
		"          getPosition: d => [d.lon, d.lat],",
		"          getElevation: d => Math.max(1, d.height || 1),",
		"          getFillColor: d => colorFor(d.top10_pct_pct, 205),",
		"          material: {ambient: 0.42, diffuse: 0.62, shininess: 24, specularColor: [255, 255, 255]}",
		"        })",
		"      ];",
		"    }",
		"    const deckgl = new deck.DeckGL({",
		"      container: 'map',",
		"      views: new deck.MapView({repeat: false}),",
		"      initialViewState: {",
		"        longitude: CENTER_LON,",
		"        latitude: CENTER_LAT,",
		"        zoom: 10.35,",
		"        pitch: 60,",
		"        bearing: -28,",
		"        minZoom: 8,",
		"        maxZoom: 15",
		"      },",
		"      controller: true,",
		"      layers: layers(),",
		"      parameters: {depthTest: true},",
		"      getTooltip: ({object}) => object && {",
		"        html: `<b>Housing price</b>: ${Math.round(object.housing_price).toLocaleString()}<br/><b>Height</b>: ${object.height}<br/><b>Top-decile phone luck</b>: ${object.top10_pct_pct}%<br/><b>N</b>: ${object.n}`,",
		"        style: {backgroundColor: 'rgba(0,0,0,.72)', color: 'white'}",
		"      },",
		"      onLoad: () => { const el = document.getElementById('loading'); if (el) el.style.display = 'none'; }",
		"    });",
		"  </script>",
		"</body>",
		"</html>"
	)
	writeLines(html, file, useBytes = TRUE)
	invisible(file)
}
write_fig2d_deckgl_hex_html <- function(dat, file, width = 1600L, height = 900L,
	tile_cache_dir = NULL, mosaic_stem = "Fig2d", tile_pad = 1L) {
	if (!requireNamespace("jsonlite", quietly = TRUE)) pacman::p_load(jsonlite)
	if (!requireNamespace("base64enc", quietly = TRUE)) pacman::p_load(base64enc)
	html_escape <- function(x) {
		x <- gsub("&", "&amp;", x, fixed = TRUE)
		x <- gsub("<", "&lt;", x, fixed = TRUE)
		x <- gsub(">", "&gt;", x, fixed = TRUE)
		x <- gsub('"', "&quot;", x, fixed = TRUE)
		x
	}
	tile_x_to_lon <- function(x, z) x / 2^z * 360 - 180
	tile_y_to_lat <- function(y, z) atan(sinh(pi * (1 - 2 * y / 2^z))) * 180 / pi
	satellite_src <- ""
	satellite_bounds <- c(NA_real_, NA_real_, NA_real_, NA_real_)
	if (!is.null(tile_cache_dir) && dir.exists(tile_cache_dir)) {
		z <- 10L
		tile_size <- 256L
		world <- tile_size * 2^z
		mx_r <- function(lon) (lon + 180) / 360 * world
		my_r <- function(lat) {
			s <- sin(lat * pi / 180)
			(0.5 - log((1 + s) / (1 - s)) / (4 * pi)) * world
		}
		ok <- is.finite(dat$lon) & is.finite(dat$lat)
		if (any(ok)) {
			xs <- mx_r(dat$lon[ok])
			ys <- my_r(dat$lat[ok])
			tile_pad <- max(0L, as.integer(tile_pad)[1])
			tx_range <- seq.int(max(0L, floor(min(xs, na.rm = TRUE) / tile_size) - tile_pad), min(2^z - 1L, floor(max(xs, na.rm = TRUE) / tile_size) + tile_pad))
			ty_range <- seq.int(max(0L, floor(min(ys, na.rm = TRUE) / tile_size) - tile_pad), min(2^z - 1L, floor(max(ys, na.rm = TRUE) / tile_size) + tile_pad))
			mosaic_file <- file.path(tile_cache_dir, sprintf("%s_satellite_mosaic_z%d.jpg", mosaic_stem, z))
			if (file.exists(mosaic_file) && !is.na(file.info(mosaic_file)$size) && file.info(mosaic_file)$size > 1024) {
				satellite_src <- paste0("data:image/jpeg;base64,", base64enc::base64encode(mosaic_file))
				west <- tile_x_to_lon(min(tx_range), z)
				east <- tile_x_to_lon(max(tx_range) + 1L, z)
				north <- tile_y_to_lat(min(ty_range), z)
				south <- tile_y_to_lat(max(ty_range) + 1L, z)
				satellite_bounds <- c(west, south, east, north)
			}
		}
	}
	dat_out <- data.frame(
		lon = suppressWarnings(as.numeric(dat$lon)),
		lat = suppressWarnings(as.numeric(dat$lat)),
		n = suppressWarnings(as.numeric(dat$n)),
		housing_price = suppressWarnings(as.numeric(dat$housing_price)),
		housing_price_bar = suppressWarnings(as.numeric(dat$housing_price_bar)),
		height_metric = as.character(dat$height_metric),
		top10_pct_pct = suppressWarnings(as.numeric(dat$top10_pct_pct)),
		height = suppressWarnings(as.numeric(dat$height)),
		stringsAsFactors = FALSE
	)
	dat_out <- dat_out[is.finite(dat_out$lon) & is.finite(dat_out$lat) & is.finite(dat_out$housing_price_bar) & is.finite(dat_out$top10_pct_pct), , drop = FALSE]
	if (!nrow(dat_out)) stop("Fig2d has no finite rows to render.", call. = FALSE)
	satellite_json <- jsonlite::toJSON(list(
		has = nzchar(satellite_src),
		src = satellite_src,
		bounds = unname(satellite_bounds)
	), auto_unbox = TRUE, digits = 10, na = "null")
	dat_json <- jsonlite::toJSON(dat_out, dataframe = "rows", auto_unbox = TRUE, digits = 7, na = "null")
	center_lon <- mean(dat_out$lon, na.rm = TRUE)
	center_lat <- mean(dat_out$lat, na.rm = TRUE)
	html <- c(
		"<!doctype html>",
		"<html>",
		"<head>",
		"  <meta charset=\"utf-8\"/>",
		"  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"/>",
		"  <title>Fig2d interactive 3D hexagon satellite map</title>",
		"  <script src=\"https://unpkg.com/deck.gl@8.9.36/dist.min.js\"></script>",
		"  <style>",
		"    html, body { margin: 0; width: 100%; height: 100%; overflow: hidden; background: #07090b; font-family: Arial, sans-serif; }",
		"    #map { position: fixed; inset: 0; width: 100vw; height: 100vh; }",
		"    #loading { position: fixed; left: 16px; bottom: 14px; color: white; background: rgba(0,0,0,.58); padding: 8px 11px; font-size: 13px; border-radius: 3px; pointer-events: none; }",
		"    #legend { position: fixed; right: 16px; top: 16px; width: min(360px, calc(100vw - 32px)); color: white; background: rgba(0,0,0,.58); padding: 12px 14px 10px; border-radius: 3px; font-size: 12px; line-height: 1.35; }",
		"    #legend-title { font-weight: 700; margin-bottom: 8px; }",
		"    #legend-bar { height: 14px; border: 1px solid rgba(255,255,255,.72); background: linear-gradient(90deg, #009392, #39B185, #9CCB86, #E9E29C, #EEB479, #E88471, #CF597E); }",
		"    #legend-ticks { display: flex; justify-content: space-between; margin-top: 5px; color: rgba(255,255,255,.86); }",
		"    #attrib { position: fixed; right: 10px; bottom: 8px; color: rgba(255,255,255,.78); background: rgba(0,0,0,.38); padding: 3px 6px; font-size: 11px; border-radius: 2px; pointer-events: none; }",
		"    .deck-tooltip { font-family: Arial, sans-serif; font-size: 12px; line-height: 1.35; }",
		"  </style>",
		"</head>",
		"<body>",
		"  <div id=\"map\"></div>",
		"  <div id=\"loading\">Fig2d: interactive 3D hexagons on Esri satellite imagery</div>",
		"  <div id=\"legend\"><div id=\"legend-title\">% top-decile phone luck</div><div id=\"legend-bar\"></div><div id=\"legend-ticks\"><span>0%</span><span>10%</span><span>20%</span><span>30%+</span></div></div>",
		"  <div id=\"attrib\">Esri satellite imagery</div>",
		"  <script>",
		paste0("    const DATA = ", dat_json, ";"),
		paste0("    const SATELLITE = ", satellite_json, ";"),
		sprintf("    const CENTER_LON = %.8f;", center_lon),
		sprintf("    const CENTER_LAT = %.8f;", center_lat),
		"    const TILE_URL = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}';",
		"    const colorRange = [[0,147,146], [57,177,133], [156,203,134], [233,226,156], [238,180,121], [232,132,113], [207,89,126]];",
		"    function embeddedSatelliteLayer() {",
		"      if (!SATELLITE.has) return null;",
		"      return new deck.BitmapLayer({",
		"        id: 'embedded-esri-satellite-mosaic',",
		"        image: SATELLITE.src,",
		"        bounds: SATELLITE.bounds",
		"      });",
		"    }",
		"    function satelliteLayer() {",
		"      return new deck.TileLayer({",
		"        id: 'esri-world-imagery',",
		"        data: TILE_URL,",
		"        minZoom: 0,",
		"        maxZoom: 19,",
		"        tileSize: 256,",
		"        refinementStrategy: 'best-available',",
		"        getTileData: ({index}) => {",
		"          const {x, y, z} = index;",
		"          return TILE_URL.replace('{z}', z).replace('{y}', y).replace('{x}', x);",
		"        },",
		"        renderSubLayers: props => {",
		"          if (!props.data) return null;",
		"          const bbox = props.tile.bbox;",
		"          const bounds = [bbox.west, bbox.south, bbox.east, bbox.north];",
		"          return new deck.BitmapLayer(props, {",
		"            id: `${props.id}-bitmap`,",
		"            image: props.data,",
		"            bounds",
		"          });",
		"        }",
		"      });",
		"    }",
		"    function hexRows(object) {",
		"      return (object && object.points ? object.points : []).map(p => p.source || p);",
		"    }",
		"    function mean(rows, key) {",
		"      const vals = rows.map(d => Number(d[key])).filter(Number.isFinite);",
		"      return vals.length ? vals.reduce((a, b) => a + b, 0) / vals.length : null;",
		"    }",
		"    function sum(rows, key) {",
		"      return rows.reduce((a, d) => a + (Number(d[key]) || 0), 0);",
		"    }",
		"    function layers() {",
		"      return [embeddedSatelliteLayer(), satelliteLayer(),",
		"        new deck.HexagonLayer({",
		"          id: 'housing-price-hexagons',",
		"          data: DATA,",
		"          pickable: true,",
		"          autoHighlight: true,",
		"          extruded: true,",
		"          radius: 150,",
		"          coverage: 0.86,",
		"          opacity: 0.82,",
		"          colorRange,",
		"          colorScaleType: 'quantize',",
		"          colorAggregation: 'MEAN',",
		"          elevationAggregation: 'MEAN',",
		"          elevationScale: 16,",
		"          elevationRange: [0, 5200],",
		"          getPosition: d => [d.lon, d.lat],",
		"          getColorWeight: d => Number(d.top10_pct_pct) || 0,",
		"          getElevationWeight: d => Number(d.housing_price_bar) || 0,",
		"          material: {ambient: 0.42, diffuse: 0.62, shininess: 24, specularColor: [255, 255, 255]}",
		"        })",
		"      ];",
		"    }",
		"    const deckgl = new deck.DeckGL({",
		"      container: 'map',",
		"      views: new deck.MapView({repeat: false}),",
		"      initialViewState: {",
		"        longitude: CENTER_LON,",
		"        latitude: CENTER_LAT,",
		"        zoom: 10.2,",
		"        pitch: 55,",
		"        bearing: -28,",
		"        minZoom: 8,",
		"        maxZoom: 15",
		"      },",
		"      controller: true,",
		"      layers: layers(),",
		"      parameters: {depthTest: true},",
		"      getTooltip: ({object}) => {",
		"        if (!object) return null;",
		"        const rows = hexRows(object);",
		"        const avgPrice = mean(rows, 'housing_price');",
		"        const avgTop10 = mean(rows, 'top10_pct_pct');",
		"        const totalN = sum(rows, 'n');",
		"        return {",
		"          html: `<b>Hexagons</b>: ${rows.length}<br/><b>Mean housing price</b>: ${avgPrice == null ? 'NA' : Math.round(avgPrice).toLocaleString()}<br/><b>Mean top-decile phone luck</b>: ${avgTop10 == null ? 'NA' : avgTop10.toFixed(2) + '%'}<br/><b>Total N</b>: ${Math.round(totalN).toLocaleString()}`,",
		"          style: {backgroundColor: 'rgba(0,0,0,.74)', color: 'white'}",
		"        };",
		"      },",
		"      onLoad: () => { const el = document.getElementById('loading'); if (el) el.style.display = 'none'; }",
		"    });",
		"    window.addEventListener('resize', () => deckgl.redraw(true));",
		"  </script>",
		"</body>",
		"</html>"
	)
	dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
	writeLines(html, file, useBytes = TRUE)
	invisible(file)
}
fig2_find_chrome <- function() {
	candidates <- c(
		Sys.getenv("CHROME_BIN", unset = ""),
		Sys.getenv("CHROME_PATH", unset = ""),
		Sys.getenv("EDGE_BIN", unset = ""),
		unname(Sys.which(c(
			"google-chrome", "google-chrome-stable", "chromium-browser", "chromium",
			"microsoft-edge", "microsoft-edge-stable", "chrome", "msedge",
			"chrome.exe", "chromium.exe", "msedge.exe"
		))),
		file.path(Sys.getenv("ProgramFiles", unset = ""), "Google", "Chrome", "Application", "chrome.exe"),
		file.path(Sys.getenv("ProgramFiles(x86)", unset = ""), "Google", "Chrome", "Application", "chrome.exe"),
		file.path(Sys.getenv("ProgramFiles", unset = ""), "Microsoft", "Edge", "Application", "msedge.exe"),
		file.path(Sys.getenv("ProgramFiles(x86)", unset = ""), "Microsoft", "Edge", "Application", "msedge.exe"),
		file.path(Sys.getenv("LocalAppData", unset = ""), "Google", "Chrome", "Application", "chrome.exe"),
		file.path(Sys.getenv("LocalAppData", unset = ""), "Microsoft", "Edge", "Application", "msedge.exe"),
		"/usr/bin/google-chrome",
		"/usr/bin/google-chrome-stable",
		"/usr/bin/chromium-browser",
		"/usr/bin/chromium",
		"/opt/google/chrome/chrome",
		"/opt/microsoft/msedge/msedge",
		"/snap/bin/chromium",
		"/mnt/c/Program Files/Google/Chrome/Application/chrome.exe",
		"/mnt/c/Program Files (x86)/Google/Chrome/Application/chrome.exe",
		"/mnt/c/Program Files/Microsoft/Edge/Application/msedge.exe",
		"/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
	)
	candidates <- unique(candidates[nzchar(candidates)])
	candidates[file.exists(candidates)][1]
}
fig2_wsl_to_windows_path <- function(path) {
	path <- normalizePath(path, winslash = "/", mustWork = FALSE)
	m <- regexec("^/mnt/([A-Za-z])/(.*)$", path)
	hit <- regmatches(path, m)[[1]]
	if (length(hit) == 3L) return(paste0(toupper(hit[2]), ":/", hit[3]))
	path
}
fig2_browser_is_windows <- function(browser) {
	browser <- normalizePath(browser, winslash = "/", mustWork = FALSE)
	grepl("\\.exe$", browser, ignore.case = TRUE) || grepl("^/mnt/[A-Za-z]/", browser)
}
fig2_browser_path <- function(path, browser) {
	path <- normalizePath(path, winslash = "/", mustWork = FALSE)
	if (fig2_browser_is_windows(browser)) path <- fig2_wsl_to_windows_path(path)
	path
}
fig2_file_url <- function(path, browser = NULL) {
	path <- normalizePath(path, winslash = "/", mustWork = TRUE)
	if (!is.null(browser) && fig2_browser_is_windows(browser)) path <- fig2_wsl_to_windows_path(path)
	path <- utils::URLencode(path, reserved = FALSE)
	if (grepl("^[A-Za-z]:/", path)) paste0("file:///", path) else paste0("file://", path)
}
render_fig2_html_png <- function(html_file, png_file, width, height, wait_ms = 5000L) {
	dir.create(dirname(png_file), showWarnings = FALSE, recursive = TRUE)
	if (file.exists(png_file)) unlink(png_file, force = TRUE)
	chrome <- fig2_find_chrome()
	if (is.na(chrome) || !nzchar(chrome)) stop("Chrome/Edge executable not found for rendering Fig2c HTML. Install chromium/google-chrome, or set CHROME_BIN to the browser executable.", call. = FALSE)
	profile_tmp <- if (fig2_browser_is_windows(chrome)) dirname(png_file) else tempdir()
	profile_dir <- tempfile("fig2-chrome-profile-", tmpdir = profile_tmp)
	dir.create(profile_dir, showWarnings = FALSE, recursive = TRUE)
	on.exit(unlink(profile_dir, recursive = TRUE, force = TRUE), add = TRUE)
	png_file2 <- fig2_browser_path(png_file, chrome)
	args <- c(
		"--headless=new",
		"--disable-gpu",
		"--force-device-scale-factor=1",
		"--force-color-profile=srgb",
		"--hide-scrollbars",
		"--run-all-compositor-stages-before-draw",
		"--allow-file-access-from-files",
		"--no-first-run",
		"--no-default-browser-check",
		"--disable-extensions",
		paste0("--user-data-dir=", fig2_browser_path(profile_dir, chrome)),
		paste0("--window-size=", as.integer(width), ",", as.integer(height)),
		paste0("--virtual-time-budget=", as.integer(wait_ms)),
		paste0("--screenshot=", png_file2),
		fig2_file_url(html_file, chrome)
	)
	out <- suppressWarnings(system2(chrome, args, stdout = TRUE, stderr = TRUE))
	if (!file.exists(png_file) || is.na(file.info(png_file)$size) || file.info(png_file)$size < 50000) {
		stop("Chrome did not create a usable Fig2c PNG: ", paste(out, collapse = " | "), call. = FALSE)
	}
	invisible(png_file)
}
render_fig2_satellite_png <- function(dat, png_file, width, height,
	tile_cache_dir, mosaic_stem = "Fig2c", show_legend = TRUE, show_note = TRUE,
	map_alpha = 0.30, bar_height_scale = 1.05, bar_width_scale = 1.70) {
	if (!requireNamespace("jpeg", quietly = TRUE)) pacman::p_load(jpeg)
	z <- 10L
	tile_size <- 256L
	world <- tile_size * 2^z
	mx_r <- function(lon) (lon + 180) / 360 * world
	my_r <- function(lat) {
		s <- sin(lat * pi / 180)
		(0.5 - log((1 + s) / (1 - s)) / (4 * pi)) * world
	}
	color_for <- function(v, alpha = 220) {
		stops_x <- c(0, 5, 10, 15, 20, 25, 30)
		stops_col <- grDevices::col2rgb(c("#009392", "#39B185", "#9CCB86", "#E9E29C", "#EEB479", "#E88471", "#CF597E"))
		v <- pmin(pmax(v, min(stops_x, na.rm = TRUE)), max(stops_x, na.rm = TRUE))
		rgb <- vapply(seq_len(3), function(i) stats::approx(stops_x, stops_col[i, ], xout = v, rule = 2)$y, numeric(length(v)))
		grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], alpha = alpha, maxColorValue = 255)
	}
	xs <- mx_r(dat$lon)
	ys <- my_r(dat$lat)
	tx_range <- seq.int(max(0L, floor(min(xs, na.rm = TRUE) / tile_size) - 1L), min(2^z - 1L, floor(max(xs, na.rm = TRUE) / tile_size) + 1L))
	ty_range <- seq.int(max(0L, floor(min(ys, na.rm = TRUE) / tile_size) - 1L), min(2^z - 1L, floor(max(ys, na.rm = TRUE) / tile_size) + 1L))
	mosaic_file <- file.path(tile_cache_dir, sprintf("%s_satellite_mosaic_z%d.jpg", mosaic_stem, z))
	if (!file.exists(mosaic_file)) stop("Fig2c satellite mosaic is missing: ", mosaic_file, call. = FALSE)
	mosaic <- jpeg::readJPEG(mosaic_file)
	if (length(dim(mosaic)) == 2L) mosaic <- array(rep(mosaic, 3), dim = c(dim(mosaic), 3))
	if (dim(mosaic)[3] > 3L) mosaic <- mosaic[, , 1:3, drop = FALSE]
	mosaic <- mosaic[dim(mosaic)[1]:1, , , drop = FALSE]
	map_alpha <- as.numeric(map_alpha)[1]
	if (!is.finite(map_alpha)) map_alpha <- 0.30
	map_alpha <- clamp(map_alpha, 0, 1)
	bar_height_scale <- as.numeric(bar_height_scale)[1]
	if (!is.finite(bar_height_scale) || bar_height_scale <= 0) bar_height_scale <- 1.05
	bar_width_scale <- as.numeric(bar_width_scale)[1]
	if (!is.finite(bar_width_scale) || bar_width_scale <= 0) bar_width_scale <- 1.70
	mosaic_rgba <- array(1, dim = c(dim(mosaic)[1], dim(mosaic)[2], 4L))
	mosaic_rgba[, , 1:3] <- mosaic
	mosaic_rgba[, , 4] <- map_alpha
	rect_width <- as.numeric(width)
	rect_height <- as.numeric(height)
	margin <- 8
	scale <- min((rect_width - 2 * margin) / (max(xs, na.rm = TRUE) - min(xs, na.rm = TRUE)), (rect_height - 2 * margin) / (max(ys, na.rm = TRUE) - min(ys, na.rm = TRUE))) * 0.96
	ox <- (rect_width - (max(xs, na.rm = TRUE) - min(xs, na.rm = TRUE)) * scale) / 2
	oy <- (rect_height - (max(ys, na.rm = TRUE) - min(ys, na.rm = TRUE)) * scale) / 2
	dat2 <- dat %>%
		mutate(
			wx = mx_r(lon),
			wy = my_r(lat),
			x = ox + (wx - min(xs, na.rm = TRUE)) * scale,
			y = oy + (wy - min(ys, na.rm = TRUE)) * scale,
			h = pmax(18, pmin(320, .data$height * bar_height_scale)),
			sx = x + h * 0.10,
			sy = y - h
		) %>%
		arrange(y)
	sat_x <- ox + (min(tx_range) * tile_size - min(xs, na.rm = TRUE)) * scale
	sat_y <- oy + (min(ty_range) * tile_size - min(ys, na.rm = TRUE)) * scale
	sat_w <- length(tx_range) * tile_size * scale
	sat_h <- length(ty_range) * tile_size * scale
	dir.create(dirname(png_file), showWarnings = FALSE, recursive = TRUE)
	grDevices::png(png_file, width = as.integer(width), height = as.integer(height), units = "px", bg = "white", type = if (capabilities("cairo")) "cairo" else getOption("bitmapType", "cairo"))
	on.exit(invisible(grDevices::dev.off()), add = TRUE)
	grid::grid.newpage()
	grid::pushViewport(grid::viewport(xscale = c(0, rect_width), yscale = c(rect_height, 0)))
	grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
	grid::grid.raster(mosaic_rgba, x = grid::unit(sat_x + sat_w / 2, "native"), y = grid::unit(sat_y + sat_h / 2, "native"), width = grid::unit(sat_w, "native"), height = grid::unit(sat_h, "native"), interpolate = TRUE)
	grid::grid.segments(grid::unit(dat2$x + 3, "native"), grid::unit(dat2$y + 3, "native"), grid::unit(dat2$sx + 3, "native"), grid::unit(dat2$sy + 3, "native"), gp = grid::gpar(col = grDevices::rgb(20, 25, 35, alpha = 118, maxColorValue = 255), lwd = 4.3 * bar_width_scale, lineend = "round"))
	grid::grid.segments(grid::unit(dat2$x, "native"), grid::unit(dat2$y, "native"), grid::unit(dat2$sx, "native"), grid::unit(dat2$sy, "native"), gp = grid::gpar(col = color_for(dat2$top10_pct_pct, 250), lwd = 2.6 * bar_width_scale, lineend = "round"))
	grid::grid.points(grid::unit(dat2$sx, "native"), grid::unit(dat2$sy, "native"), pch = 16, size = grid::unit(2.7 * bar_width_scale, "native"), gp = grid::gpar(col = color_for(dat2$top10_pct_pct, 255)))
	if (isTRUE(show_legend)) {
		legend_x <- rect_width - 18 - 430
		legend_y <- 18
		grid::grid.rect(x = grid::unit(legend_x, "native"), y = grid::unit(legend_y, "native"), width = grid::unit(430, "native"), height = grid::unit(125, "native"), just = c("left", "top"), gp = grid::gpar(fill = grDevices::rgb(0, 0, 0, alpha = 128, maxColorValue = 255), col = NA))
		grid::grid.text("% top-decile phone luck", x = grid::unit(legend_x + 14, "native"), y = grid::unit(legend_y + 24, "native"), just = c("left", "center"), gp = grid::gpar(col = "white", fontsize = 24, fontface = "bold"))
		legend_cols <- matrix(grDevices::colorRampPalette(c("#009392", "#39B185", "#9CCB86", "#E9E29C", "#EEB479", "#E88471", "#CF597E"))(256), nrow = 1)
		grid::grid.raster(legend_cols, x = grid::unit(legend_x + 14 + 402 / 2, "native"), y = grid::unit(legend_y + 55, "native"), width = grid::unit(402, "native"), height = grid::unit(18, "native"), interpolate = TRUE)
		grid::grid.rect(x = grid::unit(legend_x + 14, "native"), y = grid::unit(legend_y + 46, "native"), width = grid::unit(402, "native"), height = grid::unit(18, "native"), just = c("left", "top"), gp = grid::gpar(fill = NA, col = grDevices::rgb(255, 255, 255, alpha = 166, maxColorValue = 255), lwd = 2))
		tick_x <- legend_x + 14 + c(0, 134, 268, 402)
		grid::grid.text(c("0%", "10%", "20%", "30%"), x = grid::unit(tick_x, "native"), y = grid::unit(legend_y + 86, "native"), just = c("center", "center"), gp = grid::gpar(col = "white", fontsize = 20))
		grid::grid.text("Height: housing price / 1000", x = grid::unit(legend_x + 14, "native"), y = grid::unit(legend_y + 112, "native"), just = c("left", "center"), gp = grid::gpar(col = "white", fontsize = 24))
	}
	if (isTRUE(show_note)) {
		grid::grid.rect(x = grid::unit(18, "native"), y = grid::unit(rect_height - 16 - 34, "native"), width = grid::unit(235, "native"), height = grid::unit(34, "native"), just = c("left", "top"), gp = grid::gpar(fill = grDevices::rgb(0, 0, 0, alpha = 128, maxColorValue = 255), col = NA))
		grid::grid.text("Esri satellite imagery", x = grid::unit(31, "native"), y = grid::unit(rect_height - 33, "native"), just = c("left", "center"), gp = grid::gpar(col = "white", fontsize = 20))
	}
	grid::popViewport()
	invisible(png_file)
}
make_fig2c_panel <- function(img = NULL) {
	legend_cols <- matrix(grDevices::colorRampPalette(c("#009392", "#39B185", "#9CCB86", "#E9E29C", "#EEB479", "#E88471", "#CF597E"))(256), nrow = 1)
	p <- ggplot()
	if (!is.null(img)) {
		p <- p + annotation_custom(grid::rasterGrob(img, interpolate = TRUE), xmin = 0, xmax = 1, ymin = 0, ymax = 1.000)
	} else {
		p <- p + annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1.000, fill = "white", color = NA)
	}
	p +
		annotate("text", x = .40, y = 1.043, label = "% top-decile phone luck", fontface = "bold", size = 2.85, color = "grey15", hjust = 1) +
		annotation_custom(grid::rasterGrob(legend_cols, interpolate = TRUE), xmin = .415, xmax = .735, ymin = 1.033, ymax = 1.052) +
		annotate("rect", xmin = .415, xmax = .735, ymin = 1.033, ymax = 1.052, fill = NA, color = "grey35", linewidth = .22) +
		annotate("text", x = seq(.415, .735, length.out = 6), y = 1.012, label = paste0(seq(0, 25, 5), "%"), size = 2.4, color = "grey15") +
		coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.060), expand = FALSE, clip = "off") +
		labs(title = "c. Shenzhen map: housing price and top-decile phone luck", x = NULL, y = NULL) +
		theme_void(base_size = 8.2) +
		theme(
			plot.title = element_text(size = 9.8, face = "bold", hjust = 0.5, margin = margin(t = 4, b = 0)),
			plot.background = element_rect(fill = "white", color = NA),
			panel.border = element_blank(),
			plot.margin = margin(0, 0, 0, 0)
		)
}
make_fig2c_image_panel <- function(png_file) {
	if (!requireNamespace("png", quietly = TRUE)) pacman::p_load(png)
	img <- png::readPNG(png_file, native = FALSE)
	make_fig2c_panel(img)
}
make_fig2c_skip_panel <- function() make_fig2c_panel(NULL)
if (isTRUE(draw_2D_map) && nrow(fig2c_dat) >= 20) {
	fig2c_tile_cache_dir <- file.path(fig2_map_tile_root, "2D", "Fig2c_satellite_tiles")
	fig2c_result <- tryCatch({
		write_fig2_satellite_canvas(
			dat = fig2c_dat,
			file = fig2c_html_file,
			width = fig2c_html_width,
			height = fig2c_html_height,
			tile_cache_dir = fig2c_tile_cache_dir,
			mosaic_stem = "Fig2c",
			html_title = "Fig2c satellite 3D housing price",
			note_text = "Esri satellite imagery",
			show_legend = FALSE,
			show_note = FALSE,
			map_alpha = fig2c_map_alpha,
			bar_height_scale = fig2c_bar_height_scale,
			bar_width_scale = fig2c_bar_width_scale
		)
		fig2c_shot <- tryCatch({
			render_fig2_html_png(
				html_file = fig2c_html_file,
				png_file = fig2c_png_file,
				width = fig2c_html_width,
				height = fig2c_html_height
			)
			list(ok = TRUE, error = "")
		}, error = function(e) {
			if (file.exists(fig2c_png_file)) unlink(fig2c_png_file, force = TRUE)
			list(ok = FALSE, error = conditionMessage(e))
		})
		if (isTRUE(fig2c_shot$ok)) {
			list(plot = make_fig2c_image_panel(fig2c_png_file), status = "ok: Fig2c.html screenshot rendered to Fig2c.png by Chrome", mode = "html_chrome_png")
		} else {
			list(plot = make_fig2c_skip_panel(), status = paste0("skipped: Fig2c.html written; Chrome screenshot skipped; title and legend kept for manual splice. Error: ", fig2c_shot$error), mode = "html_chrome_skipped")
		}
	}, error = function(e) {
		stop("Fig2c 2D map rendering failed. Check internet access to Esri tiles and the cache folder ",
			fig2c_tile_cache_dir,
			". Error: ", conditionMessage(e), call. = FALSE)
	})
	p2C <- fig2c_result$plot
	fig2c_status <- fig2c_result$status
	fig2_tile_status <- fig2c_result$status
	fig2c_render_mode <- fig2c_result$mode
} else if (isTRUE(draw_2D_map)) {
	stop("Fig2c 2D map requested but insufficient coordinates after filtering. Check dat.list.rds: 接车地址经度, 接车地址纬度, 房价指数, and phone.sco.", call. = FALSE)
}
if (isTRUE(draw_3D_map)) {
	fig2_map_tile_3d_dir <- file.path(fig2_map_tile_root, "3D")
	dir.create(fig2_map_tile_3d_dir, showWarnings = FALSE, recursive = TRUE)
	fig2d_data_file <- file.path(fig2_map_tile_3d_dir, "Fig2d_house_data.tsv")
	fig2d_py_file <- file.path(fig2_map_tile_3d_dir, "Fig2d.keplergl.py")
	fig2d_config_file <- file.path(fig2_map_tile_3d_dir, "Fig2d.keplergl.config.json")
	fig2d_html_file <- file.path(dir.analysis, "Fig2d.html")
	fig2d_map_raw <- fig2_base %>%
		filter(is.finite(.data[[fig2_coord_vars[1]]]), is.finite(.data[[fig2_coord_vars[2]]]), is.finite(.data[[fig2_house_plot_var]]), is.finite(phone.sco))
	if (nrow(fig2d_map_raw) >= 100 && is.finite(fig2_top10_cut)) {
		fig2d_map_dat <- fig2d_map_raw %>%
			mutate(phone_top10 = phone.sco >= fig2_top10_cut, lon = .data[[fig2_coord_vars[1]]], lat = .data[[fig2_coord_vars[2]]], lon_bin = round(lon, 3), lat_bin = round(lat, 3)) %>%
			filter(between(lon, 113.6, 114.8), between(lat, 22.35, 22.95)) %>%
			group_by(lon = lon_bin, lat = lat_bin) %>%
			summarise(n = n(), housing_price = median(.data[[fig2_house_plot_var]], na.rm = TRUE), top10_pct = mean(phone_top10, na.rm = TRUE), .groups = "drop") %>%
			mutate(
				housing_price_sqrt = if_else(housing_price > 0, sqrt(housing_price), NA_real_),
				housing_price_log = if_else(housing_price > 0, log(housing_price), NA_real_),
				housing_price_log10 = if_else(housing_price > 0, log10(housing_price), NA_real_),
				housing_price_bar = if (identical(house_price_bar, "log")) {
					housing_price_log
				} else if (identical(house_price_bar, "log10")) {
					housing_price_log10
				} else {
					housing_price_sqrt
				}
			) %>%
			filter(n >= 10, is.finite(housing_price), is.finite(housing_price_bar), is.finite(top10_pct)) %>%
			slice_max(n, n = 1500, with_ties = FALSE)
	}
}
if (isTRUE(draw_3D_map) && nrow(fig2d_map_dat) < 20) {
	stop("geo.map.use requires Fig2d rendering, but Fig2d has insufficient map data after filtering. Check coordinates, 房价指数, and phone.sco in dat.list.rds.", call. = FALSE)
}
if (isTRUE(draw_3D_map) && nrow(fig2d_map_dat) >= 20) {
	dir.create(dirname(fig2d_data_file), showWarnings = FALSE, recursive = TRUE)
	fig2d_dat <- fig2d_map_dat %>%
		transmute(
			lon,
			lat,
			n,
			housing_price = round(housing_price, 2),
			housing_price_sqrt = round(housing_price_sqrt, 4),
			housing_price_log = round(housing_price_log, 4),
			housing_price_log10 = round(housing_price_log10, 4),
			housing_price_bar = round(housing_price_bar, 4),
			height_metric = fig2_house_bar_label,
			top10_pct = round(top10_pct, 5),
			top10_pct_pct = round(100 * top10_pct, 2),
			height = round(pmin(housing_price, fig2_house_price_cap["high"]) / 1000, 2)
		)
	data.table::fwrite(fig2d_dat, fig2d_data_file, sep = "\t", na = "")
	fig2d_status <- tryCatch({
		if (identical(geo.map.use, "deckgl_satellite")) {
			fig2d_py_file <<- NA_character_
			fig2d_config_file <<- NA_character_
			fig2d_python <<- NA_character_
			fig2d_tile_cache_dir <- file.path(fig2_map_tile_root, "3D", "Fig2d_satellite_tiles")
			fig2d_mosaic_seed_file <- file.path(fig2_map_tile_3d_dir, "Fig2d.satellite_mosaic_seed.html")
			write_fig2_satellite_canvas(
				dat = fig2d_dat,
				file = fig2d_mosaic_seed_file,
				width = fig2d_html_width,
				height = fig2d_html_height,
				tile_cache_dir = fig2d_tile_cache_dir,
				mosaic_stem = "Fig2d",
				html_title = "Fig2d satellite mosaic seed",
				note_text = "",
				show_legend = FALSE,
				show_note = FALSE,
				tile_pad = 8L
			)
			write_fig2d_deckgl_html(
				fig2d_dat,
				fig2d_html_file,
				width = fig2d_html_width,
				height = fig2d_html_height,
				tile_cache_dir = fig2d_tile_cache_dir,
				mosaic_stem = "Fig2d",
				tile_pad = 8L
			)
			fig2d_render_mode <<- "deckgl_satellite"
		} else if (identical(geo.map.use, "deckgl_hex_satellite")) {
			fig2d_py_file <<- NA_character_
			fig2d_config_file <<- NA_character_
			fig2d_python <<- NA_character_
			fig2d_tile_cache_dir <- file.path(fig2_map_tile_root, "3D", "Fig2d_satellite_tiles")
			fig2d_mosaic_seed_file <- file.path(fig2_map_tile_3d_dir, "Fig2d.satellite_mosaic_seed.html")
			write_fig2_satellite_canvas(
				dat = fig2d_dat,
				file = fig2d_mosaic_seed_file,
				width = fig2d_html_width,
				height = fig2d_html_height,
				tile_cache_dir = fig2d_tile_cache_dir,
				mosaic_stem = "Fig2d",
				html_title = "Fig2d satellite mosaic seed",
				note_text = "",
				show_legend = FALSE,
				show_note = FALSE,
				tile_pad = 8L
			)
			write_fig2d_deckgl_hex_html(
				fig2d_dat,
				fig2d_html_file,
				width = fig2d_html_width,
				height = fig2d_html_height,
				tile_cache_dir = fig2d_tile_cache_dir,
				mosaic_stem = "Fig2d",
				tile_pad = 8L
			)
			fig2d_render_mode <<- "deckgl_hex_satellite"
		} else if (identical(geo.map.use, "keplergl_3d")) {
			fig2d_python <<- find_fig2_kepler_python()
			write_fig2d_keplergl_files(
				dat = fig2d_dat,
				data_file = fig2d_data_file,
				py_file = fig2d_py_file,
				config_file = fig2d_config_file,
				html_file = fig2d_html_file,
				width = fig2d_html_width,
				height = fig2d_html_height
			)
			out <- suppressWarnings(system2(fig2d_python, shQuote(fig2d_py_file), stdout = TRUE, stderr = TRUE))
			status <- attr(out, "status")
			if (is.null(status)) status <- 0L
			if (!identical(status, 0L)) stop(paste(out, collapse = " | "), call. = FALSE)
			fig2d_render_mode <<- "keplergl_3d"
		} else {
			stop("Unsupported geo.map.use: ", geo.map.use, call. = FALSE)
		}
		if (!file.exists(fig2d_html_file) || is.na(file.info(fig2d_html_file)$size) || file.info(fig2d_html_file)$size < 1024) {
			stop("Fig2d.html was not created or is too small after HTML rendering.", call. = FALSE)
		}
		sprintf("ok: Fig2d.html rendered with geo.map.use=%s", geo.map.use)
	}, error = function(e) {
		fig2d_render_mode <<- paste0(geo.map.use, "_failed")
		stop(paste0("Fig2d 3D HTML generation failed for geo.map.use=", geo.map.use, ". Error: ", conditionMessage(e)), call. = FALSE)
	})
}

fig2_top_core <- (p2A | plot_spacer() | p2B) + plot_layout(widths = c(1.04, .065, .6))
fig2_outer_widths <- c(.04, 1, .04)
fig2_top <- (plot_spacer() | fig2_top_core | plot_spacer()) + plot_layout(widths = fig2_outer_widths)
if (isTRUE(draw_2D_map) && !is.null(p2C)) {
	fig2_bottom <- (plot_spacer() | p2C | plot_spacer()) + plot_layout(widths = fig2_outer_widths)
	Fig2 <- cowplot::plot_grid(fig2_top, NULL, fig2_bottom, ncol = 1, rel_heights = c(1.30, .10, 2.22), align = "none")
	fig2_height <- 12.4
} else {
	Fig2 <- fig2_top
	fig2_height <- 4.9
}
save_plot(Fig2, "Fig3.png", width = 13.2, height = fig2_height, dpi = 600)
fig2_density_n <- if (fig2_density_var %in% names(fig2_base)) sum(is.finite(fig2_base[[fig2_density_var]])) else 0L

writexl::write_xlsx(list(
	merge_log = geo_xia_merge_log,
	fig2_summary = fig2_base %>% summarise(n = n(), years = n_distinct(year), geo_type_use = geo.type.use, geo_type_column = fig2_geo_type_var, house_price_use = fig2_house_price_filter_label, house_price_min = house_price_min, n_before_house_price_min = fig2_n_before_house_price_min, n_after_house_price_min = fig2_n_after_house_price_min, housing_capped_n = sum(is.finite(.data[[fig2_house_plot_var]])), building_density_n = fig2_density_n, phone_score_n = sum(is.finite(phone.sco)), top10_cutoff = fig2_top10_cut, housing_cap_low = fig2_house_price_cap["low"], housing_cap_high = fig2_house_price_cap["high"], house_price_scale = house_price_scale, house_price_bar = house_price_bar, map_tile_status = fig2_tile_status),
	panelA_top_luck_cutoffs = fig2_top_luck_cuts,
	panelA_house_binned = fig2A_bin,
	panelA_endpoint_values = if ("is_endpoint" %in% names(fig2A_bin)) fig2A_bin %>% filter(.data$is_endpoint) else tibble(),
	panelA_house_trend = fig2A_trend,
	panelB_dx_house_quintile = fig2B_heat,
	panelC_map_grid = fig2_map_dat,
	panelC_render_config = tibble(html_file = fig2c_html_file, png_file = fig2c_png_file, html_width = fig2c_html_width, html_height = fig2c_html_height, status = fig2c_status, render_mode = fig2c_render_mode, height_metric = fig2_house_bar_label, map_transparent = geo.map.transparent, map_alpha = fig2c_map_alpha, bar_height_scale = fig2c_bar_height_scale, bar_width_scale = fig2c_bar_width_scale),
	panelD_kepler_data = fig2d_dat,
	panelD_kepler_config = tibble(data_file = fig2d_data_file, py_file = fig2d_py_file, config_file = fig2d_config_file, html_file = fig2d_html_file, html_width = fig2d_html_width, html_height = fig2d_html_height, python = fig2d_python, geo_map_use = geo.map.use, status = fig2d_status, render_mode = fig2d_render_mode, height_metric = fig2_house_bar_label),
	phone_environment_smd = fig2B_dat
), "Fig3.out.xlsx")

.ems120_remove_temp_rds("Fig2.rds")
ems120_maybe_exit_after("fig3")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig4. Phone-score gradient in EMS phenotype composition across 12 years
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig3_dxs_plot <- dxs.all
fig3_dxs_out <- dxs.all
if (!ems120_should_run("fig4")) {
	cap2 <- 1.10
	fig3_dat <- bind_rows(lapply(years, function(y) {
		d0 <- dat1.list[[as.character(y)]]
		if (is.null(d0) || !nrow(d0)) return(tibble())
		d0 <- d0 %>%
			filter(phone.luck %in% phone.grp.use, !is.na(dx_grp), dx_grp %in% fig3_dxs_out) %>%
			transmute(group = factor(as.character(phone.luck), levels = phone.grp.use), disease = factor(as.character(dx_grp), levels = fig3_dxs_out))
		tg <- d0 %>% count(disease, group, name = "n") %>% complete(disease = fig3_dxs_out, group = phone.grp.use, fill = list(n = 0))
		ta <- d0 %>% count(disease, name = "n_all") %>% complete(disease = fig3_dxs_out, fill = list(n_all = 0))
		Ns <- d0 %>% count(group, name = "N") %>% complete(group = phone.grp.use, fill = list(N = 0))
		tg %>%
			left_join(Ns, by = "group") %>%
			mutate(pct_group = n / N) %>%
			dplyr::select(disease, group, n, N, pct_group) %>%
			pivot_wider(names_from = group, values_from = c(n, N, pct_group), names_sep = "_") %>%
			left_join(ta, by = "disease") %>%
			mutate(
				year = y,
				N_all = sum(ta$n_all),
				pct_all = n_all / N_all,
				enrich_low = pct_group_low / pct_all,
				enrich_high = pct_group_high / pct_all,
				RR = pct_group_high / pct_group_low,
				se_logRR = sqrt(1 / pmax(n_high, 1) - 1 / pmax(N_high, 1) + 1 / pmax(n_low, 1) - 1 / pmax(N_low, 1)),
				RR_lo = exp(log(RR) - 1.96 * se_logRR),
				RR_hi = exp(log(RR) + 1.96 * se_logRR),
				p = purrr::pmap_dbl(
					list(n_high, N_high, n_low, N_low),
					function(a, A, b, B) suppressWarnings(chisq.test(matrix(c(a, A - a, b, B - b), nrow = 2))$p.value)
				),
				lo = clamp(enrich_low, 1 / cap2, cap2),
				hi = clamp(enrich_high, 1 / cap2, cap2),
				lf = enrich_low < 1 / cap2,
				hf = enrich_high > cap2
			)
	})) %>%
		group_by(disease) %>%
		mutate(
			n_total = n_low + n_high,
			N_total = N_low + N_high,
			p_adj = p.adjust(p, method = "BH"),
			sig = !is.na(p_adj) & p_adj < .01
		) %>%
		ungroup()
}
if (ems120_should_run("fig4")) {
cap2 <- 1.10
fig3_count_x_high <- 1.136
fig3_count_x_low <- 1.160

fig3_dat <- bind_rows(lapply(years, function(y) {
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())
	d0 <- d0 %>%
		filter(phone.luck %in% phone.grp.use, !is.na(dx_grp), dx_grp %in% fig3_dxs_out) %>%
		transmute(group = factor(as.character(phone.luck), levels = phone.grp.use), disease = factor(as.character(dx_grp), levels = fig3_dxs_out))
	tg <- d0 %>% count(disease, group, name = "n") %>% complete(disease = fig3_dxs_out, group = phone.grp.use, fill = list(n = 0))
	ta <- d0 %>% count(disease, name = "n_all") %>% complete(disease = fig3_dxs_out, fill = list(n_all = 0))
	Ns <- d0 %>% count(group, name = "N") %>% complete(group = phone.grp.use, fill = list(N = 0))
	out <- tg %>%
		left_join(Ns, by = "group") %>%
		mutate(pct_group = n / N) %>%
		dplyr::select(disease, group, n, N, pct_group) %>%
		pivot_wider(names_from = group, values_from = c(n, N, pct_group), names_sep = "_") %>%
		left_join(ta, by = "disease") %>%
		mutate(
			year = y,
			N_all = sum(ta$n_all),
			pct_all = n_all / N_all,
			enrich_low = pct_group_low / pct_all,
			enrich_high = pct_group_high / pct_all,
			RR = pct_group_high / pct_group_low,
			se_logRR = sqrt(1 / pmax(n_high, 1) - 1 / pmax(N_high, 1) + 1 / pmax(n_low, 1) - 1 / pmax(N_low, 1)),
			RR_lo = exp(log(RR) - 1.96 * se_logRR),
			RR_hi = exp(log(RR) + 1.96 * se_logRR),
			p = purrr::pmap_dbl(
				list(n_high, N_high, n_low, N_low),
				function(a, A, b, B) suppressWarnings(chisq.test(matrix(c(a, A - a, b, B - b), nrow = 2))$p.value)
			),
			lo = clamp(enrich_low, 1 / cap2, cap2),
			hi = clamp(enrich_high, 1 / cap2, cap2),
			lf = enrich_low < 1 / cap2,
			hf = enrich_high > cap2
		)
	out
})) %>%
	group_by(disease) %>%
	mutate(
		n_total = n_low + n_high,
		N_total = N_low + N_high,
		p_adj = p.adjust(p, method = "BH"),
		sig = !is.na(p_adj) & p_adj < .01
	) %>%
	ungroup()

fit_fig3_trend <- function(d) {
	d <- d %>% filter(is.finite(year), n_low >= 0, n_high >= 0, N_low > 0, N_high > 0) %>% arrange(year)
	if (nrow(d) < 4) return(tibble(beta_year = NA_real_, RR_change_per_year = NA_real_, p_trend = NA_real_, p_trend_high_low_interaction = NA_real_, n_years = nrow(d)))
	d <- d %>% mutate(year_c = year - min(year, na.rm = TRUE), n_total = n_low + n_high, N_total = N_low + N_high)
	fit_overall <- tryCatch(glm(cbind(n_total, N_total - n_total) ~ year_c, family = binomial(), data = d), error = function(e) NULL)
	beta <- p_overall <- NA_real_
	if (!is.null(fit_overall)) {
		td <- broom::tidy(fit_overall) %>% filter(term == "year_c")
		if (nrow(td)) { beta <- td$estimate[1]; p_overall <- td$p.value[1] }
	}
	bin <- bind_rows(
		d %>% transmute(year_c, group = factor("low", levels = c("low", "high")), n = n_low, N = N_low),
		d %>% transmute(year_c, group = factor("high", levels = c("low", "high")), n = n_high, N = N_high)
	)
	fit_int <- tryCatch(glm(cbind(n, N - n) ~ group * year_c, family = binomial(), data = bin), error = function(e) NULL)
	p_int <- NA_real_
	if (!is.null(fit_int)) {
		td_int <- broom::tidy(fit_int) %>% filter(term == "grouphigh:year_c")
		if (nrow(td_int)) p_int <- td_int$p.value[1]
	}
	tibble(beta_year = beta, RR_change_per_year = exp(beta), p_trend = p_overall, p_trend_high_low_interaction = p_int, n_years = nrow(d))
}
fit_fig3_pair <- function(d) {
	d <- d %>% filter(is.finite(year), is.finite(enrich_low), is.finite(enrich_high), n_low >= 0, n_high >= 0, N_low > 0, N_high > 0) %>% arrange(year)
	if (nrow(d) < 2) return(tibble(mean_enrich_low = mean(d$enrich_low, na.rm = TRUE), mean_enrich_high = mean(d$enrich_high, na.rm = TRUE), mean_pair_diff = NA_real_, p_pair = NA_real_, p_pair_paired_t = NA_real_, n_years_pair = nrow(d)))
	tt <- tryCatch(stats::t.test(d$enrich_high, d$enrich_low, paired = TRUE), error = function(e) NULL)
	bin <- bind_rows(
		d %>% transmute(year = factor(year), group = factor("low", levels = c("low", "high")), n = n_low, N = N_low),
		d %>% transmute(year = factor(year), group = factor("high", levels = c("low", "high")), n = n_high, N = N_high)
	)
	fit_pair <- tryCatch(glm(cbind(n, N - n) ~ group + year, family = binomial(), data = bin), error = function(e) NULL)
	p_pair <- NA_real_
	if (!is.null(fit_pair)) {
		td_pair <- broom::tidy(fit_pair) %>% filter(term == "grouphigh")
		if (nrow(td_pair)) p_pair <- td_pair$p.value[1]
	}
	tibble(
		mean_enrich_low = mean(d$enrich_low, na.rm = TRUE),
		mean_enrich_high = mean(d$enrich_high, na.rm = TRUE),
		mean_pair_diff = mean(d$enrich_high - d$enrich_low, na.rm = TRUE),
		p_pair = p_pair,
		p_pair_paired_t = ifelse(is.null(tt), NA_real_, tt$p.value),
		n_years_pair = nrow(d)
	)
}
fig3_fmt_p_compact <- function(p) {
	p <- suppressWarnings(as.numeric(p))[1]
	if (is.finite(p) && p < 1e-320) return("<1x10^-320")
	fmt_p_compact(p)
}
fig3_p_math <- function(p, suffix, digits = 2) {
	p <- suppressWarnings(as.numeric(p))[1]
	lhs <- sprintf("italic(P)[%s]", suffix)
	if (!is.finite(p)) return(paste0(lhs, "*' = NA'"))
	if (p < 1e-320) return(sprintf("%s*' < '*1~'×'~10^{-320}", lhs))
	val <- fmt_p_x_math(p, digits = digits)[1]
	val <- gsub("'x'", "'×'", val, fixed = TRUE)
	paste0(lhs, "*' = '*", val)
}
plotmath_quote <- function(x) gsub("'", "\\'", as.character(x), fixed = TRUE)
fig3_trend <- fig3_dat %>%
	group_by(disease) %>% group_modify(~ fit_fig3_trend(.x)) %>% ungroup() %>%
	left_join(fig3_dat %>% group_by(disease) %>% group_modify(~ fit_fig3_pair(.x)) %>% ungroup(), by = "disease") %>%
	mutate(
		p_trend_label = vapply(p_trend, fig3_fmt_p_compact, character(1)),
		p_pair_label = vapply(p_pair, fig3_fmt_p_compact, character(1)),
		p_trend_math = vapply(p_trend, function(x) fig3_p_math(x, "t"), character(1)),
		p_pair_math = vapply(p_pair, function(x) fig3_p_math(x, "p"), character(1)),
		trend_label = paste0(as.character(disease), " (Pt=", p_trend_label, "; Pp=", p_pair_label, ")"),
		trend_label_math = paste0("bold('", plotmath_quote(as.character(disease)), "')*' ('*", p_trend_math, "*'; '*", p_pair_math, "*')'")
	)
fig3_title_label <- setNames(fig3_trend$trend_label_math, as.character(fig3_trend$disease))

Fig3 <- wrap_plots(lapply(fig3_dxs_plot, function(dx) {
	col <- dxs.all.color[dx]
	d <- fig3_dat %>% filter(disease == dx)
	title_label <- if (dx %in% names(fig3_title_label)) parse(text = fig3_title_label[[dx]])[[1]] else dx_to_eng(dx)
	ggplot(d, aes(y = year)) +
		geom_vline(xintercept = 1, linewidth = .45) +
		geom_segment(aes(x = 1, xend = lo, yend = year), linetype = "dashed", color = "grey70", linewidth = .85) +
		geom_segment(aes(x = 1, xend = hi, yend = year), linetype = "dashed", color = col, linewidth = .85) +
		geom_point(aes(x = lo), color = "grey50", size = 2.4) +
		geom_point(aes(x = hi), color = col, size = 2.4) +
		geom_text(data = d %>% filter(sig), aes(x = hi, label = "*"), hjust = -.2, vjust = .3, size = 3.8, fontface = "bold") +
		geom_text(data = d %>% filter(lf), aes(x = lo, label = "<"), hjust = 1.15, size = 3.2) +
		geom_text(data = d %>% filter(hf), aes(x = hi, label = ">"), hjust = -.15, size = 3.2) +
		scale_x_continuous(limits = c(.90, 1.11), breaks = c(.9, 1.0, 1.1), expand = expansion(mult = c(0, 0))) +
		scale_y_continuous(breaks = years, labels = years) +
		labs(title = title_label, x = NULL, y = NULL) +
		theme_minimal(base_size = 11) +
		theme(axis.text = element_text(face = "bold", size = 7.5), axis.text.y = element_text(margin = margin(r = 1)), axis.title = element_text(size = 8.5, face = "bold"), axis.line = element_line(linewidth = .35), plot.title = element_text(face = "plain", size = 9.2, hjust = 0.5, margin = margin(b = 2)), panel.grid.minor = element_blank(), plot.margin = margin(5, 2, 5, 1))
}), ncol = 3) + plot_annotation(
	title = "Annual phenotype enrichment by phone-score group",
	theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = .5, margin = margin(b = 7)))
)
save_plot(Fig3, "Fig4.png", width = 10.2, height = 10.4, dpi = 600)

writexl::write_xlsx(list(
	panel_data_all_dxs = fig3_dat %>% dplyr::select(year, disease, n_low, n_high, n_total, N_low, N_high, N_total, pct_all, pct_group_low, pct_group_high, enrich_low, enrich_high, RR, RR_lo, RR_hi, p, p_adj, sig, lo, hi),
	panel_data_plotted = fig3_dat %>% filter(disease %in% fig3_dxs_plot),
	trend_tests = fig3_trend,
	summary_all_dxs = fig3_dat %>% group_by(disease) %>%
		summarise(
			mean_RR = round(mean(RR, na.rm = TRUE), 3),
			min_RR = round(min(RR, na.rm = TRUE), 3),
			max_RR = round(max(RR, na.rm = TRUE), 3),
			n_sig_years = sum(sig, na.rm = TRUE),
			.groups = "drop"
		)
), "Fig4.out.xlsx")
ems120_maybe_exit_after("fig4")
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig5. COVID-19 policy-transition analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Design note: keep the original 2 x 2 main-figure layout, but make the
# statistical target explicit. Panels A/C test the overall high-vs-low
# change in total EMS demand. Panels B/D test whether the disease mix changed
# differently between phone-score groups, using a global disease-mix
# interaction test plus disease-specific estimates. Low-count disease cells
# are flagged rather than suppressed.
if (ems120_should_run("fig5")) {
fig5_min_cell_n <- 20
fig5_note_size <- 2.35
fig5_note_x_frac <- .004
fig5_note_plain <- function(rr, policy_p, did_p) {
	sprintf("RR=%.2f\nPolicy P=%s\nH/L DiD P=%s", rr, fmt_p_x_math(policy_p), fmt_p_x_math(did_p))
}
fig5_note_plot_dat <- function(x, rr, policy_p, did_p, y_source) {
	y_rng <- range(suppressWarnings(as.numeric(y_source)), na.rm = TRUE)
	y_span <- diff(y_rng)
	if (!all(is.finite(y_rng)) || !is.finite(y_span) || y_span <= 0) {
		y_rng <- c(0, 1)
		y_span <- 1
	}
	tibble(
		x = x,
		y = y_rng[2] + y_span * .175 - (0:2) * y_span * .082,
		label = c(
			sprintf("bold('RR = ')*bold(%.2f)", rr),
			sprintf("bold('Policy P = ')*bold(%s)", fmt_p_x_math(policy_p)),
			sprintf("bold('H/L DiD P = ')*bold(%s)", fmt_p_x_math(did_p))
		)
	)
}
fig5_format_effect <- function(x, sig = "", low_count = FALSE) {
	sig <- replace_na(as.character(sig), "")
	low_count <- replace_na(as.logical(low_count), FALSE)
	ifelse(is.finite(x), paste0(sprintf("%.2f", x), ifelse(nzchar(sig), paste0(" ", sig), ""), ifelse(low_count, "\u2020", "")), "")
}

fit_total_policy_period <- function(dat, period_var = "period", ref_period, contrast_periods) {
	dat <- dat %>% filter(phone.luck %in% c("low", "high"), total_group_day >= 0) %>%
		mutate(high = as.integer(phone.luck == "high"), periodF = stats::relevel(factor(.data[[period_var]]), ref = ref_period))
	if (!"dow" %in% names(dat)) dat <- dat %>% mutate(dow = factor(lubridate::wday(date, label = TRUE, week_start = 1)))
	empty <- tibble(period = contrast_periods, RR = NA_real_, lo = NA_real_, hi = NA_real_, p_poisson = NA_real_, overdispersion = NA_real_, p = NA_real_, p_adj = NA_real_, sig05 = "", p_label = "NA")
	if (nrow(dat) < 10 || length(unique(dat$periodF)) < 2) return(empty)
	fit <- tryCatch(glm(total_group_day ~ periodF + high + dow, family = poisson(), data = dat), error = function(e) NULL)
	if (is.null(fit)) return(empty)
	td <- broom::tidy(fit)
	disp <- suppressWarnings(sum(stats::residuals(fit, type = "pearson")^2, na.rm = TRUE) / stats::df.residual(fit))
	if (!is.finite(disp) || disp < 1) disp <- 1
	purrr::map_dfr(contrast_periods, function(pp) {
		term <- paste0("periodF", pp)
		x <- td %>% filter(term == .env$term)
		if (!nrow(x)) return(tibble(period = pp, RR = NA_real_, lo = NA_real_, hi = NA_real_, p_poisson = NA_real_, overdispersion = disp, p = NA_real_))
		se <- x$std.error * sqrt(disp)
		tibble(period = pp, RR = exp(x$estimate), lo = exp(x$estimate - 1.96 * se), hi = exp(x$estimate + 1.96 * se), p_poisson = x$p.value, overdispersion = disp, p = 2 * stats::pnorm(abs(x$estimate / se), lower.tail = FALSE))
	}) %>% mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj), p_label = fmt_p(p_adj))
}

fit_policy_period_counts <- function(dat, period_var = "period", ref_period, contrast_periods, min_cell_n = 20) {
	purrr::map_dfr(dxs.all, function(dx) {
		d <- dat %>% filter(dx_grp == dx, phone.luck %in% c("low", "high"), total_group_day > 0) %>%
			mutate(high = as.integer(phone.luck == "high"), periodF = stats::relevel(factor(.data[[period_var]]), ref = ref_period))
		if (!"dow" %in% names(d)) d <- d %>% mutate(dow = factor(lubridate::wday(date, label = TRUE, week_start = 1)))
		if (nrow(d) < 10 || sum(d$count, na.rm = TRUE) < 20) return(tibble(dx_grp = dx, period = contrast_periods, RR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
		fit <- tryCatch(glm(count ~ periodF + high + dow + offset(log(total_group_day)), family = poisson(), data = d), error = function(e) NULL)
		if (is.null(fit)) return(tibble(dx_grp = dx, period = contrast_periods, RR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
		td <- broom::tidy(fit)
		purrr::map_dfr(contrast_periods, function(pp) {
			term <- paste0("periodF", pp)
			x <- td %>% filter(term == .env$term)
			if (!nrow(x)) return(tibble(dx_grp = dx, period = pp, RR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
			tibble(dx_grp = dx, period = pp, RR = exp(x$estimate), lo = exp(x$estimate - 1.96 * x$std.error), hi = exp(x$estimate + 1.96 * x$std.error), p = x$p.value)
		})
	}) %>%
		group_by(period) %>% mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj), p_label = fmt_p(p_adj)) %>% ungroup() %>%
		left_join(dat %>% filter(dx_grp %in% dxs.all, phone.luck %in% c("low", "high")) %>% mutate(period_chr = as.character(.data[[period_var]])) %>% filter(period_chr %in% c(ref_period, contrast_periods)) %>% group_by(dx_grp, period_chr) %>% summarise(days = n_distinct(date), calls = sum(count, na.rm = TRUE), group_calls = sum(total_group_day, na.rm = TRUE), rate = ifelse(group_calls > 0, calls / group_calls, NA_real_), .groups = "drop") %>% rename(period = period_chr), by = c("dx_grp", "period")) %>%
		mutate(low_count = is.finite(calls) & calls < min_cell_n, label = ifelse(is.finite(RR), sprintf("%.2f%s%s\n(%.2f–%.2f)", RR, sig05, ifelse(low_count, "\u2020", ""), lo, hi), ""), short_label = ifelse(is.finite(RR), sprintf("%.2f%s%s", RR, sig05, ifelse(low_count, "\u2020", "")), ""))
}

fit_policy_mix_global_lrt <- function(dat, period_var = "period", ref_period, contrast_periods) {
	purrr::map_dfr(contrast_periods, function(pp) {
		d <- dat %>%
			filter(dx_grp %in% dxs.all, phone.luck %in% c("low", "high"), total_group_day > 0) %>%
			mutate(period_chr = as.character(.data[[period_var]])) %>%
			filter(period_chr %in% c(ref_period, pp)) %>%
			mutate(high = as.integer(phone.luck == "high"), periodF = stats::relevel(factor(period_chr), ref = ref_period), dxF = stats::relevel(factor(as.character(dx_grp), levels = dxs.all), ref = dxs.all[1]))
		if (!"dow" %in% names(d)) d <- d %>% mutate(dow = factor(lubridate::wday(date, label = TRUE, week_start = 1)))
		d <- d %>% filter(!is.na(dxF), !is.na(periodF), !is.na(dow))
		if (nrow(d) < 24 || length(unique(d$periodF)) < 2) return(tibble(period = pp, df = NA_real_, deviance = NA_real_, p_chisq = NA_real_, overdispersion = NA_real_, p_overdisp = NA_real_, p_label = "NA"))
		f0 <- tryCatch(glm(count ~ dxF + periodF + high + dow + offset(log(total_group_day)), family = poisson(), data = d), error = function(e) NULL)
		f1 <- tryCatch(glm(count ~ dxF * periodF + high + dow + offset(log(total_group_day)), family = poisson(), data = d), error = function(e) NULL)
		if (is.null(f0) || is.null(f1)) return(tibble(period = pp, df = NA_real_, deviance = NA_real_, p_chisq = NA_real_, overdispersion = NA_real_, p_overdisp = NA_real_, p_label = "NA"))
		a <- tryCatch(anova(f0, f1, test = "Chisq"), error = function(e) NULL)
		if (is.null(a) || nrow(a) < 2) return(tibble(period = pp, df = NA_real_, deviance = NA_real_, p_chisq = NA_real_, overdispersion = NA_real_, p_overdisp = NA_real_, p_label = "NA"))
		dev <- suppressWarnings(as.numeric(a$Deviance[2])); df <- suppressWarnings(as.numeric(a$Df[2])); p0 <- suppressWarnings(as.numeric(a$`Pr(>Chi)`[2]))
		disp <- suppressWarnings(sum(stats::residuals(f1, type = "pearson")^2, na.rm = TRUE) / stats::df.residual(f1))
		p_over <- ifelse(is.finite(dev) & is.finite(df) & df > 0, stats::pchisq(dev / max(1, disp), df = df, lower.tail = FALSE), NA_real_)
		tibble(period = pp, df = df, deviance = dev, p_chisq = p0, overdispersion = disp, p_overdisp = p_over, p_label = fmt_p(p_over))
	}) %>% mutate(test = "global disease-specific policy period effect; overdispersion-adjusted P used for figure title")
}

# 5A/5B: March 2022 PHSM shock
phsm_start <- as.Date("2022-03-14")
phsm_end <- as.Date("2022-03-20")
phsm_pre_days <- 14
phsm_post_days <- 14
phsm_obj <- build_daily_dx(phsm_start - phsm_pre_days, phsm_end + phsm_post_days)

phsm_total <- phsm_obj$daily_total %>%
	mutate(
		period = case_when(
			date < phsm_start ~ "Pre",
			date <= phsm_end ~ "PHSM",
			TRUE ~ "Post"
		),
		period = factor(period, levels = c("Pre", "PHSM", "Post")),
		dow = factor(lubridate::wday(date, label = TRUE, week_start = 1))
	) %>%
	group_by(phone.luck) %>%
	arrange(date) %>%
	mutate(
		count7 = roll7(total_group_day),
		baseline = mean(total_group_day[period == "Pre"], na.rm = TRUE),
		index = total_group_day / baseline,
		index7 = count7 / baseline
	) %>%
	ungroup()

fig5A_test <- fit_total_period_did(
	phsm_total %>% filter(period %in% c("Pre", "PHSM")),
	period_var = "period",
	ref_period = "Pre",
	contrast_periods = "PHSM"
)
fig5A_policy <- fit_total_policy_period(
	phsm_total %>% filter(period %in% c("Pre", "PHSM")),
	period_var = "period",
	ref_period = "Pre",
	contrast_periods = "PHSM"
)
fig5A_title <- "a. PHSM: total EMS demand"
fig5A_note <- fig5_note_plain(fig5A_policy$RR[1], fig5A_policy$p_adj[1], fig5A_test$p_adj[1])
fig5A_y_label <- "7-day rolling index\n(vs pre-PHSM baseline)"
fig5A_note_x <- min(phsm_total$date, na.rm = TRUE) + as.numeric(diff(range(phsm_total$date, na.rm = TRUE))) * fig5_note_x_frac
fig5A_note_dat <- fig5_note_plot_dat(fig5A_note_x, fig5A_policy$RR[1], fig5A_policy$p_adj[1], fig5A_test$p_adj[1], phsm_total$index7)

p5A <- ggplot(phsm_total, aes(date, index7, color = phone.luck, group = phone.luck)) +
	annotate("rect", xmin = phsm_start, xmax = phsm_end, ymin = -Inf, ymax = Inf, alpha = .12, fill = "orange") +
	geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
	geom_line(linewidth = 1.05, na.rm = TRUE) +
	geom_vline(xintercept = c(phsm_start, phsm_end), linetype = "dashed", color = "orange", linewidth = .8) +
	geom_text(data = fig5A_note_dat, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = fig5_note_size, fontface = "bold", parse = TRUE) +
	scale_color_manual(values = phone_cols, name = NULL, labels = c(low = "Low score", high = "High score")) +
	guides(color = guide_legend(nrow = 1)) +
	scale_y_continuous(expand = expansion(mult = c(.05, .24))) +
	scale_x_date(labels = date_format("%b %d", locale = "en")) +
	labs(title = fig5A_title, x = NULL, y = fig5A_y_label) +
	fig_theme(base_size = 9) +
	theme(legend.position = "top", legend.direction = "horizontal", legend.margin = margin(0, 2, 0, 2), legend.text = element_text(size = 8), axis.text = element_text(size = 9, face = "bold"), axis.title = element_text(size = 9, face = "bold"), axis.title.y = element_text(margin = margin(r = -4)), plot.title = element_text(size = 10, face = "bold", hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

phsm_dx <- phsm_obj$daily_dx %>%
	mutate(
		period = case_when(date < phsm_start ~ "Pre", date <= phsm_end ~ "PHSM", TRUE ~ "Post"),
		period = factor(period, levels = c("Pre", "PHSM", "Post"))
	)

fig5B_global <- fit_mix_global_lrt(
	phsm_dx %>% filter(period %in% c("Pre", "PHSM")),
	period_var = "period",
	ref_period = "Pre",
	contrast_periods = "PHSM"
)
fig5B_policy_global <- fit_policy_mix_global_lrt(
	phsm_dx %>% filter(period %in% c("Pre", "PHSM")),
	period_var = "period",
	ref_period = "Pre",
	contrast_periods = "PHSM"
)
fig5B_dat <- fit_period_did_counts(
	phsm_dx %>% filter(period %in% c("Pre", "PHSM")),
	period_var = "period",
	ref_period = "Pre",
	contrast_periods = "PHSM",
	min_cell_n = fig5_min_cell_n
) %>%
	mutate(dx_grp = factor(dx_grp, levels = dxs.all))
fig5B_policy_dat <- fit_policy_period_counts(
	phsm_dx %>% filter(period %in% c("Pre", "PHSM")),
	period_var = "period",
	ref_period = "Pre",
	contrast_periods = "PHSM",
	min_cell_n = fig5_min_cell_n
) %>% mutate(dx_grp = factor(dx_grp, levels = dxs.all))
fig5_response_x_limits <- c(.30, 10)
fig5_response_x_breaks <- c(.3, .4, .6, .8, 1, 1.5, 2, 3, 4, 6, 8, 10)
fig5_response_text_x_cap <- 4.6
fig5B_show <- fig5B_policy_dat %>%
	left_join(fig5B_dat %>% transmute(dx_grp, period, did_RR = RR, did_lo = lo, did_hi = hi, did_p_adj = p_adj, did_sig05 = sig05, did_short_label = short_label, did_low_count = low_count), by = c("dx_grp", "period")) %>%
	mutate(
		low_count = replace_na(low_count, TRUE),
		did_low_count = replace_na(did_low_count, TRUE),
		did_sig05 = replace_na(did_sig05, ""),
		count_flag = factor(ifelse(low_count | did_low_count, paste0("<", fig5_min_cell_n, "/cell"), paste0(">=", fig5_min_cell_n, "/cell")), levels = c(paste0(">=", fig5_min_cell_n, "/cell"), paste0("<", fig5_min_cell_n, "/cell"))),
		text_x = pmin(pmax(hi, RR) * 1.10, fig5_response_text_x_cap),
		policy_label = fig5_format_effect(RR, sig05, low_count),
		label2 = ifelse(is.finite(did_RR), paste0(policy_label, " (H/L ", sprintf("%.2f", did_RR), ")"), policy_label)
	)
fig5B_title <- "c. PHSM: disease response"

p5B <- ggplot(fig5B_show, aes(x = RR, y = fct_rev(dx_grp), color = dx_grp)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
	geom_errorbar(aes(xmin = lo, xmax = hi), width = .2, orientation = "y", linewidth = .75) +
	geom_point(aes(shape = count_flag), size = 3) +
	geom_text(aes(x = text_x, label = label2), hjust = -.05, fontface = "bold", size = 2.45, show.legend = FALSE) +
	scale_x_log10(breaks = fig5_response_x_breaks, labels = function(x) sprintf("%g", x), limits = fig5_response_x_limits) +
	scale_y_discrete(expand = expansion(add = c(.45, .55))) +
	scale_shape_manual(values = setNames(c(16, 1), c(paste0(">=", fig5_min_cell_n, "/cell"), paste0("<", fig5_min_cell_n, "/cell"))), name = "Cell count") +
	scale_color_manual(values = dxs.all.color[dxs.all], guide = "none") +
	labs(title = fig5B_title, x = "Policy-period RR vs reference (log10 scale)", y = NULL) +
	fig_theme(base_size = 9) +
	theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), axis.text = element_text(size = 9, face = "bold"), axis.title = element_text(size = 9, face = "bold"), axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 6)), legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 7.8), legend.title = element_text(size = 7.8, face = "bold"), legend.margin = margin(0, 0, 0, 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# 5C/5D: Reopening event study
open_half <- as.Date("2022-11-11")
open_full <- as.Date("2022-12-07")
open_end <- open_full + 35
open_obj <- build_daily_dx(open_half - 14, open_end)

open_total <- open_obj$daily_total %>%
	mutate(
		period = case_when(
			date < open_half ~ "Pre-open",
			date < open_full ~ "Initial relaxation",
			TRUE ~ "Full opening"
		),
		period = factor(period, levels = c("Pre-open", "Initial relaxation", "Full opening")),
		dow = factor(lubridate::wday(date, label = TRUE, week_start = 1))
	) %>%
	group_by(phone.luck) %>%
	arrange(date) %>%
	mutate(
		count7 = roll7(total_group_day),
		baseline = mean(total_group_day[period == "Pre-open"], na.rm = TRUE),
		index = total_group_day / baseline,
		index7 = count7 / baseline
	) %>%
	ungroup()

fig5C_test <- fit_total_period_did(
	open_total,
	period_var = "period",
	ref_period = "Pre-open",
	contrast_periods = c("Initial relaxation", "Full opening")
)
fig5C_policy <- fit_total_policy_period(
	open_total,
	period_var = "period",
	ref_period = "Pre-open",
	contrast_periods = c("Initial relaxation", "Full opening")
)
fig5C_title <- "b. Reopening: total EMS demand"
fig5C_main_period <- "Full opening"
fig5C_policy_main <- fig5C_policy %>% filter(period == fig5C_main_period) %>% slice_head(n = 1)
fig5C_test_main <- fig5C_test %>% filter(period == fig5C_main_period) %>% slice_head(n = 1)
if (!nrow(fig5C_policy_main)) fig5C_policy_main <- tibble(period = fig5C_main_period, RR = NA_real_)
if (!nrow(fig5C_test_main)) fig5C_test_main <- tibble(period = fig5C_main_period, p_adj = NA_real_)
fig5C_note <- fig5_note_plain(fig5C_policy_main$RR[1], fig5C_policy_main$p_adj[1], fig5C_test_main$p_adj[1])
fig5C_y_label <- "7-day rolling index\n(vs pre-open baseline)"
fig5C_note_x <- min(open_total$date, na.rm = TRUE) + as.numeric(diff(range(open_total$date, na.rm = TRUE))) * fig5_note_x_frac
fig5C_note_dat <- fig5_note_plot_dat(fig5C_note_x, fig5C_policy_main$RR[1], fig5C_policy_main$p_adj[1], fig5C_test_main$p_adj[1], open_total$index7)

p5C <- ggplot(open_total, aes(date, index7, color = phone.luck, group = phone.luck)) +
	annotate("rect", xmin = open_half, xmax = open_full, ymin = -Inf, ymax = Inf, alpha = .10, fill = "orange") +
	annotate("rect", xmin = open_full, xmax = open_end, ymin = -Inf, ymax = Inf, alpha = .08, fill = "red") +
	geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
	geom_line(linewidth = 1.05, na.rm = TRUE) +
	geom_vline(xintercept = c(open_half, open_full), linetype = "dashed", color = "orange", linewidth = .8) +
	geom_label(data = fig5C_note_dat, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = fig5_note_size + .25, fontface = "bold", parse = TRUE, fill = scales::alpha("white", .86), linewidth = .25, label.padding = grid::unit(.13, "lines")) +
	scale_color_manual(values = phone_cols, name = NULL, labels = c(low = "Low score", high = "High score")) +
	guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
	scale_y_continuous(expand = expansion(mult = c(.05, .24))) +
	scale_x_date(labels = date_format("%b %d", locale = "en")) +
	labs(title = fig5C_title, x = NULL, y = fig5C_y_label) +
	fig_theme(base_size = 9) +
	theme(legend.position = "inside", legend.position.inside = c(.025, .975), legend.justification = c(0, 1), legend.direction = "vertical", legend.background = element_rect(fill = scales::alpha("white", .86), color = "grey70", linewidth = .25), legend.margin = margin(2, 4, 2, 4), legend.text = element_text(size = 8, face = "bold"), axis.text = element_text(size = 9, face = "bold"), axis.title = element_text(size = 9, face = "bold"), axis.title.y = element_text(margin = margin(r = 2)), plot.title = element_text(size = 10, face = "bold", hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

open_dx <- open_obj$daily_dx %>%
	mutate(
		period = case_when(date < open_half ~ "Pre-open", date < open_full ~ "Initial relaxation", TRUE ~ "Full opening"),
		period = factor(period, levels = c("Pre-open", "Initial relaxation", "Full opening"))
	)

fig5D_global <- fit_mix_global_lrt(
	open_dx,
	period_var = "period",
	ref_period = "Pre-open",
	contrast_periods = c("Initial relaxation", "Full opening")
)
fig5D_policy_global <- fit_policy_mix_global_lrt(
	open_dx,
	period_var = "period",
	ref_period = "Pre-open",
	contrast_periods = c("Initial relaxation", "Full opening")
)
fig5D_dat <- fit_period_did_counts(
	open_dx,
	period_var = "period",
	ref_period = "Pre-open",
	contrast_periods = c("Initial relaxation", "Full opening"),
	min_cell_n = fig5_min_cell_n
) %>%
	mutate(dx_grp = factor(dx_grp, levels = dxs.all))
fig5D_policy_dat <- fit_policy_period_counts(
	open_dx,
	period_var = "period",
	ref_period = "Pre-open",
	contrast_periods = c("Initial relaxation", "Full opening"),
	min_cell_n = fig5_min_cell_n
) %>% mutate(dx_grp = factor(dx_grp, levels = dxs.all))
fig5D_show <- fig5D_policy_dat %>%
	left_join(fig5D_dat %>% transmute(dx_grp, period, did_RR = RR, did_lo = lo, did_hi = hi, did_p_adj = p_adj, did_sig05 = sig05, did_short_label = short_label, did_low_count = low_count), by = c("dx_grp", "period")) %>%
	mutate(
		did_sig05 = replace_na(did_sig05, ""),
		dx_grp = factor(dx_grp, levels = dxs.all),
		period = factor(period, levels = c("Full opening", "Initial relaxation")),
		label = ifelse(is.finite(RR), sprintf("%.2f%s%s；H/L %.2f%s", RR, sig05, ifelse(low_count, "\u2020", ""), did_RR, did_sig05), "")
	)
fig5D_show_all <- fig5D_show %>%
	mutate(
		low_count = replace_na(low_count, TRUE),
		did_low_count = replace_na(did_low_count, TRUE),
		did_sig05 = replace_na(did_sig05, ""),
		count_flag = factor(ifelse(low_count | did_low_count, paste0("<", fig5_min_cell_n, "/cell"), paste0(">=", fig5_min_cell_n, "/cell")), levels = c(paste0(">=", fig5_min_cell_n, "/cell"), paste0("<", fig5_min_cell_n, "/cell"))),
		text_x = pmin(pmax(hi, RR) * 1.10, fig5_response_text_x_cap),
		policy_label = fig5_format_effect(RR, sig05, low_count),
		label2 = ifelse(is.finite(did_RR), paste0(policy_label, " (H/L ", sprintf("%.2f", did_RR), ")"), policy_label)
	)
fig5D_show <- fig5D_show_all %>% filter(as.character(period) == "Full opening")
fig5D_title <- "d. Reopening: disease response"

p5D <- ggplot(fig5D_show, aes(x = RR, y = fct_rev(dx_grp), color = dx_grp)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
	geom_errorbar(aes(xmin = lo, xmax = hi), width = .2, orientation = "y", linewidth = .75) +
	geom_point(aes(shape = count_flag), size = 3) +
	geom_text(aes(x = text_x, label = label2), hjust = -.05, fontface = "bold", size = 2.45, show.legend = FALSE) +
	scale_x_log10(breaks = fig5_response_x_breaks, labels = function(x) sprintf("%g", x), limits = fig5_response_x_limits) +
	scale_y_discrete(expand = expansion(add = c(.45, .55))) +
	scale_shape_manual(values = setNames(c(16, 1), c(paste0(">=", fig5_min_cell_n, "/cell"), paste0("<", fig5_min_cell_n, "/cell"))), name = "Cell count") +
	scale_color_manual(values = dxs.all.color[dxs.all], guide = "none") +
	labs(title = fig5D_title, x = "Policy-period RR vs reference (log10 scale)", y = NULL) +
	fig_theme(base_size = 9) +
	theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), axis.text = element_text(size = 9, face = "bold"), axis.title = element_text(size = 9, face = "bold"), axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 6)), legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 7.8), legend.title = element_text(size = 7.8, face = "bold"), legend.margin = margin(0, 0, 0, 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Main Fig5 now focuses on the strong reopening signal.  PHSM panels are retained
# for the supplementary policy figure and in the workbook below.
p5A_main <- p5C + labs(title = "a. Reopening: total EMS demand")
# Rebuild the full-reopening forest without repeating the H/L DiD value in each label;
# differential H/L effects are isolated in panel d.
p5B_main <- ggplot(fig5D_show, aes(x = RR, y = fct_rev(dx_grp), color = dx_grp)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
	geom_errorbar(aes(xmin = lo, xmax = hi), width = .2, orientation = "y", linewidth = .75) +
	geom_point(aes(shape = count_flag), size = 3) +
	geom_text(aes(x = text_x, label = policy_label), hjust = -.05, fontface = "bold", size = 2.5, show.legend = FALSE) +
	scale_x_log10(breaks = fig5_response_x_breaks, labels = function(x) sprintf("%g", x), limits = fig5_response_x_limits) +
	scale_y_discrete(labels = function(x) dx_to_eng(x), expand = expansion(add = c(.45, .55))) +
	scale_shape_manual(values = setNames(c(16, 1), c(paste0(">=", fig5_min_cell_n, "/cell"), paste0("<", fig5_min_cell_n, "/cell"))), name = "Cell count") +
	scale_color_manual(values = dxs.all.color[dxs.all], guide = "none") +
	labs(title = "b. Full reopening: phenotype-specific response", x = "Policy-period RR vs pre-open reference (log10 scale)", y = NULL) +
	fig_theme(base_size = 9) +
	theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5), axis.text = element_text(size = 8.4, face = "bold"), axis.title.x = element_text(size = 8.5, face = "bold", margin = margin(t = 8)), legend.position = "inside", legend.position.inside = c(.78, .94), legend.justification = c(.5, 1), legend.background = element_rect(fill = scales::alpha("white", .82), color = NA), legend.text = element_text(size = 7.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# c. Stage-specific policy effects, using the same overdispersion-adjusted policy models.
fig5_stage_heat <- fig5D_policy_dat %>%
	filter(as.character(period) %in% c("Initial relaxation", "Full opening"), dx_grp %in% dxs.all) %>%
	mutate(
		period_label = factor(as.character(period), levels = c("Initial relaxation", "Full opening"), labels = c("Initial relaxation", "Full reopening")),
		dx_label = factor(dx_to_eng(as.character(dx_grp)), levels = rev(dx_to_eng(dxs.all))),
		log2_RR = ifelse(is.finite(RR) & RR > 0, log2(RR), NA_real_),
		cell_label = ifelse(is.finite(RR), paste0(sprintf("%.2f", RR), sig05), "")
	)
p5C_main <- ggplot(fig5_stage_heat, aes(period_label, dx_label, fill = log2_RR)) +
	geom_tile(color = "white", linewidth = .35) +
	geom_text(aes(label = cell_label), size = 2.45, fontface = "bold") +
	scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "log2 RR", na.value = "grey95") +
	labs(title = "c. Phenotype shifts across reopening stages", x = NULL, y = NULL) +
	fig_theme(base_size = 8.4) +
	theme(plot.title = element_text(size = 10, hjust = .5), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 7.4), legend.position = "right")

# d. Differential high-vs-low response is shown explicitly as secondary heterogeneity.
fig5_full_global_p <- fig5D_global %>% filter(as.character(period) == "Full opening") %>% pull(p_overdisp) %>% dplyr::first(default = NA_real_)
fig5_did_full <- fig5D_dat %>%
	filter(as.character(period) == "Full opening", dx_grp %in% dxs.all) %>%
	mutate(
		dx_label = factor(dx_to_eng(as.character(dx_grp)), levels = rev(dx_to_eng(dxs.all))),
		low_count = replace_na(as.logical(low_count), TRUE),
		label_did = ifelse(is.finite(RR), paste0(sprintf("%.2f", RR), sig05, ifelse(low_count, "†", "")), "")
	)
p5D_main <- ggplot(fig5_did_full, aes(RR, dx_label, color = dx_grp)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
	geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .16, linewidth = .65) +
	geom_point(aes(shape = low_count), size = 2.45) +
	geom_text(aes(x = pmin(pmax(hi, RR) * 1.06, 3.5), label = label_did), hjust = 0, size = 2.25, fontface = "bold", show.legend = FALSE) +
	scale_x_log10(breaks = c(.3, .5, .8, 1, 1.5, 2, 3), limits = c(.25, 4)) +
	scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c(`FALSE` = paste0(">=", fig5_min_cell_n, "/cell"), `TRUE` = paste0("<", fig5_min_cell_n, "/cell")), name = "Cell count") +
	scale_color_manual(values = dxs.all.color[dxs.all], guide = "none") +
	labs(
		title = "d. High-vs-low differential response after full reopening",
		x = "DiD rate ratio: high vs low (log10 scale)", y = NULL
	) +
	fig_theme(base_size = 8.4) +
	theme(plot.title = element_text(size = 10, hjust = .5), plot.subtitle = element_text(size = 7.4), axis.text.y = element_text(size = 7.4), axis.title.x = element_text(margin = margin(t = 8)), legend.position = "none")

# Assemble rows independently.  The long phenotype labels in the second row must
# not determine the left grob width of panel a; otherwise its y title is pushed
# away from the axis.  align_panel_rows() aligns only within each row and adds
# real grid gutters/spacer rows instead of compensating with plot.margin.
Fig5 <- align_panel_rows(
	rows = list(
		list(p5A_main, p5C_main),
		list(p5B_main, p5D_main)
	),
	rel_widths = c(1.02, .98),
	rel_heights = c(.82, 1.06),
	row_gap = .08,
	side_pad = .018,
	axis_text_size = 9,
	axis_title_size = 9
)
save_plot(Fig5, "Fig5.png", width = 12.8, height = 9.6, dpi = 600, bg = "white")

phsm_period_dx_summary <- phsm_dx %>%
	group_by(period, dx_grp, phone.luck) %>%
	summarise(
		days = n_distinct(date),
		calls = sum(count, na.rm = TRUE),
		group_calls = sum(total_group_day, na.rm = TRUE),
		rate = ifelse(group_calls > 0, calls / group_calls, NA_real_),
		.groups = "drop"
	)

open_period_dx_summary <- open_dx %>%
	group_by(period, dx_grp, phone.luck) %>%
	summarise(
		days = n_distinct(date),
		calls = sum(count, na.rm = TRUE),
		group_calls = sum(total_group_day, na.rm = TRUE),
		rate = ifelse(group_calls > 0, calls / group_calls, NA_real_),
		.groups = "drop"
	)

fig5_event_windows <- tibble(
	event = c("PHSM", "Reopening initial relaxation", "Reopening full opening"),
	start = c(phsm_start, open_half, open_full),
	end = c(phsm_end, open_full - 1, open_end),
	reference = c(sprintf("%s to %s", phsm_start - phsm_pre_days, phsm_start - 1), sprintf("%s to %s", open_half - 14, open_half - 1), sprintf("%s to %s", open_half - 14, open_half - 1))
)
fig5_figure_annotations <- tibble(
	panel = c("a", "b", "c", "d"),
	contrast = c("PHSM vs Pre", "Full opening vs Pre-open", "PHSM disease response", "Full opening disease response"),
	figure_label = c(fig5A_note, fig5C_note, NA_character_, NA_character_),
	policy_global_p = c(NA_character_, NA_character_, fig5B_policy_global$p_label[1], fig5D_policy_global %>% filter(period == "Full opening") %>% pull(p_label) %>% dplyr::first(default = NA_character_)),
	high_low_global_p = c(NA_character_, NA_character_, fig5B_global$p_label[1], fig5D_global %>% filter(period == "Full opening") %>% pull(p_label) %>% dplyr::first(default = NA_character_))
)

writexl::write_xlsx(list(
	event_windows = fig5_event_windows,
	figure_annotation_tests = fig5_figure_annotations,
	low_count_rule = tibble(min_cell_n = fig5_min_cell_n, note = "Dagger/empty point flags any high/low × reference/event disease cell below this count; estimates are retained but should be interpreted cautiously."),
	PHSM_total_event_index = phsm_total,
	PHSM_period_dx_summary = phsm_period_dx_summary,
	PHSM_total_policy_models = fig5A_policy,
	PHSM_total_DID_models = fig5A_test,
	PHSM_policy_global_mix_test = fig5B_policy_global,
	PHSM_global_mix_test = fig5B_global,
	PHSM_policy_models = fig5B_policy_dat,
	PHSM_DID_models = fig5B_dat,
	PHSM_combined_policy_DID = fig5B_show,
	reopening_total_event_index = open_total,
	reopening_period_dx_summary = open_period_dx_summary,
	reopen_total_policy_models = fig5C_policy,
	reopening_total_DID_models = fig5C_test,
	reopen_policy_global_mix = fig5D_policy_global,
	reopening_global_mix_tests = fig5D_global,
	reopening_policy_models = fig5D_policy_dat,
	reopening_DID_models = fig5D_dat,
	reopen_combined_all_periods = fig5D_show_all,
	reopen_combined_plotted = fig5D_show,
	main_reopening_stage_heatmap = fig5_stage_heat,
	main_full_reopening_HL_DID = fig5_did_full,
	main_figure_note = tibble(note = "Main Fig5 focuses on reopening. PHSM is retained in FigS5 and in the PHSM workbook sheets.")
), "Fig5.out.xlsx")
ems120_maybe_exit_after("fig5")
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS2. EMS on-scene time and circular hourly mechanism (former Fig6)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (ems120_should_run("fig6")) {
fig4_base_all <- bind_rows(lapply(years, function(y) {
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())
	d0 %>%
		filter(dx_grp %in% dxs.all) %>%
		transmute(
			Year = as.integer(y),
			dx_grp = factor(as.character(dx_grp), levels = dxs.all),
			phone.luck = factor(as.character(phone.luck), levels = c("low", "middle", "high")),
			hour = suppressWarnings(as.integer(hour)),
			dispatch = suppressWarnings(as.numeric(.data[[vars.basic.ems[5]]]) / 60),
			onsite = suppressWarnings(as.numeric(现场时间) / 60)
		)
}))
fig4_base <- fig4_base_all %>% filter(dx_grp %in% dxs.vip) %>% mutate(dx_grp = factor(as.character(dx_grp), levels = dxs.vip))
fig4_all_dx_onsite_summary <- fig4_base_all %>%
	filter(is.finite(onsite)) %>%
	group_by(dx_grp, Year) %>%
	summarise(mean = mean(onsite, na.rm = TRUE), median = median(onsite, na.rm = TRUE), sd = sd(onsite, na.rm = TRUE), n = n(), .groups = "drop")

fig_hourly_all <- bind_rows(lapply(years, function(y) {
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())

	d0 %>%
		filter(
			dx_grp %in% dxs.all,
			phone.luck %in% c("low", "middle", "high"),
			is.finite(hour),
			between(as.integer(hour), 0, 23)
		) %>%
		transmute(
			year = as.integer(y),
			dx_grp = factor(as.character(dx_grp), levels = dxs.all),
			phone.luck = factor(as.character(phone.luck), levels = c("low", "middle", "high")),
			hour = as.integer(hour)
		) %>%
		count(year, dx_grp, phone.luck, hour, name = "call_n") %>%
		complete(
			year = as.integer(y),
			dx_grp = dxs.all,
			phone.luck = c("low", "middle", "high"),
			hour = 0:23,
			fill = list(call_n = 0)
		) %>%
		group_by(year, dx_grp, phone.luck) %>%
		mutate(
			total_n = sum(call_n, na.rm = TRUE),
			pct = ifelse(total_n > 0, call_n / total_n, NA_real_),
			pct_smooth = (call_n + 0.5) / (total_n + 24 * 0.5)
		) %>%
		ungroup()
}))

fig4A_dat <- fig4_base %>%
	filter(is.finite(onsite)) %>%
	group_by(dx_grp, Year) %>%
	summarise(
		mean = mean(onsite, na.rm = TRUE),
		median = median(onsite, na.rm = TRUE),
		n = n(),
		.groups = "drop"
	)

fig4_luck_summary <- function(value_col) {
	fig4_base %>%
		filter(phone.luck %in% c("low", "high"), is.finite(.data[[value_col]])) %>%
		group_by(Year, phone.luck) %>%
		summarise(
			mean = mean(.data[[value_col]], na.rm = TRUE),
			median = median(.data[[value_col]], na.rm = TRUE),
			n = n(),
			.groups = "drop"
		) %>%
		pivot_wider(names_from = phone.luck, values_from = c(mean, median, n), names_sep = "_") %>%
		mutate(
			diff = mean_high - mean_low,
			p = purrr::map_dbl(Year, function(yy) {
				d0 <- fig4_base %>% filter(Year == yy, phone.luck %in% c("low", "high"), is.finite(.data[[value_col]]))
				suppressWarnings(tryCatch(
					wilcox.test(d0[[value_col]][d0$phone.luck == "high"], d0[[value_col]][d0$phone.luck == "low"])$p.value,
					error = function(e) NA_real_
				))
			}),
			p_adj = p.adjust(p, "BH"),
			sig = !is.na(p_adj) & p_adj < .01
		)
}
fig4_add_star_x <- function(d) {
	rng <- range(c(d$mean_low, d$mean_high), na.rm = TRUE)
	eps <- if (all(is.finite(rng)) && diff(rng) > 0) 0.04 * diff(rng) else 0.1
	d %>% mutate(star_x = pmax(mean_low, mean_high, na.rm = TRUE) + eps)
}
fig4B_dat <- fig4_luck_summary("onsite") %>% fig4_add_star_x()
fig4C_dat <- fig4_luck_summary("dispatch") %>% fig4_add_star_x()
gm2_high <- mean(fig4B_dat$mean_high, na.rm = TRUE)
gm2_low <- mean(fig4B_dat$mean_low, na.rm = TRUE)
gm3_high <- mean(fig4C_dat$mean_high, na.rm = TRUE)
gm3_low <- mean(fig4C_dat$mean_low, na.rm = TRUE)
fig4_dx_cols <- setNames(grDevices::adjustcolor(unname(dxs.vip.color[dxs.vip]), alpha.f = 0.72), dxs.vip)
fig4_title_size <- 10.5
fig4_axis_title_size <- 8.5
fig4_axis_text_size <- 8.0
fig4_legend_text_size <- 7.2
fig4_width <- 11.2

p4A <- ggplot(fig4A_dat, aes(x = mean, y = Year, color = dx_grp)) +
	geom_point(size = 2.6) +
	scale_color_manual(values = fig4_dx_cols[dxs.vip], name = NULL, drop = FALSE) +
	guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
	scale_y_continuous(breaks = years, labels = years) +
	labs(title = "a. On-scene time by phenotype", x = "Time (mins)", y = "Year") +
	fig_theme(base_size = 9.5) +
	theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = fig4_legend_text_size, color = "grey20", lineheight = .75), legend.spacing.y = grid::unit(0, "pt"), legend.key.height = grid::unit(7, "pt"), legend.margin = margin(0, 0, 0, 0), plot.title = element_text(size = fig4_title_size, face = "bold", hjust = 0.5), axis.title = element_text(size = fig4_axis_title_size, face = "bold"), axis.text = element_text(size = fig4_axis_text_size, face = "bold"))

fig4_luck_plot <- function(dat, title, xintercept_high, xintercept_low, accent_color = "blue") {
	ggplot(dat, aes(y = Year)) +
		geom_segment(aes(x = mean_low, xend = mean_high, yend = Year), color = "black", linewidth = .45) +
		geom_point(aes(x = mean_low, color = "Low luck", shape = "Low luck"), size = 2.8) +
		geom_point(aes(x = mean_high, color = "High luck", shape = "High luck"), size = 3.0) +
		geom_vline(xintercept = xintercept_low, color = unname(phone_cols["low"]), linetype = "dashed", linewidth = .8) +
			geom_vline(xintercept = xintercept_high, color = accent_color, linetype = "dashed", linewidth = .8) +
		geom_text(data = dplyr::filter(dat, sig), aes(x = star_x, label = "*"), fontface = "bold", size = 5, vjust = .35) +
		labs(title = title, x = "Time (mins)", y = "Year") +
		scale_y_continuous(breaks = years, labels = years) +
		scale_color_manual(values = c("High luck" = accent_color, "Low luck" = unname(phone_cols["low"])), breaks = c("High luck", "Low luck"), name = NULL) +
		scale_shape_manual(values = c("High luck" = 17, "Low luck" = 16), breaks = c("High luck", "Low luck"), name = NULL) +
		guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
		fig_theme(base_size = 9.5) +
		theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = fig4_legend_text_size), legend.margin = margin(0, 0, 0, 0), plot.title = element_text(size = fig4_title_size, face = "bold", hjust = 0.5), axis.title = element_text(size = fig4_axis_title_size, face = "bold"), axis.text = element_text(size = fig4_axis_text_size, face = "bold"))
}
p4B <- fig4_luck_plot(fig4B_dat, "b. On-scene time by luck group", gm2_high, gm2_low)
p4C <- fig4_luck_plot(fig4C_dat, "c. Dispatch time by luck group", gm3_high, gm3_low, accent_color = "#D55E00")

# Selected circular panels can be changed here after inspecting Fig4.out.xlsx / hourly_summary_by_year.
dxs2 <- dxs.vip
years2 <- 2018:2023
fig4_circle_nrow <- 2L
fig4_circle_ncol <- ceiling(length(years2) / fig4_circle_nrow)

fig4d_smooth_circ <- function(x, k = 3L) {
	x <- as.numeric(x)
	n <- length(x)
	if (!n) return(x)
	if (n == 1) return(x)
	h <- floor(k / 2)
	sapply(seq_len(n), function(i) {
		idx <- ((i - h - 1):(i + h - 1)) %% n + 1
		mean(x[idx], na.rm = TRUE)
	})
}

yin_path_candidates <- c(file.path(dir0, "files", "YinYang.png"), "YinYang.png", file.path(dir0, "files", "YinYang_clean.png"), "YinYang_clean.png")
yin_path <- yin_path_candidates[file.exists(yin_path_candidates)][1]
yin_img <- NULL
yin_center_radius <- 0.315
if (!is.na(yin_path) && length(yin_path) == 1 && file.exists(yin_path)) {
	yin_img <- read_center_image(yin_path, circle_alpha = TRUE)
} else {
	warning("YinYang image not found. Put YinYang_clean.png in the working directory or D:/files/.")
}

fig4D_circle_dat <- fig_hourly_all %>%
	filter(
		year %in% years2,
		dx_grp %in% dxs2,
		phone.luck %in% c("low", "high")
	) %>%
	mutate(
		dx_grp = factor(as.character(dx_grp), levels = dxs2),
		phone.luck = factor(as.character(phone.luck), levels = c("low", "high"))
	)

fig4D_hourly_contrast <- fig4D_circle_dat %>%
	dplyr::select(year, dx_grp, phone.luck, hour, call_n, total_n, pct) %>%
	pivot_wider(names_from = phone.luck, values_from = c(call_n, total_n, pct), names_sep = "_") %>%
	mutate(
		diff_pp = 100 * (pct_high - pct_low)
	) %>%
	group_by(year, dx_grp) %>%
	arrange(hour, .by_group = TRUE) %>%
	mutate(diff_pp_smooth = fig4d_smooth_circ(diff_pp, k = 3L)) %>%
	ungroup()

fig4D_scale <- fig4D_hourly_contrast %>%
	group_by(year) %>%
	summarise(max_abs = max(abs(diff_pp), abs(diff_pp_smooth), na.rm = TRUE), .groups = "drop") %>%
	mutate(max_abs = pmax(max_abs, 0.45))
fig4D_scale_map <- setNames(fig4D_scale$max_abs, fig4D_scale$year)

fig4D_circle_file <- tempfile(fileext = ".png")
png(
	filename = fig4D_circle_file,
	width = fig4_width + 1.2,
	height = 8.65,
	units = "in",
	res = 420,
	bg = "transparent"
)

par(
	mfrow = c(fig4_circle_nrow, fig4_circle_ncol),
	mar = c(0.72, 0.00, 1.72, 0.00),
	oma = c(0.60, 0.00, 0.10, 0.00),
	xaxs = "i",
	yaxs = "i"
)

for (yy in years2) {
	dat2 <- fig4D_hourly_contrast %>% filter(year == yy) %>% arrange(dx_grp, hour)
	bg_col <- setNames(grDevices::adjustcolor(unname(fig4_dx_cols[dxs2]), alpha.f = 0.17), dxs2)
	y_abs <- unname(fig4D_scale_map[as.character(yy)])
	if (!is.finite(y_abs) || y_abs <= 0) y_abs <- 0.5
	circlize::circos.clear()
	circlize::circos.par(start.degree = 90, gap.degree = 2, cell.padding = c(0, 0, 0, 0), track.margin = c(0.0015, 0.0015), canvas.xlim = c(-1.10, 1.10), canvas.ylim = c(-1.13, 1.13))
	circlize::circos.initialize(factors = "all", xlim = c(0, 24))

	for (i in seq_along(dxs2)) {
		dx <- dxs2[i]
		dx_dat <- dat2 %>% filter(dx_grp == dx) %>% arrange(hour)
		circlize::circos.trackPlotRegion(
			factors = "all",
			track.index = i,
			ylim = c(-1.20 * y_abs, 1.20 * y_abs),
			bg.col = bg_col[dx],
			bg.border = NA,
			track.height = min(0.12, 0.72 / length(dxs2)),
			panel.fun = function(...) {
				circlize::circos.lines(c(0, 24), c(0, 0), col = "grey62", lwd = 0.75)
				for (jj in seq_len(nrow(dx_dat))) {
					val <- dx_dat$diff_pp[jj]
					if (!is.finite(val)) next
					fill_col <- if (val >= 0) grDevices::adjustcolor(fig4_dx_cols[dx], alpha.f = 0.80) else grDevices::adjustcolor("#4C78A8", alpha.f = 0.82)
					circlize::circos.rect(dx_dat$hour[jj], min(0, val), dx_dat$hour[jj] + 1, max(0, val), col = fill_col, border = NA)
				}
				circlize::circos.lines(dx_dat$hour + 0.5, clamp(dx_dat$diff_pp_smooth, -1.20 * y_abs, 1.20 * y_abs), col = "black", lwd = 1.00)
			}
		)
	}

	circlize::circos.axis(h = "top", major.at = 0:23, labels = sprintf("%02d", 0:23), labels.cex = 0.86, labels.font = 2, minor.ticks = 0, sector.index = "all", track.index = 1)
	if (!is.null(yin_img)) graphics::rasterImage(image = yin_img, xleft = -yin_center_radius, ybottom = -yin_center_radius, xright = yin_center_radius, ytop = yin_center_radius, interpolate = TRUE)
	title(main = yy, font.main = 2, cex.main = 1.25, col.main = "#009E73", line = 0.10)
}

dev.off()
circlize::circos.clear()

fig4D_raster <- png::readPNG(fig4D_circle_file, native = FALSE)
if (file.exists(fig4D_circle_file)) file.remove(fig4D_circle_file)
fig4D_legend_dat <- tibble(
	dx_grp = factor(dxs2, levels = dxs2),
	label = dx_to_eng(dxs2),
	x = seq(0.035, 0.875, length.out = length(dxs2)),
	x_text = x + 0.012,
	y = 0.982
)
p4D <- ggplot() +
	annotation_custom(grid::rasterGrob(fig4D_raster, interpolate = TRUE), xmin = -0.325, xmax = 1.280, ymin = -0.392, ymax = 0.952) +
	annotate("rect", xmin = 0.145, xmax = 0.165, ymin = 1.020, ymax = 1.045, fill = grDevices::adjustcolor(fig4_dx_cols[dxs2[1]], alpha.f = 0.80), color = NA) +
	annotate("text", x = 0.173, y = 1.033, label = "High > low hourly share", hjust = 0, size = 2.35, fontface = "bold") +
	annotate("rect", xmin = 0.395, xmax = 0.415, ymin = 1.020, ymax = 1.045, fill = grDevices::adjustcolor("#4C78A8", alpha.f = 0.82), color = NA) +
	annotate("text", x = 0.423, y = 1.033, label = "Low > high hourly share", hjust = 0, size = 2.35, fontface = "bold") +
	annotate("segment", x = 0.645, xend = 0.675, y = 1.033, yend = 1.033, linewidth = 0.85, color = "black") +
	annotate("text", x = 0.685, y = 1.033, label = "3-hour circular moving average", hjust = 0, size = 2.35, fontface = "bold") +
	geom_point(data = fig4D_legend_dat, aes(x = x, y = y, color = dx_grp), inherit.aes = FALSE, size = 1.85) +
	geom_text(data = fig4D_legend_dat, aes(x = x_text, y = y, label = label), inherit.aes = FALSE, fontface = "bold", size = 2.38, color = "grey20", hjust = 0) +
	scale_color_manual(values = fig4_dx_cols[dxs2], guide = "none", drop = FALSE) +
	coord_cartesian(xlim = c(0, 1), ylim = c(-0.275, 1.085), expand = FALSE, clip = "off") +
	labs(title = "c. Circadian EMS phenotype by luck group") +
	theme_void(base_size = 11) +
	theme(plot.title = element_text(face = "bold", hjust = 0.5, size = fig4_title_size, margin = margin(t = 4, b = 5)), plot.margin = margin(t = 5, r = 0, b = 24, l = 0))

# Main Fig4 keeps only the two most informative upper panels.  The former panel c
# (dispatch time) remains in Fig4.out.xlsx but is removed from the main figure;
# the circular panel is therefore relabelled from d to c.
Fig4 <- (p4A | p4B) / p4D + plot_layout(heights = c(1.20, 2.78), widths = c(1, 1))
save_plot(Fig4, "FigS2.png", width = fig4_width, height = 12.6, dpi = 600, bg = "transparent")

fig4_hourly_diag <- fig_hourly_all %>%
	filter(phone.luck %in% c("low", "high")) %>%
	group_by(dx_grp, phone.luck) %>%
	summarise(
		total_calls = sum(call_n, na.rm = TRUE),
		mean_hour = ifelse(total_calls > 0, round(sum(hour * call_n, na.rm = TRUE) / total_calls, 2), NA_real_),
		night_share = ifelse(total_calls > 0, round(sum(call_n[hour <= 5], na.rm = TRUE) / total_calls, 4), NA_real_),
		morning_share = ifelse(total_calls > 0, round(sum(call_n[hour %in% 6:11], na.rm = TRUE) / total_calls, 4), NA_real_),
		afternoon_share = ifelse(total_calls > 0, round(sum(call_n[hour %in% 12:17], na.rm = TRUE) / total_calls, 4), NA_real_),
		evening_share = ifelse(total_calls > 0, round(sum(call_n[hour %in% 18:23], na.rm = TRUE) / total_calls, 4), NA_real_),
		.groups = "drop"
	)

fig4_workflow_by_year <- setNames(lapply(years, function(y) {
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())
	d0 %>% filter(phone.luck %in% c("low", "high")) %>%
		transmute(
			Year = y,
			phone.luck = factor(as.character(phone.luck), levels = c("low", "high")),
			Dispatch = suppressWarnings(as.numeric(.data[[vars.basic.ems[5]]]) / 60),
			Driving = suppressWarnings(as.numeric(.data[[vars.basic.ems[6]]]) / 60),
			Return = suppressWarnings(as.numeric(.data[[vars.basic.ems[8]]]) / 60)
		)
}), as.character(years))

fig4_workflow_long <- bind_rows(fig4_workflow_by_year) %>%
	pivot_longer(cols = c(Dispatch, Driving, Return), names_to = "interval", values_to = "time_min") %>%
	filter(is.finite(time_min))

fig4_workflow_summary <- fig4_workflow_long %>%
	group_by(interval, Year, phone.luck) %>%
	summarise(mean = mean(time_min, na.rm = TRUE), median = median(time_min, na.rm = TRUE), n = n(), .groups = "drop") %>%
	pivot_wider(names_from = phone.luck, values_from = c(mean, median, n), names_sep = "_") %>%
	mutate(diff = mean_high - mean_low)
saveRDS(fig4_workflow_long, .ems120_register_temp_rds("Fig4.rds"))
.ems120_remove_temp_rds("Fig4.rds")

fig4_hourly_summary <- fig_hourly_all %>%
	filter(phone.luck %in% c("low", "high")) %>%
	group_by(year, dx_grp, phone.luck) %>%
	summarise(
		total_calls = sum(call_n, na.rm = TRUE),
		peak_hour = hour[which.max(pct_smooth)],
		peak_pct = round(max(pct_smooth, na.rm = TRUE), 4),
		night_share = ifelse(total_calls > 0, round(sum(call_n[hour <= 5], na.rm = TRUE) / total_calls, 4), NA_real_),
		morning_share = ifelse(total_calls > 0, round(sum(call_n[hour %in% 6:11], na.rm = TRUE) / total_calls, 4), NA_real_),
		afternoon_share = ifelse(total_calls > 0, round(sum(call_n[hour %in% 12:17], na.rm = TRUE) / total_calls, 4), NA_real_),
		evening_share = ifelse(total_calls > 0, round(sum(call_n[hour %in% 18:23], na.rm = TRUE) / total_calls, 4), NA_real_),
		.groups = "drop"
	)

fig4D_circle_summary <- fig4D_hourly_contrast %>%
	group_by(year, dx_grp) %>%
	summarise(
		peak_enrichment_hour = hour[which.max(diff_pp_smooth)],
		trough_enrichment_hour = hour[which.min(diff_pp_smooth)],
		max_diff_pp = round(max(diff_pp_smooth, na.rm = TRUE), 3),
		min_diff_pp = round(min(diff_pp_smooth, na.rm = TRUE), 3),
		night_diff_pp = round(mean(diff_pp[hour <= 5], na.rm = TRUE), 3),
		morning_diff_pp = round(mean(diff_pp[hour %in% 6:11], na.rm = TRUE), 3),
		afternoon_diff_pp = round(mean(diff_pp[hour %in% 12:17], na.rm = TRUE), 3),
		evening_diff_pp = round(mean(diff_pp[hour %in% 18:23], na.rm = TRUE), 3),
		.groups = "drop"
	)

figS2_main_wb <- list(
	panelA_data = fig4A_dat,
		all_dx_onsite_summary = fig4_all_dx_onsite_summary,
	panelB_data = fig4B_dat,
	panelC_dispatch_data = fig4C_dat,
	hourly_summary_by_year = fig4_hourly_summary,
	panelD_selected_circle_summary = fig4D_circle_summary,
	panelD_hourly_contrast = fig4D_hourly_contrast,
	hourly_diagnostics = fig4_hourly_diag,
	workflow_intervals_summary = fig4_workflow_summary,
	by_dx_high_low = fig4_base %>%
		filter(phone.luck %in% c("low", "high"), is.finite(onsite)) %>%
		group_by(dx_grp, phone.luck) %>%
		summarise(
			mean = round(mean(onsite), 2),
			median = round(median(onsite), 2),
			sd = round(sd(onsite), 2),
			n = n(),
			.groups = "drop"
		),
	panelD_config = tibble(
		years2 = paste(years2, collapse = ", "),
		dxs2 = paste(dxs2, collapse = ", "),
		layout_rows = fig4_circle_nrow,
		layout_cols = fig4_circle_ncol,
		yin_yang_file = ifelse(is.null(yin_img), NA_character_, yin_path),
		yin_center_radius = yin_center_radius
	)
)
writexl::write_xlsx(figS2_main_wb, "FigS2.out.xlsx")
ems120_maybe_exit_after("fig6")
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS1. Yearly composition of all EMS disease types
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS1_dat <- bind_rows(lapply(years, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(!is.na(dx_grp)) %>%
		transmute(year = as.integer(y), dx_grp = factor(as.character(dx_grp), levels = dxs.all0))
})) %>%
	count(year, dx_grp, name = "n") %>%
	complete(year = years, dx_grp = dxs.all0, fill = list(n = 0)) %>%
	group_by(year) %>%
	mutate(total = sum(n), pct = ifelse(total > 0, n / total, NA_real_), label = ifelse(is.finite(pct) & pct >= .018, scales::percent(pct, accuracy = .1), "")) %>%
	ungroup() %>%
	mutate(dx_grp = factor(dx_grp, levels = rev(dxs.all0)), disease_label = factor(dx_to_eng(as.character(dx_grp)), levels = dx_to_eng(rev(dxs.all0))))

pS1 <- ggplot(figS1_dat, aes(x = factor(year), y = pct, fill = disease_label)) +
	geom_col(width = .86, color = "white", linewidth = .25) +
	geom_text(aes(label = label), position = position_stack(vjust = .5), size = 2.35, color = "black") +
	scale_fill_manual(values = setNames(dxs.all.color[rev(dxs.all0)], dx_to_eng(rev(dxs.all0))), name = "Category:", drop = FALSE) +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, .02))) +
	labs(title = NULL, x = "Year", y = "Percentage") +
	fig_theme(base_size = 9.0) +
	theme(axis.text.x = element_text(face = "bold"), axis.title = element_text(face = "bold"), legend.position = "right", legend.title = element_text(face = "bold"), legend.text = element_text(size = 7.6), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

save_plot(pS1, "FigS1.png", width = 8.6, height = 8.0, dpi = 600)
writexl::write_xlsx(list(yearly_composition = figS1_dat %>% arrange(year, desc(pct)), yearly_totals = figS1_dat %>% group_by(year) %>% summarise(total = max(total), check_pct_sum = sum(pct, na.rm = TRUE), .groups = "drop"), mapping = tibble(type = dxs.all0, english = dx_to_eng(dxs.all0), chinese = vapply(dxs.type.list[dxs.all0], function(x) paste(x, collapse = "; "), character(1)))), "FigS1.out.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS2. Phone score distributions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS2_score_dat <- phone_score_bin_summary %>%
	filter(!is.na(score_bin), score_bin != "NA", n > 0) %>%
	group_by(score_name, score_bin) %>%
	summarise(n = sum(n), .groups = "drop") %>%
	left_join(phone_quarter_summary, by = "score_name") %>%
	mutate(
		score_name = factor(score_name, levels = phone_score_vars),
		score_bin = factor(score_bin, levels = c("0", sprintf("(%d,%d]", 0:9, 1:10))),
		threshold_label = sprintf("low <= %s; high >= %s", low_cutoff, high_cutoff)
	) %>% filter(!is.na(score_bin))

plot_score_hist <- function(v, ttl, col) {
	d <- figS2_score_dat %>% filter(score_name == v) %>% mutate(pct = n / sum(n))
	ggplot(d, aes(score_bin, n)) +
		geom_col(width = 0.82, fill = col, color = "white", linewidth = .25) +
		geom_text(aes(label = ifelse(pct >= .035, scales::percent(pct, accuracy = .1), "")), vjust = -.25, size = 2.35, fontface = "bold") +
		geom_text(data = d %>% distinct(threshold_label), aes(x = 5.5, y = Inf, label = threshold_label), inherit.aes = FALSE, hjust = 0.5, vjust = 1.30, size = 2.7, fontface = "bold") +
		scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, .22))) +
		labs(title = ttl, x = NULL, y = "Calls, N") +
		fig_theme(base_size = 8.6) +
		theme(plot.title = element_text(face = "bold", size = 9.3, hjust = 0.5), axis.text.x = element_text(angle = 35, hjust = 1, size = 6.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
figS2_sco3_tail <- phone_score_summary %>%
	filter(score_name == "phone.sco3", !is.na(score), score != "NA") %>%
	mutate(score_num = suppressWarnings(as.numeric(score))) %>%
	filter(is.finite(score_num)) %>%
	group_by(score_num) %>% summarise(n = sum(n), .groups = "drop") %>%
	arrange(desc(score_num)) %>%
	mutate(cum_n = cumsum(n), total_n = sum(n), pct_at_or_above = cum_n / total_n) %>%
	arrange(score_num) %>%
	mutate(se = sqrt(pct_at_or_above * (1 - pct_at_or_above) / pmax(total_n, 1)), ci_lo = pmax(0, pct_at_or_above - 1.96 * se), ci_hi = pmin(1, pct_at_or_above + 1.96 * se))
figS2_tail_compare <- phone_score_summary %>%
	filter(score_name %in% c("phone.sco1", "phone.sco2", "phone.sco3"), !is.na(score), score != "NA") %>%
	mutate(score_num = suppressWarnings(as.numeric(score))) %>%
	filter(is.finite(score_num)) %>%
	group_by(score_name, score_num) %>% summarise(n = sum(n), .groups = "drop") %>%
	group_by(score_name) %>%
	arrange(desc(score_num), .by_group = TRUE) %>%
	mutate(cum_n = cumsum(n), total_n = sum(n), pct_at_or_above = cum_n / total_n) %>%
	ungroup() %>%
	mutate(score_method = factor(dplyr::recode(score_name, `phone.sco1` = "sco1: simple rule", `phone.sco2` = "sco2: advanced rule", `phone.sco3` = "sco3: ML score"), levels = c("sco1: simple rule", "sco2: advanced rule", "sco3: ML score"))) %>%
	arrange(score_method, score_num)
pS2A <- plot_score_hist("phone.sco0", "a. sco0 (pattern rule)", "#0072B2")
pS2B <- plot_score_hist("phone.sco1", "b. sco1 (simple rule)", "#D55E00")
pS2C <- plot_score_hist("phone.sco2", "c. sco2 (advanced rule)", "#009E73")
pS2D <- ggplot(figS2_tail_compare, aes(score_num, pct_at_or_above, color = score_method, group = score_method)) +
	geom_step(linewidth = .88) +
	geom_point(size = 1.35, alpha = .82) +
	scale_color_manual(values = c("sco1: simple rule" = "#9ECAE1", "sco2: advanced rule" = "#3182BD", "sco3: ML score" = "#08519C"), name = NULL) +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, .08))) +
	scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
	labs(title = "d. Top-tail curves by score construction", x = "Score cutoff", y = "% at or above cutoff") +
	fig_theme(base_size = 8.6) +
	theme(plot.title = element_text(face = "bold", size = 9.3, hjust = 0.5), axis.title.x = element_text(margin = margin(t = 1)), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigS2 <- (pS2A | pS2B) / (pS2C | pS2D)
save_plot(FigS2, "FigS2.png", width = 8.6, height = 6.4, dpi = 600)
writexl::write_xlsx(list(score_distribution = phone_score_summary, score_distribution_integer_bins = phone_score_bin_summary, quasi_quarter_cutoffs = phone_quarter_summary, plot_data = figS2_score_dat, sco3_top_tail_curve = figS2_sco3_tail, three_score_top_tail_curves = figS2_tail_compare), "FigS2.out.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS3. Address-type vs geo.type1 concordance and address keywords
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS3_geo_base <- bind_rows(lapply(years, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		transmute(
			year = as.integer(y),
			address = as.character(.data[[vars.basic.ems[2]]]),
			addr_type = as.character(.data[["地址类型"]]),
			geo_type1 = as.character(.data[["geo.type1"]])
		)
})) %>%
	mutate(
		address = replace_na(address, ""),
		addr_type = ifelse(is.na(addr_type) | trimws(addr_type) == "", "Missing", trimws(addr_type)),
		geo_type1 = ifelse(is.na(geo_type1) | trimws(geo_type1) == "", "Missing", trimws(geo_type1))
	)
figS3_geo_type_levels <- figS3_geo_base %>%
	filter(geo_type1 != "Missing") %>%
	count(geo_type1, sort = TRUE) %>%
	slice_head(n = 8) %>%
	pull(geo_type1)
figS3_addr_type_levels <- figS3_geo_base %>%
	filter(addr_type != "Missing") %>%
	count(addr_type, sort = TRUE) %>%
	slice_head(n = 8) %>%
	pull(addr_type)
figS3_geo_type_levels <- unique(c(intersect(c("住宅区", "工作区"), figS3_geo_base$geo_type1), figS3_geo_type_levels))
figS3_addr_type_levels <- unique(c(intersect(c("住宅区", "工作区"), figS3_geo_base$addr_type), figS3_addr_type_levels))
figS3_geo_base2 <- figS3_geo_base %>%
	mutate(
		geo_type_plot = ifelse(geo_type1 %in% figS3_geo_type_levels, geo_type1, "Other"),
		addr_type_plot = ifelse(addr_type %in% figS3_addr_type_levels, addr_type, "Other")
	)
figS3_geo_trend <- figS3_geo_base2 %>%
	count(year, geo_type_plot, name = "n") %>%
	group_by(year) %>%
	mutate(total = sum(n), pct = n / total) %>%
	ungroup() %>%
	mutate(geo_type_plot = factor(geo_type_plot, levels = unique(c(figS3_geo_type_levels, "Other"))))
pS3GeoA <- ggplot(figS3_geo_trend, aes(year, pct, color = geo_type_plot, group = geo_type_plot)) +
	geom_line(linewidth = .78) +
	geom_point(size = 1.6) +
	scale_x_continuous(breaks = years) +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
	labs(title = "a. geo.type1 composition by year", x = NULL, y = "Percentage", color = NULL) +
	fig_theme(base_size = 8.6) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", legend.text = element_text(size = 7.2), plot.title = element_text(size = 10.2, face = "bold", hjust = .5))
figS3_geo_agree <- figS3_geo_base %>% filter(addr_type != "Missing", geo_type1 != "Missing") %>% summarise(n = n(), exact_agreement = mean(addr_type == geo_type1), .groups = "drop")
figS3_geo_heat <- figS3_geo_base2 %>%
	filter(addr_type_plot != "Missing", geo_type_plot != "Missing") %>%
	count(geo_type_plot, addr_type_plot, name = "n") %>%
	complete(geo_type_plot = unique(c(figS3_geo_type_levels, "Other")), addr_type_plot = unique(c(figS3_addr_type_levels, "Other")), fill = list(n = 0)) %>%
	group_by(geo_type_plot) %>%
	mutate(row_total = sum(n), pct = ifelse(row_total > 0, n / row_total, NA_real_), is_diag = as.character(geo_type_plot) == as.character(addr_type_plot), label = ifelse(n > 0, sprintf("%s
%.1f%%", scales::comma(n), 100 * pct), "")) %>%
	ungroup() %>%
	mutate(
		addr_type_plot = factor(addr_type_plot, levels = unique(c(figS3_addr_type_levels, "Other"))),
		geo_type_plot = factor(geo_type_plot, levels = rev(unique(c(figS3_geo_type_levels, "Other"))))
	)
figS3_geo_heat_title <- sprintf("b. Recorded address type vs predicted geographic type (%.1f%% exact)", 100 * figS3_geo_agree$exact_agreement[1])
pS3GeoB <- ggplot(figS3_geo_heat, aes(addr_type_plot, geo_type_plot, fill = pct)) +
	geom_tile(color = "white", linewidth = .25) +
	geom_text(data = dplyr::filter(figS3_geo_heat, !is_diag), aes(label = label), size = 2.15) +
	geom_text(data = dplyr::filter(figS3_geo_heat, is_diag), aes(label = label), size = 2.15, fontface = "bold") +
	scale_fill_gradient(low = "grey95", high = "#2166AC", labels = scales::percent_format(accuracy = 1), name = NULL, na.value = "grey98") +
	labs(title = figS3_geo_heat_title, x = "地址类型", y = "geo.type1") +
	fig_theme(base_size = 8.0) +
	theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(size = 10.2, face = "bold", hjust = .5), legend.position = "right")
geo_word_terms <- c(
	"小区", "花园", "家园", "公寓", "宿舍", "社区", "住宅", "村", "新村", "城中村",
	"公司", "工厂", "工业园", "科技园", "产业园", "写字楼", "大厦", "办公", "市场", "商场",
	"医院", "学校", "酒店", "宾馆", "地铁", "车站", "路口", "公园", "广场", "机场", "港口", "餐厅"
)
figS3_cloud_levels <- unique(c(intersect(c("住宅区", "工作区"), figS3_geo_base$geo_type1), figS3_geo_type_levels))
figS3_cloud_levels <- figS3_cloud_levels[figS3_cloud_levels != "Missing"]
figS3_cloud_levels <- head(figS3_cloud_levels, 6)
figS3_geo_terms <- purrr::map_dfr(figS3_cloud_levels, function(g) {
	dg <- figS3_geo_base %>% filter(geo_type1 == g)
	term_counts <- tibble(term = geo_word_terms, n = vapply(geo_word_terms, function(tt) sum(stringr::str_detect(dg$address, stringr::fixed(tt)), na.rm = TRUE), integer(1))) %>% filter(n > 0)
	if (nrow(term_counts) < 8) {
		chunks <- stringr::str_extract_all(paste(dg$address, collapse = " "), "[\\p{Han}]{2,6}")[[1]]
		fallback <- tibble(term = chunks) %>%
			mutate(term = stringr::str_replace_all(term, "深圳市|广东省|深圳|广东|宝安区|南山区|福田区|罗湖区|龙岗区|龙华区|坪山区|光明区|盐田区", "")) %>%
			filter(nchar(term) >= 2, !term %in% c("附近", "门口", "对面", "旁边", "地址", "未知", "深圳市")) %>%
			count(term, name = "n", sort = TRUE) %>%
			slice_head(n = 18)
		term_counts <- bind_rows(term_counts, fallback) %>% group_by(term) %>% summarise(n = sum(n), .groups = "drop")
	}
	term_counts %>% arrange(desc(n), term) %>% slice_head(n = 16) %>% mutate(geo_type1 = g)
}) %>%
	group_by(geo_type1) %>%
	mutate(
		rank = row_number(),
		angle = rank * pi * (3 - sqrt(5)),
		r = sqrt(rank) / sqrt(max(rank)),
		x = 0.50 + 0.48 * r * cos(angle),
		y = 0.50 + 0.40 * r * sin(angle),
		size2 = scales::rescale(sqrt(n), to = c(3.2, 9.2)),
		geo_label = factor(geo_type1, levels = figS3_cloud_levels)
	) %>%
	ungroup()
figS3_geo_label_dat <- figS3_geo_terms %>% distinct(geo_type1, geo_label) %>% mutate(x = .5, y = -.08)
pS3GeoC <- ggplot(figS3_geo_terms, aes(x, y, label = term, size = size2, color = geo_label)) +
	geom_text(fontface = "bold", alpha = .88, check_overlap = TRUE) +
	geom_text(data = figS3_geo_label_dat, aes(x = x, y = y, label = geo_label), inherit.aes = FALSE, fontface = "bold", size = 3.1, color = "black") +
	facet_wrap(~geo_label, ncol = 3) +
	scale_size_identity() +
	coord_cartesian(xlim = c(0, 1), ylim = c(-0.14, 1), expand = FALSE, clip = "off") +
	labs(title = "c. Address-keyword word clouds by geo.type1") +
	theme_void(base_size = 8.5) +
	theme(strip.text = element_blank(), legend.position = "none", plot.title = element_text(size = 10.2, face = "bold", hjust = .5), panel.spacing = grid::unit(.6, "lines"), plot.margin = margin(5, 0, 5, 0))
cv_summary_file <- file.path(dir_dx_cv, "cv_summary.csv")
cv_perclass_file <- file.path(dir_dx_cv, "cv_per_class_grouped.csv")
cv_confusion_file <- file.path(dir_dx_cv, "cv_confusion_grouped.csv")
cv_fold_file <- file.path(dir_dx_cv, "cv_fold_metrics.csv")
cv_bootstrap_file <- file.path(dir_dx_cv, "cv_bootstrap_metrics.csv")
cv_calibration_file <- file.path(dir_dx_cv, "cv_confidence_calibration.csv")
cv_sheets <- list()
# Keep FigS3 exactly as in the OLD pipeline.  CV diagnostics are a separate
# companion figure and never change the legacy figure's layout or workbook.
FigS3 <- (pS3GeoA | pS3GeoB) / pS3GeoC + plot_layout(heights = c(1, 1.85), widths = c(1.05, 1.15))

if (all(file.exists(c(cv_summary_file, cv_perclass_file, cv_confusion_file, cv_fold_file)))) {
	cv_summary <- read.csv(cv_summary_file, check.names = FALSE, stringsAsFactors = FALSE)
	cv_perclass <- read.csv(cv_perclass_file, check.names = FALSE, stringsAsFactors = FALSE)
	cv_folds <- read.csv(cv_fold_file, check.names = FALSE, stringsAsFactors = FALSE)
	cv_bootstrap <- if (file.exists(cv_bootstrap_file)) read.csv(cv_bootstrap_file, check.names = FALSE, stringsAsFactors = FALSE) else tibble()
	cv_calibration <- if (file.exists(cv_calibration_file)) read.csv(cv_calibration_file, check.names = FALSE, stringsAsFactors = FALSE) else tibble()
	cv_cm0 <- read.csv(cv_confusion_file, check.names = FALSE, stringsAsFactors = FALSE)
	names(cv_cm0)[1] <- "true_group"
	cv_cm <- cv_cm0 %>%
		pivot_longer(-true_group, names_to = "pred_group", values_to = "n") %>%
		group_by(true_group) %>%
		mutate(total = sum(n), row_pct = ifelse(total > 0, n / total, NA_real_)) %>%
		ungroup() %>%
		mutate(
			true_group = factor(true_group, levels = rev(dxs.all0)),
			pred_group = factor(pred_group, levels = dxs.all0)
		)

	cv_metric_map <- c(accuracy = "Accuracy", macro_f1 = "Macro-F1", weighted_f1 = "Weighted-F1")
	cv_fold_long <- cv_folds %>%
		filter(level == "grouped", system %in% c("Keyword", "MacBERT")) %>%
		select(any_of(c("fold", "system", names(cv_metric_map)))) %>%
		pivot_longer(cols = all_of(names(cv_metric_map)), names_to = "metric", values_to = "value") %>%
		mutate(metric = factor(cv_metric_map[metric], levels = unname(cv_metric_map)), system = factor(system, levels = c("Keyword", "MacBERT")))

	pCV_A <- ggplot(cv_fold_long, aes(metric, value, color = system)) +
		geom_boxplot(aes(group = interaction(metric, system)), position = position_dodge(width = .62), width = .46, outlier.shape = NA, alpha = .12) +
		geom_point(position = position_jitterdodge(jitter.width = .08, dodge.width = .62), size = 1.8, alpha = .85) +
		scale_color_manual(values = c(Keyword = "grey50", MacBERT = "#2166AC")) +
		scale_y_continuous(limits = c(.60, .94), breaks = seq(.6, .9, .05), expand = expansion(mult = c(.02, .04))) +
		labs(title = "a. Fold-to-fold stability", subtitle = "Five held-out folds for each grouped metric", x = NULL, y = "Performance", color = NULL) +
		fig2_pub_theme(base_size = 8.3) +
		theme(legend.position = "top")

	cv_pc <- cv_perclass %>%
		filter(level == "grouped", system == "MacBERT") %>%
		mutate(class = as.character(class))
	pCV_B <- ggplot(cv_pc, aes(recall, precision, size = support, fill = f1, label = class)) +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey68", linewidth = .55) +
		geom_point(shape = 21, color = "grey25", alpha = .90, stroke = .35) +
		geom_text(nudge_y = .015, size = 2.0, check_overlap = TRUE) +
		scale_fill_gradient(low = "grey90", high = "#2166AC", limits = c(0.5, 1), name = "F1") +
		scale_size_continuous(range = c(2.2, 6.0), labels = scales::comma, name = "Expert N") +
		coord_fixed(xlim = c(.50, 1.01), ylim = c(.50, 1.01)) +
		labs(title = "b. Per-phenotype precision-recall profile", subtitle = "Bubble size reflects the expert-labelled sample size", x = "Recall", y = "Precision") +
		fig2_pub_theme(base_size = 7.8) +
		theme(legend.position = "right")

	cv_top_errors <- cv_cm %>%
		filter(as.character(true_group) != as.character(pred_group), is.finite(row_pct), row_pct > 0) %>%
		mutate(path = sprintf("%s → %s", as.character(true_group), as.character(pred_group))) %>%
		arrange(desc(row_pct), desc(n)) %>%
		slice_head(n = 15) %>%
		mutate(path = factor(path, levels = rev(path)))
	pCV_C <- ggplot(cv_top_errors, aes(row_pct, path)) +
		geom_col(fill = "#8DB3E2", width = .68) +
		geom_text(aes(label = sprintf("%.1f%%", 100 * row_pct)), hjust = -0.15, size = 2.2, fontface = "bold", color = "grey20") +
		scale_x_continuous(limits = c(0, max(.08, max(cv_top_errors$row_pct, na.rm = TRUE) * 1.18)), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, .02))) +
		labs(title = "c. Largest off-diagonal pathways", subtitle = "Most residual errors are concentrated in a small number of phenotype pairs", x = "Row-normalised error share", y = NULL) +
		fig2_pub_theme(base_size = 7.8) +
		theme(axis.text.y = element_text(size = 6.6))

	cv_utility <- if (exists("fig2nlp_utility")) fig2nlp_utility else tibble()
	pCV_D1 <- if (nrow(cv_calibration)) {
		ggplot(cv_calibration, aes(mean_confidence, observed_accuracy)) +
			geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = .55) +
			geom_line(linewidth = .70, color = "grey35") +
			geom_point(aes(size = n), color = "#2166AC", alpha = .92) +
			coord_fixed(xlim = c(.50, 1.00), ylim = c(.50, 1.00), expand = FALSE) +
			labs(title = "d. Confidence calibration", x = "Mean confidence", y = "Observed accuracy", size = "N") +
			fig2_pub_theme(base_size = 7.6) +
			theme(legend.position = "right")
	} else empty_panel("d. Confidence calibration", "Calibration output unavailable")
	pCV_D2 <- if (nrow(cv_utility)) {
		cv_utility_long <- cv_utility %>%
			select(threshold, Coverage = coverage, Accuracy = accuracy, `Macro-F1` = macro_f1) %>%
			pivot_longer(-threshold, names_to = "metric", values_to = "value")
		ggplot(cv_utility_long, aes(threshold, value, color = metric, linetype = metric)) +
			geom_vline(xintercept = .80, linetype = "dotted", color = "grey42", linewidth = .55) +
			geom_line(linewidth = .80) +
			scale_color_manual(values = c(Coverage = "grey45", Accuracy = "#2166AC", `Macro-F1` = "#008B8B")) +
			scale_linetype_manual(values = c(Coverage = "dashed", Accuracy = "solid", `Macro-F1` = "solid")) +
			scale_y_continuous(limits = c(.55, 1.02), labels = scales::percent_format(accuracy = 1)) +
			labs(title = "e. Confidence-threshold utility", x = "Minimum confidence", y = "Retained fraction / performance", color = NULL, linetype = NULL) +
			fig2_pub_theme(base_size = 7.4) +
			theme(legend.position = "top", legend.text = element_text(size = 6.6))
	} else empty_panel("e. Confidence-threshold utility", "Utility output unavailable")
	cv_sheets <- list(
		CV_summary = cv_summary,
		CV_bootstrap_95CI = cv_bootstrap,
		CV_fold_metrics = cv_folds,
		CV_per_class_grouped = cv_perclass,
		CV_confusion_grouped = cv_cm0,
		CV_calibration = cv_calibration,
		CV_top_errors = cv_top_errors
	)
} else {
	cat("Cross-validation diagnostic panels are unavailable because the required CV files are missing.\n")
}

save_plot(FigS3, "FigS3.png", width = 11.2, height = 10.8, dpi = 600)
writexl::write_xlsx(list(
	geo_concordance_summary = figS3_geo_agree,
	geo_type1_yearly_trend = figS3_geo_trend,
	address_type_geo_type1_heatmap = figS3_geo_heat,
	geo_type1_word_cloud_terms = figS3_geo_terms,
	configuration = tibble(note = "FigS3 compares raw 地址类型 with ML-derived geo.type1. Low exact agreement means the two fields should not be treated as interchangeable without validation.")
), "FigS3.out.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS4a-b. Sensitivity: phone-score definitions and tail cutoffs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phone_score_labels <- c(
	phone.sco0 = "Pattern rule\n(sco0)",
	phone.sco1 = "Simple rule\n(sco1)",
	phone.sco2 = "Advanced rule\n(sco2)",
	phone.sco3 = "ML learning\n(sco3)"
)

calc_tail_rr <- function(dat, low_col = "low", high_col = "high") {
	if (!nrow(dat)) return(tibble())
	tg <- dat %>% count(score_name, contrast, year, disease, group, name = "n") %>% complete(score_name, contrast, year, disease, group = c(low_col, high_col), fill = list(n = 0))
	Ns <- dat %>% count(score_name, contrast, year, group, name = "N") %>% complete(score_name, contrast, year, group = c(low_col, high_col), fill = list(N = 0))
	tg %>%
		left_join(Ns, by = c("score_name", "contrast", "year", "group")) %>%
		pivot_wider(names_from = group, values_from = c(n, N), names_sep = "_") %>%
		mutate(
			RR = (n_high / pmax(N_high, 1)) / (n_low / pmax(N_low, 1)),
			se_logRR = sqrt(1 / pmax(n_high, 1) - 1 / pmax(N_high, 1) + 1 / pmax(n_low, 1) - 1 / pmax(N_low, 1)),
			RR_lo = exp(log(RR) - 1.96 * se_logRR),
			RR_hi = exp(log(RR) + 1.96 * se_logRR),
			p = purrr::pmap_dbl(list(n_high, N_high, n_low, N_low), function(a, A, b, B) suppressWarnings(tryCatch(chisq.test(matrix(c(a, A - a, b, B - b), nrow = 2))$p.value, error = function(e) NA_real_)))
		)
}

fit_binom_score <- function(d) {
	if (nrow(d) < 3 || length(unique(d$score_z[is.finite(d$score_z)])) < 2) return(tibble(OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_))
	fit <- tryCatch(glm(cbind(case, total - case) ~ score_z, family = binomial(), data = d), error = function(e) NULL)
	if (is.null(fit) || !"score_z" %in% rownames(coef(summary(fit)))) return(tibble(OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_))
	b <- coef(summary(fit))["score_z", ]
	tibble(OR = exp(b["Estimate"]), OR_lo = exp(b["Estimate"] - 1.96 * b["Std. Error"]), OR_hi = exp(b["Estimate"] + 1.96 * b["Std. Error"]), p = b["Pr(>|z|)"])
}

sens_score_base <- bind_rows(lapply(years, function(y) {
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())
	d0 %>%
		filter(!is.na(dx_grp), dx_grp %in% fig3_dxs_out) %>%
		transmute(year = as.integer(y), disease = factor(as.character(dx_grp), levels = fig3_dxs_out), phone.sco0, phone.sco1, phone.sco2, phone.sco3)
}))

sens_score_long <- sens_score_base %>%
	pivot_longer(cols = all_of(phone_score_vars), names_to = "score_name", values_to = "score") %>%
	mutate(score = suppressWarnings(as.numeric(score))) %>%
	filter(is.finite(score), !is.na(disease))

sens_tail_cutoffs <- sens_score_long %>%
	group_by(score_name) %>%
	summarise(q10 = quantile(score, .10, na.rm = TRUE), q25 = quantile(score, .25, na.rm = TRUE), q75 = quantile(score, .75, na.rm = TRUE), q90 = quantile(score, .90, na.rm = TRUE), .groups = "drop")

sens_tail_dat <- sens_score_long %>%
	left_join(sens_tail_cutoffs, by = "score_name") %>%
	mutate(
		`Top/bottom 10%` = case_when(score <= q10 ~ "low", score >= q90 ~ "high", TRUE ~ NA_character_),
		`Top/bottom quartile` = case_when(score <= q25 ~ "low", score >= q75 ~ "high", TRUE ~ NA_character_)
	) %>%
	pivot_longer(cols = c(`Top/bottom 10%`, `Top/bottom quartile`), names_to = "contrast", values_to = "group") %>%
	filter(group %in% c("low", "high")) %>%
	mutate(group = factor(group, levels = c("low", "high")), contrast = factor(contrast, levels = c("Top/bottom 10%", "Top/bottom quartile")))

sens_tail_rr <- calc_tail_rr(sens_tail_dat) %>%
	group_by(score_name, contrast, disease) %>%
	mutate(p_adj = p.adjust(p, "BH"), sig = !is.na(p_adj) & p_adj < .01) %>%
	ungroup()

meta_log_effect <- function(est, lo, hi) {
	ok <- is.finite(est) & is.finite(lo) & is.finite(hi) & est > 0 & lo > 0 & hi > 0
	if (sum(ok) < 1) return(tibble(meta = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
	y <- log(est[ok]); se <- (log(hi[ok]) - log(lo[ok])) / (2 * 1.96); ok2 <- is.finite(y) & is.finite(se) & se > 0
	if (!any(ok2)) return(tibble(meta = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
	w <- 1 / se[ok2]^2; b <- sum(w * y[ok2]) / sum(w); seb <- sqrt(1 / sum(w)); z <- b / seb
	tibble(meta = exp(b), lo = exp(b - 1.96 * seb), hi = exp(b + 1.96 * seb), p = 2 * stats::pnorm(abs(z), lower.tail = FALSE))
}

sens_tail_summary <- sens_tail_rr %>%
	group_by(score_name, contrast, disease) %>%
	group_modify(~ meta_log_effect(.x$RR, .x$RR_lo, .x$RR_hi) %>% mutate(n_sig_years = sum(.x$sig, na.rm = TRUE), n_years = sum(is.finite(.x$RR)))) %>%
	ungroup() %>%
	rename(mean_RR = meta, mean_RR_lo = lo, mean_RR_hi = hi, p_meta = p) %>%
	group_by(contrast) %>% mutate(p_adj_meta = p.adjust(p_meta, "BH"), sig05_meta = sig_star05(p_adj_meta)) %>% ungroup() %>%
	mutate(score_label = factor(phone_score_labels[score_name], levels = phone_score_labels), disease_label = factor(dx_to_eng(disease), levels = rev(dx_to_eng(dxs.all))), ci_label = ifelse(is.finite(mean_RR), sprintf("%.2f (%.2f–%.2f)%s", mean_RR, mean_RR_lo, mean_RR_hi, sig05_meta), NA_character_))

sens_score_tot <- sens_score_long %>% count(score_name, year, score, name = "total")
sens_score_moments <- sens_score_tot %>%
	group_by(score_name) %>%
	summarise(mu = weighted.mean(score, w = total, na.rm = TRUE), sd = sqrt(weighted.mean((score - weighted.mean(score, w = total, na.rm = TRUE))^2, w = total, na.rm = TRUE)), .groups = "drop")

sens_cont_agg <- sens_score_long %>%
	count(score_name, year, disease, score, name = "case") %>%
	complete(score_name, year, disease, score, fill = list(case = 0)) %>%
	left_join(sens_score_tot, by = c("score_name", "year", "score")) %>%
	left_join(sens_score_moments, by = "score_name") %>%
	mutate(score_z = ifelse(is.finite(sd) & sd > 0, (score - mu) / sd, NA_real_)) %>%
	filter(is.finite(score_z), is.finite(total), total > 0)

sens_cont_rr <- sens_cont_agg %>%
	group_by(score_name, year, disease) %>%
	group_modify(~ fit_binom_score(.x)) %>%
	ungroup() %>%
	group_by(score_name, disease) %>%
	mutate(p_adj = p.adjust(p, "BH"), sig = !is.na(p_adj) & p_adj < .01) %>%
	ungroup()

sens_cont_summary <- sens_cont_rr %>%
	group_by(score_name, disease) %>%
	group_modify(~ meta_log_effect(.x$OR, .x$OR_lo, .x$OR_hi) %>% mutate(n_sig_years = sum(.x$sig, na.rm = TRUE), n_years = sum(is.finite(.x$OR)))) %>%
	ungroup() %>%
	rename(mean_OR = meta, mean_OR_lo = lo, mean_OR_hi = hi, p_meta = p) %>%
	mutate(p_adj_meta = p.adjust(p_meta, "BH"), sig05_meta = sig_star05(p_adj_meta), score_label = factor(phone_score_labels[score_name], levels = phone_score_labels), disease_label = factor(dx_to_eng(disease), levels = rev(dx_to_eng(dxs.all))), ci_label = ifelse(is.finite(mean_OR), sprintf("%.2f (%.2f–%.2f)%s", mean_OR, mean_OR_lo, mean_OR_hi, sig05_meta), NA_character_))

pS3A <- sens_tail_summary %>%
	filter(disease %in% dxs.all) %>%
	ggplot(aes(x = mean_RR, y = disease_label, xmin = mean_RR_lo, xmax = mean_RR_hi, color = score_label)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
	geom_errorbar(position = position_dodge(width = .62), width = .18, orientation = "y", linewidth = .55) +
	geom_point(position = position_dodge(width = .62), size = 1.8) +
	geom_text(aes(label = sig05_meta, x = mean_RR_hi), position = position_dodge(width = .62), hjust = -.25, size = 3.6, fontface = "bold", show.legend = FALSE) +
	facet_wrap(~contrast, nrow = 1) +
	scale_x_log10(breaks = c(.85, .9, 1, 1.1, 1.2, 1.3), limits = c(.84, 1.34)) +
	labs(title = "a. Tail-cutoff robustness: pooled RR", x = "RR: high vs low phone-score group", y = NULL, color = NULL) +
	fig_theme(base_size = 9.0) +
	theme(strip.text = element_text(face = "bold"), legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 8), plot.title = element_text(size = 11, face = "bold"))

pS3B <- sens_cont_summary %>%
	filter(disease %in% dxs.all) %>%
	ggplot(aes(x = mean_OR, y = disease_label, xmin = mean_OR_lo, xmax = mean_OR_hi, color = score_label)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
	geom_errorbar(position = position_dodge(width = .62), width = .18, orientation = "y", linewidth = .55) +
	geom_point(position = position_dodge(width = .62), size = 1.8) +
	geom_text(aes(label = sig05_meta, x = mean_OR_hi), position = position_dodge(width = .62), hjust = -.25, size = 3.6, fontface = "bold", show.legend = FALSE) +
	scale_x_log10(breaks = c(.96, .98, 1, 1.02, 1.04, 1.06), limits = c(.955, 1.065)) +
	labs(title = "b. Continuous-score robustness: pooled per-SD OR", x = "OR per 1-SD higher phone score", y = NULL, color = NULL) +
	fig_theme(base_size = 9.0) +
	theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 8), plot.title = element_text(size = 11, face = "bold"))

FigS4_phone_robustness_panel <- pS3A | pS3B
figS4_phone_robustness_workbook <- list(
	phone_score_cutoffs = sens_tail_cutoffs,
	phone_tail_yearly = sens_tail_rr,
	phone_tail_summary = sens_tail_summary,
	phone_continuous_yearly = sens_cont_rr,
	phone_continuous_summary = sens_cont_summary
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS4c-d. Sensitivity: adjusted disease mix and permutation negative control
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sens_env_cols <- unique(fig2_num_vars)
sens_env_cols <- sens_env_cols[nzchar(sens_env_cols)]
sens_env_cols <- sens_env_cols[vapply(sens_env_cols, function(v) {
	any(vapply(dat1.list[as.character(years)], function(d0) !is.null(d0) && v %in% names(d0), logical(1)))
}, logical(1))]
if (!exists("get_col0", mode = "function")) {
	get_col0 <- function(d, v) {
		if (v %in% names(d)) return(d[[v]])
		rep(NA, nrow(d))
	}
}
if (!exists("make_qcat", mode = "function")) {
	make_qcat <- function(x, n = 5) {
		x <- suppressWarnings(as.numeric(x))
		out <- rep("Missing", length(x))
		ok <- is.finite(x)
		if (!any(ok)) return(factor(out))
		qs <- unique(as.numeric(stats::quantile(x[ok], probs = seq(0, 1, length.out = n + 1), na.rm = TRUE, type = 7)))
		if (length(qs) < 2) {
			out[ok] <- "Q1"
		} else {
			out[ok] <- paste0("Q", cut(x[ok], breaks = qs, include.lowest = TRUE, labels = FALSE))
		}
		factor(out)
	}
}
if (!exists("fit_binom_or", mode = "function")) {
	fit_binom_or <- function(d, adjusted = FALSE) {
		d <- d %>% filter(is.finite(case), is.finite(total), total > 0, case >= 0, case <= total)
		if (nrow(d) < 2 || n_distinct(d$phone_high) < 2) return(tibble(OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_, n_strata = nrow(d)))
		covars <- character(0)
		if (isTRUE(adjusted)) {
			covars <- setdiff(names(d), c("phone_high", "case", "total"))
			covars <- covars[vapply(d[covars], function(x) n_distinct(x, na.rm = FALSE) > 1, logical(1))]
		}
		form <- stats::as.formula(paste("cbind(case, total - case) ~", paste(c("phone_high", covars), collapse = " + ")))
		fit <- tryCatch(suppressWarnings(stats::glm(form, family = stats::binomial(), data = d)), error = function(e) NULL)
		if (is.null(fit)) return(tibble(OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_, n_strata = nrow(d)))
		co <- summary(fit)$coefficients
		if (!"phone_high" %in% rownames(co)) return(tibble(OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_, n_strata = nrow(d)))
		b <- co["phone_high", "Estimate"]; se <- co["phone_high", "Std. Error"]; p <- co["phone_high", "Pr(>|z|)"]
		tibble(OR = exp(b), OR_lo = exp(b - 1.96 * se), OR_hi = exp(b + 1.96 * se), p = p, n_strata = nrow(d))
	}
}
sens_adj_base <- bind_rows(lapply(years, function(y) {
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())
	env_df <- if (length(sens_env_cols)) as_tibble(setNames(lapply(sens_env_cols, function(v) suppressWarnings(as.numeric(get_col0(d0, v)))), paste0("env", seq_along(sens_env_cols)))) else tibble(.rows = nrow(d0))
	bind_cols(
		tibble(
			year = as.integer(y),
			date = as.Date(get_col0(d0, "日期")),
			disease = factor(as.character(get_col0(d0, "dx_grp")), levels = fig3_dxs_out),
			phone_group = as.character(get_col0(d0, "phone.luck")),
			lon = suppressWarnings(as.numeric(get_col0(d0, vars.basic.ems[11]))),
			lat = suppressWarnings(as.numeric(get_col0(d0, vars.basic.ems[12]))),
			geo_type = as.character(get_col0(d0, "geo.type1"))
		),
		env_df
	)
})) %>%
	filter(phone_group %in% c("low", "high"), !is.na(disease), disease %in% fig3_dxs_out) %>%
	mutate(
		phone_high = as.integer(phone_group == "high"),
		month = ifelse(is.na(date), 0L, lubridate::month(date)),
		dow = ifelse(is.na(date), "Missing", as.character(lubridate::wday(date, label = TRUE, week_start = 1))),
		geo_cell_raw = ifelse(is.finite(lon) & is.finite(lat), paste0(round(lon / 0.02), "_", round(lat / 0.02)), "Missing"),
		geo_type = ifelse(is.na(geo_type) | geo_type == "", "Missing", geo_type)
	)

if (nrow(sens_adj_base)) {
	top_cells <- sens_adj_base %>% count(geo_cell_raw, sort = TRUE) %>% slice_head(n = 150) %>% pull(geo_cell_raw)
	sens_adj_base <- sens_adj_base %>% mutate(geo_cell = ifelse(geo_cell_raw %in% top_cells, geo_cell_raw, "Sparse"))
	for (i in seq_along(sens_env_cols)) sens_adj_base[[paste0("env", i, "_q")]] <- make_qcat(sens_adj_base[[paste0("env", i)]])
}

sens_adj_res <- bind_rows(lapply(dxs.all, function(dx) {
	d <- sens_adj_base %>% mutate(case0 = as.integer(disease == dx))
	d_crude <- d %>% group_by(phone_high) %>% summarise(case = sum(case0), total = n(), .groups = "drop")
	d_adj <- d %>% group_by(phone_high, year, month, dow, geo_type, across(any_of(paste0("env", seq_along(sens_env_cols), "_q")))) %>% summarise(case = sum(case0), total = n(), .groups = "drop")
	bind_rows(
		fit_binom_or(d_crude, adjusted = FALSE) %>% mutate(model = "Crude"),
		fit_binom_or(d_adj, adjusted = TRUE) %>% mutate(model = "Calendar + geo/environment adjusted")
	) %>% mutate(disease = dx)
})) %>%
	group_by(model) %>% mutate(p_adj = p.adjust(p, "BH"), sig = !is.na(p_adj) & p_adj < .05, sig05 = sig_star05(p_adj), p_label = fmt_p(p_adj), ci_label = ifelse(is.finite(OR), sprintf("%.2f (%.2f–%.2f)", OR, OR_lo, OR_hi), "")) %>% ungroup() %>%
	mutate(disease_label = factor(dx_to_eng(disease), levels = rev(dx_to_eng(dxs.all))), model = factor(model, levels = c("Crude", "Calendar + geo/environment adjusted")))

figS4_or_rng <- range(c(sens_adj_res$OR_lo, sens_adj_res$OR_hi), na.rm = TRUE)
figS4_or_lim <- c(floor(figS4_or_rng[1] * 20) / 20, ceiling(figS4_or_rng[2] * 20) / 20)
figS4_or_lim <- range(c(figS4_or_lim, 1), na.rm = TRUE)
figS4_or_breaks <- c(.5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.3, 1.5, 2)
figS4_or_breaks <- figS4_or_breaks[figS4_or_breaks >= figS4_or_lim[1] & figS4_or_breaks <= figS4_or_lim[2]]

pS4 <- ggplot(sens_adj_res, aes(y = disease_label, x = OR, xmin = OR_lo, xmax = OR_hi, color = model)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
	geom_segment(aes(x = OR_lo, xend = OR_hi, yend = disease_label), position = position_dodge(width = .55), linewidth = .7) +
	geom_point(position = position_dodge(width = .55), size = 2.4) +
	geom_text(data = sens_adj_res %>% filter(sig), aes(label = sig05, x = OR_hi), position = position_dodge(width = .55), hjust = -.15, size = 4, show.legend = FALSE) +
	scale_x_log10(breaks = figS4_or_breaks, labels = function(x) sprintf("%g", x)) +
	coord_cartesian(xlim = figS4_or_lim) +
	labs(title = NULL, x = "Odds ratio for phenotype membership: high vs low phone-luck group", y = NULL, color = NULL) +
	fig_theme(base_size = 9.2) +
	theme(legend.position = "bottom")
figS4b_shift <- sens_adj_res %>%
	select(disease, model, OR) %>%
	pivot_wider(names_from = model, values_from = OR) %>%
	mutate(adj_vs_crude_log_ratio = log(`Calendar + geo/environment adjusted`) - log(Crude), disease_label = factor(dx_to_eng(disease), levels = rev(dx_to_eng(dxs.all))))
pS4b <- ggplot(figS4b_shift, aes(adj_vs_crude_log_ratio, disease_label)) +
	geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
	geom_segment(aes(x = 0, xend = adj_vs_crude_log_ratio, yend = disease_label), linewidth = .75, color = "grey65") +
	geom_point(size = 2.6, color = "black") +
	scale_x_continuous(labels = function(x) sprintf("%+.2f", x)) +
	labs(title = NULL, x = "log(OR adjusted / OR crude)", y = NULL) +
	fig_theme(base_size = 9.2)
# FigS4 is assembled after the permutation negative-control panel below; pS4 is retained as the left panel.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FigS4 right panel: negative-control permutation of phone-luck labels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sensitivity_n_perm <- as.integer(Sys.getenv("EMS120_SENS_N_PERM", unset = "1000"))
if (!is.finite(sensitivity_n_perm) || sensitivity_n_perm < 100) sensitivity_n_perm <- 1000
set.seed(12014)

perm_base <- fig3_dat %>%
	filter(disease %in% dxs.all) %>%
	transmute(year, disease, n_high = as.integer(n_high), n_low = as.integer(n_low), N_high = as.integer(N_high), N_low = as.integer(N_low), n_disease = as.integer(n_high + n_low)) %>%
	filter(n_disease > 0, N_high > 0, N_low > 0)

perm_obs <- perm_base %>%
	group_by(disease) %>%
	summarise(obs_logRR = weighted.mean(log(((n_high + 0.5) / (N_high + 1)) / ((n_low + 0.5) / (N_low + 1))), w = n_disease, na.rm = TRUE), obs_RR = exp(obs_logRR), .groups = "drop")

perm_null <- bind_rows(lapply(split(perm_base, perm_base$disease), function(d) {
	vals <- replicate(sensitivity_n_perm, {
		x <- stats::rhyper(nrow(d), d$N_high, d$N_low, d$n_disease)
		weighted.mean(log(((x + 0.5) / (d$N_high + 1)) / (((d$n_disease - x) + 0.5) / (d$N_low + 1))), w = d$n_disease, na.rm = TRUE)
	})
	tibble(disease = unique(d$disease), perm_id = seq_along(vals), null_logRR = as.numeric(vals), null_RR = exp(null_logRR))
}))

perm_summary <- perm_null %>%
	group_by(disease) %>%
	summarise(null_lo = quantile(null_logRR, .025, na.rm = TRUE), null_hi = quantile(null_logRR, .975, na.rm = TRUE), null_mean = mean(null_logRR, na.rm = TRUE), .groups = "drop") %>%
	left_join(perm_obs, by = "disease") %>%
	left_join(perm_null %>% group_by(disease) %>% summarise(p_perm_two_sided = (1 + sum(abs(null_logRR) >= abs(perm_obs$obs_logRR[match(disease[1], perm_obs$disease)]), na.rm = TRUE)) / (n() + 1), .groups = "drop"), by = "disease") %>%
	mutate(disease_label = factor(dx_to_eng(disease), levels = dx_to_eng(dxs.all)))

pS5 <- perm_null %>%
	mutate(disease_label = factor(dx_to_eng(disease), levels = dx_to_eng(dxs.all))) %>%
	ggplot(aes(null_logRR)) +
	geom_histogram(bins = 40, fill = "grey80", color = "white") +
	geom_vline(data = perm_summary, aes(xintercept = obs_logRR), color = "#B2182B", linewidth = .85, linetype = "dashed") +
	geom_text(data = perm_summary, aes(x = obs_logRR, y = Inf, label = paste0("Pperm=", fmt_p(p_perm_two_sided))), color = "#B2182B", vjust = 1.35, hjust = -.05, size = 2.7, fontface = "bold", inherit.aes = FALSE) +
	facet_wrap(~disease_label, ncol = 3, scales = "free_y", axes = "all_x") +
	scale_x_continuous(labels = function(x) sprintf("%.2f", exp(x))) +
	labs(title = NULL, x = "Weighted mean RR under permutation; red line = observed", y = "Permutation count") +
	fig_theme(base_size = 8.6) +
	theme(strip.text = element_text(face = "bold"))
pS5b <- perm_summary %>%
	mutate(disease_label = factor(dx_to_eng(disease), levels = rev(dx_to_eng(dxs.all)))) %>%
	ggplot(aes(y = disease_label)) +
	geom_violin(
		data = perm_null %>% mutate(disease_label = factor(dx_to_eng(disease), levels = rev(dx_to_eng(dxs.all)))),
		aes(x = null_logRR, y = disease_label),
		inherit.aes = FALSE,
		fill = "grey82",
		color = "grey55",
		linewidth = .35,
		scale = "width",
		trim = TRUE
	) +
	geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
	geom_point(aes(x = obs_logRR), size = 2.5, color = "#B2182B") +
	geom_text(aes(x = null_mean, label = paste0("Pperm=", fmt_p(p_perm_two_sided))), hjust = 0.5, size = 2.7, fontface = "bold", color = "#B2182B") +
	scale_x_continuous(labels = function(x) sprintf("%.2f", exp(x))) +
	labs(title = NULL, x = "Observed RR and permutation density", y = NULL) +
	fig_theme(base_size = 9.0) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigS4 <- (pS3A | pS3B) +
	plot_layout(widths = c(1.18, 1))
save_plot(FigS4, "FigS4.png", width = 11.8, height = 6.8, dpi = 600)
unlink(c(Sys.glob("FigS4b*"), Sys.glob("FigS5b*")), force = TRUE)
writexl::write_xlsx(c(figS4_phone_robustness_workbook, list(
	adjusted_model_results = sens_adj_res,
	adjusted_model_shift = figS4b_shift,
	adjusted_model_base_summary = sens_adj_base %>% summarise(n = n(), years = n_distinct(year), geo_cells = n_distinct(geo_cell), env_covariates = paste(sens_env_cols, collapse = "; ")),
	permutation_observed_vs_null = perm_summary,
	permutation_null = perm_null,
	permutation_input = perm_base,
	configuration = tibble(n_permutations = sensitivity_n_perm, permutation_unit = "within-year phone-luck labels", statistic = "disease-count-weighted mean log(RR) across years", note = "Current FigS4 keeps only panels a-b; adjustment/permutation diagnostics remain workbook-only.")
)), "FigS4.out.xlsx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS5. ml_dx keyword and transformer summary
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS5_death_terms <- c("无生命体征", "心跳呼吸停止", "呼吸心跳停止", "心跳停止", "呼吸停止", "已死", "死亡", "死者", "尸体", "心肺复苏", "猝死", "瞳孔散大", "无反应", "意识丧失")
figS5_death_regex <- paste(figS5_death_terms, collapse = "|")
figS5_base <- bind_rows(lapply(years, function(y) {
	d <- dat1.list[[as.character(y)]]; if (is.null(d) || !nrow(d)) return(tibble())
	txt_cols <- intersect(vars.dxs, names(d))
	raw_text <- if (length(txt_cols)) do.call(paste, c(lapply(d[txt_cols], function(z) replace_na(as.character(z), "")), sep = " ")) else rep("", nrow(d))
	d %>% mutate(.raw_text = raw_text) %>%
		transmute(year = as.integer(y), dx0_raw = trimws(as.character(dx.type0)), dx1_raw = trimws(as.character(dx.type1)), dx1_grp = factor(as.character(dx_grp), levels = dxs.all), reason0 = as.character(dx.type0.reason), raw_text = .raw_text)
}))
figS5_base <- figS5_base %>%
	left_join(map_grp %>% rename(dx0_raw = dx_raw, dx0_grp = dx_grp), by = "dx0_raw") %>%
	mutate(
		dx0_grp = factor(dx0_grp, levels = dxs.all),
		agree_grp = !is.na(dx0_grp) & !is.na(dx1_grp) & dx0_grp == dx1_grp,
		death_explicit_text = stringr::str_detect(raw_text, figS5_death_regex),
		death_explicit_reason = stringr::str_detect(replace_na(reason0, ""), figS5_death_regex),
		death_explicit_any = death_explicit_text | death_explicit_reason,
		kw_rule_family = case_when(
			is.na(reason0) | reason0 == "" ~ "No keyword reason",
			stringr::str_detect(reason0, figS5_death_regex) ~ "Death-specific terms",
			stringr::str_detect(reason0, "胸痛|心梗|心肌梗死|心悸|心律|心衰|中风|脑卒中|意识障碍|抽搐") ~ "CVD/neurologic terms",
			stringr::str_detect(reason0, "呼吸困难|呼吸|气促|气喘|窒息|缺氧|肺炎") ~ "Respiratory terms",
			stringr::str_detect(reason0, "中毒|酒精|醉酒|服药|农药|一氧化碳|吸毒") ~ "Intoxication terms",
			stringr::str_detect(reason0, "精神|自杀|自伤|跳楼|情绪") ~ "Psychiatric terms",
			TRUE ~ "Other keyword terms"
		)
	)

tmpS5 <- plot_dx_raw_trend(dat1.list, years, "dx.type1", dxs.type.list, level = "raw", show = "percent")
raw_col <- setNames(dxs.all.color[map_grp$dx_grp], map_grp$dx_raw)
pS5A <- tmpS5$plot +
	scale_color_manual(values = raw_col, labels = dx_raw_to_type_label, name = NULL) +
	scale_fill_manual(values = raw_col, labels = dx_raw_to_type_label, name = NULL) +
	labs(title = "a. Transformer raw disease-category trends", x = NULL) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.title = element_blank(), plot.title = element_text(size = 10, face = "bold"), plot.margin = margin(5, 1, 5, 5))

figS5_concordance_summary <- figS5_base %>% summarise(n = n(), compared = sum(!is.na(dx0_grp) & !is.na(dx1_grp)), overall_agreement = mean(agree_grp, na.rm = TRUE))
figS5_concordance_summary_vip <- figS5_base %>% filter(dx0_grp %in% dxs.vip, dx1_grp %in% dxs.vip) %>% summarise(n = n(), compared = sum(!is.na(dx0_grp) & !is.na(dx1_grp)), agreement = mean(agree_grp, na.rm = TRUE))
figS5B_title <- sprintf("b. Keyword vs transformer group concordance (%.1f%%)", 100 * figS5_concordance_summary_vip$agreement[1])

figS5B_dat <- figS5_base %>%
	filter(!is.na(dx0_grp), !is.na(dx1_grp), dx0_grp %in% dxs.vip, dx1_grp %in% dxs.vip) %>%
	count(dx0_grp, dx1_grp, name = "n") %>%
	complete(dx0_grp = dxs.vip, dx1_grp = dxs.vip, fill = list(n = 0)) %>%
	group_by(dx1_grp) %>%
	mutate(row_total = sum(n), pct = ifelse(row_total > 0, n / row_total, NA_real_), is_diag = as.character(dx0_grp) == as.character(dx1_grp), label = ifelse(n > 0, sprintf("%s\n%.1f%%", scales::comma(n), 100 * pct), "")) %>%
	ungroup() %>%
	mutate(dx0_label = factor(dx_to_eng(dx0_grp), levels = dx_to_eng(dxs.vip)), dx1_label = factor(dx_to_eng(dx1_grp), levels = dx_to_eng(dxs.vip)))
pS5B <- ggplot(figS5B_dat, aes(dx0_label, fct_rev(dx1_label), fill = pct)) +
	geom_tile(color = "white", linewidth = .3) +
	geom_text(data = dplyr::filter(figS5B_dat, !is_diag), aes(label = label), size = 2.25) +
	geom_text(data = dplyr::filter(figS5B_dat, is_diag), aes(label = label), size = 2.25, fontface = "bold") +
	scale_fill_gradient(low = "grey95", high = "#2166AC", labels = scales::percent_format(accuracy = 1), name = NULL, na.value = "grey98") +
	labs(title = figS5B_title, x = "Keyword classification", y = "Transformer classification") +
	fig_theme(base_size = 8.0) +
	theme(
		axis.text.x = element_text(angle = 35, hjust = 1),
		axis.title.x = element_text(margin = margin(t = 14)),
		axis.title.y = element_text(margin = margin(r = 14)),
		plot.title = element_text(size = 10, face = "bold"),
		plot.margin = margin(5, 8, 12, 8)
	)

set.seed(12011)
figS5C_terms <- figS5_base %>%
	filter(dx1_grp %in% dxs.all, !is.na(reason0), reason0 != "") %>%
	mutate(
		reason0 = stringr::str_remove(reason0, "^(?:RuleV[0-9]*|Rule|FieldKW):\\s*"),
		reason0 = stringr::str_remove(reason0, ";score=.*$")
	) %>%
	tidyr::separate_rows(reason0, sep = "[,;，； ]+") %>%
	mutate(
		term = stringr::str_trim(reason0),
		# New FieldKW reasons encode field=value or field~value.  Panel c is a
		# keyword cloud, so retain only the value just as in the OLD output.
		term = stringr::str_remove(term, "^[^=~]+[=~]")
	) %>%
	filter(nchar(term) >= 2, !stringr::str_detect(term, "score="), !term %in% c("Rule:", "Rule", "RuleV3", "RuleV4", "FieldKW", "NoMatch", "EmptyText", "NoKeywordModel")) %>%
	count(dx1_grp, term, name = "n") %>%
	group_by(dx1_grp) %>%
	slice_max(n, n = 16, with_ties = FALSE) %>%
	arrange(dx1_grp, desc(n), term) %>%
	mutate(
		rank = row_number(),
		angle = rank * pi * (3 - sqrt(5)),
		r = sqrt(rank) / sqrt(max(rank)),
		x = 0.50 + 0.48 * r * cos(angle),
		y = 0.50 + 0.41 * r * sin(angle),
		size2 = scales::rescale(sqrt(n), to = c(3.2, 9.3)),
		dx1_label = factor(dx_to_eng(dx1_grp), levels = dx_to_eng(dxs.all))
	) %>%
	ungroup()
figS5C_labels <- figS5C_terms %>% distinct(dx1_grp, dx1_label) %>% mutate(x = 0.5, y = -0.08)
pS5C <- ggplot(figS5C_terms, aes(x, y, label = term, size = size2, color = dx1_grp)) +
	geom_text(fontface = "bold", alpha = .88, check_overlap = TRUE) +
	geom_text(data = figS5C_labels, aes(x = x, y = y, label = dx1_label), inherit.aes = FALSE, fontface = "bold", size = 3.0, color = "black") +
	facet_wrap(~dx1_label, ncol = 3) +
	scale_color_manual(values = dxs.all.color[dxs.all], guide = "none", drop = FALSE) +
	scale_size_identity() +
	coord_cartesian(xlim = c(0, 1), ylim = c(-0.14, 1), expand = FALSE, clip = "off") +
	labs(title = "c. Keyword-token word clouds by EMS phenotype") +
	theme_void(base_size = 8.6) +
	theme(strip.text = element_blank(), strip.background = element_blank(), plot.title = element_text(size = 10, face = "bold"), panel.spacing = grid::unit(.6, "lines"), plot.margin = margin(5, 0, 5, 0))

# FigS6. Death-specific validation
figS6A_dat <- figS5_base %>%
	filter(dx1_grp == "Death", !is.na(dx0_grp), dx0_grp %in% dxs.all) %>%
	count(year, dx0_grp, name = "n") %>%
	complete(year = years, dx0_grp = dxs.all, fill = list(n = 0)) %>%
	group_by(year) %>% mutate(total = sum(n), pct = ifelse(total > 0, n / total, NA_real_)) %>% ungroup() %>%
	mutate(dx0_grp = factor(dx0_grp, levels = dxs.all))
pS6A <- ggplot(figS6A_dat, aes(year, pct, fill = dx0_grp)) +
	geom_col(width = .82) +
	scale_fill_manual(values = dxs.all.color[dxs.all], labels = dx_to_eng, name = "Keyword group") +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
	scale_x_continuous(breaks = years) +
	labs(title = "a. Death-call keyword mix", x = NULL, y = "% of transformer-Death") +
	fig_theme(base_size = 8.0) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", legend.text = element_text(size = 6.5), legend.key.size = grid::unit(.28, "cm"), plot.title = element_text(size = 10, face = "bold"))

figS6B_dat <- figS5_base %>%
	filter(dx1_grp == "Death") %>%
	group_by(year) %>%
	summarise(n = n(), keyword_confirmed_death = mean(dx0_grp == "Death", na.rm = TRUE), explicit_death_terms_in_text = mean(death_explicit_text, na.rm = TRUE), explicit_death_terms_anywhere = mean(death_explicit_any, na.rm = TRUE), .groups = "drop") %>%
	pivot_longer(cols = c(keyword_confirmed_death, explicit_death_terms_in_text, explicit_death_terms_anywhere), names_to = "metric", values_to = "pct") %>%
	mutate(metric = factor(metric, levels = c("keyword_confirmed_death", "explicit_death_terms_in_text", "explicit_death_terms_anywhere"), labels = c("Keyword confirms Death", "Death terms in raw text", "Death terms in text/reason")))
figS6B_duplicate_check <- figS6B_dat %>%
	select(year, metric, pct) %>%
	pivot_wider(names_from = metric, values_from = pct) %>%
	summarise(
		raw_text_equals_text_or_reason = all(abs(`Death terms in raw text` - `Death terms in text/reason`) < 1e-12, na.rm = TRUE),
		.groups = "drop"
	)
figS6B_plot_dat <- figS6B_dat
if (isTRUE(figS6B_duplicate_check$raw_text_equals_text_or_reason[1])) {
	figS6B_plot_dat <- figS6B_plot_dat %>% filter(metric != "Death terms in text/reason")
}
pS6B <- ggplot(figS6B_plot_dat, aes(year, pct, color = metric, group = metric)) +
	geom_hline(yintercept = .50, linetype = "dotted", color = "grey55") +
	geom_line(linewidth = .9) +
	geom_point(size = 1.8) +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
	scale_x_continuous(breaks = years) +
	labs(title = "b. Death-term support", x = NULL, y = "% of transformer-Death", color = NULL) +
	fig_theme(base_size = 8.0) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 7), plot.title = element_text(size = 10, face = "bold"))

figS5_yearly_agreement <- figS5_base %>% group_by(year) %>% summarise(n = n(), compared = sum(!is.na(dx0_grp) & !is.na(dx1_grp)), agreement = mean(agree_grp, na.rm = TRUE), death_transformer_n = sum(dx1_grp == "Death", na.rm = TRUE), death_keyword_confirmed_pct = mean(dx0_grp[dx1_grp == "Death"] == "Death", na.rm = TRUE), death_explicit_text_pct = mean(death_explicit_text[dx1_grp == "Death"], na.rm = TRUE), .groups = "drop")
figS6_death_rule_family <- figS5_base %>% filter(dx1_grp == "Death") %>% count(dx0_grp, kw_rule_family, name = "n") %>% group_by(dx0_grp) %>% mutate(pct = n / sum(n)) %>% ungroup()
figS6C_dat <- figS6_death_rule_family %>%
	filter(!is.na(dx0_grp)) %>%
	mutate(dx0_label = factor(dx_to_eng(dx0_grp), levels = rev(dx_to_eng(dxs.all))))
pS6C <- ggplot(figS6C_dat, aes(kw_rule_family, dx0_label, fill = pct)) +
	geom_tile(color = "white", linewidth = .25) +
	geom_text(aes(label = ifelse(n > 0, sprintf("%s (%.1f%%)", scales::comma(n), 100 * pct), "")), size = 2.15, fontface = "bold") +
	scale_fill_gradientn(colors = c("#F7FCF0", "#C7E9B4", "#41B6C4", "#225EA8"), labels = scales::percent_format(accuracy = 1), name = NULL) +
	labs(title = "c. Keyword-rule families", x = NULL, y = NULL) +
	fig_theme(base_size = 8.0) +
	theme(axis.text.x = element_text(angle = 25, hjust = 1), plot.title = element_text(size = 10, face = "bold"), legend.position = "top")
figS5_top <- cowplot::plot_grid(pS5A, NULL, pS5B, nrow = 1, rel_widths = c(1.2, .08, 1), align = "h", axis = "tb")
figS5_death_panel <- align_panel_rows(
	rows = list(
		list(pS6A, pS6B),
		list(pS6C)
	),
	rel_widths = list(c(1.2, 1), 1),
	rel_heights = c(1, 1.18),
	row_gap = .08,
	side_pad = .018,
	axis_text_size = 8,
	axis_title_size = 8
)
pS5C <- pS5C + theme(plot.title = element_text(size = 10, face = "bold", hjust = .5))
FigS5 <- cowplot::plot_grid(figS5_top, NULL, pS5C, ncol = 1, rel_heights = c(1, .08, 2.20), align = "none")
save_plot(FigS5, "FigS5.png", width = 12.4, height = 13.0, dpi = 600)
save_plot(figS5_death_panel, "FigS6.png", width = 9.2, height = 9.2, dpi = 600)
unlink("FigS5_death_validation.png", force = TRUE)
writexl::write_xlsx(list(
	raw_category_trends = tmpS5$data,
	concordance_summary_all = figS5_concordance_summary,
	concordance_summary_vip = figS5_concordance_summary_vip,
	keyword_transformer_concordance_vip = figS5B_dat,
	word_cloud_tokens_all_dxs = figS5C_terms,
	configuration = tibble(note = "FigS5 contains the ml_dx keyword and transformer summary panels a-c; death-specific validation is FigS6.")
), "FigS5.out.xlsx")
writexl::write_xlsx(list(
	yearly_agreement = figS5_yearly_agreement,
	transformer_death_by_keyword_group = figS6A_dat,
	transformer_death_explicit_terms = figS6B_dat,
	transformer_death_explicit_plot = figS6B_plot_dat,
	transformer_death_duplicate_check = figS6B_duplicate_check,
	death_rule_family = figS6_death_rule_family,
	death_rule_family_plot = figS6C_dat,
	configuration = tibble(note = "FigS6 is the death-specific validation formerly written as FigS5_death_validation.png.")
), "FigS6.out.xlsx")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS7. Monthly composition and six-year hourly distribution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS7_month <- bind_rows(lapply(yrs_trend, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(!is.na(dx_grp), !is.na(.data[["日期"]]), dx_grp %in% dxs.all) %>%
		transmute(year = y, date = as.Date(.data[["日期"]]), dx = factor(as.character(dx_grp), levels = dxs.all)) %>%
		mutate(month_start = floor_date(date, "month")) %>%
		count(year, month_start, dx, name = "n") %>%
		group_by(year, month_start) %>%
		mutate(total_calls_month = sum(n), pct = n / total_calls_month) %>%
		ungroup()
}))
figS7_hourly_all <- bind_rows(lapply(years, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(dx_grp %in% dxs.all, is.finite(hour), between(as.integer(hour), 0, 23)) %>%
		transmute(year = as.integer(y), dx_grp = factor(as.character(dx_grp), levels = dxs.all), hour = as.integer(hour)) %>%
		count(year, dx_grp, hour, name = "call_n") %>%
		complete(year = as.integer(y), dx_grp = dxs.all, hour = 0:23, fill = list(call_n = 0))
}))
figS7_years <- 2018:2023
figS7_hourly <- figS7_hourly_all %>%
	filter(year %in% figS7_years) %>%
	group_by(dx_grp, year, hour) %>%
	summarise(call_n = sum(call_n, na.rm = TRUE), .groups = "drop") %>%
	group_by(dx_grp, year) %>%
	mutate(total_n = sum(call_n, na.rm = TRUE), pct = ifelse(total_n > 0, call_n / total_n, NA_real_)) %>%
	ungroup() %>%
	mutate(year = factor(year, levels = figS7_years))
figS7_year_cols <- setNames(rainbow(length(figS7_years)), as.character(figS7_years))
pFigS7A <- ggplot(figS7_month, aes(month_start, fct_rev(dx), fill = pct)) +
	geom_tile() +
	scale_fill_viridis_c(name = NULL, labels = scales::percent_format(accuracy = 1), guide = guide_colorbar(title = NULL, barwidth = grid::unit(12, "cm"), barheight = grid::unit(.45, "cm"))) +
	scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
	labs(title = "a. Monthly composition of the EMS spectrum", x = NULL, y = NULL) +
	fig_theme(base_size = 8.8) +
	theme(legend.position = "bottom", legend.direction = "horizontal", plot.title = element_text(size = 10.5, face = "bold"))
pFigS7B <- ggplot(figS7_hourly, aes(hour, pct, color = year, group = year)) +
	geom_line(linewidth = .85) +
	geom_point(size = 1.15) +
	facet_wrap(~dx_grp, scales = "free_y", ncol = 3, labeller = as_labeller(setNames(dx_to_eng(dxs.all), dxs.all))) +
	scale_color_manual(values = figS7_year_cols, name = NULL, guide = guide_legend(nrow = 1, byrow = TRUE)) +
	scale_x_continuous(breaks = seq(0, 24, 4), labels = sprintf("%02d", seq(0, 24, 4))) +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
	labs(title = "b. Hourly distribution of EMS phenotypes during six COVID-era years", x = "Dispatch hour", y = NULL) +
	fig_theme(base_size = 8.2) +
	theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 7), legend.margin = margin(t = 0, r = 0, b = 2, l = 0), plot.title = element_text(size = 10, face = "bold"), axis.title = element_text(size = 8.2, face = "bold"), strip.text = element_text(size = 8.2, face = "bold"))
FigS7 <- pFigS7A / pFigS7B + plot_layout(heights = c(1.0, 1.65))
save_plot(FigS7, "FigS7.png", width = 10.5, height = 11.0, dpi = 600)
writexl::write_xlsx(list(
	monthly_composition = figS7_month,
	hourly_six_years = figS7_hourly,
	configuration = tibble(years = paste(figS7_years, collapse = ", "), n_phenotypes = dplyr::n_distinct(figS7_hourly$dx_grp), panel = "FigS7", moved_panel = "Panel a monthly composition plus panel b hourly distribution")
), "FigS7.out.xlsx")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS8. Quantitative check of high-vs-low circadian differences
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS8_base <- bind_rows(lapply(years, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(dx_grp %in% dxs.all, phone.luck %in% c("low", "high"), is.finite(hour), between(as.integer(hour), 0, 23)) %>%
		transmute(
			year = as.integer(y),
			dx_grp = factor(as.character(dx_grp), levels = dxs.all),
			phone.luck = factor(as.character(phone.luck), levels = c("low", "high")),
			hour = as.integer(hour)
		)
}))

figS8_hour <- figS8_base %>%
	count(year, dx_grp, phone.luck, hour, name = "call_n") %>%
	complete(year = years, dx_grp = dxs.all, phone.luck = c("low", "high"), hour = 0:23, fill = list(call_n = 0)) %>%
	group_by(year, dx_grp, phone.luck) %>%
	mutate(total_n = sum(call_n), pct = ifelse(total_n > 0, call_n / total_n, NA_real_)) %>%
	ungroup()

figS8_dist_stats <- function(high_n, low_n) {
	high_n <- as.numeric(high_n); low_n <- as.numeric(low_n)
	Nh <- sum(high_n, na.rm = TRUE); Nl <- sum(low_n, na.rm = TRUE)
	if (!is.finite(Nh) || !is.finite(Nl) || Nh <= 0 || Nl <= 0) {
		return(tibble(total_high = Nh, total_low = Nl, total_variation = NA_real_, max_abs_diff_pp = NA_real_, peak_hour_high = NA_integer_, peak_hour_low = NA_integer_, chi_square = NA_real_, p_chisq = NA_real_, cramer_v = NA_real_))
	}
	ph <- high_n / Nh; pl <- low_n / Nl
	mat <- rbind(high = high_n, low = low_n)
	chi <- suppressWarnings(tryCatch(chisq.test(mat), error = function(e) NULL))
	chi_stat <- if (is.null(chi)) NA_real_ else unname(chi$statistic)
	p <- if (is.null(chi)) NA_real_ else unname(chi$p.value)
	N <- sum(mat)
	v <- ifelse(is.finite(chi_stat) && N > 0, sqrt(chi_stat / N), NA_real_)
	tibble(
		total_high = Nh,
		total_low = Nl,
		total_variation = 0.5 * sum(abs(ph - pl), na.rm = TRUE),
		max_abs_diff_pp = 100 * max(abs(ph - pl), na.rm = TRUE),
		peak_hour_high = which.max(ph) - 1L,
		peak_hour_low = which.max(pl) - 1L,
		chi_square = chi_stat,
		p_chisq = p,
		cramer_v = v
	)
}

figS8_tmp <- figS8_hour %>%
	arrange(year, dx_grp, phone.luck, hour) %>%
	group_by(year, dx_grp) %>%
	summarise(
		high_n = list(call_n[phone.luck == "high"][order(hour[phone.luck == "high"])]),
		low_n = list(call_n[phone.luck == "low"][order(hour[phone.luck == "low"])]),
		.groups = "drop"
	)

figS8_stats <- bind_cols(
	figS8_tmp %>% select(year, dx_grp),
	purrr::map2_dfr(figS8_tmp$high_n, figS8_tmp$low_n, figS8_dist_stats)
) %>%
	group_by(year) %>%
	mutate(p_adj = p.adjust(p_chisq, method = "BH"), sig05 = sig_star05(p_adj)) %>%
	ungroup() %>%
	mutate(
		dx_grp = factor(as.character(dx_grp), levels = rev(dxs.all)),
		label = ifelse(is.finite(max_abs_diff_pp), sprintf("%.1f%s", max_abs_diff_pp, sig05), "")
	)

# Corrected all-year block summary. The denominator is the all-hour, all-year total within each phenotype and phone group.
figS8_block <- figS8_hour %>%
	mutate(block = factor(case_when(
		hour <= 5 ~ "00-05 night",
		hour <= 11 ~ "06-11 morning",
		hour <= 17 ~ "12-17 afternoon",
		TRUE ~ "18-23 evening"
	), levels = c("00-05 night", "06-11 morning", "12-17 afternoon", "18-23 evening"))) %>%
	group_by(dx_grp, phone.luck, block) %>%
	summarise(call_n = sum(call_n), .groups = "drop") %>%
	group_by(dx_grp, phone.luck) %>%
	mutate(total_n = sum(call_n), pct = ifelse(total_n > 0, call_n / total_n, NA_real_)) %>%
	ungroup() %>%
	select(dx_grp, phone.luck, block, call_n, total_n, pct) %>%
	pivot_wider(names_from = phone.luck, values_from = c(call_n, total_n, pct), names_sep = "_") %>%
	mutate(
		diff_pp = 100 * (pct_high - pct_low),
		dx_grp = factor(as.character(dx_grp), levels = rev(dxs.all))
	)

pS8A <- ggplot(figS8_stats, aes(x = factor(year), y = dx_grp, fill = max_abs_diff_pp)) +
	geom_tile(color = "white", linewidth = .25) +
	geom_text(aes(label = label), size = 2.35, fontface = "bold") +
	scale_fill_viridis_c(name = "Max hourly H-L difference", option = "C", na.value = "grey95", guide = guide_colorbar(direction = "horizontal", title.position = "left", barwidth = grid::unit(4.6, "cm"), barheight = grid::unit(.35, "cm"))) +
	labs(title = "a. High- and low-score circadian profiles", x = NULL, y = NULL) +
	fig_theme(base_size = 8.6) +
	theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_text(size = 7.5, face = "bold"), legend.text = element_text(size = 7.5), legend.margin = margin(0, 0, 2, 0), plot.title = element_text(size = 10.5, face = "bold", hjust = .5), axis.text.x = element_text(angle = 45, hjust = 1))

pS8B <- ggplot(figS8_block, aes(x = diff_pp, y = dx_grp, fill = block)) +
	geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
	geom_col(position = position_dodge(width = .72), width = .62) +
	labs(title = "b. Corrected time-block difference: high minus low", x = "High-low difference in within-phenotype hourly share", y = NULL, fill = NULL) +
	fig_theme(base_size = 8.6) +
	theme(plot.title = element_text(size = 10.5, face = "bold", hjust = .5), legend.position = "top", legend.text = element_text(size = 7.5))

FigS8 <- pS8A | pS8B
save_plot(FigS8, "FigS8.png", width = 12.8, height = 7.6, dpi = 600)
writexl::write_xlsx(list(
	hourly_counts = figS8_hour,
	circadian_difference_tests = figS8_stats,
	time_block_high_low_difference = figS8_block,
	configuration = tibble(note = "Corrected FigS8. max_abs_diff_pp is the largest absolute hourly percentage-point difference between high and low groups within phenotype and year. time-block pct uses all-year totals within phenotype and phone group as denominators.")
), "FigS8.out.xlsx")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS9. Six-panel circular on-scene-duration contrast by phone-score group
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS9c_dxs <- dxs.vip
figS9c_years <- 2018:2023
figS9_yin_path_candidates <- c(
	file.path(dir0, "files", "YinYang.png"),
	"YinYang.png",
	file.path(dir0, "files", "YinYang_clean.png"),
	"YinYang_clean.png"
)
figS9_yin_path <- figS9_yin_path_candidates[file.exists(figS9_yin_path_candidates)][1]
figS9_yin_img <- NULL
figS9_yin_center_radius <- 0.315
if (!is.na(figS9_yin_path)) {
	figS9_yin_img <- read_center_image(figS9_yin_path, circle_alpha = TRUE)
} else {
	warning("YinYang image not found. Put YinYang_clean.png in the working directory or D:/files/.")
}
figS9c_smooth_circ <- function(x, k = 3L) {
	x <- as.numeric(x)
	n <- length(x)
	if (!n) return(x)
	if (n == 1) return(x)
	h <- floor(k / 2)
	sapply(seq_len(n), function(i) {
		idx <- ((i - h - 1):(i + h - 1)) %% n + 1
		mean(x[idx], na.rm = TRUE)
	})
}
figS9c_base <- bind_rows(lapply(figS9c_years, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(dx_grp %in% figS9c_dxs, phone.luck %in% c("low", "high"), is.finite(hour), between(as.integer(hour), 0, 23)) %>%
		transmute(
			year = as.integer(y),
			dx_grp = factor(as.character(dx_grp), levels = figS9c_dxs),
			phone.luck = factor(as.character(phone.luck), levels = c("low", "high")),
			hour = as.integer(hour),
			onsite = suppressWarnings(as.numeric(现场时间) / 60)
		) %>%
		filter(is.finite(onsite), onsite >= 0, onsite <= 180)
}))
figS9c_hour <- figS9c_base %>%
	group_by(year, dx_grp, phone.luck, hour) %>%
	summarise(mean_onsite = mean(onsite, na.rm = TRUE), median_onsite = median(onsite, na.rm = TRUE), n = n(), .groups = "drop") %>%
	complete(year = figS9c_years, dx_grp = figS9c_dxs, phone.luck = c("low", "high"), hour = 0:23, fill = list(mean_onsite = NA_real_, median_onsite = NA_real_, n = 0))
figS9c_stats <- figS9c_base %>%
	group_by(year, dx_grp, hour) %>%
	group_modify(~ {
		d <- .x
		if (sum(d$phone.luck == "high") < 15 || sum(d$phone.luck == "low") < 15) return(tibble(mean_high = NA_real_, mean_low = NA_real_, diff_min = NA_real_, p_lm = NA_real_))
		d <- d %>% mutate(phone_high = as.integer(phone.luck == "high"))
		fit <- tryCatch(lm(log1p(onsite) ~ phone_high, data = d), error = function(e) NULL)
		p_lm <- NA_real_
		if (!is.null(fit)) {
			td <- broom::tidy(fit) %>% filter(term == "phone_high")
			if (nrow(td)) p_lm <- td$p.value[1]
		}
		tibble(mean_high = mean(d$onsite[d$phone.luck == "high"], na.rm = TRUE), mean_low = mean(d$onsite[d$phone.luck == "low"], na.rm = TRUE), diff_min = mean_high - mean_low, p_lm = p_lm)
	}) %>%
	ungroup() %>%
	group_by(year, dx_grp) %>%
	mutate(p_adj = p.adjust(p_lm, "BH"), sig05 = sig_star05(p_adj)) %>%
	ungroup()
figS9c_contrast <- figS9c_hour %>%
	select(year, dx_grp, phone.luck, hour, mean_onsite, median_onsite, n) %>%
	pivot_wider(names_from = phone.luck, values_from = c(mean_onsite, median_onsite, n), names_sep = "_") %>%
	mutate(
		diff_mean_min = mean_onsite_high - mean_onsite_low,
		diff_median_min = median_onsite_high - median_onsite_low
	) %>%
	group_by(year, dx_grp) %>%
	arrange(hour, .by_group = TRUE) %>%
	mutate(diff_mean_smooth = figS9c_smooth_circ(diff_mean_min, k = 3L)) %>%
	ungroup() %>%
	left_join(figS9c_stats %>% select(year, dx_grp, hour, p_lm, p_adj, sig05), by = c("year", "dx_grp", "hour"))
figS9c_scale <- figS9c_contrast %>%
	group_by(year) %>%
	summarise(max_abs = max(abs(diff_mean_min), abs(diff_mean_smooth), na.rm = TRUE), .groups = "drop") %>%
	mutate(max_abs = pmax(max_abs, 0.35))
figS9c_scale_map <- setNames(figS9c_scale$max_abs, figS9c_scale$year)
figS9c_cols <- setNames(grDevices::adjustcolor(unname(dxs.vip.color[figS9c_dxs]), alpha.f = 0.74), figS9c_dxs)
figS9c_circle_file <- tempfile(fileext = ".png")
png(figS9c_circle_file, width = 12.4, height = 8.65, units = "in", res = 420, bg = "transparent")
par(mfrow = c(2, 3), mar = c(0.72, 0.00, 1.72, 0.00), oma = c(0.60, 0.00, 0.10, 0.00), xaxs = "i", yaxs = "i")
for (yy in figS9c_years) {
	d0 <- figS9c_contrast %>% filter(year == yy) %>% arrange(dx_grp, hour)
	y_abs <- unname(figS9c_scale_map[as.character(yy)])
	if (!is.finite(y_abs) || y_abs <= 0) y_abs <- 0.5
	circlize::circos.clear()
	circlize::circos.par(start.degree = 90, gap.degree = 2, cell.padding = c(0, 0, 0, 0), track.margin = c(0.0015, 0.0015), canvas.xlim = c(-1.10, 1.10), canvas.ylim = c(-1.13, 1.13))
	circlize::circos.initialize(factors = "all", xlim = c(0, 24))
	for (i in seq_along(figS9c_dxs)) {
		dx <- figS9c_dxs[i]
		dx_dat <- d0 %>% filter(dx_grp == dx) %>% arrange(hour)
		circlize::circos.trackPlotRegion(
			factors = "all",
			track.index = i,
			ylim = c(-1.20 * y_abs, 1.20 * y_abs),
			bg.col = grDevices::adjustcolor(figS9c_cols[dx], alpha.f = .18),
			bg.border = NA,
			track.height = min(0.12, 0.72 / length(figS9c_dxs)),
			panel.fun = function(...) {
				circlize::circos.lines(c(0, 24), c(0, 0), col = "grey62", lwd = 0.75)
				for (jj in seq_len(nrow(dx_dat))) {
					val <- dx_dat$diff_mean_min[jj]
					if (!is.finite(val)) next
					fill_col <- if (val >= 0) grDevices::adjustcolor(figS9c_cols[dx], alpha.f = .80) else grDevices::adjustcolor("#4C78A8", alpha.f = .82)
					circlize::circos.rect(dx_dat$hour[jj], min(0, val), dx_dat$hour[jj] + 1, max(0, val), col = fill_col, border = NA)
				}
				circlize::circos.lines(dx_dat$hour + 0.5, clamp(dx_dat$diff_mean_smooth, -1.20 * y_abs, 1.20 * y_abs), col = "black", lwd = 1.00)
				d_star <- dx_dat %>% filter(nzchar(sig05), is.finite(p_adj))
				if (nrow(d_star)) {
					for (kk in seq_len(nrow(d_star))) {
						v0 <- d_star$diff_mean_smooth[kk]
						if (!is.finite(v0)) v0 <- 0
						y_star <- ifelse(v0 >= 0, 1.06 * y_abs, -1.06 * y_abs)
						circlize::circos.text(d_star$hour[kk] + .5, y_star, labels = d_star$sig05[kk], cex = .42, font = 2, col = ifelse(v0 >= 0, figS9c_cols[dx], "#4C78A8"), facing = "inside", niceFacing = TRUE)
					}
				}
			}
		)
	}
	circlize::circos.axis(h = "top", major.at = 0:23, labels = sprintf("%02d", 0:23), labels.cex = .86, labels.font = 2, minor.ticks = 0, sector.index = "all", track.index = 1)
	if (!is.null(figS9_yin_img)) graphics::rasterImage(figS9_yin_img, -figS9_yin_center_radius, -figS9_yin_center_radius, figS9_yin_center_radius, figS9_yin_center_radius, interpolate = TRUE)
	title(main = yy, font.main = 2, cex.main = 1.25, col.main = "#009E73", line = 0.10)
}
dev.off()
circlize::circos.clear()
figS9c_raster <- png::readPNG(figS9c_circle_file, native = FALSE)
unlink(figS9c_circle_file, force = TRUE)
figS9c_legend <- tibble(
	dx_grp = factor(figS9c_dxs, levels = figS9c_dxs),
	label = dx_to_eng(figS9c_dxs),
	x = seq(0.130, 0.770, length.out = length(figS9c_dxs)),
	x_text = x + .009,
	y = .985
)
pS9 <- ggplot() +
	annotation_custom(grid::rasterGrob(figS9c_raster, interpolate = TRUE), xmin = -0.255, xmax = 1.255, ymin = -0.275, ymax = 0.920) +
	annotate("rect", xmin = .145, xmax = .165, ymin = 1.008, ymax = 1.033, fill = grDevices::adjustcolor(figS9c_cols[1], alpha.f = .80), color = NA) +
	annotate("text", x = .173, y = 1.021, label = "High > low on-scene time", hjust = 0, size = 2.35, fontface = "bold") +
	annotate("rect", xmin = .395, xmax = .415, ymin = 1.008, ymax = 1.033, fill = grDevices::adjustcolor("#4C78A8", alpha.f = .82), color = NA) +
	annotate("text", x = .423, y = 1.021, label = "Low > high on-scene time", hjust = 0, size = 2.35, fontface = "bold") +
	annotate("segment", x = .645, xend = .675, y = 1.021, yend = 1.021, linewidth = .85, color = "black") +
	annotate("text", x = .685, y = 1.021, label = "3-hour circular moving average", hjust = 0, size = 2.35, fontface = "bold") +
	geom_point(data = figS9c_legend, aes(x = x, y = y, color = dx_grp), size = 1.85) +
	geom_text(data = figS9c_legend, aes(x = x_text, y = y, label = label), hjust = 0, size = 2.55, fontface = "bold", color = "grey20") +
	coord_cartesian(xlim = c(0, 1), ylim = c(-0.275, 1.075), expand = FALSE, clip = "off") +
	scale_color_manual(values = figS9c_cols, guide = "none") +
	labs(title = "Circular on-scene-duration contrast by phone-score group") +
	theme_void(base_size = 10) +
	theme(
		plot.title = element_text(size = 11, face = "bold", hjust = .5, margin = margin(t = 3, b = 7)),
		plot.background = element_rect(fill = "white", color = NA),
		panel.background = element_rect(fill = "white", color = NA),
		plot.margin = margin(4, 2, 7, 2)
	)
save_plot(pS9, "FigS9.png", width = 12.4, height = 8.95, dpi = 600, bg = "white")
writexl::write_xlsx(list(
	hourly_onscene_summary = figS9c_hour,
	hourly_onscene_contrast = figS9c_contrast,
	hourly_high_low_tests = figS9c_stats,
	configuration = tibble(selected_phenotypes = paste(figS9c_dxs, collapse = ", "), selected_years = paste(figS9c_years, collapse = ", "), note = "Six yearly circular panels show the high-minus-low difference in mean on-scene duration by hour for the eight VIP phenotypes. Bars show raw differences in minutes, black lines show 3-hour circular averages, and stars use hour-specific BH-adjusted P values from log1p(on-scene duration) ~ high within each year and phenotype.")
), "FigS9.out.xlsx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS10. Long-tail phone-score gradient heatmaps
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS10_period_levels <- c("2013-2016 baseline", "2017-2019 pre-COVID", "2020 first wave", "2021-2022 pre-reopening", "2023 post-reopening", "2024 sustained tail")
figS10_ref_period <- "2017-2019 pre-COVID"
figS10_base <- bind_rows(lapply(years, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(dx_grp %in% dxs.all, phone.luck %in% c("low", "high"), !is.na(.data[["日期"]])) %>%
		transmute(year = as.integer(y), dx_grp = factor(as.character(dx_grp), levels = dxs.all), phone.luck = factor(as.character(phone.luck), levels = c("low", "high")))
})) %>%
	mutate(period = case_when(
		year <= 2016 ~ "2013-2016 baseline",
		year <= 2019 ~ "2017-2019 pre-COVID",
		year == 2020 ~ "2020 first wave",
		year %in% 2021:2022 ~ "2021-2022 pre-reopening",
		year == 2023 ~ "2023 post-reopening",
		year == 2024 ~ "2024 sustained tail",
		TRUE ~ NA_character_
	)) %>%
	filter(!is.na(period)) %>%
	mutate(period = factor(period, levels = figS10_period_levels))
figS10_counts <- figS10_base %>% count(period, dx_grp, phone.luck, name = "n") %>% complete(period = figS10_period_levels, dx_grp = dxs.all, phone.luck = c("low", "high"), fill = list(n = 0))
figS10_den <- figS10_base %>% count(period, phone.luck, name = "N")
figS10_rr <- figS10_counts %>%
	left_join(figS10_den, by = c("period", "phone.luck")) %>%
	pivot_wider(names_from = phone.luck, values_from = c(n, N), names_sep = "_") %>%
	mutate(
		pct_low = n_low / N_low,
		pct_high = n_high / N_high,
		RR = pct_high / pct_low,
		se_logRR = sqrt(1 / pmax(n_high, 1) - 1 / pmax(N_high, 1) + 1 / pmax(n_low, 1) - 1 / pmax(N_low, 1)),
		RR_lo = exp(log(RR) - 1.96 * se_logRR),
		RR_hi = exp(log(RR) + 1.96 * se_logRR),
		p = 2 * pnorm(abs(log(RR) / se_logRR), lower.tail = FALSE)
	) %>%
	group_by(period) %>%
	mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj)) %>%
	ungroup() %>%
	mutate(dx_grp = factor(as.character(dx_grp), levels = rev(dxs.all)), label = ifelse(is.finite(RR), sprintf("%.2f%s", RR, sig05), ""))
figS10_ref <- figS10_rr %>% filter(period == figS10_ref_period) %>% select(dx_grp, RR_ref = RR, se_ref = se_logRR)
figS10_shift <- figS10_rr %>%
	left_join(figS10_ref, by = "dx_grp") %>%
	filter(period != figS10_ref_period) %>%
	mutate(
		RRR_vs_ref = RR / RR_ref,
		se_logRRR = sqrt(se_logRR^2 + se_ref^2),
		RRR_lo = exp(log(RRR_vs_ref) - 1.96 * se_logRRR),
		RRR_hi = exp(log(RRR_vs_ref) + 1.96 * se_logRRR),
		p_shift = 2 * pnorm(abs(log(RRR_vs_ref) / se_logRRR), lower.tail = FALSE)
	) %>%
	group_by(period) %>%
	mutate(p_adj_shift = p.adjust(p_shift, "BH"), sig05_shift = sig_star05(p_adj_shift)) %>%
	ungroup() %>%
	mutate(dx_grp = factor(as.character(dx_grp), levels = rev(dxs.all)), label_shift = ifelse(is.finite(RRR_vs_ref), sprintf("%.2f%s", RRR_vs_ref, sig05_shift), ""))
pS10A <- ggplot(figS10_rr, aes(x = period, y = dx_grp, fill = log2(RR))) +
	geom_tile(color = "white", linewidth = .25) +
	geom_text(aes(label = label), size = 2.05, fontface = "bold") +
	scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "log2 RR\nhigh vs low", na.value = "grey95") +
	labs(title = "a. Period-level high-vs-low RR", x = NULL, y = NULL) +
	fig_theme(base_size = 8.2) +
	theme(plot.title = element_text(size = 10.4, face = "bold", hjust = .5), axis.text.x = element_text(angle = 35, hjust = 1))
pS10B <- ggplot(figS10_shift, aes(x = period, y = dx_grp, fill = log2(RRR_vs_ref))) +
	geom_tile(color = "white", linewidth = .25) +
	geom_text(aes(label = label_shift), size = 2.05, fontface = "bold") +
	scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "log2 ratio\nvs 2017-2019", na.value = "grey95") +
	labs(title = "b. Change relative to 2017-2019", x = NULL, y = NULL) +
	fig_theme(base_size = 8.2) +
	theme(plot.title = element_text(size = 10.4, face = "bold", hjust = .5), axis.text.x = element_text(angle = 35, hjust = 1))
FigS10 <- pS10A / pS10B + plot_layout(heights = c(1, 1))
save_plot(FigS10, "FigS10.png", width = 11.8, height = 13.0, dpi = 600)
writexl::write_xlsx(list(
	period_high_low_RR = figS10_rr,
	gradient_shift_vs_preCOVID = figS10_shift,
	configuration = tibble(reference_period = figS10_ref_period, note = "Panel a shows absolute high-vs-low phenotype RR within each period. Panel b shows each period-specific high-vs-low RR divided by the 2017-2019 high-vs-low RR.")
), "FigS10.out.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS10b. Quarterly dynamics of high-vs-low phone-score gradients
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS10_quarter_focus <- intersect(c("Violence", "Intoxication", "Psychiatric", "Traffic", "CVD", "Ob/Gyn", "Pediatric", "Digestive"), dxs.all)
figS10_quarter_base <- bind_rows(lapply(years, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(dx_grp %in% dxs.all, phone.luck %in% c("low", "high"), !is.na(.data[["日期"]])) %>%
		transmute(
			year = as.integer(y),
			date = as.Date(.data[["日期"]]),
			quarter = lubridate::floor_date(as.Date(.data[["日期"]]), "quarter"),
			dx_grp = factor(as.character(dx_grp), levels = dxs.all),
			phone.luck = factor(as.character(phone.luck), levels = c("low", "high"))
		)
})) %>%
	filter(year >= 2017)

figS10_quarter_counts <- figS10_quarter_base %>%
	count(quarter, dx_grp, phone.luck, name = "n") %>%
	complete(quarter = sort(unique(figS10_quarter_base$quarter)), dx_grp = dxs.all, phone.luck = c("low", "high"), fill = list(n = 0))
figS10_quarter_den <- figS10_quarter_base %>% count(quarter, phone.luck, name = "N")
figS10_quarter_rr <- figS10_quarter_counts %>%
	left_join(figS10_quarter_den, by = c("quarter", "phone.luck")) %>%
	pivot_wider(names_from = phone.luck, values_from = c(n, N), names_sep = "_") %>%
	mutate(
		RR = (n_high / N_high) / (n_low / N_low),
		se_logRR = sqrt(1 / pmax(n_high, 1) - 1 / pmax(N_high, 1) + 1 / pmax(n_low, 1) - 1 / pmax(N_low, 1)),
		RR_lo = exp(log(RR) - 1.96 * se_logRR),
		RR_hi = exp(log(RR) + 1.96 * se_logRR),
		p = purrr::pmap_dbl(list(n_high, N_high, n_low, N_low), function(a, A, b, B) suppressWarnings(tryCatch(chisq.test(matrix(c(a, A - a, b, B - b), nrow = 2))$p.value, error = function(e) NA_real_)))
	) %>%
	group_by(quarter) %>%
	mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj)) %>%
	ungroup() %>%
	mutate(dx_label = factor(dx_to_eng(dx_grp), levels = dx_to_eng(figS10_quarter_focus)))

figS10_quarter_event_lines <- tibble(
	date = as.Date(c("2020-01-23", "2022-03-14", "2022-12-07")),
	event = c("COVID-19 PHSM", "2022 PHSM", "Full reopening")
)
figS10_quarter_plot_dat <- figS10_quarter_rr %>% filter(dx_grp %in% figS10_quarter_focus)
pS10C <- ggplot(figS10_quarter_plot_dat, aes(quarter, RR, color = dx_grp, fill = dx_grp)) +
	geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
	geom_vline(data = figS10_quarter_event_lines, aes(xintercept = date), inherit.aes = FALSE, linetype = "dotted", color = "grey35", linewidth = .45) +
	geom_ribbon(aes(ymin = RR_lo, ymax = RR_hi), alpha = .12, color = NA) +
	geom_line(linewidth = .75, na.rm = TRUE) +
	geom_point(size = 1.4, na.rm = TRUE) +
	facet_wrap(~dx_label, ncol = 2, scales = "free_y") +
	scale_color_manual(values = dxs.all.color[dxs.all], guide = "none") +
	scale_fill_manual(values = dxs.all.color[dxs.all], guide = "none") +
	scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(.01, .02))) +
	labs(title = "c. Quarterly high-vs-low phone-score gradients", x = NULL, y = "RR: high vs low phone-score group") +
	fig_theme(base_size = 8.8) +
	theme(plot.title = element_text(size = 11, face = "bold", hjust = .5), axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
FigS10 <- ((pS10A / pS10B) | pS10C) + plot_layout(widths = c(1.05, 1.0))
save_plot(FigS10, "FigS10.png", width = 15.2, height = 13.0, dpi = 600)
writexl::write_xlsx(list(
	period_high_low_RR = figS10_rr,
	gradient_shift_vs_preCOVID = figS10_shift,
	quarterly_RR_all_dxs = figS10_quarter_rr,
	quarterly_RR_plotted = figS10_quarter_plot_dat,
	event_lines = figS10_quarter_event_lines,
	configuration = tibble(reference_period = figS10_ref_period, note = "Period heatmaps and quarterly gradients are merged into the current FigS10.")
), "FigS10.out.xlsx")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS11a. Phone-score decile dose-response inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS11_decile_base <- bind_rows(lapply(years, function(y) {
	d <- dat1.list[[as.character(y)]]
	if (is.null(d) || !nrow(d)) return(tibble())
	d %>%
		filter(dx_grp %in% dxs.all, is.finite(phone.sco)) %>%
		transmute(
			year = as.integer(y),
			dx_grp = factor(as.character(dx_grp), levels = dxs.all),
			phone.sco = as.numeric(phone.sco)
		)
}))
figS11_decile_base$score_decile <- score_decile_by_year(figS11_decile_base$phone.sco, figS11_decile_base$year)
figS11_decile_base <- figS11_decile_base %>% filter(is.finite(score_decile), score_decile >= 1, score_decile <= 10)

figS11_decile_counts <- figS11_decile_base %>%
	count(score_decile, dx_grp, name = "n") %>%
	complete(score_decile = 1:10, dx_grp = dxs.all, fill = list(n = 0))
figS11_decile_den <- figS11_decile_base %>% count(score_decile, name = "N")
figS11_decile <- figS11_decile_counts %>%
	left_join(figS11_decile_den, by = "score_decile") %>%
	mutate(pct = n / N) %>%
	group_by(dx_grp) %>%
	mutate(pct_D1 = pct[score_decile == 1][1], RR_vs_D1 = pct / pct_D1) %>%
	ungroup() %>%
	mutate(dx_grp = factor(as.character(dx_grp), levels = rev(dxs.all)), label = ifelse(is.finite(RR_vs_D1), sprintf("%.2f", RR_vs_D1), ""))

figS11_fit_decile_trend <- function(dx) {
	d <- figS11_decile_base %>%
		mutate(case = as.integer(dx_grp == dx)) %>%
		count(year, score_decile, case, name = "n") %>%
		pivot_wider(names_from = case, values_from = n, values_fill = 0) %>%
		rename(control = `0`, event = `1`) %>%
		mutate(total = event + control)
	if (!nrow(d) || sum(d$event, na.rm = TRUE) < 20) return(tibble(dx_grp = dx, OR_per_decile = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_))
	fit <- tryCatch(glm(cbind(event, total - event) ~ score_decile + factor(year), family = binomial(), data = d), error = function(e) NULL)
	if (is.null(fit)) return(tibble(dx_grp = dx, OR_per_decile = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_))
	td <- broom::tidy(fit) %>% filter(term == "score_decile")
	if (!nrow(td)) return(tibble(dx_grp = dx, OR_per_decile = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_))
	tibble(dx_grp = dx, OR_per_decile = exp(td$estimate[1]), OR_lo = exp(td$estimate[1] - 1.96 * td$std.error[1]), OR_hi = exp(td$estimate[1] + 1.96 * td$std.error[1]), p = td$p.value[1])
}
figS11_decile_trend <- bind_rows(lapply(dxs.all, figS11_fit_decile_trend)) %>%
	mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj), dx_label = factor(dx_to_eng(dx_grp), levels = rev(dx_to_eng(dxs.all))))
figS11_panel_titles <- c(
	A = "a. Disease membership across phone-score deciles",
	B = "b. Linear decile trend",
	C = "c. Sequential adjustment",
	D = "d. Full-model attenuation"
)

pS11A <- ggplot(figS11_decile, aes(x = factor(score_decile), y = dx_grp, fill = log2(RR_vs_D1))) +
	geom_tile(color = "white", linewidth = .25) +
	geom_text(aes(label = label), size = 2.25, fontface = "bold") +
	scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "log2 RR\nvs decile 1", na.value = "grey95") +
	labs(title = figS11_panel_titles[["A"]], x = "Phone-score decile within year", y = NULL) +
	fig_theme(base_size = 8.4) +
	theme(plot.title = element_text(size = 10.5, face = "bold", hjust = .5))
pS11B <- ggplot(figS11_decile_trend, aes(x = OR_per_decile, y = dx_label, xmin = OR_lo, xmax = OR_hi)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
	geom_errorbar(width = .18, orientation = "y", linewidth = .65, color = "grey35") +
	geom_point(size = 2.3, color = "black") +
	geom_text(aes(label = sig05, x = OR_hi), hjust = -.20, size = 3.4, fontface = "bold") +
	scale_x_log10(breaks = c(.96, .98, 1, 1.02, 1.04, 1.06)) +
	labs(title = figS11_panel_titles[["B"]], x = "OR per one higher phone-score decile", y = NULL) +
	fig_theme(base_size = 8.4) +
	theme(plot.title = element_text(size = 10.5, face = "bold", hjust = .5))
FigS11_decile_panel <- (pS11A | pS11B) + plot_layout(widths = c(1.25, .95))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS11b. Sequential adjustment and attenuation inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS11_adj_env_cols <- unique(fig2_num_vars)
figS11_adj_env_cols <- figS11_adj_env_cols[nzchar(figS11_adj_env_cols)]
figS11_adj_env_cols <- figS11_adj_env_cols[vapply(figS11_adj_env_cols, function(v) {
	any(vapply(dat1.list[as.character(years)], function(d0) !is.null(d0) && v %in% names(d0), logical(1)))
}, logical(1))]

figS11_adj_base <- bind_rows(lapply(years, function(y) {
	d0 <- dat1.list[[as.character(y)]]
	if (is.null(d0) || !nrow(d0)) return(tibble())
	env_df <- if (length(figS11_adj_env_cols)) as_tibble(setNames(lapply(figS11_adj_env_cols, function(v) suppressWarnings(as.numeric(get_col0(d0, v)))), paste0("env", seq_along(figS11_adj_env_cols)))) else tibble(.rows = nrow(d0))
	bind_cols(
		tibble(
			year = as.integer(y),
			date = as.Date(get_col0(d0, "日期")),
			hour = suppressWarnings(as.integer(get_col0(d0, "hour"))),
			disease = factor(as.character(get_col0(d0, "dx_grp")), levels = dxs.all),
			phone_group = as.character(get_col0(d0, "phone.luck")),
			age = suppressWarnings(as.numeric(get_col0(d0, vars.dxs[2]))),
			sex = as.character(get_col0(d0, vars.dxs[1])),
			lon = suppressWarnings(as.numeric(get_col0(d0, vars.basic.ems[11]))),
			lat = suppressWarnings(as.numeric(get_col0(d0, vars.basic.ems[12]))),
			geo_type = as.character(get_col0(d0, "geo.type1"))
		),
		env_df
	)
})) %>%
	filter(phone_group %in% c("low", "high"), !is.na(disease), disease %in% dxs.all) %>%
	mutate(
		phone_high = as.integer(phone_group == "high"),
		month = ifelse(is.na(date), 0L, lubridate::month(date)),
		dow = ifelse(is.na(date), "Missing", as.character(lubridate::wday(date, label = TRUE, week_start = 1))),
		hour = ifelse(is.finite(hour), hour, -1L),
		sex = ifelse(is.na(sex) | sex == "", "Missing", sex),
		age_band = cut(age, breaks = c(-Inf, 18, 35, 50, 65, 80, Inf), labels = c("<=18", "19-35", "36-50", "51-65", "66-80", ">80"), right = TRUE),
		age_band = forcats::fct_na_value_to_level(age_band, level = "Missing"),
		geo_type = ifelse(is.na(geo_type) | geo_type == "", "Missing", geo_type),
		geo_cell_raw = ifelse(is.finite(lon) & is.finite(lat), paste0(round(lon / 0.03), "_", round(lat / 0.03)), "Missing")
	)
if (nrow(figS11_adj_base)) {
	figS11_adj_top_cells <- figS11_adj_base %>% count(geo_cell_raw, sort = TRUE) %>% slice_head(n = 120) %>% pull(geo_cell_raw)
	figS11_adj_base <- figS11_adj_base %>% mutate(geo_cell = ifelse(geo_cell_raw %in% figS11_adj_top_cells, geo_cell_raw, "Sparse"))
	for (i in seq_along(figS11_adj_env_cols)) figS11_adj_base[[paste0("env", i, "_q")]] <- make_qcat(figS11_adj_base[[paste0("env", i)]])
}

figS11_adj_models <- list(
	`Crude` = character(0),
	`+ age + sex` = c("age_band", "sex"),
	`+ calendar + hour` = c("age_band", "sex", "year", "month", "dow", "hour"),
	`+ geo/environment` = c("age_band", "sex", "year", "month", "dow", "hour", "geo_type", "geo_cell", paste0("env", seq_along(figS11_adj_env_cols), "_q"))
)
figS11_fit_adjustment <- function(dx, model_name, covars) {
	d <- figS11_adj_base %>% mutate(case0 = as.integer(disease == dx))
	covars <- covars[covars %in% names(d)]
	if (length(covars)) {
		d_agg <- d %>% group_by(phone_high, across(all_of(covars))) %>% summarise(case = sum(case0), total = n(), .groups = "drop")
		covars2 <- covars[vapply(d_agg[covars], function(x) n_distinct(x, na.rm = FALSE) > 1, logical(1))]
	} else {
		d_agg <- d %>% group_by(phone_high) %>% summarise(case = sum(case0), total = n(), .groups = "drop")
		covars2 <- character(0)
	}
	if (nrow(d_agg) < 2 || n_distinct(d_agg$phone_high) < 2 || sum(d_agg$case, na.rm = TRUE) < 20) return(tibble(dx_grp = dx, model = model_name, OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_, n_strata = nrow(d_agg)))
	form <- as.formula(paste("cbind(case, total - case) ~", paste(c("phone_high", covars2), collapse = " + ")))
	fit <- tryCatch(suppressWarnings(glm(form, family = binomial(), data = d_agg)), error = function(e) NULL)
	if (is.null(fit)) return(tibble(dx_grp = dx, model = model_name, OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_, n_strata = nrow(d_agg)))
	td <- broom::tidy(fit) %>% filter(term == "phone_high")
	if (!nrow(td)) return(tibble(dx_grp = dx, model = model_name, OR = NA_real_, OR_lo = NA_real_, OR_hi = NA_real_, p = NA_real_, n_strata = nrow(d_agg)))
	tibble(dx_grp = dx, model = model_name, OR = exp(td$estimate[1]), OR_lo = exp(td$estimate[1] - 1.96 * td$std.error[1]), OR_hi = exp(td$estimate[1] + 1.96 * td$std.error[1]), p = td$p.value[1], n_strata = nrow(d_agg))
}
figS11_adj_res <- bind_rows(lapply(dxs.all, function(dx) {
	bind_rows(lapply(names(figS11_adj_models), function(m) figS11_fit_adjustment(dx, m, figS11_adj_models[[m]])))
})) %>%
	group_by(model) %>% mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj)) %>% ungroup() %>%
	mutate(model = factor(model, levels = names(figS11_adj_models)), dx_label = factor(dx_to_eng(dx_grp), levels = rev(dx_to_eng(dxs.all))))
figS11_adj_shift <- figS11_adj_res %>%
	select(dx_grp, model, OR) %>%
	pivot_wider(names_from = model, values_from = OR) %>%
	mutate(
		full_vs_crude_log_ratio = log(`+ geo/environment`) - log(Crude),
		dx_label = factor(dx_to_eng(dx_grp), levels = rev(dx_to_eng(dxs.all)))
	)

pS11C <- ggplot(figS11_adj_res, aes(x = OR, y = dx_label, xmin = OR_lo, xmax = OR_hi, color = model)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
	geom_errorbar(position = position_dodge(width = .70), width = .14, orientation = "y", linewidth = .55) +
	geom_point(position = position_dodge(width = .70), size = 1.8) +
	scale_x_log10(breaks = c(.75, .85, .9, 1, 1.1, 1.2, 1.35)) +
	labs(title = figS11_panel_titles[["C"]], x = "OR: high vs low phone-score group", y = NULL, color = NULL) +
	fig_theme(base_size = 8.8) +
	theme(plot.title = element_text(size = 10.5, face = "bold", hjust = .5), legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 7.2))
pS11D <- ggplot(figS11_adj_shift, aes(x = full_vs_crude_log_ratio, y = dx_label)) +
	geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
	geom_segment(aes(x = 0, xend = full_vs_crude_log_ratio, yend = dx_label), linewidth = .65, color = "grey60") +
	geom_point(size = 2.2, color = "black") +
	scale_x_continuous(labels = function(x) sprintf("%+.2f", x)) +
	labs(title = figS11_panel_titles[["D"]], x = "log(fully adjusted OR / crude OR)", y = NULL) +
	fig_theme(base_size = 8.8) +
	theme(plot.title = element_text(size = 10.5, face = "bold", hjust = .5))
FigS11_adjustment_panel <- (pS11C | pS11D) + plot_layout(widths = c(1.35, .85))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS11. Dose-response and sequential adjustment
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FigS11 <- (FigS11_decile_panel / FigS11_adjustment_panel) + plot_layout(heights = c(1, 1))
save_plot(FigS11, "FigS11.png", width = 15.2, height = 15.2, dpi = 600)
writexl::write_xlsx(list(
	decile_plot_data = figS11_decile,
	decile_trend_models = figS11_decile_trend,
	sequential_model_results = figS11_adj_res,
	attenuation_summary = figS11_adj_shift,
	model_covariates = tibble(model = names(figS11_adj_models), covariates = vapply(figS11_adj_models, paste, character(1), collapse = "; ")),
	base_summary = figS11_adj_base %>% summarise(n = n(), years = n_distinct(year), geo_cells = n_distinct(geo_cell), env_covariates = paste(figS11_adj_env_cols, collapse = "; ")),
	configuration = tibble(note = "This integrated figure combines phone-score decile dose-response and sequential adjustment/attenuation analyses.")
), "FigS11.out.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS12. Combined PHSM and reopening event-window analyses
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figS12_collect_daily <- function(start_date, end_date, dx_set) {
	yrs <- sort(unique(lubridate::year(seq.Date(start_date, end_date, by = "day"))))
	raw <- bind_rows(lapply(yrs, function(y) {
		d <- dat1.list[[as.character(y)]]
		if (is.null(d) || !nrow(d)) return(tibble())
		d %>% transmute(date = as.Date(.data[["日期"]]), dx_grp = factor(as.character(dx_grp), levels = dxs.all), phone.luck = factor(as.character(phone.luck), levels = c("low", "middle", "high"))) %>% filter(date >= start_date, date <= end_date, dx_grp %in% dx_set, phone.luck %in% c("low", "high"))
	}))
	raw %>% count(date, dx_grp, phone.luck, name = "count") %>% complete(date = seq.Date(start_date, end_date, by = "day"), dx_grp = dx_set, phone.luck = c("low", "high"), fill = list(count = 0))
}
figS12_fit_did <- function(daily, dx, pre_start, pre_end, post_start, post_end) {
	d <- daily %>% filter(dx_grp == dx, (date >= pre_start & date <= pre_end) | (date >= post_start & date <= post_end)) %>% mutate(post = as.integer(date >= post_start & date <= post_end), high = as.integer(phone.luck == "high"))
	if (nrow(d) < 4 || sum(d$count, na.rm = TRUE) < 10) return(tibble(dx_grp = dx, RR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
	fit <- tryCatch(glm(count ~ post * high, family = poisson(), data = d), error = function(e) NULL)
	if (is.null(fit)) return(tibble(dx_grp = dx, RR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
	td <- broom::tidy(fit) %>% filter(term == "post:high")
	if (!nrow(td)) return(tibble(dx_grp = dx, RR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
	tibble(dx_grp = dx, RR = exp(td$estimate[1]), lo = exp(td$estimate[1] - 1.96 * td$std.error[1]), hi = exp(td$estimate[1] + 1.96 * td$std.error[1]), p = td$p.value[1])
}
figS12_plot_did <- function(dat, ttl) {
	ggplot(dat, aes(x = RR, y = dx_label, color = dx_grp)) +
		geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
		geom_errorbar(aes(xmin = lo, xmax = hi), width = .18, orientation = "y", linewidth = .70) +
		geom_point(size = 2.35) +
		geom_text(aes(label = label, x = hi), hjust = -.18, size = 2.85, fontface = "bold", show.legend = FALSE) +
		scale_x_log10(breaks = c(.25, .4, .6, .8, 1, 1.5, 2, 3), limits = c(.25, 3.2)) +
		scale_color_manual(values = dxs.all.color[dxs.all], guide = "none") +
		labs(title = ttl, x = "DID rate ratio: high vs low", y = NULL) +
		fig_theme(base_size = 8.0) +
		theme(plot.title = element_text(size = 10.2, face = "bold", hjust = .5))
}
# PHSM window
figS12_phsm_dxs <- intersect(c("CVD", "Respiratory", "Endocrine", "Psychiatric"), dxs.all)
figS12_phsm_label <- setNames(dx_to_eng(figS12_phsm_dxs), figS12_phsm_dxs)
figS12_phsm_daily <- figS12_collect_daily(as.Date("2022-03-07"), as.Date("2022-03-27"), figS12_phsm_dxs) %>% mutate(dx_label = factor(figS12_phsm_label[as.character(dx_grp)], levels = figS12_phsm_label[figS12_phsm_dxs]))
figS12_phsm_did1 <- bind_rows(lapply(figS12_phsm_dxs, function(dx) figS12_fit_did(figS12_phsm_daily, dx, as.Date("2022-03-07"), as.Date("2022-03-13"), as.Date("2022-03-14"), as.Date("2022-03-20")))) %>% mutate(comparison = "Pre vs PHSM")
figS12_phsm_did2 <- bind_rows(lapply(figS12_phsm_dxs, function(dx) figS12_fit_did(figS12_phsm_daily, dx, as.Date("2022-03-14"), as.Date("2022-03-20"), as.Date("2022-03-21"), as.Date("2022-03-27")))) %>% mutate(comparison = "PHSM vs post")
figS12_phsm_did <- bind_rows(figS12_phsm_did1, figS12_phsm_did2) %>% group_by(comparison) %>% mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj), dx_label = factor(figS12_phsm_label[as.character(dx_grp)], levels = rev(figS12_phsm_label[figS12_phsm_dxs])), label = ifelse(is.finite(RR), sprintf("%.2f%s", RR, sig05), "")) %>% ungroup()
pS12A <- ggplot(figS12_phsm_daily, aes(date, count, color = phone.luck, group = phone.luck)) +
	geom_line(linewidth = .82) +
	geom_vline(xintercept = as.Date(c("2022-03-14", "2022-03-20")), linetype = "dashed", color = "orange", linewidth = .7) +
	facet_wrap(~dx_label, ncol = 1, scales = "free_y") +
	scale_color_manual(values = c(low = "grey60", high = "#D55E00"), labels = c(low = "Low score", high = "High score"), name = NULL) +
	scale_x_date(date_breaks = "4 days", date_labels = "%b %d") +
	labs(title = "a. March 2022 PHSM", x = NULL, y = "Daily calls") +
	fig_theme(base_size = 8.0) +
	theme(plot.title = element_text(size = 10.6, face = "bold", hjust = .5), legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
pS12B <- figS12_plot_did(figS12_phsm_did %>% filter(comparison == "Pre vs PHSM"), "b. PHSM DID: pre vs during")
pS12C <- figS12_plot_did(figS12_phsm_did %>% filter(comparison == "PHSM vs post"), "c. PHSM DID: during vs post")
# Reopening window
figS12_open_dxs <- intersect(c("Violence", "Traffic", "Fall", "Intoxication", "CVD", "Respiratory", "Endocrine", "Psychiatric"), dxs.all)
figS12_open_label <- c(Violence = "Violence", Traffic = "Accident", Fall = "Fall", Intoxication = "Poisoning", CVD = "CVD", Respiratory = "Respiratory", Endocrine = "Endocrine", Psychiatric = "Psychiatric")
figS12_open_daily <- figS12_collect_daily(as.Date("2022-11-01"), as.Date("2023-01-01"), figS12_open_dxs) %>% mutate(dx_label = factor(figS12_open_label[as.character(dx_grp)], levels = figS12_open_label[figS12_open_dxs]))
figS12_open_did1 <- bind_rows(lapply(figS12_open_dxs, function(dx) figS12_fit_did(figS12_open_daily, dx, as.Date("2022-11-04"), as.Date("2022-11-10"), as.Date("2022-11-11"), as.Date("2022-11-17")))) %>% mutate(comparison = "Semi-removal")
figS12_open_did2 <- bind_rows(lapply(figS12_open_dxs, function(dx) figS12_fit_did(figS12_open_daily, dx, as.Date("2022-11-30"), as.Date("2022-12-06"), as.Date("2022-12-07"), as.Date("2022-12-13")))) %>% mutate(comparison = "Full reopening")
figS12_open_did <- bind_rows(figS12_open_did1, figS12_open_did2) %>% group_by(comparison) %>% mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj), dx_label = factor(figS12_open_label[as.character(dx_grp)], levels = rev(figS12_open_label[figS12_open_dxs])), label = ifelse(is.finite(RR), sprintf("%.2f%s", RR, sig05), "")) %>% ungroup()
pS12D <- ggplot(figS12_open_daily, aes(date, count, color = phone.luck, group = phone.luck)) +
	geom_line(linewidth = .70) +
	geom_vline(xintercept = as.Date(c("2022-11-11", "2022-12-07")), linetype = "dashed", color = "orange", linewidth = .65) +
	facet_wrap(~dx_label, ncol = 1, scales = "free_y") +
	scale_color_manual(values = c(low = "grey60", high = "#D55E00"), labels = c(low = "Low score", high = "High score"), name = NULL) +
	scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
	labs(title = "d. Late-2022 reopening", x = NULL, y = "Daily calls") +
	fig_theme(base_size = 7.4) +
	theme(plot.title = element_text(size = 10.6, face = "bold", hjust = .5), legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
pS12E <- figS12_plot_did(figS12_open_did %>% filter(comparison == "Semi-removal"), "e. Reopening DID: Nov 11 semi-removal")
pS12F <- figS12_plot_did(figS12_open_did %>% filter(comparison == "Full reopening"), "f. Reopening DID: Dec 7 full reopening")
FigS12 <- ((pS12A | (pS12B / pS12C)) / (pS12D | (pS12E / pS12F))) + plot_layout(heights = c(1, 1), widths = c(1.12, .88))
save_plot(FigS12, "FigS12.png", width = 12.0, height = 13.2, dpi = 600)
writexl::write_xlsx(list(
	phsm_daily_counts = figS12_phsm_daily,
	phsm_DID_models = figS12_phsm_did,
	reopening_daily_counts = figS12_open_daily,
	reopening_DID_models = figS12_open_did,
	configuration = tibble(note = "This combined supplementary figure keeps PHSM and reopening within one figure. DID uses Poisson interaction models count ~ post * high for paired 7-day windows.")
), "FigS12.out.xlsx")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FINAL supplementary-figure reassembly (2026-08-12)
#
# IMPORTANT:
#   * No analysis or ML code before the figure section is changed.
#   * Fig1-Fig5 are the main figures; the former Fig6 is exported as FigS1.
#   * Existing supplementary analyses/panels are reused; only the FINAL
#     supplementary assembly, numbering, and workbook packaging are changed.
#
# Working supplementary structure after visual review:
#   FigS0  Annual AI-reconstructed phenotype composition
#   FigS1  On-scene time and circular hourly mechanism
#   FigS2  Phone-score construction
#   FigS3  removed (MacBERT CV integrated into main Fig2)
#   FigS4  Geographic phenotype validation
#   FigS5  PHSM companion to main Fig5
#   FigS6-FigS12 retain current working numbers pending final renumbering
#
# TEST figures are no longer emitted; the former TEST.Fig1 content has been promoted and expanded into main Fig2.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Preserve all legacy supplementary workbooks in memory BEFORE official
# Final FigS1-FigS12 workbooks overwrite their assigned output file names.
ems120_read_xlsx_all <- function(path) {
	if (!file.exists(path)) return(list())
	sh <- readxl::excel_sheets(path)
	out <- lapply(sh, function(s) suppressMessages(readxl::read_excel(path, sheet = s)))
	names(out) <- sh
	out
}

ems120_prefix_xlsx_sheets <- function(x, prefix) {
	if (!length(x)) return(list())
	out <- list()
	used <- character()
	for (i in seq_along(x)) {
		nm0 <- paste0(prefix, "_", names(x)[i])
		nm0 <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", nm0)
		nm0 <- gsub("\\s+", "_", nm0)
		nm0 <- substr(nm0, 1, 31)
		nm <- nm0
		k <- 1L
		while (nm %in% used) {
			k <- k + 1L
			suf <- paste0("_", k)
			nm <- paste0(substr(nm0, 1, max(1, 31 - nchar(suf))), suf)
		}
		used <- c(used, nm)
		out[[nm]] <- x[[i]]
	}
	out
}

ems120_merge_xlsx_sources <- function(source_list, note) {
	out <- list()
	for (nm in names(source_list)) {
		out <- c(out, ems120_prefix_xlsx_sheets(source_list[[nm]], nm))
	}
	out$configuration <- tibble(
		final_figure_reassembly = TRUE,
		note = note
	)
	out
}

legacy_wb <- list(
	S1 = ems120_read_xlsx_all("FigS1.out.xlsx"),
	S2 = ems120_read_xlsx_all("FigS2.out.xlsx"),
	S3 = ems120_read_xlsx_all("FigS3.out.xlsx"),
	S4 = ems120_read_xlsx_all("FigS4.out.xlsx"),
	S5 = ems120_read_xlsx_all("FigS5.out.xlsx"),
	S6 = ems120_read_xlsx_all("FigS6.out.xlsx"),
	S7 = ems120_read_xlsx_all("FigS7.out.xlsx"),
	S8 = ems120_read_xlsx_all("FigS8.out.xlsx"),
	S9 = ems120_read_xlsx_all("FigS9.out.xlsx"),
	S10 = ems120_read_xlsx_all("FigS10.out.xlsx"),
	S11 = ems120_read_xlsx_all("FigS11.out.xlsx"),
	S12 = ems120_read_xlsx_all("FigS12.out.xlsx")
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Build the TEST-only post-2020 robustness panel used by official FigS8
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fig8_robust_available <- FALSE
fig8_robust_counts <- tibble()
fig8_robust_res <- tibble()
pS8D_final <- NULL

if (exists("figS11_adj_res") &&
	all(vapply(dat1.list, function(d) is.null(d) || !nrow(d) ||
		all(c("dx.type1.confidence", "phone.repeated_gt5") %in% names(d)), logical(1)))) {

	fig8_focus <- intersect(
		c("Violence", "Intoxication", "Psychiatric", "Traffic",
		  "CVD", "Urinary", "Respiratory", "Ob/Gyn"),
		dxs.all
	)
	fig8_focus_labels <- rev(dx_to_eng(fig8_focus))

	fig8_robust_count_one <- function(d0, y, analysis) {
		if (is.null(d0) || !nrow(d0) || !y %in% 2020:2024) return(tibble())
		kw_map <- setNames(map_grp$dx_grp, map_grp$dx_raw)
		d <- tibble(
			phone_group = as.character(d0$phone.luck),
			main_group = as.character(d0$dx_grp),
			keyword_group = unname(kw_map[trimws(as.character(d0$dx.type0))]),
			confidence = suppressWarnings(as.numeric(d0$dx.type1.confidence)),
			repeated = dplyr::coalesce(as.logical(d0$phone.repeated_gt5), FALSE)
		) %>%
			filter(phone_group %in% c("low", "high"))

		if (analysis == "Exclude repeated phone >5/year") d <- d %>% filter(!repeated)
		if (analysis == "MacBERT confidence >=0.80") {
			d <- d %>% filter(is.finite(confidence), confidence >= .80)
		}
		d$phenotype <- if (analysis == "Keyword phenotype") d$keyword_group else d$main_group
		d <- d %>% filter(phenotype %in% dxs.all)
		if (!nrow(d)) return(tibble())

		tot <- d %>% count(phone_group, name = "total")
		cas <- d %>% filter(phenotype %in% fig8_focus) %>% count(phone_group, phenotype, name = "case")

		tidyr::crossing(phone_group = c("low", "high"), phenotype = fig8_focus) %>%
			left_join(tot, by = "phone_group") %>%
			left_join(cas, by = c("phone_group", "phenotype")) %>%
			mutate(
				case = replace_na(case, 0L),
				total = replace_na(total, 0L),
				year = as.integer(y),
				analysis = analysis
			)
	}

	fig8_robust_analyses <- c(
		"Primary MacBERT",
		"Exclude repeated phone >5/year",
		"MacBERT confidence >=0.80",
		"Keyword phenotype"
	)

	fig8_robust_counts <- bind_rows(lapply(fig8_robust_analyses, function(a) {
		bind_rows(lapply(2020:2024, function(y) {
			fig8_robust_count_one(dat1.list[[as.character(y)]], y, a)
		}))
	}))

	fig8_fit_robust <- function(dd) {
		dd <- dd %>%
			mutate(high = as.integer(phone_group == "high")) %>%
			filter(total > 0)

		if (nrow(dd) < 4 || sum(dd$case) < 20 || n_distinct(dd$high) < 2) {
			return(tibble(
				OR = NA_real_, lo = NA_real_, hi = NA_real_,
				p = NA_real_, n = sum(dd$total)
			))
		}

		fit <- tryCatch(
			glm(cbind(case, total - case) ~ high + factor(year),
				family = binomial(), data = dd),
			error = function(e) NULL
		)
		if (is.null(fit)) {
			return(tibble(
				OR = NA_real_, lo = NA_real_, hi = NA_real_,
				p = NA_real_, n = sum(dd$total)
			))
		}

		td <- broom::tidy(fit) %>% filter(term == "high")
		if (!nrow(td)) {
			return(tibble(
				OR = NA_real_, lo = NA_real_, hi = NA_real_,
				p = NA_real_, n = sum(dd$total)
			))
		}

		tibble(
			OR = exp(td$estimate[1]),
			lo = exp(td$estimate[1] - 1.96 * td$std.error[1]),
			hi = exp(td$estimate[1] + 1.96 * td$std.error[1]),
			p = td$p.value[1],
			n = sum(dd$total)
		)
	}

	fig8_robust_res <- fig8_robust_counts %>%
		group_by(analysis, phenotype) %>%
		group_modify(~ fig8_fit_robust(.x)) %>%
		ungroup() %>%
		group_by(analysis) %>%
		mutate(
			p_adj = p.adjust(p, "BH"),
			sig05 = sig_star05(p_adj)
		) %>%
		ungroup() %>%
		mutate(
			analysis = factor(analysis, levels = fig8_robust_analyses),
			phenotype_label = factor(dx_to_eng(phenotype), levels = fig8_focus_labels)
		)

	pS8D_final <- ggplot(
		fig8_robust_res,
		aes(OR, phenotype_label, xmin = lo, xmax = hi, shape = analysis)
	) +
		geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
		geom_errorbar(
			position = position_dodge(width = .64),
			width = .14, orientation = "y", linewidth = .5
		) +
		geom_point(position = position_dodge(width = .64), size = 2.0) +
		scale_x_log10(breaks = c(.8, .9, 1, 1.1, 1.2, 1.3)) +
		labs(
			title = "d. Post-2020 robustness to data and phenotype definitions",
			x = "Year-adjusted OR: high vs low",
			y = NULL,
			shape = NULL
		) +
		fig_theme(base_size = 8.5) +
		theme(
			legend.position = "top",
			legend.text = element_text(size = 6.8),
			plot.title = element_text(face = "bold", hjust = .5, size = 10)
		)

	fig8_robust_available <- TRUE
}

if (!isTRUE(fig8_robust_available)) {
	pS8D_final <- empty_panel(
		"d. Post-2020 robustness to data and phenotype definitions",
		"Required confidence/repeated-phone fields unavailable"
	)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Working supplementary-figure assembly after visual review
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Requested working numbering:
#   FigS0  Annual AI-reconstructed phenotype composition (former FigS5)
#   FigS1  On-scene time and circular hourly mechanism (former FigS2)
#   FigS2  Phone-score construction (former FigS1)
#   FigS3  removed: its MacBERT OOF information is now integrated into main Fig2
#   FigS4  Geographic phenotype validation
#   FigS5  PHSM companion to main Fig5 (former annual composition moved to S0)
#   FigS6-FigS12 retain their current numbering so manuscript cross-checking is easy.
# Only plot assembly / styling changes below; analytical objects are reused.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Clear stale staging files whose final number has changed.  legacy_wb has
# already captured their workbook content above.
unlink(c(paste0("FigS", 0:13, ".png"), paste0("FigS", 0:13, ".out.xlsx")), force = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS0. Annual AI-reconstructed EMS phenotype composition (former FigS5)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pS0_final <- pS1 +
	labs(title = "Annual AI-reconstructed EMS phenotype composition") +
	theme(
		panel.grid = element_blank(),
		plot.title = element_text(size = 11, face = "bold", hjust = .5)
	)

save_plot(pS0_final, "FigS1.png", width = 8.1, height = 6.8, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(
		list(annual = legacy_wb$S1),
		"FigS1 is the annual AI-reconstructed EMS phenotype composition."
	),
	"FigS1.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS1. EMS on-scene time and circular hourly mechanism (former FigS2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This figure contains the headline high-vs-low on-scene-time result in panel b.
save_plot(Fig4, "FigS3.png", width = fig4_width, height = 12.6, dpi = 600, bg = "white")
writexl::write_xlsx(figS2_main_wb, "FigS3.out.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS2. Phone-score construction (former FigS1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pS2A_final <- pS2A + labs(title = "a. sco0 (pattern rule)")
pS2B_final <- pS2B + labs(title = "b. sco1 (simple rule)")
pS2C_final <- pS2C + labs(title = "c. sco2 (advanced rule)")
pS2D_final <- pS2D +
	labs(title = "d. Top-tail curves by score construction") +
	guides(color = guide_legend(nrow = 1)) +
	theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.key.width = grid::unit(16, "pt"))

FigS2_final <- (
	(pS2A_final | pS2B_final) /
	(pS2C_final | pS2D_final)
) + plot_layout(heights = c(1, 1))

save_plot(FigS2_final, "FigS4.png", width = 9.44, height = 8.4, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(
		list(score = legacy_wb$S2),
		"FigS4 contains phone-score construction and distribution diagnostics."
	),
	"FigS4.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS4. Confusion structure and largest error pathways
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pS4A_final <- p2NLP_mainC +
	labs(title = "a. OOF confusion matrix", subtitle = NULL) +
	scale_fill_gradientn(colors = c("#2CA25F", "#2B8CBE", "#F28E2B"), values = scales::rescale(c(0, .12, 1)), limits = c(0, 1), labels = scales::percent_format(accuracy = 1), name = "Row %", na.value = "grey96") +
	theme(axis.text.x = element_text(angle = 52, hjust = 1, vjust = 1, size = 6.7, face = "bold"), axis.text.y = element_text(size = 6.8, face = "bold"), legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"))
pS4B_final <- ggplot(cv_top_errors, aes(row_pct, path, fill = row_pct)) +
	geom_col(width = .68) +
	geom_text(aes(label = sprintf("%.1f%% (n=%s)", 100 * row_pct, scales::comma(n))), hjust = -.10, size = 2.55, fontface = "bold", color = "grey18") +
	scale_fill_gradient(low = "#D7ECFF", high = "#005A9C", guide = "none") +
	scale_x_continuous(limits = c(0, max(.08, max(cv_top_errors$row_pct, na.rm = TRUE) * 1.30)), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, .02))) +
	labs(title = "b. Largest off-diagonal pathways", x = "Row-normalised error share", y = NULL) +
	fig2_pub_theme(base_size = 7.8) +
	theme(axis.text = element_text(face = "bold"), axis.text.y = element_text(size = 6.8, face = "bold"), axis.title = element_text(face = "bold"), plot.subtitle = element_blank())
FigS4_final <- pS4A_final | pS4B_final

save_plot(FigS4_final, "FigS5.png", width = 13.2, height = 7.0, dpi = 600, bg = "white")
writexl::write_xlsx(
	list(
		panel_a_confusion_long = fig2nlp_cm_long,
		panel_a_confusion_matrix = fig2nlp_cm0,
		panel_b_error_pathways = cv_top_errors,
		configuration = tibble(
			figure = "FigS4. OOF confusion structure and largest error pathways",
			layout = "1 row x 2 equal-width columns",
			panel_order = "a OOF confusion matrix; b largest off-diagonal pathways",
			note = "Only data used by the two plotted panels are included."
		)
	),
	"FigS5.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS5. Geographic phenotype validation (former FigS4)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pS4A_geo <- pS3GeoA + labs(title = "a. geo.type1 composition by year")
pS4B_geo <- pS3GeoB +
	labs(title = sub("^b\\.\\s*", "b. ", figS3_geo_heat_title), x = "Recorded address type", y = "Predicted geographic type") +
	theme(axis.text = element_text(face = "bold"), axis.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))
figS4_geo_name_map <- tibble(
	geo_type1 = figS3_cloud_levels,
	geo_label_en = dplyr::case_when(
		geo_type1 == "\u4f4f\u5b85\u533a" ~ "Residential area",
		geo_type1 == "\u5de5\u4f5c\u533a" ~ "Workplace",
		geo_type1 == "\u516c\u5171\u573a\u6240" ~ "Public place",
		geo_type1 == "\u8857\u9053" ~ "Street",
		geo_type1 == "\u5bbe\u9986(\u9152\u5e97)" ~ "Hotel",
		geo_type1 == "\u8bca\u6240" ~ "Clinic",
		geo_type1 == "\u533b\u7597\u673a\u6784" ~ "Healthcare facility",
		geo_type1 == "\u4ea4\u901a\u573a\u6240" ~ "Transportation area",
		geo_type1 == "\u6559\u80b2\u673a\u6784" ~ "Educational facility",
		geo_type1 == "\u5176\u4ed6" ~ "Other",
		TRUE ~ paste("Geographic type", seq_along(figS3_cloud_levels))
	)
)
figS4_geo_terms_en <- figS3_geo_terms %>%
	left_join(figS4_geo_name_map, by = "geo_type1") %>%
	mutate(geo_label_en = factor(geo_label_en, levels = figS4_geo_name_map$geo_label_en))
figS4_geo_labels_en <- figS4_geo_terms_en %>% distinct(geo_label_en) %>% mutate(x = .5, y = -.035)
pS4C_geo <- ggplot(figS4_geo_terms_en, aes(x, y, label = term, size = size2, color = geo_label_en)) +
	geom_text(fontface = "bold", alpha = .90, check_overlap = TRUE) +
	geom_text(data = figS4_geo_labels_en, aes(x = x, y = y, label = geo_label_en), inherit.aes = FALSE, fontface = "bold", size = 3.35, color = "black") +
	facet_wrap(~geo_label_en, ncol = 3) +
	scale_size_identity() +
	coord_cartesian(xlim = c(.01, .99), ylim = c(-.07, .98), expand = FALSE, clip = "off") +
	labs(title = "c. Address-keyword word clouds by geographic type") +
	theme_void(base_size = 8.8) +
	theme(strip.text = element_blank(), legend.position = "none", plot.title = element_text(size = 10.2, face = "bold", hjust = .5), panel.spacing = grid::unit(.18, "lines"), plot.margin = margin(3, 0, 3, 0))

FigS4_final <- (
	(pS4A_geo | pS4B_geo) /
	pS4C_geo
) + plot_layout(heights = c(.88, 1.35), widths = c(1.0, 1.08))

save_plot(FigS4_final, "FigS6.png", width = 11.8, height = 11.4, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(
		list(geo = legacy_wb$S3),
		"FigS6 contains geo.type1 temporal composition, raw-address concordance, and address-keyword validation."
	),
	"FigS6.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS6. March 2022 PHSM companion to main Fig5 (former FigS5)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main Fig5 is intentionally reopening-focused.  Keep the earlier PHSM results
# visible in the supplement, but separate the overall policy-period response from
# the high-vs-low interaction rather than combining both numbers in one label.
pS5A_final <- p5A +
	labs(title = "a. March 2022 PHSM: total EMS demand", y = "7-day rolling index") +
	theme(plot.title = element_text(face = "bold", hjust = 0, size = 10.5))

figS5_phsm_plot <- fig5B_show %>%
	mutate(
		dx_label = factor(dx_to_eng(as.character(dx_grp)), levels = rev(dx_to_eng(dxs.all))),
		policy_lab = ifelse(is.finite(RR), sprintf("%.2f%s", RR, sig05), ""),
		did_lab = ifelse(is.finite(did_RR), sprintf("%.2f%s", did_RR, did_sig05), "")
	)

pS5B_final <- ggplot(figS5_phsm_plot, aes(RR, dx_label, xmin = lo, xmax = hi)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = .5) +
	geom_errorbar(width = .14, orientation = "y", linewidth = .58, color = "grey35") +
	geom_point(size = 2.15, color = "#2166AC") +
	geom_text(aes(label = policy_lab, x = pmin(hi * 1.06, 8.6)), hjust = 0, size = 2.0, fontface = "bold") +
	scale_x_log10(breaks = c(.3, .4, .6, .8, 1, 1.5, 2, 3, 4, 6, 8, 10), limits = c(.30, 10)) +
	labs(title = "b. Phenotype response during PHSM", x = "RR vs pre-PHSM reference", y = NULL) +
	theme_classic(base_size = 8.3) +
	theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7), plot.title = element_text(face = "bold", hjust = 0, size = 10.4))

pS5C_final <- ggplot(figS5_phsm_plot, aes(did_RR, dx_label, xmin = did_lo, xmax = did_hi)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = .5) +
	geom_errorbar(width = .14, orientation = "y", linewidth = .58, color = "grey35", na.rm = TRUE) +
	geom_point(size = 2.15, color = "#B2182B", na.rm = TRUE) +
	geom_text(aes(label = did_lab, x = pmin(did_hi * 1.05, 3.0)), hjust = 0, size = 2.0, fontface = "bold", na.rm = TRUE) +
	scale_x_log10(breaks = c(.25, .4, .6, .8, 1, 1.5, 2, 3), limits = c(.25, 3.2)) +
	labs(title = "c. Differential response: high vs low score", x = "H/L x PHSM interaction RR", y = NULL) +
	theme_classic(base_size = 8.3) +
	theme(panel.grid = element_blank(), axis.text.y = element_text(size = 7), plot.title = element_text(face = "bold", hjust = 0, size = 10.4))

FigS5_final <- align_panel_rows(
	rows = list(
		list(pS5A_final),
		list(pS5B_final, pS5C_final)
	),
	rel_widths = list(1, c(1, 1)),
	rel_heights = c(.78, 1.25),
	row_gap = .08,
	side_pad = .018,
	axis_text_size = 8.3,
	axis_title_size = 8.3
)
save_plot(FigS5_final, "FigS12.png", width = 12.2, height = 9.2, dpi = 600, bg = "white")
writexl::write_xlsx(list(
	PHSM_total_event_index = phsm_total,
	PHSM_total_policy_models = fig5A_policy,
	PHSM_total_DID_models = fig5A_test,
	PHSM_policy_global_mix_test = fig5B_policy_global,
	PHSM_global_mix_test = fig5B_global,
	PHSM_policy_models = fig5B_policy_dat,
	PHSM_DID_models = fig5B_dat,
	PHSM_event_window_daily = figS12_phsm_daily,
	PHSM_event_window_DID = figS12_phsm_did,
	configuration = tibble(note = "FigS12 is the PHSM companion to reopening-focused main Fig5 and is placed last to match manuscript order. No PHSM model is refit in this final plotting block.")
), "FigS12.out.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS7. Transformer phenotype overview (former FigS6)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pS6A_final <- pS5A + labs(title = "a. Transformer raw disease-category trends")
pS6A_final <- pS6A_final + theme(legend.spacing.y = grid::unit(7, "pt"), legend.key.height = grid::unit(13, "pt"), legend.text = element_text(face = "bold"))
pS6B_final <- pS5C +
	labs(title = "b. Keyword-token word clouds by EMS phenotype") +
	coord_cartesian(xlim = c(.01, .99), ylim = c(-.08, .99), expand = FALSE, clip = "off") +
	theme(panel.spacing = grid::unit(.20, "lines"), plot.margin = margin(2, 0, 2, 0))
FigS6_final <- pS6A_final / patchwork::plot_spacer() / pS6B_final + plot_layout(heights = c(.72, .10, 1.58))

save_plot(FigS6_final, "FigS7.png", width = 11.4, height = 12.8, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(
		list(transformer = legacy_wb$S5),
		"FigS7 contains the transformer phenotype overview; it was formerly FigS6. Keyword-vs-transformer concordance is workbook-only."
	),
	"FigS7.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS8. Independent validation of MacBERT-identified Death calls (former FigS7)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The former three-panel figure was difficult to interpret.  The rule-family
# matrix remains in the workbook, while the plotted figure now focuses on the
# two questions that matter: (1) independent textual support over time and
# (2) what the keyword comparator calls the same MacBERT-Death records.
figS7_support_dat <- figS6B_plot_dat %>%
	filter(metric %in% c("Keyword confirms Death", "Death terms in raw text"))

pS7A_final <- ggplot(figS7_support_dat, aes(year, pct, color = metric, group = metric)) +
	geom_hline(yintercept = .50, linetype = "dotted", color = "grey65", linewidth = .45) +
	geom_line(linewidth = .90) +
	geom_point(size = 1.9) +
	scale_color_manual(values = c("Keyword confirms Death" = "#2166AC", "Death terms in raw text" = "#B2182B")) +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, .02))) +
	scale_x_continuous(breaks = years) +
	labs(title = "a. Independent support for MacBERT-identified Death calls", x = NULL, y = "% of MacBERT-Death calls", color = NULL) +
	theme_classic(base_size = 9) +
	theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", legend.justification = "left", plot.title = element_text(face = "bold", hjust = 0, size = 10.4))

figS7_top_codes <- figS6A_dat %>%
	group_by(dx0_grp) %>% summarise(N = sum(n), .groups = "drop") %>%
	arrange(desc(N)) %>% slice_head(n = 5) %>% pull(dx0_grp) %>% as.character()
figS7_mix_dat <- figS6A_dat %>%
	mutate(
		kw_code = as.character(dx0_grp),
		keyword_group = ifelse(kw_code %in% figS7_top_codes, dx_to_eng(kw_code), "All remaining")
	) %>%
	group_by(year, keyword_group) %>% summarise(n = sum(n), .groups = "drop") %>%
	group_by(year) %>% mutate(pct = n / sum(n)) %>% ungroup()
figS7_group_levels <- c(dx_to_eng(figS7_top_codes), "All remaining")
figS7_mix_dat <- figS7_mix_dat %>% mutate(keyword_group = factor(keyword_group, levels = rev(figS7_group_levels)))
figS7_cols <- c(setNames(unname(dxs.all.color[figS7_top_codes]), dx_to_eng(figS7_top_codes)), "All remaining" = "grey82")

pS7B_final <- ggplot(figS7_mix_dat, aes(year, pct, fill = keyword_group)) +
	geom_col(width = .82, color = "white", linewidth = .15) +
	geom_text(aes(label = ifelse(as.character(keyword_group) == "Death" & pct >= .04, scales::percent(pct, accuracy = 1), "")), position = position_stack(vjust = .5), size = 2.25, fontface = "bold") +
	scale_fill_manual(values = figS7_cols, drop = FALSE, name = "Keyword comparator") +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) +
	scale_x_continuous(breaks = years) +
	labs(title = "b. Keyword assignments within MacBERT-Death calls", x = NULL, y = "% of MacBERT-Death calls") +
	theme_classic(base_size = 9) +
	theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right", legend.text = element_text(size = 7), plot.title = element_text(face = "bold", hjust = 0, size = 10.4))

# Remove the former left-hand line-only support panel; retain the composition panel.
FigS7_final <- pS7B_final + labs(title = "Keyword assignments within MacBERT-identified Death calls")
save_plot(FigS7_final, "FigS8.png", width = 7.2, height = 6.2, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(
		list(death = legacy_wb$S6),
		"FigS8 retains only the keyword-assignment composition panel from former FigS7; the left-hand line-only panel was removed."
	),
	"FigS8.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS9. Descriptive temporal EMS spectrum (former FigS8)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pS8A_final <- pFigS7A + labs(title = "a. Monthly composition of the EMS spectrum")
pS8B_final <- pFigS7B + labs(title = "b. Hourly distribution of EMS phenotypes during six COVID-era years")
FigS8_final <- pS8A_final / pS8B_final + plot_layout(heights = c(1.0, 1.65))

save_plot(FigS8_final, "FigS2.png", width = 10.5, height = 11.0, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(list(temporal = legacy_wb$S7), "FigS2 contains the descriptive monthly and hourly temporal EMS spectrum."),
	"FigS2.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS10. High-vs-low circadian and on-scene-time contrasts (former FigS9)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (FALSE) {
pS9A_final <- pS8A + labs(title = "a. High- and low-score circadian profiles")
pS9B_final <- pS8B + labs(title = "b. Corrected time-block difference: high minus low")
pS9C_final <- pS9 +
	labs(title = "c. Circular on-scene-duration contrast by phone-score group") +
	theme(plot.title = element_text(size = 11.2, face = "bold", hjust = .5), plot.margin = margin(0, 0, 0, 0))
FigS9_final <- ((pS9A_final | pS9B_final) / pS9C_final) + plot_layout(heights = c(.62, 1.78))

save_plot(FigS9_final, "FigS10.png", width = 13.4, height = 14.8, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(
		list(circadian = legacy_wb$S8, onsite = legacy_wb$S9),
		"FigS10 combines circadian contrasts with the circular on-scene-duration analysis; it was formerly FigS9."
	),
	"FigS10.out.xlsx"
)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS11. Robustness to phone-score definition (former FigS10)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Robustness is a matrix question: do phenotype effects retain their direction
# across alternative score definitions? Heatmaps show that pattern more directly
# and with much less visual clutter than four over-plotted forest series.
figS10_tail_plot <- sens_tail_summary %>%
	filter(disease %in% dxs.all, is.finite(mean_RR), mean_RR > 0) %>%
	mutate(
		effect = log2(mean_RR),
		cell = sprintf("%.2f%s", mean_RR, sig05_meta),
		score_label = factor(score_label, levels = phone_score_labels),
		disease_label = factor(dx_to_eng(disease), levels = rev(dx_to_eng(dxs.all)))
	)
figS10_cont_plot <- sens_cont_summary %>%
	filter(disease %in% dxs.all, is.finite(mean_OR), mean_OR > 0) %>%
	mutate(
		effect = log2(mean_OR),
		cell = sprintf("%.2f%s", mean_OR, sig05_meta),
		score_label = factor(score_label, levels = phone_score_labels),
		disease_label = factor(dx_to_eng(disease), levels = rev(dx_to_eng(dxs.all)))
	)
figS10_tail_lim <- max(abs(figS10_tail_plot$effect), na.rm = TRUE); if (!is.finite(figS10_tail_lim) || figS10_tail_lim < 1e-6) figS10_tail_lim <- .25
figS10_cont_lim <- max(abs(figS10_cont_plot$effect), na.rm = TRUE); if (!is.finite(figS10_cont_lim) || figS10_cont_lim < 1e-6) figS10_cont_lim <- .08

pS10A_final <- ggplot(figS10_tail_plot, aes(score_label, disease_label, fill = effect)) +
	geom_tile(color = "white", linewidth = .28) +
	geom_text(aes(label = cell), size = 2.05, fontface = "bold") +
	facet_grid(. ~ contrast) +
	scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-figS10_tail_lim, figS10_tail_lim), oob = scales::squish, name = "log2 RR") +
	labs(title = "a. Tail-cutoff robustness", x = NULL, y = NULL) +
	theme_minimal(base_size = 8.2) +
	theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 28, hjust = 1, size = 7, face = "bold"), axis.text.y = element_text(size = 7, face = "bold"), axis.title = element_text(face = "bold"), legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0, size = 10.5))

pS10B_final <- ggplot(figS10_cont_plot, aes(score_label, disease_label, fill = effect)) +
	geom_tile(color = "white", linewidth = .28) +
	geom_text(aes(label = cell), size = 2.10, fontface = "bold") +
	scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-figS10_cont_lim, figS10_cont_lim), oob = scales::squish, name = "log2 OR") +
	labs(title = "b. Continuous-score robustness (per 1-SD)", x = NULL, y = NULL) +
	theme_minimal(base_size = 8.2) +
	theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 28, hjust = 1, size = 7, face = "bold"), axis.text.y = element_text(size = 7, face = "bold"), axis.title = element_text(face = "bold"), legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0, size = 10.5))

FigS10_final <- (pS10A_final | pS10B_final) + plot_layout(widths = c(1.55, 1.0))
save_plot(FigS10_final, "FigS9.png", width = 11.4, height = 7.5, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(
		list(score_robustness = legacy_wb$S4),
		"FigS9 displays phone-score robustness analyses as effect heatmaps."
	),
	"FigS9.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS12. Period-level and quarterly gradient dynamics (former FigS11)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pS11A_final <- pS10A + labs(title = "a. Period-level high-vs-low RR")
pS11B_final <- pS10B + labs(title = "b. Change relative to 2017-2019")
pS11C_final <- pS10C + labs(title = "c. Quarterly high-vs-low phone-score gradients")
FigS11_final <- ((pS11A_final / pS11B_final) | pS11C_final) + plot_layout(widths = c(1.04, 1.0))

save_plot(FigS11_final, "FigS10.png", width = 15.2, height = 13.0, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(list(period_gradient = legacy_wb$S10), "FigS10 contains period-level heatmaps and quarterly gradient dynamics."),
	"FigS10.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS13. Dose-response, adjustment, and post-2020 robustness (former FigS12)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Panels c-d are converted from multi-series forest plots to effect matrices.
# The analyses are unchanged; the new layout makes the stability question visible.
pS12A_final <- pS11A +
	labs(title = "a. Disease membership across phone-score deciles") +
	theme(panel.grid = element_blank(), plot.title = element_text(size = 10.5, face = "bold", hjust = 0))
pS12B_final <- pS11B +
	labs(title = "b. Linear decile trend") +
	theme(panel.grid = element_blank(), plot.title = element_text(size = 10.5, face = "bold", hjust = 0))

figS12_adj_plot <- figS11_adj_res %>%
	filter(is.finite(OR), OR > 0) %>%
	mutate(
		effect = log2(OR),
		cell = sprintf("%.2f%s", OR, sig05),
		model_short = factor(as.character(model), levels = names(figS11_adj_models)),
		dx_label = factor(dx_to_eng(dx_grp), levels = rev(dx_to_eng(dxs.all)))
	)
figS12_adj_lim <- max(abs(figS12_adj_plot$effect), na.rm = TRUE); if (!is.finite(figS12_adj_lim) || figS12_adj_lim < 1e-6) figS12_adj_lim <- .25
pS12C_final <- ggplot(figS12_adj_plot, aes(model_short, dx_label, fill = effect)) +
	geom_tile(color = "white", linewidth = .28) +
	geom_text(aes(label = cell), size = 2.05, fontface = "bold") +
	scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-figS12_adj_lim, figS12_adj_lim), oob = scales::squish, name = "log2 OR") +
	labs(title = "c. Sequential adjustment", x = NULL, y = NULL) +
	theme_minimal(base_size = 8.2) +
	theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 24, hjust = 1, size = 7, face = "bold"), axis.text.y = element_text(size = 7, face = "bold"), axis.title = element_text(face = "bold"), legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0, size = 10.5))

figS12_robust_plot <- fig8_robust_res %>%
	filter(is.finite(OR), OR > 0) %>%
	mutate(
		effect = log2(OR),
		cell = sprintf("%.2f%s", OR, sig05),
		analysis_short = dplyr::recode(as.character(analysis),
			"Primary MacBERT" = "Primary",
			"Exclude repeated phone >5/year" = "No repeated\nphone >5/y",
			"MacBERT confidence >=0.80" = "Confidence\n>=0.80",
			"Keyword phenotype" = "Keyword\nphenotype"
		),
		analysis_short = factor(analysis_short, levels = c("Primary", "No repeated\nphone >5/y", "Confidence\n>=0.80", "Keyword\nphenotype"))
	)
figS12_robust_lim <- max(abs(figS12_robust_plot$effect), na.rm = TRUE); if (!is.finite(figS12_robust_lim) || figS12_robust_lim < 1e-6) figS12_robust_lim <- .25
pS12D_final <- if (nrow(figS12_robust_plot)) {
	ggplot(figS12_robust_plot, aes(analysis_short, phenotype_label, fill = effect)) +
		geom_tile(color = "white", linewidth = .28) +
		geom_text(aes(label = cell), size = 2.15, fontface = "bold") +
		scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-figS12_robust_lim, figS12_robust_lim), oob = scales::squish, name = "log2 OR") +
		labs(title = "d. Post-2020 data/phenotype robustness", x = NULL, y = NULL) +
		theme_minimal(base_size = 8.2) +
		theme(panel.grid = element_blank(), axis.text.x = element_text(size = 7, face = "bold"), axis.text.y = element_text(size = 7, face = "bold"), axis.title = element_text(face = "bold"), legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0, size = 10.5))
} else empty_panel("d. Post-2020 data/phenotype robustness", "Required sensitivity fields unavailable")

FigS12_final <- (
	(pS12A_final | pS12B_final) /
	(pS12C_final | pS12D_final)
) + plot_layout(heights = c(1.02, .98), widths = c(1.18, 1.0))

save_plot(FigS12_final, "FigS11.png", width = 13.2, height = 10.8, dpi = 600, bg = "white")
figS12_wb <- ems120_merge_xlsx_sources(
	list(dose_adjustment = legacy_wb$S11),
	"FigS11 combines dose-response, sequential adjustment, and post-2020 robustness. Panels c-d are visualized as effect matrices."
)
figS12_wb$post2020_sens_counts <- fig8_robust_counts
figS12_wb$post2020_sens_models <- fig8_robust_res
writexl::write_xlsx(figS12_wb, "FigS11.out.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Legacy PHSM/reopening export layouts are retained below for reference only.
# They are disabled because FigS11 and FigS12 are now occupied by the shifted
# former FigS9 and FigS10 outputs requested in the final numbering.
if (FALSE) {
# 🚩 Legacy FigS11. PHSM event-window analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reopening is no longer stacked underneath PHSM; it is moved to FigS12.
pS11A_final <- p5A + labs(title = "a. PHSM: total EMS demand")
pS11B_final <- p5B + labs(title = "b. PHSM: phenotype-specific response")
pS11C_final <- pS12A + labs(title = "c. March 2022 PHSM daily counts")
pS11D_final <- pS12B + labs(title = "d. PHSM DID: pre vs during")
pS11E_final <- pS12C + labs(title = "e. PHSM DID: during vs post")

FigS11_final <- (
	(pS11A_final | pS11B_final) /
	pS11C_final /
	(pS11D_final | pS11E_final)
) + plot_layout(heights = c(.78, 1.28, .78))

save_plot(FigS11_final, "FigS11.png", width = 15.2, height = 13.8, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(
		list(policy = legacy_wb$S12),
		"FigS11 contains PHSM only. The reopening analyses formerly shown as FigS9f-h are separated into FigS12."
	),
	"FigS11.out.xlsx"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 FigS12. Late-2022 reopening analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Former FigS9f-h, now a dedicated reopening figure.
pS12A_final <- pS12D + labs(title = "a. Late-2022 reopening daily counts")
pS12B_final <- pS12E + labs(title = "b. Reopening DID: Nov 11 semi-removal")
pS12C_final <- pS12F + labs(title = "c. Reopening DID: Dec 7 full reopening")

FigS12_final <- (
	pS12A_final /
	(pS12B_final | pS12C_final)
) + plot_layout(heights = c(1.45, .75))

save_plot(FigS12_final, "FigS12.png", width = 12.2, height = 10.6, dpi = 600, bg = "white")
writexl::write_xlsx(
	ems120_merge_xlsx_sources(
		list(reopening = legacy_wb$S12),
		"FigS12 is the dedicated late-2022 reopening analysis formerly shown as FigS9f-h."
	),
	"FigS12.out.xlsx"
)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 TEST / obsolete companion outputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The former TEST.Fig1 content is integrated into main Fig2.
# TEST.Fig5 robustness is incorporated into official FigS10.
unlink(c("TEST.Fig1.png", "TEST.Fig1.xlsx", "TEST.Fig5.png", "TEST.Fig5.xlsx"), force = TRUE)

cat("EMS120 pipeline finished successfully.\n")
