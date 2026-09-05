# 🚩 CHNS: MAHA replication / validation pipeline
# Single-cohort replication with explicit CHNS inputs and CFCT FOODCODE maps.

dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
project_dir <- Sys.getenv("CHNS_PROJECT_DIR", unset = file.path(dir0, "data", "chns"))
raw_dir <- file.path(project_dir, "raw")
taobao_dir <- file.path(project_dir, "CHNS2.0-TaoB", "CHNS原始数据")
id2019_dir <- file.path(taobao_dir, "Master_ID_201908")
foodcode_dir <- file.path(project_dir, "MAHA_foodcode")
food_map_default <- file.path(foodcode_dir, "chns_foodcode_map.tsv")
pexam_public_default <- file.path(raw_dir, "individual", "pexampub12.sav")
nutrient_map_default <- file.path(foodcode_dir, "chns_foodcode_nutrients.tsv")

outdir <- Sys.getenv("MAHA_OUTDIR", unset = file.path(dir0, "analysis", "maha"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

pacman::p_load(tidyverse, haven, survival, broom, writexl, patchwork, cowplot, readxl, scales)
.this_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.this_file <- if (length(.this_arg)) sub("^--file=", "", .this_arg[[1]]) else file.path(dir0, "scripts", "maha", "f", "MAHA_chns.R")
.this_dir <- dirname(normalizePath(.this_file, winslash = "/", mustWork = FALSE))
source(file.path(.this_dir, "comm.f.R"))

cohort_prefix <- "chns"
maha_auxdir <- Sys.getenv("MAHA_AUXDIR", unset = file.path(outdir, cohort_prefix))
dir.create(maha_auxdir, recursive = TRUE, showWarnings = FALSE)
cohort_file <- function(path) file.path(maha_auxdir, basename(path))
pub_file <- function(path) file.path(outdir, basename(path))
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 360, bg = "white") {
  ggplot2::ggsave(cohort_file(filename), plot = plot, width = width, height = height, dpi = dpi, bg = bg)
}
save_pub_plot <- function(plot, filename, width = 8, height = 6, dpi = 360, bg = "white") {
  ggplot2::ggsave(pub_file(filename), plot = plot, width = width, height = height, dpi = dpi, bg = bg)
}
write_xlsx <- function(x, path, ...) writexl::write_xlsx(x, cohort_file(path), ...)
write_pub_xlsx <- function(x, path, ...) writexl::write_xlsx(x, pub_file(path), ...)
step_header <- function(x) cat("\n", strrep("=", 88), "\n", x, "\n", strrep("=", 88), "\n", sep = "")

analysis_spec <- "primary"
if (!identical(analysis_spec, "primary")) stop("CHNS primary model must remain primary in this publication pipeline.", call. = FALSE)

baseline_req <- trimws(Sys.getenv("CHNS_BASELINE_WAVE", unset = "2011"))
followup_req <- trimws(Sys.getenv("CHNS_FOLLOWUP_END", unset = "2015"))
min_followup_years <- suppressWarnings(as.integer(Sys.getenv("CHNS_MIN_FOLLOWUP_YEARS", unset = "4")))
min_age <- suppressWarnings(as.integer(Sys.getenv("CHNS_MIN_AGE", unset = "40")))
min_diet_days <- suppressWarnings(as.integer(Sys.getenv("CHNS_MIN_DIET_DAYS", unset = "2")))
maha_jobs <- suppressWarnings(as.integer(Sys.getenv("MAHA_JOBS", unset = ifelse(.Platform$OS.type == "windows", "1", "4"))))
null_max <- suppressWarnings(as.integer(Sys.getenv("MAHA_NULL_MAX", unset = "0")))
full_map_min_cov <- suppressWarnings(as.numeric(Sys.getenv("CHNS_MIN_MAP_COVERAGE", unset = "0.90")))
require_nutrient_map <- suppressWarnings(as.integer(Sys.getenv("CHNS_REQUIRE_NUTRIENT_MAP", unset = "1")))
require_enhanced_map <- suppressWarnings(as.integer(Sys.getenv("CHNS_REQUIRE_ENHANCED_MAP", unset = "1")))
run_mortality_sweep <- suppressWarnings(as.integer(Sys.getenv("CHNS_RUN_MORTALITY_SWEEP", unset = "1")))
mortality_sweep_child <- identical(Sys.getenv("CHNS_MORTALITY_SWEEP_CHILD", unset = "0"), "1")
mortality_sweep_wave_req <- trimws(Sys.getenv("CHNS_MORTALITY_SWEEP_WAVES", unset = "2004,2006,2009,2011"))
mortality_model_age <- suppressWarnings(as.integer(Sys.getenv("CHNS_MORTALITY_MODEL_AGE", unset = "40")))
mortality_event_target <- suppressWarnings(as.integer(Sys.getenv("CHNS_MORTALITY_EVENT_TARGET", unset = "500")))
if (!require_nutrient_map %in% c(0L,1L)) require_nutrient_map <- 1L
if (!require_enhanced_map %in% c(0L,1L)) require_enhanced_map <- 1L
if (!run_mortality_sweep %in% c(0L,1L)) run_mortality_sweep <- 1L
if (!is.finite(mortality_model_age) || mortality_model_age < 18L) mortality_model_age <- 40L
if (!is.finite(mortality_event_target) || mortality_event_target < 1L) mortality_event_target <- 500L

if (!is.finite(min_followup_years) || min_followup_years < 1) min_followup_years <- 4L
if (!is.finite(min_age) || min_age < 18) min_age <- 40L
if (!is.finite(min_diet_days) || min_diet_days < 2) min_diet_days <- 2L
if (!is.finite(maha_jobs) || maha_jobs < 1) maha_jobs <- 1L
if (!is.finite(null_max) || null_max < 0) null_max <- 0L
if (!is.finite(full_map_min_cov) || full_map_min_cov <= 0 || full_map_min_cov > 1) full_map_min_cov <- .90

log_file <- cohort_file("maha.log")
if (file.exists(log_file)) file.remove(log_file)
log_con <- file(log_file, "wt")
sink(log_con, split = TRUE)
on.exit({ while (sink.number() > 0) sink(); close(log_con) }, add = TRUE)
options(width = 220, warn = 1)


# 🚩 Utilities
num <- function(x) {
  if (inherits(x, "haven_labelled") || inherits(x, "labelled")) x <- haven::zap_labels(x)
  suppressWarnings(as.numeric(x))
}
recode_chns_school_years <- function(x) {
  z <- num(x)
  # Estimated years: 6 primary + 3 lower secondary + 3 upper secondary;
  # technical school follows lower secondary; college 6+ years is capped at 18.
  dplyr::case_when(
    z == 0       ~ 0,
    z %in% 11:16 ~ z - 10,
    z %in% 21:26 ~ z - 14,
    z %in% 27:29 ~ z - 17,
    z %in% 31:36 ~ z - 18,
    TRUE         ~ NA_real_
  )
}
cache_identical_inputs <- function(fn) {
  # Run-local only: reuse intervals only when all input data and the definition
  # are identical. Downstream models still fit independently on these intervals.
  entries <- list()
  function(d,mortality_definition) {
    for(entry in entries) {
      if(identical(entry$definition,mortality_definition) && identical(entry$data,d))
        return(entry$result)
    }
    message("[TIME UPDATED] Building intervals: ",mortality_definition,"; input rows=",nrow(d))
    result <- fn(d,mortality_definition)
    entries[[length(entries)+1L]] <<- list(data=d,definition=mortality_definition,result=result)
    message("[TIME UPDATED] Cached intervals: ",nrow(result)," rows")
    result
  }
}
zstd <- function(x) {
  x <- num(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
first_num <- function(x) {
  x <- num(x); x <- x[is.finite(x)]
  if (length(x)) x[1] else NA_real_
}
qscore <- function(x, reverse = FALSE) {
  x <- num(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (sum(ok) < 5 || length(unique(x[ok])) < 2) return(out)
  r <- dplyr::percent_rank(x[ok])
  out[ok] <- as.numeric(as.character(cut(r, c(-Inf,.2,.4,.6,.8,Inf),
                                            labels = c(0,2.5,5,7.5,10), right = TRUE)))
  if (reverse) out[ok] <- 10 - out[ok]
  out
}
qscore5 <- function(x, reverse = FALSE) {
  x <- num(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (sum(ok) < 5 || length(unique(x[ok])) < 2) return(out)
  r <- dplyr::percent_rank(x[ok])
  out[ok] <- as.numeric(as.character(cut(r, c(-Inf,.2,.4,.6,.8,Inf),
                                            labels = 1:5, right = TRUE)))
  if (reverse) out[ok] <- 6 - out[ok]
  out
}
qscore5_grouped <- function(x, group, reverse = FALSE) {
  x <- num(x); g <- as.character(group); out <- rep(NA_real_, length(x))
  for (lev in unique(g[!is.na(g)])) {
    ii <- which(g == lev & is.finite(x))
    if (length(ii)) out[ii] <- qscore5(x[ii], reverse = reverse)
  }
  out
}
energy_residual_grouped <- function(x, kcal, group) {
  x <- num(x); e <- num(kcal); g <- as.character(group); out <- rep(NA_real_, length(x))
  for (lev in unique(g[!is.na(g)])) {
    ii <- which(g == lev & is.finite(x) & is.finite(e))
    if (length(ii) < 20L || stats::sd(e[ii]) == 0) next
    fit <- stats::lm(x[ii] ~ e[ii])
    out[ii] <- stats::residuals(fit) + mean(x[ii], na.rm = TRUE)
  }
  out
}
energy_density <- function(x, kcal, multiplier = 1000) {
  x <- num(x); e <- num(kcal)
  ifelse(is.finite(x) & is.finite(e) & e > 0, x / e * multiplier, NA_real_)
}
score100 <- function(x) {
  x <- num(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (sum(ok) < 2 || diff(range(x[ok])) == 0) return(out)
  out[ok] <- (x[ok] - min(x[ok])) / diff(range(x[ok])) * 100
  out
}
f3c <- function(x) {
  x <- num(x); out <- rep(NA_character_, length(x)); ok <- is.finite(x)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2) return(factor(out, levels = c("low","middle","high")))
  r <- dplyr::percent_rank(x[ok])
  out[ok] <- as.character(cut(r, c(-Inf,1/3,2/3,Inf), labels = c("low","middle","high")))
  factor(out, levels = c("low","middle","high"))
}
rowmean_min <- function(d, vars, min_n = 1L) {
  vars <- vars[vars %in% names(d)]
  if (!length(vars)) return(rep(NA_real_, nrow(d)))
  m <- as.matrix(d[, vars, drop = FALSE])
  nok <- rowSums(is.finite(m))
  o <- rowMeans(m, na.rm = TRUE)
  o[nok < min_n] <- NA_real_
  o
}
rowmean_vecs <- function(..., min_n = 1L) {
  m <- cbind(...)
  nok <- rowSums(is.finite(m))
  o <- rowMeans(m, na.rm = TRUE)
  o[nok < min_n] <- NA_real_
  o
}
sum_known <- function(x) {
  x <- num(x)
  if (!any(is.finite(x))) return(NA_real_)
  sum(x, na.rm = TRUE)
}
weighted_sum_known <- function(g, flag) {
  g <- num(g); flag <- num(flag)
  ok <- is.finite(g) & is.finite(flag)
  if (!any(ok)) return(NA_real_)
  sum(g[ok] * flag[ok], na.rm = TRUE)
}
safe_cor <- function(x,y) {
  x <- num(x); y <- num(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}
yes1 <- function(x) {
  z <- num(x)
  ifelse(z == 1, 1L, ifelse(z %in% c(0,2), 0L, NA_integer_))
}
pref15 <- function(x) {
  # Prefer the embedded CHNS value labels rather than assuming numeric orientation.
  # This protects against releases in which the raw 1-5 coding is reversed.
  labs <- attr(x,"labels")
  z <- num(x); out <- rep(NA_real_,length(z))
  if(!is.null(labs) && length(labs)) {
    n <- tolower(names(labs)); val <- suppressWarnings(as.numeric(labs))
    rec <- rep(NA_real_,length(val))
    rec[grepl("dislike.*very|very.*dislike",n)] <- 1
    rec[grepl("dislike.*some|some.*dislike",n)] <- 2
    rec[grepl("neutral",n)] <- 3
    rec[grepl("like.*some|some.*like",n)] <- 4
    rec[grepl("like.*very|very.*like",n)] <- 5
    for(i in which(is.finite(rec))) out[z==val[i]] <- rec[i]
    # Explicit 'not eat' / non-consumption categories remain missing.
    if(any(is.finite(out))) return(out)
  }
  # Fallback consistent with the published CHNS food-preference analysis.
  out[z %in% 1:5] <- z[z %in% 1:5]
  out
}
moderate_alcohol_score <- function(g, female) {
  g <- num(g); female <- as.integer(female)
  out <- rep(NA_real_, length(g)); ok <- is.finite(g) & female %in% c(0L,1L); out[ok] <- 0
  out[ok & female == 1L & g >= 5 & g <= 15] <- 10
  out[ok & female == 0L & g >= 10 & g <= 25] <- 10
  out
}
good_covariates <- function(d, covs, min_complete = .50) {
  covs <- intersect(covs, names(d))
  covs[vapply(covs, function(v) {
    x <- d[[v]]
    mean(!is.na(x)) >= min_complete && length(unique(x[!is.na(x)])) > 1
  }, logical(1))]
}
parse_year_request <- function(x, label) {
  if (!nzchar(x) || tolower(x) == "auto") return(NA_integer_)
  z <- suppressWarnings(as.integer(x))
  if (!is.finite(z)) stop(label, " must be an integer year or 'auto': ", x, call. = FALSE)
  z
}
valid_years <- function(x) {
  z <- num(x)
  sort(unique(z[z >= 1989 & z <= 2035]))
}
pick <- function(d, cands, required = TRUE) {
  cands <- toupper(cands)
  nm <- cands[cands %in% names(d)]
  if (length(nm)) return(nm[1])
  if (required) stop("Required variable missing: ", paste(cands, collapse = " | "), call. = FALSE)
  NA_character_
}
colnum <- function(d, nm) {
  if (is.na(nm) || !nm %in% names(d)) rep(NA_real_, nrow(d)) else num(d[[nm]])
}
read_any <- function(f, col_select = NULL) {
  ext <- tolower(tools::file_ext(f))
  if (ext == "sas7bdat") {
    # CHNS SAS names use mixed case (for example, IDind), while downstream
    # processing standardizes all names to upper case. Select requested fields
    # case-insensitively so the pre-read projection does not silently drop them.
    d <- if (is.null(col_select)) haven::read_sas(f) else haven::read_sas(
      f, col_select = matches(paste0("^(?:", paste(col_select, collapse = "|"), ")$"), ignore.case = TRUE)
    )
  } else if (ext == "sav") {
    d <- haven::read_sav(f)
    if (!is.null(col_select)) d <- d %>% dplyr::select(any_of(col_select))
  } else if (ext == "xpt") {
    d <- haven::read_xpt(f)
    if (!is.null(col_select)) d <- d %>% dplyr::select(any_of(col_select))
  } else if (ext == "csv") {
    d <- readr::read_csv(f, show_col_types = FALSE)
  } else if (ext %in% c("tsv","txt")) {
    d <- readr::read_tsv(f, show_col_types = FALSE)
  } else if (ext %in% c("xlsx","xls")) {
    d <- readxl::read_excel(f)
  } else stop("Unsupported file type: ", f)
  d <- as_tibble(d); names(d) <- toupper(names(d)); d
}
print_source <- function(nm, x, required = TRUE) {
  ok <- is.character(x[["path"]]) && length(x[["path"]]) == 1 && !is.na(x[["path"]]) && file.exists(x[["path"]])
  cat(sprintf("  %-12s %-7s %-14s : %s%s\n", nm, if(ok)"OK" else if(required)"MISSING" else "OPTIONAL",
              x[["source"]], ifelse(is.na(x[["path"]]),"<not found>",x[["path"]]),
              if(ok)"" else if(required)" 文件不存在" else ""))
  invisible(ok)
}
linear_adequacy <- function(x, target, max_score) {
  x <- num(x); out <- pmin(pmax(x / target, 0), 1) * max_score
  out[!is.finite(x)] <- NA_real_; out
}
linear_types_score <- function(n) {
  n <- num(n)
  out <- ifelse(!is.finite(n), NA_real_, ifelse(n <= 5, 0, ifelse(n >= 12, 10, (n - 5) / 7 * 10)))
  out
}


# 🚩 1) Input-file QC and source provenance
step_header("CHNS input-file QC: single primary cohort")

req_path <- function(env_name, default_path) {
  x <- Sys.getenv(env_name, unset = default_path)
  normalizePath(x, winslash = "/", mustWork = FALSE)
}
opt_path <- function(env_name, default_path) {
  x <- Sys.getenv(env_name, unset = default_path)
  if (file.exists(x)) normalizePath(x, winslash = "/", mustWork = FALSE) else NA_character_
}

src <- list(
  NUTR3 = c(path=req_path("CHNS_NUTR3", file.path(raw_dir,"diet","nutr3_00.sas7bdat")), source="CHNS_public_use_plus_user_download"),
  C12DIET = c(path=req_path("CHNS_C12DIET", file.path(raw_dir,"diet","c12diet.sas7bdat")), source="CHNS_public_use"),
  MASTER = c(path=req_path("CHNS_MASTER", file.path(id2019_dir,"mast_pub_12.sas7bdat")), source="CHNS_2019_archive"),
  ROSTER = c(path=req_path("CHNS_ROSTER", file.path(id2019_dir,"rst_12.sas7bdat")), source="CHNS_2019_archive"),
  PEXAM = c(path=req_path("CHNS_PEXAM", pexam_public_default), source="CHNS_public_use_extra_copy_user_approved"),
  PACT = c(path=req_path("CHNS_PACT", file.path(raw_dir,"individual","pact_12.sas7bdat")), source="CHNS_public_use"),
  EDUC = c(path=opt_path("CHNS_EDUC", file.path(raw_dir,"individual","educ_12.sas7bdat")), source="CHNS_public_use_optional"),
  HHINC = c(path=opt_path("CHNS_HHINC", file.path(raw_dir,"constructed","hhinc_10.sas7bdat")), source="CHNS_public_use_optional"),
  JOBS = c(path=opt_path("CHNS_JOBS", file.path(raw_dir,"individual","jobs_12.sas7bdat")), source="CHNS_public_use_optional"),
  SURVEY = c(path=opt_path("CHNS_SURVEY", file.path(id2019_dir,"surveys_pub_12.sas7bdat")), source="CHNS_2019_archive_optional"),
  HEALTH = c(path=opt_path("CHNS_HEALTH", file.path(raw_dir,"individual","hlth_12.sas7bdat")), source="CHNS_public_use_optional")
)

food_map_file <- Sys.getenv("CHNS_FOOD_MAP", unset = food_map_default)
nutrient_map_file <- Sys.getenv("CHNS_NUTRIENT_MAP", unset = nutrient_map_default)
if (!file.exists(food_map_file)) food_map_file <- NA_character_
if (!file.exists(nutrient_map_file)) nutrient_map_file <- NA_character_

cat("[CHNS INPUT CHECK] project directory: ", project_dir, "\n", sep = "")
cat("[CHNS INPUT CHECK] raw directory: ", raw_dir, "\n", sep = "")
required_names <- c("NUTR3","C12DIET","MASTER","ROSTER","PEXAM","PACT")
required_ok <- logical()
for (nm in names(src)) {
  ok <- print_source(nm, src[[nm]], required = nm %in% required_names)
  if (nm %in% required_names) required_ok <- c(required_ok, ok)
}
cat("  FOOD_MAP     ", if(!is.na(food_map_file))"OK      classification : " else "MISSING required       : ",
    ifelse(is.na(food_map_file), food_map_default, food_map_file), ifelse(is.na(food_map_file)," 文件不存在",""), "\n", sep = "")
cat("  NUTRIENT_MAP ", if(!is.na(nutrient_map_file))"OK      code-level map       : " else if(require_nutrient_map==1L)"MISSING required       : " else "OPTIONAL missing       : ",
    ifelse(is.na(nutrient_map_file), nutrient_map_default, nutrient_map_file), "\n", sep = "")

if (!all(required_ok) || is.na(food_map_file)) {
  stop("CHNS required raw input/classification map missing. See ", log_file, call. = FALSE)
}
if (require_nutrient_map == 1L && is.na(nutrient_map_file)) {
  stop("CHNS_NUTRIENT_MAP is required in the publication pipeline but was not found: ", nutrient_map_default, call. = FALSE)
}

file_meta <- function(p) {
  if (is.na(p) || !file.exists(p)) return(tibble(exists=FALSE, bytes=NA_real_, mtime=NA_character_, md5=NA_character_))
  fi <- file.info(p)
  tibble(exists=TRUE, bytes=as.numeric(fi$size), mtime=format(fi$mtime, "%Y-%m-%d %H:%M:%S"), md5=unname(tools::md5sum(p)))
}
source_manifest <- tibble(
  input = names(src),
  path = vapply(src, `[[`, character(1), "path"),
  source = vapply(src, `[[`, character(1), "source")
) %>% bind_rows(tibble(input=c("FOOD_MAP","NUTRIENT_MAP"),
                       path=c(food_map_file,nutrient_map_file),
                       source=c("CFCT6_classification_user_curated","code-level_foodcode_nutrient_finegroup_user_curated"))) %>%
  rowwise() %>% mutate(meta=list(file_meta(path))) %>% unnest_wider(meta) %>% ungroup()
readr::write_tsv(source_manifest, cohort_file("input_source_manifest.tsv"))
cat("[CHNS MAP PROVENANCE] code-level files used in this run:\n")
print(source_manifest %>% filter(input %in% c("FOOD_MAP","NUTRIENT_MAP")) %>% select(input,path,bytes,mtime,md5))


# 🚩 2) Read modules and resolve waves
step_header("CHNS variable and wave QC")

# nutr3_00 is large: only keep fields required here.
nutr3 <- read_any(src$NUTR3[["path"]], c("IDIND","COMMID","WAVE","FOODCODE","VD","T1","T2","V39","V39C","V41","V42"))
c12 <- read_any(src$C12DIET[["path"]])
master <- read_any(src$MASTER[["path"]])
roster <- read_any(src$ROSTER[["path"]])
pexam <- read_any(src$PEXAM[["path"]])
pexam_long <- pexam
pact <- read_any(src$PACT[["path"]])
educ <- if (!is.na(src$EDUC[["path"]])) read_any(src$EDUC[["path"]]) else NULL
hhinc <- if (!is.na(src$HHINC[["path"]])) read_any(src$HHINC[["path"]]) else NULL
jobs_module <- if (!is.na(src$JOBS[["path"]])) read_any(src$JOBS[["path"]]) else NULL
survey <- if (!is.na(src$SURVEY[["path"]])) read_any(src$SURVEY[["path"]]) else NULL
health <- if (!is.na(src$HEALTH[["path"]])) read_any(src$HEALTH[["path"]]) else NULL


# Official longitudinal pact_12 is a large person-year raw module. Reject small/derived substitutes.
if (nrow(pact) < 100000 || ncol(pact) < 250) {
  stop("PACT raw-data sanity check failed: expected official longitudinal pact_12 (~125k rows, 276 variables), got ",
       nrow(pact), " rows and ", ncol(pact), " variables from ", src$PACT[["path"]],
       ". Do not use pact12.sav or pact12 - 1/2/3/4.sav author-analysis files.", call. = FALSE)
}
cat("[CHNS RAW SANITY] PACT rows=",nrow(pact),"; variables=",ncol(pact)," [PASS]\n",sep="")

idn <- pick(nutr3,"IDIND"); wn <- pick(nutr3,"WAVE"); fn <- pick(nutr3,"FOODCODE"); an <- pick(nutr3,"V39"); dn <- pick(nutr3,"VD")
proc_n <- pick(nutr3,"V39C",FALSE); prep_n <- pick(nutr3,"V42",FALSE); loc_n <- pick(nutr3,"V41",FALSE)
comm_n <- pick(nutr3,"COMMID",FALSE); prov_n <- pick(nutr3,"T1",FALSE); urban_n <- pick(nutr3,"T2",FALSE)

idc <- pick(c12,"IDIND"); wc <- pick(c12,"WAVE")
idm <- pick(master,"IDIND")
idr <- pick(roster,"IDIND"); wr <- pick(roster,"WAVE")
idp <- pick(pexam,"IDIND"); wp <- pick(pexam,c("WAVE","YEAR"))
idpa <- pick(pact,"IDIND"); wpa <- pick(pact,"WAVE")

nutr_waves <- valid_years(nutr3[[wn]])
c12_waves <- valid_years(c12[[wc]])
pex_waves <- valid_years(pexam[[wp]])
idpl <- idp; wpl <- wp
pex_long_waves <- valid_years(pexam_long[[wpl]])
pact_waves <- valid_years(pact[[wpa]])
roster_waves <- valid_years(roster[[wr]])
common_baseline <- Reduce(intersect, list(nutr_waves,c12_waves,pex_waves,pact_waves))

cat("  NUTR3 waves : ", paste(nutr_waves,collapse=", "), "\n", sep="")
cat("  C12DIET waves: ", paste(c12_waves,collapse=", "), "\n", sep="")
cat("  PEXAMPUB waves: ", paste(pex_waves,collapse=", "), "\n", sep="")
cat("  PACT waves  : ", paste(pact_waves,collapse=", "), "\n", sep="")
cat("  ROSTER waves: ", paste(roster_waves,collapse=", "), "\n", sep="")
cat("  Common baseline waves: ", paste(common_baseline,collapse=", "), "\n", sep="")

baseline_wave <- parse_year_request(baseline_req,"CHNS_BASELINE_WAVE")
followup_end <- parse_year_request(followup_req,"CHNS_FOLLOWUP_END")

if (is.na(baseline_wave)) {
  eligible <- common_baseline[vapply(common_baseline, function(w) any(roster_waves >= w + min_followup_years), logical(1))]
  if (!length(eligible)) eligible <- common_baseline[vapply(common_baseline, function(w) any(roster_waves > w), logical(1))]
  if (!length(eligible)) stop("No common dietary/exam/preference wave has later roster follow-up.", call.=FALSE)
  baseline_wave <- max(eligible)
}
if (!baseline_wave %in% common_baseline) {
  stop("Requested baseline ",baseline_wave," is unavailable jointly in NUTR3/C12DIET/PEXAM/PACT.",call.=FALSE)
}
if (is.na(followup_end)) {
  later <- roster_waves[roster_waves > baseline_wave]
  if (!length(later)) stop("No roster follow-up later than baseline ",baseline_wave,call.=FALSE)
  followup_end <- max(later)
}
if (!followup_end %in% roster_waves) stop("Requested follow-up ",followup_end," not present in ROSTER.",call.=FALSE)

cat("[CHNS WAVE RESOLUTION] baseline=",baseline_wave,"; follow-up=",followup_end,
    "; gap=",followup_end-baseline_wave," years\n",sep="")


# 🚩 3) CFCT6 classification map + optional code-level nutrient map
step_header("CHNS FOODCODE classification and map QC")

class_map <- read_any(food_map_file)
required_flags <- c("CEREAL","POTATO_STARCH","LEGUME","VEGETABLE","LEAF_STEM_FLOWER_VEG_PROXY",
                    "FRUIT","BERRY","NUT_SEED","RED_MEAT","POULTRY","DAIRY","CHEESE","EGG",
                    "FISH_SHELLFISH","SWEET_PASTRY","FAST_FOOD","CONVENIENCE_FOOD","SNACK_FOOD",
                    "BEVERAGE","CARBONATED_DRINK","ALCOHOL_BEVERAGE","SUGAR_CANDY_PRESERVE",
                    "ANIMAL_FAT","PLANT_OIL","CONDIMENT","SALT_MSG_OTHER","UPF_CATEGORY_PROXY",
                    "SSB_PROXY","SWEET_SNACK_PROXY")
required_class <- c("CFCT6_PREFIX3","CLASS2_TEXT","SUBCLASS1_TEXT","CLASS_NAME_CN","SUBCLASS_NAME_CN",required_flags)
miss_class <- setdiff(required_class,names(class_map))
if(length(miss_class)) stop("Classification FOOD_MAP schema incomplete. Missing: ",paste(miss_class,collapse=", "),call.=FALSE)
class_map <- class_map %>%
  mutate(CODE_PREFIX3=sprintf("%03d",as.integer(num(CFCT6_PREFIX3))),
         CLASS2=as.integer(num(CLASS2_TEXT)), SUBCLASS1=as.integer(num(SUBCLASS1_TEXT))) %>%
  distinct(CODE_PREFIX3,.keep_all=TRUE)
for(v in required_flags) class_map[[v]] <- num(class_map[[v]])
required_join <- c("CODE_PREFIX3","CLASS2","SUBCLASS1","CLASS_NAME_CN","SUBCLASS_NAME_CN",required_flags)

food0 <- nutr3 %>%
  transmute(
    IDIND=num(.data[[idn]]), WAVE=num(.data[[wn]]), FOODCODE=num(.data[[fn]]),
    GRAMS=num(.data[[an]]), DAY=num(.data[[dn]]),
    COMMID=colnum(nutr3,comm_n), PROVINCE=colnum(nutr3,prov_n), URBAN=colnum(nutr3,urban_n),
    PROCESSED_CODE=colnum(nutr3,proc_n), PREP_METHOD=colnum(nutr3,prep_n), MEAL_LOCATION=colnum(nutr3,loc_n)
  ) %>%
  filter(WAVE==baseline_wave,is.finite(IDIND),is.finite(FOODCODE),is.finite(GRAMS),GRAMS>=0) %>%
  mutate(FOODCODE6=sprintf("%06d",as.integer(FOODCODE)),CODE_PREFIX3=substr(FOODCODE6,1,3)) %>%
  left_join(class_map %>% select(all_of(required_join)),by="CODE_PREFIX3")
if (!nrow(food0)) stop("No nutr3_00 records at baseline wave ",baseline_wave,call.=FALSE)

class_cov <- sum(food0$GRAMS[!is.na(food0$CLASS2)],na.rm=TRUE)/sum(food0$GRAMS,na.rm=TRUE)
cat(sprintf("[CFCT6 classification] gram-weighted coverage %.1f%%\n",100*class_cov))
if(class_cov < .95) warning(sprintf("CFCT6 prefix classification covers only %.1f%% of baseline food grams; inspect unmatched FOODCODEs.",100*class_cov))

# Native CHNS processing/preparation indicators and class-derived proxies.
food0 <- food0 %>% mutate(
  OFFICIAL_PROCESSED=ifelse(is.finite(PROCESSED_CODE),as.integer(PROCESSED_CODE==3),NA_integer_),
  RESTAURANT_MADE=ifelse(is.finite(PROCESSED_CODE),as.integer(PROCESSED_CODE==2),NA_integer_),
  DEEPFRIED=ifelse(is.finite(PREP_METHOD),as.integer(PREP_METHOD==3),NA_integer_)
)

# Code-level audit for transparent manual completion later.
partial_map <- food0 %>% group_by(FOODCODE,FOODCODE6,CODE_PREFIX3) %>% summarise(
  CLASS2=first(na.omit(CLASS2),default=NA),
  SUBCLASS1=first(na.omit(SUBCLASS1),default=NA),
  CLASS_NAME_CN=first(na.omit(CLASS_NAME_CN),default=NA_character_),
  SUBCLASS_NAME_CN=first(na.omit(SUBCLASS_NAME_CN),default=NA_character_),
  grams=sum(GRAMS,na.rm=TRUE),records=n(),
  V39C_processed_fraction=ifelse(any(is.finite(PROCESSED_CODE)),mean(PROCESSED_CODE==3,na.rm=TRUE),NA_real_),
  V42_deepfried_fraction=ifelse(any(is.finite(PREP_METHOD)),mean(PREP_METHOD==3,na.rm=TRUE),NA_real_),
  .groups="drop"
) %>% arrange(desc(grams))
readr::write_tsv(partial_map, cohort_file("food_map.code_audit.tsv"))

# Prefer CFCT 2002/2004 numeric FOODCODE rows; evaluate coverage by field.
full_map <- NULL
code_map_available <- FALSE
full_map_mode <- FALSE
enhanced_map_mode <- FALSE
full_map_coverage <- NA_real_
map_field_coverage <- tibble(field=character(),gram_coverage=double())
map_field_ok <- setNames(logical(),character())
required_code_fields <- c("WHOLE_GRAIN","REFINED_GRAIN","LOWFAT_DAIRY","GREEN_LEAFY",
                           "PROCESSED_MEAT","SSB","SODIUM_MG_100G","SFAT_G_100G","MUFA_G_100G",
                           "PUFA_G_100G","ALCOHOL_G_100G")
fine_group_fields <- c("WHOLE_GRAIN","REFINED_GRAIN","LOWFAT_DAIRY","GREEN_LEAFY","PROCESSED_MEAT","SSB")
fat_fields <- c("SFAT_G_100G","MUFA_G_100G","PUFA_G_100G")

if(!is.na(nutrient_map_file)) {
  tmp <- read_any(nutrient_map_file)
  miss <- setdiff(required_code_fields,names(tmp))
  if(length(miss) || (!"FOODCODE6" %in% names(tmp) && !"FOODCODE" %in% names(tmp))) {
    warning("NUTRIENT_MAP schema incomplete. Missing: ",
            paste(c(miss,if(!"FOODCODE6"%in%names(tmp) && !"FOODCODE"%in%names(tmp))"FOODCODE or FOODCODE6" else character()),collapse=", "))
  } else {
    if("FOODCODE6" %in% names(tmp)) {
      tmp <- tmp %>% mutate(FOODCODE6=sprintf("%06d",as.integer(num(FOODCODE6))))
    } else {
      tmp <- tmp %>% mutate(FOODCODE6=sprintf("%06d",as.integer(num(FOODCODE))))
    }
    tmp <- tmp %>% filter(grepl("^[0-9]{6}$",FOODCODE6)) %>% distinct(FOODCODE6,.keep_all=TRUE)
    for(v in required_code_fields) tmp[[v]] <- num(tmp[[v]])
    gram_by_code <- food0 %>% group_by(FOODCODE6) %>% summarise(grams=sum(GRAMS,na.rm=TRUE),.groups="drop")
    aud <- gram_by_code %>% left_join(tmp,by="FOODCODE6")
    den <- sum(aud$grams,na.rm=TRUE)
    if(is.finite(den) && den>0) {
      full_map_coverage <- sum(aud$grams[!is.na(aud$WHOLE_GRAIN)],na.rm=TRUE)/den
      map_field_coverage <- tibble(field=required_code_fields) %>% rowwise() %>%
        mutate(gram_coverage=sum(aud$grams[is.finite(num(aud[[field]]))],na.rm=TRUE)/den) %>% ungroup()
      map_field_ok <- setNames(map_field_coverage$gram_coverage>=full_map_min_cov,map_field_coverage$field)
      full_map <- tmp
      code_map_available <- TRUE
      fine_ok <- all(map_field_ok[fine_group_fields],na.rm=TRUE)
      sodium_ok <- isTRUE(map_field_ok[["SODIUM_MG_100G"]])
      fat_ok <- all(map_field_ok[fat_fields],na.rm=TRUE)
      alcohol_ok <- isTRUE(map_field_ok[["ALCOHOL_G_100G"]])
      enhanced_map_mode <- isTRUE(fine_ok)
      full_map_mode <- isTRUE(fine_ok && sodium_ok && fat_ok && alcohol_ok)
      cat("[CODE-LEVEL FOODCODE COVERAGE] ",
          paste0(map_field_coverage$field,"=",sprintf("%.1f%%",100*map_field_coverage$gram_coverage),collapse="; "),"\n",sep="")
      if(!enhanced_map_mode)
        warning(sprintf("Code-level fine-group coverage is below %.1f%% for at least one core field; proxy mode retained.",100*full_map_min_cov))
      else if(!full_map_mode)
        warning("Code-level fine groups are sufficiently covered, but one or more nutrient domains (sodium/fat/alcohol) are below threshold; foodcode-enhanced mode will use transparent proxy fallback for those domains.")
    }
  }
}
if(!exists("sodium_ok")) sodium_ok <- FALSE
if(!exists("fat_ok")) fat_ok <- FALSE
if(!exists("alcohol_ok")) alcohol_ok <- FALSE
if (require_enhanced_map == 1L && !enhanced_map_mode) {
  stop(sprintf("Code-level FOODCODE fine-group map coverage is insufficient for publication scoring (< %.0f%% in at least one required fine-group field). Inspect %s and %s.",
               100*full_map_min_cov, food_map_file, nutrient_map_file), call. = FALSE)
}

analysis_mode <- if(full_map_mode) {
  "FOODCODE-complete mode: code-level map with all required nutrient domains above threshold"
} else if(enhanced_map_mode) {
  "FOODCODE-enhanced mode: code-level fine groups + field-specific nutrient/proxy fallback"
} else {
  "CFCT6 classification mode: MAHA-CN/DASH-CN/MIND-CN/MEDI-CN proxies"
}
cat("[CHNS MAP MODE] ",analysis_mode,"\n",sep="")


# 🚩 4) Aggregate dietary intake
step_header("CHNS dietary harmonization")

food <- food0
if(code_map_available) food <- food %>% left_join(full_map %>% select(FOODCODE6,all_of(required_code_fields)),by="FOODCODE6")
for(v in required_code_fields) if(!v %in% names(food)) food[[v]] <- NA_real_
for(v in required_code_fields) food[[v]] <- num(food[[v]])

food_day <- food %>% group_by(IDIND,DAY) %>% summarise(
  total_food_g=sum(GRAMS,na.rm=TRUE), n_food_types=n_distinct(FOODCODE),
  cereal=weighted_sum_known(GRAMS,CEREAL), potato_starch=weighted_sum_known(GRAMS,POTATO_STARCH),
  legumes=weighted_sum_known(GRAMS,LEGUME), vegetable=weighted_sum_known(GRAMS,VEGETABLE),
  leaf_stem_flower_veg=weighted_sum_known(GRAMS,LEAF_STEM_FLOWER_VEG_PROXY),
  fruit=weighted_sum_known(GRAMS,FRUIT), berry=weighted_sum_known(GRAMS,BERRY),
  nuts=weighted_sum_known(GRAMS,NUT_SEED), red_meat=weighted_sum_known(GRAMS,RED_MEAT),
  poultry=weighted_sum_known(GRAMS,POULTRY), dairy=weighted_sum_known(GRAMS,DAIRY),
  cheese=weighted_sum_known(GRAMS,CHEESE), eggs=weighted_sum_known(GRAMS,EGG),
  fish=weighted_sum_known(GRAMS,FISH_SHELLFISH),
  sweet_pastry_proxy=weighted_sum_known(GRAMS,SWEET_PASTRY),
  fast_food_proxy=weighted_sum_known(GRAMS,FAST_FOOD), convenience_food_proxy=weighted_sum_known(GRAMS,CONVENIENCE_FOOD),
  snack_food_proxy=weighted_sum_known(GRAMS,SNACK_FOOD), beverage_g=weighted_sum_known(GRAMS,BEVERAGE),
  carbonated_drink_g=weighted_sum_known(GRAMS,CARBONATED_DRINK), alcohol_beverage_g=weighted_sum_known(GRAMS,ALCOHOL_BEVERAGE),
  sugar_candy_preserve_g=weighted_sum_known(GRAMS,SUGAR_CANDY_PRESERVE),
  animal_fat_g=weighted_sum_known(GRAMS,ANIMAL_FAT), plant_oil_g=weighted_sum_known(GRAMS,PLANT_OIL),
  condiment_g=weighted_sum_known(GRAMS,CONDIMENT), salt_msg_proxy_g=weighted_sum_known(GRAMS,SALT_MSG_OTHER),
  upf_class_proxy_g=weighted_sum_known(GRAMS,UPF_CATEGORY_PROXY), ssb_proxy_g=weighted_sum_known(GRAMS,SSB_PROXY),
  sweet_snack_proxy_g=weighted_sum_known(GRAMS,SWEET_SNACK_PROXY),
  processed_food_g=weighted_sum_known(GRAMS,OFFICIAL_PROCESSED), restaurant_made_g=weighted_sum_known(GRAMS,RESTAURANT_MADE),
  deepfried_g=weighted_sum_known(GRAMS,DEEPFRIED),
  whole_grain=weighted_sum_known(GRAMS,WHOLE_GRAIN), refined_grain=weighted_sum_known(GRAMS,REFINED_GRAIN),
  lowfat_dairy=weighted_sum_known(GRAMS,LOWFAT_DAIRY), green_leafy=weighted_sum_known(GRAMS,GREEN_LEAFY),
  processed_meat=weighted_sum_known(GRAMS,PROCESSED_MEAT), ssb_code=weighted_sum_known(GRAMS,SSB),
  sodium_mg=sum_known(GRAMS/100*SODIUM_MG_100G), sfat_g=sum_known(GRAMS/100*SFAT_G_100G),
  mufa_g=sum_known(GRAMS/100*MUFA_G_100G), pufa_g=sum_known(GRAMS/100*PUFA_G_100G),
  alcohol_g_food=sum_known(GRAMS/100*ALCOHOL_G_100G), .groups="drop"
)
food_person <- food_day %>% group_by(IDIND) %>% summarise(
  n_diet_days=n_distinct(DAY), across(-DAY,~mean(.x,na.rm=TRUE)), .groups="drop"
) %>% filter(n_diet_days>=min_diet_days)

geo_person <- food0 %>% group_by(IDIND) %>% summarise(
  COMMID=first_num(COMMID),province=factor(first_num(PROVINCE)),urban=factor(first_num(URBAN)),.groups="drop")

# c12diet
kcal_c <- pick(c12,c("D3KCAL","KCAL"))
protein_c <- pick(c12,c("D3PROTN","D3PROT","PROTEIN"))
fat_c <- pick(c12,c("D3FAT","FAT"),FALSE)
carb_c <- pick(c12,c("D3CARBO","CARBO"),FALSE)
c12_base <- c12 %>% filter(num(.data[[wc]])==baseline_wave)
c12b <- c12_base %>%
  transmute(IDIND=num(.data[[idc]]),kcal=num(.data[[kcal_c]]),protein_g=num(.data[[protein_c]]),
            fat_total_g=colnum(c12_base,fat_c),carb_g=colnum(c12_base,carb_c)) %>%
  group_by(IDIND) %>% summarise(across(everything(),~mean(.x,na.rm=TRUE)),.groups="drop")

# pexampub12: anthropometrics, alcohol, smoking, disease indicators, BP
weight_p <- pick(pexam,c("WEIGHT","WT"),FALSE); height_p <- pick(pexam,c("HEIGHT","HT"),FALSE)
u25 <- pick(pexam,"U25",FALSE); u27 <- pick(pexam,"U27",FALSE); u40 <- pick(pexam,"U40",FALSE); u41 <- pick(pexam,"U41",FALSE)
u232 <- pick(pexam,"U232",FALSE); u234 <- pick(pexam,"U234",FALSE); u236 <- pick(pexam,"U236",FALSE)
u22 <- pick(pexam,"U22",FALSE); u24a <- pick(pexam,"U24A",FALSE); u24j <- pick(pexam,"U24J",FALSE)
u24l <- pick(pexam,"U24L",FALSE); u24w <- pick(pexam,"U24W",FALSE)
sys <- intersect(c("SYSTOL1","SYSTOL2","SYSTOL3"),names(pexam))
dia <- intersect(c("DIASTOL1","DIASTOL2","DIASTOL3"),names(pexam))

pex_base <- pexam %>% filter(num(.data[[wp]])==baseline_wave)
sugar_mat <- cbind(colnum(pex_base,u232),colnum(pex_base,u234),colnum(pex_base,u236))
sugar_week <- rowSums(sugar_mat,na.rm=TRUE); sugar_week[rowSums(is.finite(sugar_mat))==0] <- NA_real_
sbp_base <- if(length(sys)) rowMeans(as.data.frame(lapply(pex_base[sys],num)),na.rm=TRUE) else rep(NA_real_,nrow(pex_base))
dbp_base <- if(length(dia)) rowMeans(as.data.frame(lapply(pex_base[dia],num)),na.rm=TRUE) else rep(NA_real_,nrow(pex_base))
sbp_base[is.nan(sbp_base)] <- NA_real_; dbp_base[is.nan(dbp_base)] <- NA_real_
pexb <- pex_base %>% transmute(
    IDIND=num(.data[[idp]]),
    weight_kg=colnum(pex_base,weight_p),height_cm=colnum(pex_base,height_p),
    ever_smoke=if(!is.na(u25))yes1(.data[[u25]])else NA_integer_,
    current_smoke=if(!is.na(u27))yes1(.data[[u27]])else NA_integer_,
    alcohol_drinker=if(!is.na(u40))yes1(.data[[u40]])else NA_integer_,
    alcohol_frequency=colnum(pex_base,u41),
    sugary_drinks_week=sugar_week,
    htn_self=if(!is.na(u22))yes1(.data[[u22]])else NA_integer_,
    t2dm_self=if(!is.na(u24a))yes1(.data[[u24a]])else NA_integer_,
    mi_self=if(!is.na(u24j))yes1(.data[[u24j]])else NA_integer_,
    stroke_self=if(!is.na(u24l))yes1(.data[[u24l]])else NA_integer_,
    cancer_self=if(!is.na(u24w))yes1(.data[[u24w]])else NA_integer_,
    sbp=sbp_base,dbp=dbp_base
  ) %>% group_by(IDIND) %>% summarise(across(everything(),~first_num(.x)),.groups="drop")


# The same longitudinal pexampub12 public-use copy supplies baseline examination data
# and later-wave prospective outcomes (including 2015 when present).
pexfu <- NULL
if(followup_end %in% pex_long_waves) {
  fu22 <- pick(pexam_long,"U22",FALSE); fu24a <- pick(pexam_long,"U24A",FALSE)
  fu24j <- pick(pexam_long,"U24J",FALSE); fu24l <- pick(pexam_long,"U24L",FALSE); fu24w <- pick(pexam_long,"U24W",FALSE)
  fu_weight <- pick(pexam_long,c("WEIGHT","WT"),FALSE); fu_height <- pick(pexam_long,c("HEIGHT","HT"),FALSE)
  fu_sys <- intersect(c("SYSTOL1","SYSTOL2","SYSTOL3"),names(pexam_long))
  fu_dia <- intersect(c("DIASTOL1","DIASTOL2","DIASTOL3"),names(pexam_long))
  pex_fu_base <- pexam_long %>% filter(num(.data[[wpl]])==followup_end)
  sbp_fu_vec <- if(length(fu_sys))rowMeans(as.data.frame(lapply(pex_fu_base[fu_sys],num)),na.rm=TRUE)else rep(NA_real_,nrow(pex_fu_base))
  dbp_fu_vec <- if(length(fu_dia))rowMeans(as.data.frame(lapply(pex_fu_base[fu_dia],num)),na.rm=TRUE)else rep(NA_real_,nrow(pex_fu_base))
  sbp_fu_vec[is.nan(sbp_fu_vec)] <- NA_real_; dbp_fu_vec[is.nan(dbp_fu_vec)] <- NA_real_
  pexfu <- pex_fu_base %>% transmute(
    IDIND=num(.data[[idpl]]),
    weight_kg_fu=colnum(pex_fu_base,fu_weight),height_cm_fu=colnum(pex_fu_base,fu_height),
    sbp_fu=sbp_fu_vec,dbp_fu=dbp_fu_vec,
    htn_self_fu=if(!is.na(fu22))yes1(.data[[fu22]])else NA_integer_,
    t2dm_fu=if(!is.na(fu24a))yes1(.data[[fu24a]])else NA_integer_,
    mi_fu=if(!is.na(fu24j))yes1(.data[[fu24j]])else NA_integer_,
    stroke_fu=if(!is.na(fu24l))yes1(.data[[fu24l]])else NA_integer_,
    cancer_fu=if(!is.na(fu24w))yes1(.data[[fu24w]])else NA_integer_
  ) %>% group_by(IDIND) %>% summarise(across(everything(),~first_num(.x)),.groups="drop")
}

# PACT food preferences + sleep
pref_vars <- c(
  fast=pick(pact,c("U389","U389_FF"),FALSE),
  salty=pick(pact,c("U390","U390_SSF"),FALSE),
  fruit=pick(pact,c("U391","U391_F"),FALSE),
  veg=pick(pact,c("U392","U392_V"),FALSE),
  ssb=pick(pact,c("U393","U393_SD"),FALSE)
)
sleep_var <- pick(pact,c("U324","TIMEBED"),FALSE)

pact_base <- pact %>% filter(num(.data[[wpa]])==baseline_wave)
pactb <- pact_base %>%
  transmute(
    IDIND=num(.data[[idpa]]),
    pref_fast=if(!is.na(pref_vars["fast"]))pref15(.data[[pref_vars["fast"]]])else NA_real_,
    pref_salty=if(!is.na(pref_vars["salty"]))pref15(.data[[pref_vars["salty"]]])else NA_real_,
    pref_fruit=if(!is.na(pref_vars["fruit"]))pref15(.data[[pref_vars["fruit"]]])else NA_real_,
    pref_veg=if(!is.na(pref_vars["veg"]))pref15(.data[[pref_vars["veg"]]])else NA_real_,
    pref_ssb=if(!is.na(pref_vars["ssb"]))pref15(.data[[pref_vars["ssb"]]])else NA_real_,
    sleep_hours=colnum(pact_base,sleep_var)
  ) %>% group_by(IDIND) %>% summarise(across(everything(),~first_num(.x)),.groups="drop") %>%
  mutate(
    pref_healthy=rowMeans(cbind(pref_fruit,pref_veg),na.rm=TRUE),
    pref_unhealthy=rowMeans(cbind(pref_fast,pref_salty,pref_ssb),na.rm=TRUE),
    pref_health_index=rowMeans(cbind(zstd(pref_fruit),zstd(pref_veg),-zstd(pref_fast),-zstd(pref_salty),-zstd(pref_ssb)),na.rm=TRUE)
  )
pactb$pref_healthy[!is.finite(pactb$pref_healthy)] <- NA_real_
pactb$pref_unhealthy[!is.finite(pactb$pref_unhealthy)] <- NA_real_
pactb$pref_health_index[!is.finite(pactb$pref_health_index)] <- NA_real_

# Education
educb <- NULL
if(!is.null(educ)) {
  ide <- pick(educ,"IDIND"); we <- pick(educ,"WAVE")
  a11 <- pick(educ,c("A11","EDUYRS"),FALSE); a12 <- pick(educ,c("A12","EDUCATION"),FALSE)
  educ_base <- educ %>% filter(num(.data[[we]])==baseline_wave)
  schooling <- colnum(educ_base,a11)
  if (identical(a11,"A11")) {
    schooling <- recode_chns_school_years(schooling)
  }
  educb <- educ_base %>%
    transmute(IDIND=num(.data[[ide]]),education_years=schooling,education_level=colnum(educ_base,a12)) %>%
    group_by(IDIND) %>% summarise(across(everything(),~first_num(.x)),.groups="drop")
}

# Household income
incb <- NULL
if(!is.null(hhinc)) {
  wh <- pick(hhinc,"WAVE"); hh <- pick(hhinc,"HHID")
  incv <- pick(hhinc,c("HHINCPC","HHINC_PC","HHINCPC_CPI","HHINC_CPI"),FALSE)
  if(!is.na(incv)) {
    hhinc_base <- hhinc %>% filter(num(.data[[wh]])==baseline_wave)
    incb <- hhinc_base %>%
      transmute(HHID=num(.data[[hh]]),income_pc=num(.data[[incv]])) %>%
      group_by(HHID) %>% summarise(income_pc=first_num(income_pc),.groups="drop")
  }
}

# Master / mortality
gender_m <- pick(master,c("GENDER","SEX"))
birth_m <- pick(master,c("WEST_DOB_Y","DOB_Y","BIRTH_Y","MOON_DOB_Y"))
dody_m <- pick(master,c("DOD_Y","DODY"),FALSE)
dodrpt_m <- pick(master,c("DOD_RPT","DOD_REPORT"),FALSE)
masterb <- master %>% transmute(
  IDIND=num(.data[[idm]]),female=as.integer(num(.data[[gender_m]])==2),
  birth_year=num(.data[[birth_m]]),DOD_Y=colnum(master,dody_m),DOD_RPT=colnum(master,dodrpt_m)
) %>% group_by(IDIND) %>% summarise(across(everything(),~first_num(.x)),.groups="drop")

# Retain both end-wave status and last observed contact year.
hh_r <- pick(roster,"HHID",FALSE)
roster_base0 <- roster %>% filter(num(.data[[wr]])==baseline_wave)
roster_base <- roster_base0 %>%
  transmute(IDIND=num(.data[[idr]]),HHID=colnum(roster_base0,hh_r)) %>% distinct(IDIND,.keep_all=TRUE)
roster_follow <- roster %>%
  transmute(IDIND=num(.data[[idr]]), CONTACT_YEAR=num(.data[[wr]])) %>%
  filter(is.finite(IDIND),is.finite(CONTACT_YEAR),CONTACT_YEAR>=baseline_wave,CONTACT_YEAR<=followup_end) %>%
  group_by(IDIND) %>%
  summarise(
    last_contact_year=max(CONTACT_YEAR),
    observed_after_baseline=as.integer(any(CONTACT_YEAR>baseline_wave)),
    observed_at_followup=as.integer(any(CONTACT_YEAR==followup_end)),
    .groups="drop"
  )

dat <- food_person %>%
  left_join(geo_person,by="IDIND") %>%
  left_join(c12b,by="IDIND") %>%
  left_join(pexb,by="IDIND") %>%
  {if(!is.null(pexfu)) left_join(.,pexfu,by="IDIND") else .} %>%
  left_join(pactb,by="IDIND") %>%
  left_join(masterb,by="IDIND") %>%
  left_join(roster_base,by="IDIND") %>%
  left_join(roster_follow,by="IDIND")
if(!is.null(educb)) dat <- dat %>% left_join(educb,by="IDIND")
if(!is.null(incb) && "HHID" %in% names(dat)) dat <- dat %>% left_join(incb,by="HHID")
# Optional modules must never make the core pipeline fail simply because a covariate
# was unavailable in a particular download/format release.
for(v in c("income_pc","education_years","education_level","weight_kg_fu","height_cm_fu","sbp_fu","dbp_fu",
           "htn_self_fu","t2dm_fu","mi_fu","stroke_fu","cancer_fu")) {
  if(!v %in% names(dat)) dat[[v]] <- NA_real_
}

dat <- dat %>% mutate(
  age=baseline_wave-birth_year,
  bmi=ifelse(is.finite(weight_kg)&is.finite(height_cm)&height_cm>100,weight_kg/(height_cm/100)^2,NA_real_),
  log_income=log1p(ifelse(is.finite(income_pc)&income_pc>=0,income_pc,NA_real_)),
  htn_measured=ifelse(is.finite(sbp)|is.finite(dbp),as.integer((sbp>=140)|(dbp>=90)),NA_integer_),
  hypertension=ifelse(htn_self==1|htn_measured==1,1L,ifelse(htn_self==0&htn_measured==0,0L,NA_integer_)),
  obesity=ifelse(is.finite(bmi),as.integer(bmi>=28),NA_integer_),
  bmi_fu=ifelse(is.finite(weight_kg_fu)&is.finite(height_cm_fu)&height_cm_fu>100,weight_kg_fu/(height_cm_fu/100)^2,NA_real_),
  htn_measured_fu=ifelse(is.finite(sbp_fu)|is.finite(dbp_fu),as.integer((sbp_fu>=140)|(dbp_fu>=90)),NA_integer_),
  hypertension_fu=ifelse(htn_self_fu==1|htn_measured_fu==1,1L,ifelse(htn_self_fu==0&htn_measured_fu==0,0L,NA_integer_)),
  obesity_fu=ifelse(is.finite(bmi_fu),as.integer(bmi_fu>=28),NA_integer_),
  inc_hypertension=ifelse(hypertension==0 & hypertension_fu%in%0:1,as.integer(hypertension_fu==1),NA_integer_),
  inc_obesity=ifelse(obesity==0 & obesity_fu%in%0:1,as.integer(obesity_fu==1),NA_integer_),
  inc_t2dm=ifelse(t2dm_self==0 & t2dm_fu%in%0:1,as.integer(t2dm_fu==1),NA_integer_),
  inc_mi=ifelse(mi_self==0 & mi_fu%in%0:1,as.integer(mi_fu==1),NA_integer_),
  inc_stroke=ifelse(stroke_self==0 & stroke_fu%in%0:1,as.integer(stroke_fu==1),NA_integer_),
  inc_cancer=ifelse(cancer_self==0 & cancer_fu%in%0:1,as.integer(cancer_fu==1),NA_integer_),
  delta_bmi=ifelse(is.finite(bmi)&is.finite(bmi_fu),bmi_fu-bmi,NA_real_),
  delta_sbp=ifelse(is.finite(sbp)&is.finite(sbp_fu),sbp_fu-sbp,NA_real_),
  delta_dbp=ifelse(is.finite(dbp)&is.finite(dbp_fu),dbp_fu-dbp,NA_real_),
  observed_at_followup=replace_na(as.integer(observed_at_followup),0L),
  observed_after_baseline=replace_na(as.integer(observed_after_baseline),0L),
  death_by_year=is.finite(DOD_Y)&DOD_Y>baseline_wave&DOD_Y<=followup_end,
  death_by_report=is.finite(DOD_RPT)&DOD_RPT>baseline_wave&DOD_RPT<=followup_end,
  death_year=case_when(
    death_by_year & death_by_report ~ pmin(DOD_Y,DOD_RPT),
    death_by_year ~ DOD_Y,
    death_by_report ~ DOD_RPT,
    TRUE ~ NA_real_
  ),
  death_reported=as.integer(death_by_year|death_by_report),
  followup_known=as.integer(observed_at_followup==1|death_reported==1),
  followup_known_last_contact=as.integer(observed_after_baseline==1|death_reported==1),
  death=ifelse(followup_known==1,death_reported,NA_integer_),
  protein_g_kg=protein_g/weight_kg,
  nuts_legumes=nuts+legumes,
  cereal_potato=cereal+potato_starch,
  meat_poultry_eggs=red_meat+poultry+eggs,
  fat_proxy_ratio=ifelse(is.finite(plant_oil_g)|is.finite(animal_fat_g),plant_oil_g/(plant_oil_g+animal_fat_g+0.1),NA_real_),
  healthy_fat_ratio=ifelse(is.finite(mufa_g)|is.finite(pufa_g)|is.finite(sfat_g),(mufa_g+pufa_g)/(sfat_g+0.1),NA_real_),
  upf_native=processed_food_g,
  upf_combined_proxy=rowMeans(cbind(zstd(processed_food_g),zstd(upf_class_proxy_g)),na.rm=TRUE),
  alcohol_proxy=alcohol_drinker
) %>%
  filter(age>=min_age,age<=100,is.finite(kcal),kcal>=500,kcal<=5000)

cat("[CHNS LINKAGE] baseline N=",nrow(dat),"; known follow-up N=",sum(dat$followup_known==1,na.rm=TRUE),
    "; deaths=",sum(dat$death==1,na.rm=TRUE),"\n",sep="")


# 🚩 5) CHDI + MAHA-CN / established-score proxy construction
step_header("CHNS diet-score construction")

# CHDI present-component version (0-70), using directly observed food-group intake.
dat <- dat %>% mutate(
  food_types_per_day=n_food_types,
  chdi_types=linear_types_score(food_types_per_day),
  chdi_cereal=linear_adequacy(cereal_potato/kcal*1000,140,10),
  chdi_veg=linear_adequacy(vegetable/kcal*1000,270,10),
  chdi_fruit=linear_adequacy(fruit/kcal*1000,110,10),
  chdi_dairy=linear_adequacy(dairy/kcal*1000,100,10),
  chdi_legnuts=linear_adequacy(nuts_legumes/kcal*1000,10,10),
  chdi_meategg=linear_adequacy(meat_poultry_eggs/kcal*1000,50,5),
  chdi_aquatic=linear_adequacy(fish/kcal*1000,30,5),
  diet.chdi.sum=rowSums(cbind(chdi_types,chdi_cereal,chdi_veg,chdi_fruit,chdi_dairy,chdi_legnuts,chdi_meategg,chdi_aquatic),na.rm=FALSE)
)

# Build publication scores from fine-group maps when coverage passes the threshold.

# Harmonized food quantities used in several scores.
dat <- dat %>% mutate(
  red_processed_meat_g = ifelse(is.finite(red_meat) | is.finite(processed_meat),
                                rowSums(cbind(red_meat, processed_meat), na.rm = TRUE), NA_real_),
  fried_fast_g = ifelse(is.finite(deepfried_g) | is.finite(fast_food_proxy),
                        rowSums(cbind(deepfried_g, fast_food_proxy), na.rm = TRUE), NA_real_),
  sweets_pastries_g = ifelse(is.finite(sweet_snack_proxy_g) | is.finite(sugar_candy_preserve_g),
                             rowSums(cbind(sweet_snack_proxy_g, sugar_candy_preserve_g), na.rm = TRUE), NA_real_),
  upf_pub_g = ifelse(
    is.finite(refined_grain) | is.finite(sweets_pastries_g) | is.finite(ssb_code) |
      is.finite(processed_meat) | is.finite(deepfried_g) | is.finite(fast_food_proxy) |
      is.finite(convenience_food_proxy),
    rowSums(cbind(refined_grain, sweets_pastries_g, ssb_code, processed_meat,
                  deepfried_g, fast_food_proxy, convenience_food_proxy), na.rm = TRUE),
    NA_real_
  )
)

# Build prespecified sex-specific residual and density score inputs.
dat <- dat %>% mutate(
  maha_den_protein = energy_density(protein_g, kcal),
  maha_den_dairy = energy_density(dairy, kcal),
  maha_den_veg = energy_density(vegetable, kcal),
  maha_den_fruit = energy_density(fruit, kcal),
  maha_den_wholegrain = energy_density(whole_grain, kcal),
  maha_den_upf = energy_density(upf_pub_g, kcal),
  maha_den_alcohol = energy_density(alcohol_g_food, kcal),
  maha_den_sodium = energy_density(sodium_mg, kcal),
  maha_res_protein = energy_residual_grouped(protein_g, kcal, female),
  maha_res_dairy = energy_residual_grouped(dairy, kcal, female),
  maha_res_veg = energy_residual_grouped(vegetable, kcal, female),
  maha_res_fruit = energy_residual_grouped(fruit, kcal, female),
  maha_res_wholegrain = energy_residual_grouped(whole_grain, kcal, female),
  maha_res_upf = energy_residual_grouped(upf_pub_g, kcal, female),
  maha_res_alcohol = energy_residual_grouped(alcohol_g_food, kcal, female),
  maha_res_sodium = energy_residual_grouped(sodium_mg, kcal, female),
  dash_res_fruit = energy_residual_grouped(fruit, kcal, female),
  dash_res_veg = energy_residual_grouped(vegetable, kcal, female),
  dash_res_wholegrain = energy_residual_grouped(whole_grain, kcal, female),
  dash_res_lowfatdairy = energy_residual_grouped(lowfat_dairy, kcal, female),
  dash_res_nutslegumes = energy_residual_grouped(nuts_legumes, kcal, female),
  dash_res_sodium = energy_residual_grouped(if (sodium_ok) sodium_mg else salt_msg_proxy_g, kcal, female),
  dash_res_redmeat = energy_residual_grouped(red_processed_meat_g, kcal, female),
  dash_res_ssb = energy_residual_grouped(ssb_code, kcal, female),
  dash_den_fruit = energy_density(fruit, kcal),
  dash_den_veg = energy_density(vegetable, kcal),
  dash_den_wholegrain = energy_density(whole_grain, kcal),
  dash_den_lowfatdairy = energy_density(lowfat_dairy, kcal),
  dash_den_nutslegumes = energy_density(nuts_legumes, kcal),
  dash_den_sodium = energy_density(if (sodium_ok) sodium_mg else salt_msg_proxy_g, kcal),
  dash_den_redmeat = energy_density(red_processed_meat_g, kcal),
  dash_den_ssb = energy_density(ssb_code, kcal)
)

# Canonical MAHA components: same conceptual nine domains used in UKB/NHANES.
dat <- dat %>% mutate(
  maha_c_protein = qscore(protein_g_kg),
  maha_c_dairy = qscore(dairy),
  maha_c_veg = qscore(vegetable),
  maha_c_fruit = qscore(fruit),
  maha_c_wholegrain = qscore(whole_grain),
  maha_c_fat = if (fat_ok) qscore(healthy_fat_ratio) else qscore(fat_proxy_ratio),
  maha_c_upf = qscore(upf_pub_g, TRUE),
  maha_c_alcohol = if (alcohol_ok) qscore(alcohol_g_food, TRUE) else ifelse(alcohol_proxy == 0, 10, ifelse(alcohol_proxy == 1, 0, NA_real_)),
  maha_c_sodium = if (sodium_ok) qscore(sodium_mg, TRUE) else qscore(salt_msg_proxy_g, TRUE),
  maha_c_redmeat = qscore(red_processed_meat_g, TRUE),

  # Primary MAHA and sensitivity specifications aligned to UKB/NHANES naming.
  diet.maha.sum = rowmean_vecs(maha_c_protein, maha_c_dairy, maha_c_veg, maha_c_fruit,
                               maha_c_wholegrain, maha_c_fat, maha_c_upf, maha_c_alcohol,
                               maha_c_sodium, min_n = 6),
  # "balanced" removes the most definition-sensitive alcohol domain.
  diet.maha_bal.sum = rowmean_vecs(maha_c_protein, maha_c_dairy, maha_c_veg, maha_c_fruit,
                                   maha_c_wholegrain, maha_c_fat, maha_c_upf, maha_c_sodium,
                                   min_n = 5),
  # "strict" removes protein/dairy and adds low red/processed meat.
  diet.maha_strict.sum = rowmean_vecs(maha_c_veg, maha_c_fruit, maha_c_wholegrain, maha_c_fat,
                                      maha_c_upf, maha_c_alcohol, maha_c_sodium, maha_c_redmeat,
                                      min_n = 5),
  diet.maha_nodairy.sum = rowmean_vecs(maha_c_protein, maha_c_veg, maha_c_fruit, maha_c_wholegrain,
                                       maha_c_fat, maha_c_upf, maha_c_alcohol, maha_c_sodium,
                                       min_n = 5),
  diet.maha_noprotein.sum = rowmean_vecs(maha_c_dairy, maha_c_veg, maha_c_fruit, maha_c_wholegrain,
                                         maha_c_fat, maha_c_upf, maha_c_alcohol, maha_c_sodium,
                                         min_n = 5),

  # Energy-standardized MAHA variants. The original MAHA is retained unchanged
  # as Model 1/default reference. These alternatives test whether total intake drives the
  # association; neither weights nor cut points use mortality outcomes.
  maha_den_c_protein = qscore5_grouped(maha_den_protein, female),
  maha_den_c_dairy = qscore5_grouped(maha_den_dairy, female),
  maha_den_c_veg = qscore5_grouped(maha_den_veg, female),
  maha_den_c_fruit = qscore5_grouped(maha_den_fruit, female),
  maha_den_c_wholegrain = qscore5_grouped(maha_den_wholegrain, female),
  maha_den_c_fat = qscore5_grouped(if (fat_ok) healthy_fat_ratio else fat_proxy_ratio, female),
  maha_den_c_upf = qscore5_grouped(maha_den_upf, female, TRUE),
  maha_den_c_alcohol = if (alcohol_ok) qscore5_grouped(maha_den_alcohol, female, TRUE) else maha_c_alcohol / 2.5 + 1,
  maha_den_c_sodium = if (sodium_ok) qscore5_grouped(maha_den_sodium, female, TRUE) else qscore5_grouped(energy_density(salt_msg_proxy_g,kcal), female, TRUE),
  diet.maha_density.sum = rowmean_vecs(maha_den_c_protein, maha_den_c_dairy, maha_den_c_veg,
                                       maha_den_c_fruit, maha_den_c_wholegrain, maha_den_c_fat,
                                       maha_den_c_upf, maha_den_c_alcohol, maha_den_c_sodium,
                                       min_n = 6),
  maha_res_c_protein = qscore5_grouped(maha_res_protein, female),
  maha_res_c_dairy = qscore5_grouped(maha_res_dairy, female),
  maha_res_c_veg = qscore5_grouped(maha_res_veg, female),
  maha_res_c_fruit = qscore5_grouped(maha_res_fruit, female),
  maha_res_c_wholegrain = qscore5_grouped(maha_res_wholegrain, female),
  maha_res_c_fat = qscore5_grouped(if (fat_ok) healthy_fat_ratio else fat_proxy_ratio, female),
  maha_res_c_upf = qscore5_grouped(maha_res_upf, female, TRUE),
  maha_res_c_alcohol = if (alcohol_ok) qscore5_grouped(maha_res_alcohol, female, TRUE) else maha_c_alcohol / 2.5 + 1,
  maha_res_c_sodium = if (sodium_ok) qscore5_grouped(maha_res_sodium, female, TRUE) else qscore5_grouped(energy_residual_grouped(salt_msg_proxy_g,kcal,female), female, TRUE),
  diet.maha_residual.sum = rowmean_vecs(maha_res_c_protein, maha_res_c_dairy, maha_res_c_veg,
                                        maha_res_c_fruit, maha_res_c_wholegrain, maha_res_c_fat,
                                        maha_res_c_upf, maha_res_c_alcohol, maha_res_c_sodium,
                                        min_n = 6),

  # DASH: code-level whole grains / low-fat dairy / SSB / processed meat when mapped;
  # sodium falls back only if code-level nutrient coverage is insufficient.
  dash_c_fruit = qscore(fruit),
  dash_c_veg = qscore(vegetable),
  dash_c_wholegrain = qscore(whole_grain),
  dash_c_lowfatdairy = qscore(lowfat_dairy),
  dash_c_nutslegumes = qscore(nuts_legumes),
  dash_c_sodium = if (sodium_ok) qscore(sodium_mg, TRUE) else qscore(salt_msg_proxy_g, TRUE),
  dash_c_redmeat = qscore(red_processed_meat_g, TRUE),
  dash_c_sweets = qscore(rowSums(cbind(sweets_pastries_g, ssb_code), na.rm = TRUE), TRUE),
  diet.dash.sum = rowmean_vecs(dash_c_fruit, dash_c_veg, dash_c_wholegrain, dash_c_lowfatdairy,
                               dash_c_nutslegumes, dash_c_sodium, dash_c_redmeat, dash_c_sweets,
                               min_n = 5),

  # Literature-faithful Fung DASH: eight components, sex-specific quintiles,
  # energy-residual inputs, and sweetened beverages (not pastries) as the
  # negative beverage component. Mean and sum are linearly equivalent per SD.
  dash_fung_c_fruit = qscore5_grouped(dash_res_fruit, female),
  dash_fung_c_veg = qscore5_grouped(dash_res_veg, female),
  dash_fung_c_wholegrain = qscore5_grouped(dash_res_wholegrain, female),
  dash_fung_c_lowfatdairy = qscore5_grouped(dash_res_lowfatdairy, female),
  dash_fung_c_nutslegumes = qscore5_grouped(dash_res_nutslegumes, female),
  dash_fung_c_sodium = qscore5_grouped(dash_res_sodium, female, TRUE),
  dash_fung_c_redmeat = qscore5_grouped(dash_res_redmeat, female, TRUE),
  dash_fung_c_ssb = qscore5_grouped(dash_res_ssb, female, TRUE),
  diet.dash_fung.sum = rowmean_vecs(dash_fung_c_fruit, dash_fung_c_veg, dash_fung_c_wholegrain,
                                    dash_fung_c_lowfatdairy, dash_fung_c_nutslegumes,
                                    dash_fung_c_sodium, dash_fung_c_redmeat, dash_fung_c_ssb,
                                    min_n = 8),
  dash_den_c_fruit = qscore5_grouped(dash_den_fruit, female),
  dash_den_c_veg = qscore5_grouped(dash_den_veg, female),
  dash_den_c_wholegrain = qscore5_grouped(dash_den_wholegrain, female),
  dash_den_c_lowfatdairy = qscore5_grouped(dash_den_lowfatdairy, female),
  dash_den_c_nutslegumes = qscore5_grouped(dash_den_nutslegumes, female),
  dash_den_c_sodium = qscore5_grouped(dash_den_sodium, female, TRUE),
  dash_den_c_redmeat = qscore5_grouped(dash_den_redmeat, female, TRUE),
  dash_den_c_ssb = qscore5_grouped(dash_den_ssb, female, TRUE),
  diet.dash_density.sum = rowmean_vecs(dash_den_c_fruit, dash_den_c_veg, dash_den_c_wholegrain,
                                       dash_den_c_lowfatdairy, dash_den_c_nutslegumes,
                                       dash_den_c_sodium, dash_den_c_redmeat, dash_den_c_ssb,
                                       min_n = 8),

  # Use the harmonized 13-component MIND implementation.
  mind_c_green = qscore(green_leafy),
  mind_c_otherveg = qscore(pmax(vegetable - green_leafy, 0)),
  mind_c_berry = qscore(berry),
  mind_c_nuts = qscore(nuts),
  mind_c_wholegrain = qscore(whole_grain),
  mind_c_fish = qscore(fish),
  mind_c_poultry = qscore(poultry),
  mind_c_beans = qscore(legumes),
  mind_c_fat = if (fat_ok) qscore(healthy_fat_ratio) else qscore(fat_proxy_ratio),
  mind_c_meat = qscore(red_processed_meat_g, TRUE),
  mind_c_fried = qscore(fried_fast_g, TRUE),
  mind_c_sweets = qscore(sweets_pastries_g, TRUE),
  mind_c_cheese = qscore(cheese, TRUE),
  diet.mind.sum = rowmean_vecs(mind_c_green, mind_c_otherveg, mind_c_berry, mind_c_nuts,
                               mind_c_wholegrain, mind_c_fish, mind_c_poultry, mind_c_beans,
                               mind_c_fat, mind_c_meat, mind_c_fried, mind_c_sweets, mind_c_cheese,
                               min_n = 8),

  # Mediterranean comparator, with code-level whole grain/fat/alcohol information when available.
  medi_c_fruit = qscore(fruit),
  medi_c_veg = qscore(vegetable),
  medi_c_legumes = qscore(legumes),
  medi_c_wholegrain = qscore(whole_grain),
  medi_c_nuts = qscore(nuts),
  medi_c_fish = qscore(fish),
  medi_c_fat = if (fat_ok) qscore(healthy_fat_ratio) else qscore(fat_proxy_ratio),
  medi_c_alcohol = if (alcohol_ok) moderate_alcohol_score(alcohol_g_food, female) else NA_real_,
  medi_c_meat = qscore(red_processed_meat_g, TRUE),
  medi_c_dairy = qscore(dairy, TRUE),
  diet.medi.sum = rowmean_vecs(medi_c_fruit, medi_c_veg, medi_c_legumes, medi_c_wholegrain,
                               medi_c_nuts, medi_c_fish, medi_c_fat, medi_c_alcohol,
                               medi_c_meat, medi_c_dairy, min_n = 6)
)

active_maha <- "diet.maha.sum"
active_label <- "MAHA"
component_keys <- c("protein","dairy","veg","fruit","wholegrain","fat","upf","alcohol","sodium")
component_cols <- paste0("maha_c_", component_keys)
component_labels <- c(
  protein="Protein", dairy="Dairy", veg="Vegetables", fruit="Fruit",
  wholegrain="Whole grains", fat=if(fat_ok)"Healthy fats" else "Plant-vs-animal fat proxy",
  upf="Low refined/processed foods", alcohol=if(alcohol_ok)"Low alcohol" else "No past-year alcohol proxy",
  sodium=if(sodium_ok)"Low sodium" else "Low salt/MSG seasoning proxy"
)
comparison_cols <- c("diet.dash.sum","diet.mind.sum","diet.medi.sum","diet.chdi.sum")
comparison_labels <- c("DASH","MIND","MEDI","CHDI")
validation_score_set <- c(
  MAHA="diet.maha.sum",
  "MAHA-balanced"="diet.maha_bal.sum",
  "MAHA-strict"="diet.maha_strict.sum",
  "MAHA-no dairy"="diet.maha_nodairy.sum",
  "MAHA-no protein"="diet.maha_noprotein.sum",
  DASH="diet.dash.sum",
  MIND="diet.mind.sum",
  MEDI="diet.medi.sum",
  CHDI="diet.chdi.sum"
)
algorithm_score_set <- c(
  MAHA="diet.maha.sum",
  "MAHA-density"="diet.maha_density.sum",
  "MAHA-residual"="diet.maha_residual.sum",
  DASH="diet.dash_fung.sum",
  "DASH-reference"="diet.dash.sum",
  "DASH-density"="diet.dash_density.sum",
  MIND="diet.mind.sum",
  MEDI="diet.medi.sum"
)
score_algorithm_manifest <- tibble(
  Score=names(algorithm_score_set),Variable=unname(algorithm_score_set),
  Definition=c(
    "Original nine-domain MAHA; retained unchanged",
    "MAHA components per 1000 kcal with sex-specific quintiles",
    "MAHA residual energy adjustment within sex followed by sex-specific quintiles",
    "Fung eight-component DASH; residual energy adjustment, sex-specific quintiles, SSB-only negative domain",
    "Previous CHNS DASH proxy retained for exact Model 1 default reproduction",
    "Fung eight-component DASH using intake per 1000 kcal and sex-specific quintiles",
    "Existing harmonized 13-component MIND",
    "Existing harmonized Mediterranean score"
  ),
  Reference=c(
    "Prespecified original MAHA implementation",
    "Prespecified MAHA energy-density sensitivity",
    "Prespecified MAHA residual-energy sensitivity",
    "Fung et al. Arch Intern Med. 2008;168:713-720. doi:10.1001/archinte.168.7.713",
    "Reference CHNS implementation retained for exact reproduction",
    "Fung domains with a commonly used per-1000-kcal sensitivity implementation",
    "Existing harmonized MIND implementation",
    "Existing harmonized Mediterranean implementation"
  ),
  Outcome_informed=FALSE
)
all_score_vars <- unique(c(unname(validation_score_set),unname(algorithm_score_set)))

score_qc <- tibble(
  variable=unique(c(component_cols,active_maha,unname(algorithm_score_set))),
  finite_n=vapply(unique(c(component_cols,active_maha,unname(algorithm_score_set))),
                  function(v)sum(is.finite(num(dat[[v]]))),integer(1))
)
cat("[CHNS SCORE QC] ",paste0(score_qc$variable,"=",score_qc$finite_n,collapse="; "),"\n",sep="")
anchor_cols <- comparison_cols

# Standardized/rescaled versions used by plots and sensitivity analyses.
all_score_vars <- unique(c(unname(validation_score_set),unname(algorithm_score_set)))
for(v in all_score_vars) if(v %in% names(dat)) {
  stem <- sub("\\.sum$","",v); dat[[paste0(stem,".s100")]] <- score100(dat[[v]]); dat[[paste0(stem,".3c")]] <- f3c(dat[[v]])
}


# 🚩 Model helpers
fit_bin_score <- function(d,outcome,score,covs=c("age","female","province","urban","log_income","education_years","current_smoke","bmi","kcal")) {
  dd <- d; dd$score <- zstd(score)
  covs <- good_covariates(dd,covs)
  use <- c(outcome,"score",covs); dd <- dd %>% dplyr::select(all_of(use)) %>% drop_na()
  if(nrow(dd)<100 || sum(dd[[outcome]]==1)<10 || length(unique(dd$score))<2)
    return(tibble(N=nrow(dd),Events=sum(dd[[outcome]]==1),beta=NA_real_,se=NA_real_,OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,Covariates=paste(covs,collapse=" + ")))
  fm <- as.formula(paste(outcome,"~ score",if(length(covs))paste0(" + ",paste(covs,collapse=" + "))else""))
  fit <- tryCatch(glm(fm,data=dd,family=binomial()),error=function(e)NULL)
  if(is.null(fit)) return(tibble(N=nrow(dd),Events=sum(dd[[outcome]]==1),beta=NA_real_,se=NA_real_,OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,Covariates=paste(covs,collapse=" + ")))
  tt <- broom::tidy(fit) %>% filter(term=="score")
  tibble(N=nrow(dd),Events=sum(dd[[outcome]]==1),beta=tt$estimate[1],se=tt$std.error[1],
         OR=exp(tt$estimate[1]),CI_low=exp(tt$estimate[1]-1.96*tt$std.error[1]),
         CI_high=exp(tt$estimate[1]+1.96*tt$std.error[1]),P=tt$p.value[1],Covariates=paste(covs,collapse=" + "))
}
fit_linear_pref <- function(d,outcome,score,covs=c("age","female","province","urban","log_income","education_years")) {
  dd <- d; dd$score <- zstd(score)
  covs <- good_covariates(dd,covs); dd <- dd %>% dplyr::select(all_of(c(outcome,"score",covs))) %>% drop_na()
  if(nrow(dd)<100) return(tibble(N=nrow(dd),beta=NA_real_,se=NA_real_,P=NA_real_))
  fm <- as.formula(paste(outcome,"~ score",if(length(covs))paste0(" + ",paste(covs,collapse=" + "))else""))
  fit <- lm(fm,data=dd); tt <- broom::tidy(fit) %>% filter(term=="score")
  tibble(N=nrow(dd),beta=tt$estimate[1],se=tt$std.error[1],P=tt$p.value[1])
}
fit_linear_score <- function(d,outcome,score,covs=c("age","female","province","urban","log_income","education_years","current_smoke","bmi","kcal")) {
  dd <- d; dd$score <- zstd(score)
  covs <- good_covariates(dd,covs); dd <- dd %>% dplyr::select(all_of(c(outcome,"score",covs))) %>% drop_na()
  if(nrow(dd)<100 || length(unique(dd$score))<2) return(tibble(N=nrow(dd),beta=NA_real_,se=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,Covariates=paste(covs,collapse=" + ")))
  fm <- as.formula(paste(outcome,"~ score",if(length(covs))paste0(" + ",paste(covs,collapse=" + "))else""))
  fit <- tryCatch(lm(fm,data=dd),error=function(e)NULL)
  if(is.null(fit)) return(tibble(N=nrow(dd),beta=NA_real_,se=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,Covariates=paste(covs,collapse=" + ")))
  tt <- broom::tidy(fit) %>% filter(term=="score")
  tibble(N=nrow(dd),beta=tt$estimate[1],se=tt$std.error[1],CI_low=tt$estimate[1]-1.96*tt$std.error[1],CI_high=tt$estimate[1]+1.96*tt$std.error[1],P=tt$p.value[1],Covariates=paste(covs,collapse=" + "))
}
# Publication helper: direct MAHA-versus-established coefficient comparison in a
# mutually adjusted logistic model. The reported OR ratio is exp(beta_MAHA-beta_trad):
# values >1 indicate that MAHA is less protective than the established comparator.
fit_head2head_bin <- function(d, trad_var, trad_lab, maha_var = "diet.maha.sum",
                             outcome = "death",
                             covs = c("age","female","province","urban","log_income","education_years","current_smoke","bmi","kcal")) {
  dd <- d
  dd$maha_z <- zstd(dd[[maha_var]])
  dd$trad_z <- zstd(dd[[trad_var]])
  covs <- good_covariates(dd, covs)
  dd <- dd %>% dplyr::select(all_of(c(outcome,"maha_z","trad_z",covs))) %>% drop_na()
  if (nrow(dd) < 100 || sum(dd[[outcome]] == 1) < 10) return(tibble())
  fm <- as.formula(paste(outcome, "~ maha_z + trad_z", if(length(covs)) paste0(" + ",paste(covs,collapse=" + ")) else ""))
  fit <- tryCatch(glm(fm, data=dd, family=binomial()), error=function(e) NULL)
  if (is.null(fit)) return(tibble())
  cf <- coef(fit); vc <- vcov(fit)
  if (!all(c("maha_z","trad_z") %in% names(cf))) return(tibble())
  b1 <- unname(cf[["maha_z"]]); b2 <- unname(cf[["trad_z"]])
  dff <- b1 - b2
  se_diff <- sqrt(vc["maha_z","maha_z"] + vc["trad_z","trad_z"] - 2*vc["maha_z","trad_z"])
  tibble(
    Outcome="All-cause mortality", Comparison=paste0("MAHA vs ",trad_lab), Comp2=trad_lab,
    N_total=nrow(dd), N_event=sum(dd[[outcome]]==1),
    OR_maha_mutual=exp(b1), OR_trad_mutual=exp(b2),
    logOR_difference_maha_minus_trad=dff, OR_ratio_maha_vs_trad=exp(dff),
    SE_difference=se_diff,
    LCI_ratio=exp(dff-1.96*se_diff), UCI_ratio=exp(dff+1.96*se_diff),
    P_equal_coefficients=2*pnorm(-abs(dff/se_diff))
  )
}

make_joint_logistic_risk <- function(d, facet_var, line_var, facet_lab, line_lab,
                                     outcome="death",
                                     covs=c("age","female","province","urban","log_income","education_years","current_smoke","bmi","kcal")) {
  dd <- d
  dd$facet <- factor(dd[[facet_var]], levels=c("low","middle","high"))
  dd$line <- factor(dd[[line_var]], levels=c("low","middle","high"))
  covs <- good_covariates(dd,covs)
  keep <- c(outcome,"facet","line",covs)
  dd <- dd %>% dplyr::select(all_of(keep)) %>% drop_na()
  if (nrow(dd) < 100 || sum(dd[[outcome]]==1) < 10) return(tibble())
  fm <- as.formula(paste(outcome,"~ facet + line",if(length(covs))paste0(" + ",paste(covs,collapse=" + "))else""))
  fit <- glm(fm,data=dd,family=binomial())
  nd <- expand.grid(facet=factor(c("low","middle","high"),levels=c("low","middle","high")),
                    line=factor(c("low","middle","high"),levels=c("low","middle","high")))
  for(v in covs) {
    x <- dd[[v]]
    if(is.numeric(x)) nd[[v]] <- mean(x,na.rm=TRUE)
    else if(is.factor(x)) nd[[v]] <- factor(names(which.max(table(x))),levels=levels(x))
    else nd[[v]] <- names(which.max(table(x)))
  }
  pr <- predict(fit,newdata=nd,type="link",se.fit=TRUE)
  nd$risk <- plogis(pr$fit); nd$risk_lower <- plogis(pr$fit-1.96*pr$se.fit); nd$risk_upper <- plogis(pr$fit+1.96*pr$se.fit)
  cells <- dd %>% count(facet,line,name="N")
  nd <- left_join(nd,cells,by=c("facet","line"))
  nd %>% transmute(
    Outcome="All-cause mortality", Followup_years=followup_end-baseline_wave,
    Facet_var=facet_lab, Facet_group=as.character(facet), Line_var=line_lab, Line_group=as.character(line),
    N, risk, risk_lower, risk_upper,
    risk_pct=100*risk, risk_lower_pct=100*risk_lower, risk_upper_pct=100*risk_upper
  )
}

make_precision_table <- function(mort_tbl, margin_hr=1.05, prior_95_hr=1.10) {
  mort_tbl %>% filter(is.finite(beta),is.finite(se),se>0) %>% rowwise() %>% do({
    r <- .
    tost <- ap_tost_one(r$beta, r$se, margin_hr=margin_hr)
    pw <- ap_power_one(r$se, margin_hr=margin_hr)
    bf <- ap_bf01_one(r$beta,r$se,prior_95_hr=prior_95_hr)
    tibble(
      Diet=r$Diet,N_total=r$N,N_event=r$Events,beta=r$beta,se=r$se,OR=r$OR,LCI95=r$CI_low,UCI95=r$CI_high,P=r$P,
      OR_95CI=ap_ci_lab(r$OR,r$CI_low,r$CI_high),
      OR_90CI=ap_ci_lab(r$OR,tost$CI90_low_HR,tost$CI90_high_HR),
      margin_hr=margin_hr,CI90_low_OR=tost$CI90_low_HR,CI90_high_OR=tost$CI90_high_HR,
      TOST_p=tost$TOST_p,Equivalence_5pct=ifelse(isTRUE(tost$equivalent),"Equivalent within OR 0.95-1.05","Not equivalent/inconclusive"),
      prior_95_hr=prior_95_hr,BF01=bf,BF_summary=case_when(bf>=3~"BF favors near-zero",bf<=1/3~"BF favors non-zero",TRUE~"BF inconclusive"),
      MDE_OR_lower=pw$MDE_HR_lower,MDE_OR_upper=pw$MDE_HR_upper,
      MDE80_OR_label=sprintf("%.2f-%.2f",pw$MDE_HR_lower,pw$MDE_HR_upper),
      Powered_5pct=ifelse(isTRUE(pw$powered_for_margin),"Powered for OR 1.05","Not powered for OR 1.05")
    )
  }) %>% ungroup()
}

par_lapply <- function(X,FUN) {
  if(.Platform$OS.type!="windows" && maha_jobs>1) parallel::mclapply(X,FUN,mc.cores=maha_jobs,mc.preschedule=TRUE) else lapply(X,FUN)
}


# 🚩 Primary CHNS mortality model helpers
primary_covariates <- c("age","female","province","urban","log_income","education_years","current_smoke","bmi","kcal")

add_primary_survival_fields <- function(d, mortality_definition="end_wave") {
  d <- as_tibble(d)
  dod_y <- if ("DOD_Y" %in% names(d)) num(d$DOD_Y) else rep(NA_real_, nrow(d))
  dod_rpt <- if ("DOD_RPT" %in% names(d)) num(d$DOD_RPT) else rep(NA_real_, nrow(d))
  death_time <- ifelse(
    is.finite(dod_y) & dod_y > baseline_wave & dod_y <= followup_end,
    dod_y - baseline_wave,
    ifelse(is.finite(dod_rpt) & dod_rpt > baseline_wave & dod_rpt <= followup_end,
           dod_rpt - baseline_wave, NA_real_)
  )
  death_reported <- if ("death" %in% names(d)) num(d$death) == 1 else rep(FALSE,nrow(d))
  has_postbaseline_death_date <- (is.finite(dod_y) & dod_y > baseline_wave) |
    (is.finite(dod_rpt) & dod_rpt > baseline_wave)
  invalid_reported_death <- !is.na(death_reported) & death_reported & !has_postbaseline_death_date
  death_event <- ifelse(invalid_reported_death,NA_integer_,as.integer(is.finite(death_time)))
  if (identical(mortality_definition,"end_wave")) {
    eligible <- (if ("followup_known" %in% names(d)) d$followup_known==1 else is.finite(d$death)) & !invalid_reported_death
    event <- ifelse(eligible,death_event,NA_integer_)
    censor_time <- rep(followup_end-baseline_wave,nrow(d))
  } else if (identical(mortality_definition,"last_contact")) {
    last_contact <- if ("last_contact_year" %in% names(d)) num(d$last_contact_year) else rep(NA_real_,nrow(d))
    eligible <- !invalid_reported_death & (death_event==1 | (is.finite(last_contact)&last_contact>baseline_wave))
    event <- ifelse(eligible,death_event,NA_integer_)
    censor_time <- pmin(last_contact,followup_end,na.rm=TRUE)-baseline_wave
    censor_time[!is.finite(last_contact)] <- NA_real_
  } else stop("Unknown mortality definition: ",mortality_definition,call.=FALSE)
  d$surv_event <- event
  d$surv_time <- ifelse(event == 1, pmax(death_time, 1/365), censor_time)
  d$surv_time[!is.finite(d$surv_time) | d$surv_time <= 0] <- NA_real_
  d$mortality_definition <- mortality_definition
  d$invalid_reported_death_date <- as.integer(invalid_reported_death)
  d
}

q4safe <- function(x) {
  x <- num(x)
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (sum(ok) < 8 || length(unique(x[ok])) < 4) return(factor(out, levels=c("Q1","Q2","Q3","Q4")))
  r <- dplyr::percent_rank(x[ok])
  out[ok] <- as.character(cut(r, c(-Inf,.25,.50,.75,Inf), labels=c("Q1","Q2","Q3","Q4"), right=TRUE))
  factor(out, levels=c("Q1","Q2","Q3","Q4"))
}

primary_filter <- function(d, require_three_days=FALSE, exclude_first_year=FALSE, exclude_baseline_disease=FALSE,
                           mortality_definition="end_wave") {
  dd <- add_primary_survival_fields(d,mortality_definition=mortality_definition) %>% mutate(.row_id=seq_len(nrow(d)))
  if (require_three_days) dd <- dd %>% filter(n_diet_days >= 3)
  if (exclude_first_year) dd <- dd %>% filter(surv_event == 0 | surv_time > 1)
  if (exclude_baseline_disease) {
    disease_vars <- intersect(c("hypertension","t2dm_self","mi_self","stroke_self","cancer_self"), names(dd))
    if (length(disease_vars)) {
      known <- rowSums(is.finite(as.matrix(dd[, disease_vars, drop=FALSE])))
      cases <- rowSums(as.matrix(dd[, disease_vars, drop=FALSE]) == 1, na.rm=TRUE)
      dd$baseline_disease <- ifelse(known == 0, NA_integer_, as.integer(cases > 0))
      dd <- dd %>% filter(baseline_disease == 0)
    }
  }
  dd
}

fit_primary_cox <- function(d, score, covs=primary_covariates, require_three_days=FALSE,
                            no_bmi=FALSE, exclude_first_year=FALSE, exclude_baseline_disease=FALSE,
                            mortality_definition="end_wave") {
  dd <- primary_filter(d,require_three_days,exclude_first_year,exclude_baseline_disease,
                       mortality_definition=mortality_definition)
  score_map <- tibble(.row_id=seq_len(nrow(d)),score_raw=num(score))
  dd <- dd %>% left_join(score_map,by=".row_id") %>% mutate(score_z=zstd(score_raw))
  covs_requested <- unique(as.character(covs))
  if(no_bmi) covs_requested <- setdiff(covs_requested,"bmi")
  covs_use <- intersect(covs_requested,names(dd))
  covs_use <- good_covariates(dd,covs_use)
  covs_dropped <- setdiff(covs_requested,covs_use)
  dd <- dd %>% select(all_of(c("surv_time","surv_event","score_z",covs_use))) %>% drop_na()
  empty <- function() tibble(N=nrow(dd),Events=sum(dd$surv_event==1,na.rm=TRUE),Person_years=sum(dd$surv_time,na.rm=TRUE),
                             beta=NA_real_,se=NA_real_,HR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,
                             Covariates=paste(covs_use,collapse=" + "),
                             Requested_covariates=paste(covs_requested,collapse=" + "),
                             Dropped_covariates=paste(covs_dropped,collapse=" + "))
  if(nrow(dd)<100 || sum(dd$surv_event==1,na.rm=TRUE)<10 || length(unique(dd$score_z))<2) return(empty())
  fm <- as.formula(paste0("Surv(surv_time,surv_event) ~ score_z",
                          if(length(covs_use)) paste0(" + ",paste(covs_use,collapse=" + ")) else ""))
  fit <- tryCatch(coxph(fm,data=dd,ties="efron"),error=function(e)NULL)
  if(is.null(fit)) return(empty())
  tt <- broom::tidy(fit) %>% filter(term=="score_z")
  if(!nrow(tt)) return(empty())
  tibble(N=nrow(dd),Events=sum(dd$surv_event==1,na.rm=TRUE),Person_years=sum(dd$surv_time,na.rm=TRUE),
         beta=tt$estimate[1],se=tt$std.error[1],HR=exp(tt$estimate[1]),
         CI_low=exp(tt$estimate[1]-1.96*tt$std.error[1]),CI_high=exp(tt$estimate[1]+1.96*tt$std.error[1]),
         P=tt$p.value[1],Covariates=paste(covs_use,collapse=" + "),
         Requested_covariates=paste(covs_requested,collapse=" + "),
         Dropped_covariates=paste(covs_dropped,collapse=" + "))
}

fit_primary_quartile <- function(d, score, covs=primary_covariates) {
  dd <- primary_filter(d)
  qmap <- tibble(.row_id=seq_len(nrow(d)),q4=q4safe(score))
  dd <- dd %>% left_join(qmap,by=".row_id")
  covs_use <- good_covariates(dd,intersect(covs,names(dd)))
  dd <- dd %>% select(all_of(c("surv_time","surv_event","q4",covs_use))) %>% drop_na()
  if(nrow(dd)<100 || sum(dd$surv_event==1)<10 || nlevels(droplevels(dd$q4))<4)
    return(tibble(N=nrow(dd),Events=sum(dd$surv_event==1),Q4_vs_Q1=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_))
  fm <- as.formula(paste0("Surv(surv_time,surv_event) ~ q4",
                          if(length(covs_use)) paste0(" + ",paste(covs_use,collapse=" + ")) else ""))
  fit <- tryCatch(coxph(fm,data=dd,ties="efron"),error=function(e)NULL)
  if(is.null(fit)) return(tibble(N=nrow(dd),Events=sum(dd$surv_event==1),Q4_vs_Q1=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_))
  tt <- broom::tidy(fit) %>% filter(term=="q4Q4")
  if(!nrow(tt)) return(tibble(N=nrow(dd),Events=sum(dd$surv_event==1),Q4_vs_Q1=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_))
  tibble(N=nrow(dd),Events=sum(dd$surv_event==1),Q4_vs_Q1=exp(tt$estimate[1]),
         CI_low=exp(tt$estimate[1]-1.96*tt$std.error[1]),CI_high=exp(tt$estimate[1]+1.96*tt$std.error[1]),P=tt$p.value[1])
}

fit_primary_head2head <- function(d, trad_var, trad_lab, maha_var="diet.maha.sum", covs=primary_covariates) {
  dd <- primary_filter(d)
  dd$maha_z <- zstd(dd[[maha_var]])
  dd$trad_z <- zstd(dd[[trad_var]])
  covs_use <- good_covariates(dd,intersect(covs,names(dd)))
  dd <- dd %>% select(all_of(c("surv_time","surv_event","maha_z","trad_z",covs_use))) %>% drop_na()
  if(nrow(dd)<100 || sum(dd$surv_event==1)<10) return(tibble())
  fm <- as.formula(paste0("Surv(surv_time,surv_event) ~ maha_z + trad_z",
                          if(length(covs_use)) paste0(" + ",paste(covs_use,collapse=" + ")) else ""))
  fit <- tryCatch(coxph(fm,data=dd,ties="efron"),error=function(e)NULL)
  if(is.null(fit)) return(tibble())
  cf <- coef(fit); vc <- vcov(fit)
  if(!all(c("maha_z","trad_z")%in%names(cf))) return(tibble())
  b1 <- unname(cf[["maha_z"]]); b2 <- unname(cf[["trad_z"]]); dff <- b1-b2
  se_diff <- sqrt(vc["maha_z","maha_z"]+vc["trad_z","trad_z"]-2*vc["maha_z","trad_z"])
  tibble(Outcome="All-cause mortality",Comparison=paste0("MAHA vs ",trad_lab),Comp2=trad_lab,
         N_total=nrow(dd),N_event=sum(dd$surv_event==1),HR_maha_mutual=exp(b1),HR_trad_mutual=exp(b2),
         logHR_difference_maha_minus_trad=dff,HR_ratio_maha_vs_trad=exp(dff),Ratio=exp(dff),
         SE_difference=se_diff,LCI_ratio=exp(dff-1.96*se_diff),UCI_ratio=exp(dff+1.96*se_diff),
         P_equal_coefficients=2*pnorm(-abs(dff/se_diff)),Effect_measure="HR")
}

make_precision_table_hr <- function(mort_tbl, margin_hr=1.05, prior_95_hr=1.10) {
  mort_tbl %>% filter(is.finite(beta),is.finite(se),se>0) %>% rowwise() %>% do({
    r <- .
    tost <- ap_tost_one(r$beta,r$se,margin_hr=margin_hr)
    pw <- ap_power_one(r$se,margin_hr=margin_hr)
    bf <- ap_bf01_one(r$beta,r$se,prior_95_hr=prior_95_hr)
    tibble(Diet=r$Diet,N_total=r$N,N_event=r$Events,beta=r$beta,se=r$se,Effect=r$HR,HR=r$HR,
           LCI95=r$CI_low,UCI95=r$CI_high,P=r$P,Effect_measure="HR",
           Effect_95CI=ap_ci_lab(r$HR,r$CI_low,r$CI_high),
           Effect_90CI=ap_ci_lab(r$HR,tost$CI90_low_HR,tost$CI90_high_HR),
           margin_effect=margin_hr,CI90_low=tost$CI90_low_HR,CI90_high=tost$CI90_high_HR,
           TOST_p=tost$TOST_p,Equivalence_5pct=ifelse(isTRUE(tost$equivalent),"Equivalent within HR 0.95-1.05","Not equivalent/inconclusive"),
           prior_95_effect=prior_95_hr,BF01=bf,
           BF_summary=case_when(bf>=3~"BF favors near-zero",bf<=1/3~"BF favors non-zero",TRUE~"BF inconclusive"),
           MDE_lower=pw$MDE_HR_lower,MDE_upper=pw$MDE_HR_upper,
           MDE80_label=sprintf("%.2f-%.2f",pw$MDE_HR_lower,pw$MDE_HR_upper),
           Powered_5pct=ifelse(isTRUE(pw$powered_for_margin),"Powered for HR 1.05","Not powered for HR 1.05"))
  }) %>% ungroup()
}


# 🚩 7) FigS1 anti-random validation
step_header("CHNS FigS1: anti-random validation")

# The negative control deliberately keeps the SAME components, SAME participants,
# SAME equal weighting and SAME score scale; only component directions are changed.
# This directly addresses the criticism that a newly defined diet score could be arbitrary.
dval <- dat %>% drop_na(all_of(component_cols))
X <- as.matrix(dval[,component_cols,drop=FALSE])
k <- ncol(X)
if(nrow(dval)>=300 && k>=5) {
  grid <- expand.grid(rep(list(c(-1,1)),k),KEEP.OUT.ATTRS=FALSE)
  obs_i <- which(rowSums(grid==1)==k); nul_i <- setdiff(seq_len(nrow(grid)),obs_i)
  if(null_max>0 && length(nul_i)>null_max) nul_i <- nul_i[unique(round(seq(1,length(nul_i),length.out=null_max)))]
  sel <- c(obs_i,nul_i); grid <- grid[sel,,drop=FALSE]

  # Convergent anchor from established/China-established pattern scores.
  az <- sapply(anchor_cols,function(v) zstd(dval[[v]]))
  if(is.null(dim(az))) az <- matrix(az,ncol=1)
  anchor <- rowMeans(az,na.rm=TRUE); anchor[!is.finite(anchor)] <- NA_real_

  # Completely independent questionnaire preference criterion.
  pref_anchor <- zstd(dval$pref_health_index)

  # Independent observed-behavior criterion: healthier means less processed/deep-fried/
  # sugary-drink exposure and a higher plant-vs-animal fat-source ratio.
  bz <- cbind(-zstd(dval$processed_food_g),-zstd(dval$deepfried_g),-zstd(dval$sugary_drinks_week),zstd(dval$fat_proxy_ratio))
  behavior_anchor <- rowMeans(bz,na.rm=TRUE); behavior_anchor[!is.finite(behavior_anchor)] <- NA_real_

  pat <- bind_rows(par_lapply(seq_len(nrow(grid)),function(i){
    signs <- as.numeric(grid[i,])
    raw <- rowMeans(sweep(X,2,signs,"*"))
    m <- fit_primary_cox(dval,raw)
    z_mort <- if(nrow(m) && is.finite(m$beta[1]) && is.finite(m$se[1]) && m$se[1]>0) -m$beta[1]/m$se[1] else NA_real_
    tibble(
      pattern=paste(ifelse(signs>0,"+","-"),collapse=""),
      is_prespecified=all(signs==1),
      rho_established=safe_cor(raw,anchor),
      rho_preference=safe_cor(raw,pref_anchor),
      rho_behavior=safe_cor(raw,behavior_anchor),
      protective_Z=z_mort,
      mortality_N=if(nrow(m))m$N[1]else NA_real_,
      mortality_events=if(nrow(m))m$Events[1]else NA_real_,
      mortality_HR=if(nrow(m))m$HR[1]else NA_real_
    )
  }))
  obs <- pat %>% filter(is_prespecified); nul <- pat %>% filter(!is_prespecified)
  emp_high <- function(x,obsx) (1+sum(x>=obsx,na.rm=TRUE))/(1+sum(is.finite(x)))
  pct <- function(x,obsx) mean(x<=obsx,na.rm=TRUE)
  sumval <- tibble(
    Cohort="CHNS",Mode=analysis_mode,N=nrow(dval),Components=k,Null_scores=nrow(nul),
    Established_anchor_rho=obs$rho_established,
    Established_anchor_percentile=pct(nul$rho_established,obs$rho_established),
    Established_anchor_empirical_P=emp_high(nul$rho_established,obs$rho_established),
    Preference_rho=obs$rho_preference,
    Preference_percentile=pct(nul$rho_preference,obs$rho_preference),
    Preference_empirical_P=emp_high(nul$rho_preference,obs$rho_preference),
    Behavior_rho=obs$rho_behavior,
    Behavior_percentile=pct(nul$rho_behavior,obs$rho_behavior),
    Behavior_empirical_P=emp_high(nul$rho_behavior,obs$rho_behavior),
    Mortality_protective_Z=obs$protective_Z,
    Mortality_percentile=if(is.finite(obs$protective_Z))pct(nul$protective_Z,obs$protective_Z)else NA_real_,
    Mortality_empirical_P=if(is.finite(obs$protective_Z))emp_high(nul$protective_Z,obs$protective_Z)else NA_real_
  )

  primary_raw <- rowMeans(X)
  loo <- map_dfr(seq_along(component_keys),function(j){
    raw <- rowMeans(X[,-j,drop=FALSE])
    mm <- fit_primary_cox(dval,raw)
    tibble(
      omitted=component_labels[[component_keys[j]]],
      rho_primary=safe_cor(raw,primary_raw),
      rho_established=safe_cor(raw,anchor),
      rho_preference=safe_cor(raw,pref_anchor),
      mortality_HR=if(nrow(mm))mm$HR[1]else NA_real_,
      mortality_CI_low=if(nrow(mm))mm$CI_low[1]else NA_real_,
      mortality_CI_high=if(nrow(mm))mm$CI_high[1]else NA_real_
    )
  })

  p1 <- ggplot(nul,aes(rho_established))+geom_histogram(bins=35)+geom_vline(xintercept=obs$rho_established,linewidth=.9)+
    labs(title="a. Established-score convergence",x="Spearman rho with DASH/MIND/MEDI/CHDI anchor",y="Matched arbitrary scores")+theme_classic(11)
  p2 <- ggplot(nul,aes(rho_preference))+geom_histogram(bins=35)+geom_vline(xintercept=obs$rho_preference,linewidth=.9)+
    labs(title="b. Independent food-preference criterion",x="Spearman rho with healthy preference index",y="Matched arbitrary scores")+theme_classic(11)
  if(is.finite(obs$protective_Z) && sum(is.finite(nul$protective_Z))>=20) {
    p3 <- ggplot(nul,aes(protective_Z))+geom_histogram(bins=35)+geom_vline(xintercept=obs$protective_Z,linewidth=.9)+
      labs(title="c. Mortality criterion",x="Protective Z (-beta/SE)",y="Matched arbitrary scores")+theme_classic(11)
  } else {
    p3 <- ggplot(nul,aes(rho_behavior))+geom_histogram(bins=35)+geom_vline(xintercept=obs$rho_behavior,linewidth=.9)+
      labs(title="c. Independent observed-behavior criterion",x="Spearman rho with healthier behavior anchor",y="Matched arbitrary scores")+theme_classic(11)
  }
  lp <- loo %>% mutate(omitted=factor(omitted,levels=rev(omitted)))
  p4 <- ggplot(lp,aes(rho_primary,omitted))+geom_segment(aes(x=0,xend=rho_primary,yend=omitted))+geom_point()+
    coord_cartesian(xlim=c(0,1))+labs(title="d. Leave-one-domain-out rank stability",x="Spearman rho with primary score",y=NULL)+theme_classic(11)
  FigS1 <- (p1|p2)/(p3|p4)+plot_annotation(title=paste0("CHNS: matched negative-control validation of ",active_label))
  save_plot(FigS1,"FigS1.validate.png",14,10.5,450)
  write_xlsx(list(summary=sumval,direction_null=pat,leave_one_out=loo),"FigS1.validate.out.xlsx")
} else {
  warning("Insufficient complete cases for matched anti-random validation.")
}


# 🚩 8) FigS2 construct profile
step_header("CHNS FigS2: construct profile")
dcp <- dat %>% drop_na(all_of(component_cols),all_of(active_maha))
if(nrow(dcp)>=200) {
  XX <- as.data.frame(dcp[,component_cols]); names(XX)<-component_keys
  XX$score <- dcp[[active_maha]]; XX$decile <- ntile(XX$score,10)
  profile <- XX %>% group_by(decile) %>% summarise(across(all_of(component_keys),mean),.groups="drop") %>%
    pivot_longer(all_of(component_keys),names_to="component",values_to="mean_score") %>%
    mutate(Component=unname(component_labels[component]))
  metrics <- map_dfr(component_keys,function(kk){
    rest <- rowMeans(as.matrix(XX[,setdiff(component_keys,kk),drop=FALSE]))
    tibble(Component=component_labels[[kk]],item_rest_rho=safe_cor(XX[[kk]],rest),
           D10_minus_D1=mean(XX[[kk]][XX$decile==10])-mean(XX[[kk]][XX$decile==1]))
  })
  conv <- tibble(Comparator=comparison_labels,
                 Spearman_rho=vapply(comparison_cols,function(v)safe_cor(dcp[[active_maha]],dcp[[v]]),numeric(1)))
  p1 <- ggplot(profile,aes(decile,mean_score,color=Component,group=Component))+geom_line()+geom_point()+
    labs(title="a. Component profile",x=paste(active_label,"decile"),y="Mean component score",color=NULL)+theme_classic(10)
  mp <- metrics %>% mutate(Component=factor(Component,levels=rev(Component)))
  p2 <- ggplot(mp,aes(D10_minus_D1,Component))+geom_vline(xintercept=0,linetype=2)+geom_segment(aes(x=0,xend=D10_minus_D1,yend=Component))+geom_point()+
    labs(title="b. Decile separation",x="Decile 10 - decile 1",y=NULL)+theme_classic(11)
  p3 <- ggplot(mp,aes(item_rest_rho,Component))+geom_vline(xintercept=0,linetype=2)+geom_point()+
    labs(title="c. Corrected item-rest rho",subtitle="Descriptive only; score is formative",x="Spearman rho",y=NULL)+theme_classic(11)
  p4 <- ggplot(conv,aes(Spearman_rho,reorder(Comparator,Spearman_rho)))+geom_segment(aes(x=0,xend=Spearman_rho,yend=reorder(Comparator,Spearman_rho)))+geom_point()+
    coord_cartesian(xlim=c(-1,1))+labs(title="d. Convergent validity",x="Spearman rho",y=NULL)+theme_classic(11)
  save_plot((p1|p2)/(p3|p4)+plot_annotation(title=paste0("CHNS: ",active_label," construct profile")),
            "FigS2.construct_profile.png",14,10.5,450)
  write_xlsx(list(decile_profile=profile,component_metrics=metrics,convergent_validity=conv),"FigS2.construct_profile.out.xlsx")
}


# 🚩 9) Figprimary independent food-preference validation
step_header("CHNS Figprimary: independent food-preference validation")
pref_outcomes <- c("pref_fruit","pref_veg","pref_fast","pref_salty","pref_ssb","pref_health_index")
pref_lab <- c(pref_fruit="Fruit preference",pref_veg="Vegetable preference",pref_fast="Fast-food preference",
              pref_salty="Salty-snack preference",pref_ssb="Sugary-drink preference",pref_health_index="Healthy preference index")
pref_fit <- map_dfr(pref_outcomes,function(y){
  fit_linear_pref(dat,y,dat[[active_maha]]) %>% mutate(Outcome=pref_lab[[y]],Variable=y)
})
pref_fit$Expected <- ifelse(pref_fit$Variable %in% c("pref_fruit","pref_veg","pref_health_index"),"Positive","Negative")
ddp <- dat %>% filter(is.finite(.data[[active_maha]])) %>% mutate(score_decile=ntile(.data[[active_maha]],10))
pref_profile <- ddp %>% group_by(score_decile) %>% summarise(across(all_of(pref_outcomes),~mean(.x,na.rm=TRUE)),.groups="drop") %>%
  pivot_longer(-score_decile,names_to="Variable",values_to="Mean") %>% mutate(Outcome=unname(pref_lab[Variable]))
pref_corr <- tibble(
  Criterion=c("Healthy preference index","Fruit preference","Vegetable preference","Fast-food preference","Salty-snack preference","Sugary-drink preference"),
  Spearman_rho=c(safe_cor(dat[[active_maha]],dat$pref_health_index),safe_cor(dat[[active_maha]],dat$pref_fruit),
                 safe_cor(dat[[active_maha]],dat$pref_veg),safe_cor(dat[[active_maha]],dat$pref_fast),
                 safe_cor(dat[[active_maha]],dat$pref_salty),safe_cor(dat[[active_maha]],dat$pref_ssb))
)
empty_panel <- function(title, message = "No analyzable score-preference pairs") {
  ggplot() + annotate("text", x=0, y=0, label=message, size=4) +
    xlim(-1,1) + ylim(-1,1) + labs(title=title) + theme_void()
}
if(any(is.finite(pref_fit$beta))) {
  p3a <- pref_fit %>% filter(is.finite(beta),is.finite(se)) %>% mutate(Outcome=factor(Outcome,levels=rev(Outcome))) %>%
    ggplot(aes(beta,Outcome))+geom_vline(xintercept=0,linetype=2)+
    geom_segment(aes(x=beta-1.96*se,xend=beta+1.96*se,yend=Outcome))+geom_point()+
    labs(title="a. Adjusted preference association",x=paste("Change in preference per 1 SD",active_label),y=NULL)+theme_classic(11)
} else p3a <- empty_panel("a. Adjusted preference association")
pref_profile_plot <- pref_profile %>% filter(Variable!="pref_health_index",is.finite(Mean))
if(nrow(pref_profile_plot)) {
  p3b <- pref_profile_plot %>%
    ggplot(aes(score_decile,Mean,group=Outcome))+geom_line()+geom_point()+facet_wrap(~Outcome,scales="free_y")+
    labs(title="b. Preference gradient across score deciles",x=paste(active_label,"decile"),y="Mean preference (1-5)") + theme_classic(9)
} else p3b <- empty_panel("b. Preference gradient across score deciles")
if(any(is.finite(pref_corr$Spearman_rho))) {
  p3c <- pref_corr %>% filter(is.finite(Spearman_rho)) %>%
    ggplot(aes(Spearman_rho,reorder(Criterion,Spearman_rho)))+geom_vline(xintercept=0,linetype=2)+
    geom_segment(aes(x=0,xend=Spearman_rho,yend=reorder(Criterion,Spearman_rho)))+geom_point()+
    coord_cartesian(xlim=c(-1,1))+labs(title="c. Independent convergent validity",x="Spearman rho",y=NULL)+theme_classic(11)
} else p3c <- empty_panel("c. Independent convergent validity")
save_plot((p3a|p3c)/p3b+plot_annotation(title="CHNS: intake-based score validated against independent food preferences"),
          "Figprimary.food_preference_validation.png",14,10.5,450)
write_xlsx(list(adjusted_associations=pref_fit,decile_profile=pref_profile,correlations=pref_corr),"Figprimary.food_preference_validation.out.xlsx")


# 🚩 8) FigS4 independent processing/preparation/seasoning/fat-source validation
step_header("CHNS FigS4: processing, preparation and CFCT subclass validation")
proc_cov <- tibble(
  Variable=c("V39C processed-food classification","V42 preparation method","CFCT6 class/subclass map"),
  Record_coverage=c(mean(is.finite(food0$PROCESSED_CODE)),mean(is.finite(food0$PREP_METHOD)),mean(!is.na(food0$CLASS2))),
  Gram_coverage=c(sum(food0$GRAMS[is.finite(food0$PROCESSED_CODE)],na.rm=TRUE)/sum(food0$GRAMS,na.rm=TRUE),
                  sum(food0$GRAMS[is.finite(food0$PREP_METHOD)],na.rm=TRUE)/sum(food0$GRAMS,na.rm=TRUE),class_cov)
)
proc_vars <- c("processed_food_g","upf_class_proxy_g","deepfried_g","sugary_drinks_week","carbonated_drink_g","salt_msg_proxy_g","fat_proxy_ratio","alcohol_beverage_g")
proc_labs <- c(processed_food_g="V39C processed food (g/day)",upf_class_proxy_g="CFCT processed/refined proxy (g/day)",deepfried_g="Deep-fried food (g/day)",
               sugary_drinks_week="Questionnaire sugary drinks/week",carbonated_drink_g="CFCT carbonated drinks (g/day)",salt_msg_proxy_g="Salt/MSG/other seasoning proxy (g/day)",
               fat_proxy_ratio="Plant/(plant+animal) fat-source ratio",alcohol_beverage_g="Alcoholic beverages (g/day)")
proc_profile <- dat %>% filter(is.finite(.data[[active_maha]])) %>% mutate(decile=ntile(.data[[active_maha]],10)) %>%
  group_by(decile) %>% summarise(across(all_of(proc_vars),~mean(.x,na.rm=TRUE)),.groups="drop") %>%
  pivot_longer(-decile,names_to="Variable",values_to="Value") %>% mutate(Metric=unname(proc_labs[Variable]))
proc_corr <- tibble(Metric=unname(proc_labs[proc_vars]),Variable=proc_vars,
                    Spearman_rho=vapply(proc_vars,function(v)safe_cor(dat[[active_maha]],dat[[v]]),numeric(1)))
p4a <- proc_cov %>% pivot_longer(c(Record_coverage,Gram_coverage),names_to="Coverage",values_to="Value") %>%
  ggplot(aes(Variable,Value,fill=Coverage))+geom_col(position="dodge")+scale_y_continuous(labels=scales::percent)+
  labs(title="a. Independent field/map coverage",x=NULL,y="Coverage",fill=NULL)+theme_classic(10)+theme(axis.text.x=element_text(angle=20,hjust=1))
p4b <- ggplot(proc_profile,aes(decile,Value))+geom_line()+geom_point()+facet_wrap(~Metric,scales="free_y",ncol=4)+
  labs(title="b. External behavioral gradient across MAHA-CN deciles",x=paste(active_label,"decile"),y=NULL)+theme_classic(9)
p4c <- ggplot(proc_corr,aes(Spearman_rho,reorder(Metric,Spearman_rho)))+geom_vline(xintercept=0,linetype=2)+geom_segment(aes(x=0,xend=Spearman_rho,yend=reorder(Metric,Spearman_rho)))+geom_point()+
  coord_cartesian(xlim=c(-1,1))+labs(title="c. Correlation with independent intake/behavior proxies",x="Spearman rho",y=NULL)+theme_classic(10)
save_plot((p4a|p4c)/p4b+plot_annotation(title="CHNS: processing, cooking and food-class construct checks"),
          "FigS4.processing_preparation_validation.png",15,10.5,450)
write_xlsx(list(field_coverage=proc_cov,decile_profile=proc_profile,correlations=proc_corr),"FigS4.processing_preparation_validation.out.xlsx")


# 🚩 11) FigS5 socioeconomic / geographic gradients
step_header("CHNS FigS5: socioeconomic and geographic gradients")
socio_long <- bind_rows(
  if("log_income"%in%names(dat)) tibble(Factor="Household income",rho=safe_cor(dat[[active_maha]],dat$log_income)) else NULL,
  if("education_years"%in%names(dat)) tibble(Factor="Estimated schooling, y",rho=safe_cor(dat[[active_maha]],dat$education_years)) else NULL,
  tibble(Factor="Urban residence (1/2 code)",rho=safe_cor(dat[[active_maha]],num(dat$urban)))
)
prov_score <- dat %>% group_by(province) %>% summarise(N=n(),mean_score=mean(.data[[active_maha]],na.rm=TRUE),
                                                       mean_chdi=mean(diet.chdi.sum,na.rm=TRUE),.groups="drop") %>% filter(N>=30)
p5a <- ggplot(socio_long,aes(rho,reorder(Factor,rho)))+geom_vline(xintercept=0,linetype=2)+geom_segment(aes(x=0,xend=rho,yend=reorder(Factor,rho)))+geom_point()+
  coord_cartesian(xlim=c(-1,1))+labs(title="a. Socioeconomic gradient",x="Spearman rho",y=NULL)+theme_classic(11)
p5b <- ggplot(prov_score,aes(mean_score,reorder(as.character(province),mean_score)))+geom_point(aes(size=N))+
  labs(title="b. Provincial variation",x=paste("Mean",active_label),y="Province code",size="N")+theme_classic(11)
save_plot(p5a|p5b,"FigS5.socioeconomic_geographic.png",12,7,420)
write_xlsx(list(socioeconomic=socio_long,province=prov_score),"FigS5.socioeconomic_geographic.out.xlsx")


# 🚩 12) FigS6 cardiometabolic validation
step_header("CHNS FigS6: cardiometabolic validation")
outcomes <- c(hypertension="Hypertension",t2dm_self="Diabetes",mi_self="Myocardial infarction",
              stroke_self="Stroke",cancer_self="Cancer",obesity="Obesity")
validation_scores <- validation_score_set
phen_all <- map_dfr(names(outcomes),function(y){
  imap_dfr(validation_scores,function(v,lab){
    fit_bin_score(dat,y,dat[[v]],covs=c("age","female","province","urban","log_income","education_years","current_smoke","kcal")) %>%
      mutate(Outcome=outcomes[[y]],Score=lab,Analysis="Cross-sectional at baseline")
  })
})

incident_outcomes <- c(inc_hypertension="Incident hypertension",inc_obesity="Incident obesity",inc_t2dm="Incident diabetes",
                       inc_mi="Incident myocardial infarction",inc_stroke="Incident stroke",inc_cancer="Incident cancer")
inc_all <- map_dfr(names(incident_outcomes),function(y){
  imap_dfr(validation_scores,function(v,lab){
    fit_bin_score(dat,y,dat[[v]],covs=c("age","female","province","urban","log_income","education_years","current_smoke","bmi","kcal")) %>%
      mutate(Outcome=incident_outcomes[[y]],Score=lab,Analysis=paste0("Prospective to ",followup_end))
  })
})

change_outcomes <- c(delta_bmi="Change in BMI",delta_sbp="Change in SBP (mm Hg)",delta_dbp="Change in DBP (mm Hg)")
change_baseline <- c(delta_bmi="bmi",delta_sbp="sbp",delta_dbp="dbp")
change_all <- map_dfr(names(change_outcomes),function(y){
  imap_dfr(validation_scores,function(v,lab){
    cv <- c("age","female","province","urban","log_income","education_years","current_smoke","kcal",change_baseline[[y]])
    fit_linear_score(dat,y,dat[[v]],covs=cv) %>%
      mutate(Outcome=change_outcomes[[y]],Score=lab,Analysis=paste0("Change to ",followup_end))
  })
})

p6a <- phen_all %>% filter(is.finite(OR)) %>% mutate(Outcome=factor(Outcome,levels=rev(unique(Outcome)))) %>%
  ggplot(aes(OR,Outcome,shape=Score))+geom_vline(xintercept=1,linetype=2)+
  geom_segment(aes(x=CI_low,xend=CI_high,yend=Outcome),position=position_dodge(width=.35))+
  geom_point(position=position_dodge(width=.35),size=2.1)+
  labs(title="a. Baseline cardiometabolic face validity",x="Adjusted OR per 1 SD healthier score",y=NULL,shape=NULL)+
  theme_classic(9)+theme(legend.position="bottom")
if(any(is.finite(inc_all$OR))) {
  p6b <- inc_all %>% filter(is.finite(OR)) %>% mutate(Outcome=factor(Outcome,levels=rev(unique(Outcome)))) %>%
    ggplot(aes(OR,Outcome,shape=Score))+geom_vline(xintercept=1,linetype=2)+
    geom_segment(aes(x=CI_low,xend=CI_high,yend=Outcome),position=position_dodge(width=.35))+
    geom_point(position=position_dodge(width=.35),size=2.1)+
    labs(title=paste0("b. Prospective disease validation to ",followup_end),x="Adjusted OR per 1 SD healthier score",y=NULL,shape=NULL)+
    theme_classic(9)+theme(legend.position="bottom")
} else {
  p6b <- ggplot()+annotate("text",x=0,y=0,label=paste0("No usable prospective PEXAM disease outcomes at ",followup_end),size=4)+
    xlim(-1,1)+ylim(-1,1)+labs(title="b. Prospective disease validation")+theme_void()
}
if(any(is.finite(change_all$beta))) {
  p6c <- change_all %>% filter(is.finite(beta)) %>% mutate(Outcome=factor(Outcome,levels=rev(unique(Outcome)))) %>%
    ggplot(aes(beta,Outcome,shape=Score))+geom_vline(xintercept=0,linetype=2)+
    geom_segment(aes(x=CI_low,xend=CI_high,yend=Outcome),position=position_dodge(width=.35))+
    geom_point(position=position_dodge(width=.35),size=2.1)+
    labs(title=paste0("c. Prospective continuous change to ",followup_end),x="Adjusted change per 1 SD healthier score",y=NULL,shape=NULL)+
    theme_classic(9)+theme(legend.position="bottom")
} else {
  p6c <- ggplot()+annotate("text",x=0,y=0,label="No usable paired BMI/BP measurements",size=4)+xlim(-1,1)+ylim(-1,1)+theme_void()
}
save_plot((p6a|p6b)/p6c,"FigS6.cardiometabolic_validation.png",14,11,420)
write_xlsx(list(cross_sectional=phen_all,prospective_incident=inc_all,prospective_change=change_all),"FigS6.cardiometabolic_validation.out.xlsx")


# 🚩 13) FigS7 distributions/concordance
step_header("CHNS FigS7: score distributions and concordance")
score_vars <- validation_score_set[names(validation_score_set) != "CHDI"]
score_df <- dat %>% dplyr::select(all_of(unname(score_vars)))
names(score_df) <- names(score_vars)
score_long <- score_df %>% mutate(across(everything(),score100)) %>% pivot_longer(everything(),names_to="Diet",values_to="Score")
density_curve <- score_long %>% filter(is.finite(Score)) %>% group_by(Diet) %>% group_modify(~ {
  dn <- density(.x$Score, from=0, to=100, n=256, na.rm=TRUE)
  tibble(score=dn$x,density=dn$y)
}) %>% ungroup()
cm <- cor(score_df,use="pairwise.complete.obs",method="spearman")
p7a <- ggplot(score_long,aes(Score))+geom_density()+facet_wrap(~Diet,scales="free_y",ncol=1)+
  labs(title="a. Score distributions",x="Rescaled 0-100",y="Density")+theme_classic(10)
cd <- as.data.frame(as.table(cm))
p7b <- ggplot(cd,aes(Var2,Var1,fill=Freq))+geom_tile()+geom_text(aes(label=sprintf("%.2f",Freq)))+
  scale_fill_gradient2(midpoint=0,limits=c(-1,1))+labs(title="b. Spearman concordance",x=NULL,y=NULL,fill=NULL)+theme_minimal(10)+
  theme(axis.text.x=element_text(angle=35,hjust=1),panel.grid=element_blank())
save_plot(p7a|p7b,"FigS7.score_distribution_concordance.png",11,8,420)
write_xlsx(list(summary=score_long%>%group_by(Diet)%>%summarise(N=sum(is.finite(Score)),mean=mean(Score,na.rm=TRUE),sd=sd(Score,na.rm=TRUE),.groups="drop"), density_curve=density_curve,
                spearman=as.data.frame(cm)),"FigS7.score_distribution_concordance.out.xlsx")


# 🚩 12) FigS8 food-map QC / unresolved code-level components
step_header("CHNS FigS8: FOODCODE map QC")
class_prev <- food0 %>% group_by(CLASS2,CLASS_NAME_CN) %>% summarise(grams=sum(GRAMS,na.rm=TRUE),records=n(),.groups="drop") %>%
  mutate(gram_fraction=grams/sum(grams,na.rm=TRUE)) %>% arrange(desc(gram_fraction))
unmatched <- partial_map %>% filter(is.na(CLASS2))
missing_code <- tibble(Field=required_code_fields) %>% mutate(
  gram_coverage=ifelse(Field %in% map_field_coverage$field,map_field_coverage$gram_coverage[match(Field,map_field_coverage$field)],0),
  Status=ifelse(gram_coverage>=full_map_min_cov,"coverage meets threshold","below threshold / unavailable")
)
p8a <- class_prev %>% slice_head(n=15) %>% mutate(Group=factor(paste0(CLASS2," ",CLASS_NAME_CN),levels=rev(paste0(CLASS2," ",CLASS_NAME_CN)))) %>%
  ggplot(aes(gram_fraction,Group))+geom_col()+scale_x_continuous(labels=scales::percent)+labs(title="a. Baseline food grams by CFCT class",x="Fraction of food grams",y=NULL)+theme_classic(10)
p8b <- map_field_coverage
if(!nrow(p8b)) p8b <- missing_code %>% transmute(field=Field,gram_coverage=gram_coverage)
p8bplot <- p8b %>% mutate(field=factor(field,levels=rev(field))) %>% ggplot(aes(gram_coverage,field))+geom_col()+scale_x_continuous(labels=scales::percent,limits=c(0,1))+
  labs(title="b. Code-level nutrient/fine-group coverage",subtitle=paste0("Threshold = ",scales::percent(full_map_min_cov),"; ",analysis_mode),x="Gram-weighted coverage",y=NULL)+theme_classic(9)
save_plot(p8a|p8bplot,"FigS8.food_mapping_qc.png",13,7.5,420)
write_xlsx(list(foodcode_audit=partial_map,class_prevalence=class_prev,unmatched_codes=unmatched,
                foodcode_field_coverage=map_field_coverage,unresolved_foodcode_fields=missing_code,
                metadata=tibble(map_mode=analysis_mode,class_gram_coverage=class_cov,foodcode_map_coverage=full_map_coverage)),
           "FigS8.food_mapping_qc.out.xlsx")


# 🚩 15) Main mortality comparison
step_header("CHNS main mortality comparison")

dm <- dat %>% filter(followup_known==1) %>% drop_na(death)
mort_scores <- validation_score_set

mort <- imap_dfr(mort_scores,function(v,lab) {
  fit_primary_cox(dm,dm[[v]]) %>%
    mutate(Diet=lab,Effect=HR,Effect_measure="HR",Method="Cox",.before=1)
})

# Reproducibility audit only.  Event/N changes are now expected when the user
# changes wave, age, or follow-up definition, so no result-direction guard is
# allowed to stop the analysis.
.guard_maha <- mort %>% filter(Diet=="MAHA") %>% slice(1)
.guard_cmp <- mort %>% filter(Diet %in% c("DASH","MIND","MEDI"))
primary_reproducibility_audit <- tibble(
  Expected_configuration="2011 baseline, age >=40, end-wave complete case",
  Expected_N=7387L,Expected_deaths=219L,
  Observed_N=if(nrow(.guard_maha)).guard_maha$N[1]else NA_integer_,
  Observed_deaths=if(nrow(.guard_maha)).guard_maha$Events[1]else NA_integer_,
  Exact_match=isTRUE(nrow(.guard_maha)==1L && .guard_maha$N[1]==7387L && .guard_maha$Events[1]==219L),
  Established_all_HR_below_1=isTRUE(nrow(.guard_cmp)==3L && all(is.finite(.guard_cmp$HR)) && all(.guard_cmp$HR<1)),
  Note="Audit only; neither N/events nor association direction is used to select a model."
)

quartiles <- imap_dfr(mort_scores,function(v,lab) {
  fit_primary_quartile(dm,dm[[v]]) %>%
    mutate(Diet=lab,Effect_measure="HR",Method="Cox",.before=1)
})

pub_precision <- make_precision_table_hr(mort)
pub_head2head <- bind_rows(
  fit_primary_head2head(dm,"diet.dash.sum","DASH"),
  fit_primary_head2head(dm,"diet.mind.sum","MIND"),
  fit_primary_head2head(dm,"diet.medi.sum","MEDI")
)

reference_event_n <- if(nrow(mort)&&any(is.finite(mort$Events))) max(mort$Events,na.rm=TRUE) else NA_integer_
pM <- mort %>% mutate(Diet=factor(Diet,levels=rev(names(mort_scores)))) %>%
  ggplot(aes(HR,Diet))+geom_vline(xintercept=1,linetype=2)+
  geom_segment(aes(x=CI_low,xend=CI_high,yend=Diet))+geom_point(size=2.7)+
  labs(title="CHNS all-cause mortality — Model 1 (default)",
       subtitle=paste0("Fixed 2011 baseline; reference algorithms; ",reference_event_n," deaths in the complete-case model"),
       x="Adjusted HR per 1 SD healthier score",y=NULL)+theme_classic(12)
save_plot(pM,"Fig6.mortality_validation.png",8.5,6.5,420)
write_xlsx(list(mortality=mort,quartiles=quartiles,precision=pub_precision,head2head=pub_head2head,
                model_definition=tibble(Model="Model1_default_reference",Baseline=2011L,
                  Design="Fixed baseline",Energy_covariate=TRUE,
                  Note="Exact submitted reference benchmark; Models 4-5 are prespecified high-event cumulative-average sensitivities"),
                reproducibility_audit=primary_reproducibility_audit),
           "Fig6.mortality_validation.out.xlsx")

sensitivity_specs <- tribble(
  ~Analysis,                      ~require_three_days, ~no_bmi, ~exclude_first_year, ~exclude_baseline_disease,
  "Primary",                      FALSE,               FALSE,   FALSE,               FALSE,
  "All 3 diet-record days",       TRUE,                FALSE,   FALSE,               FALSE,
  "Exclude deaths in first year", FALSE,               FALSE,   TRUE,                FALSE,
  "Exclude baseline disease",     FALSE,               FALSE,   FALSE,               TRUE
)

mort_sens <- pmap_dfr(sensitivity_specs,function(Analysis,require_three_days,no_bmi,exclude_first_year,exclude_baseline_disease) {
  imap_dfr(mort_scores,function(v,lab) {
    fit_primary_cox(dm,dm[[v]],require_three_days=require_three_days,no_bmi=no_bmi,
                    exclude_first_year=exclude_first_year,exclude_baseline_disease=exclude_baseline_disease) %>%
      mutate(Diet=lab,Analysis=Analysis,.before=1)
  })
})

# Evaluate the prespecified mortality model grid without outcome-based selection.
grid_scores <- c(MAHA="diet.maha.sum",DASH="diet.dash.sum",MIND="diet.mind.sum",MEDI="diet.medi.sum")
# Focused mortality set: the submitted fixed-wave reference plus the selected
# high-event socioeconomic model and its BMI sensitivity.  No other alternative
# model is fitted or written to the final workbooks.
adjustment_specs <- list(
  Model4_socioeconomic=c("age","female","province","urban","log_income","education_years"),
  Model5_bmi=c("age","female","province","urban","log_income","education_years","bmi")
)
final_model_roles <- c(
  Model4_socioeconomic="Primary high-event sensitivity",
  Model5_bmi="BMI / potential-overadjustment sensitivity"
)
final_model_notes <- c(
  Model4_socioeconomic="Age >=40 cumulative-average exposure across 2004-2011; adjusts for province, urban residence, household income and education",
  Model5_bmi="Adds BMI to the socioeconomic model; interpreted as a possible mediator/overadjustment sensitivity"
)
mortality_definitions <- tibble(
  Mortality_definition=c("end_wave","last_contact"),
  Definition_label=c(
    "Observed at end wave or death",
    "Censor at last observed contact"
  ),
  Publication_eligible=c(TRUE,TRUE),
  Definition_note=c(
    "Conservative complete follow-up definition used in the submitted draft",
    "Retains participants lost before the end wave and censors them at their last observed wave"
  )
)
sample_specs <- "four_score_common"

mortality_flow <- pmap_dfr(mortality_definitions,function(Mortality_definition,Definition_label,Publication_eligible,Definition_note) {
  ff <- add_primary_survival_fields(dat,Mortality_definition) %>% filter(is.finite(surv_time),surv_event%in%c(0,1))
  tibble(
    Baseline_wave=baseline_wave,Minimum_age=min_age,
    Mortality_definition=Mortality_definition,Definition_label=Definition_label,
    Publication_eligible=Publication_eligible,Definition_note=Definition_note,
    N_before_covariates=nrow(ff),Deaths_before_covariates=sum(ff$surv_event==1),
    Person_years_before_covariates=sum(ff$surv_time),
    Events_ge_target=sum(ff$surv_event==1)>=mortality_event_target,
    Event_target=mortality_event_target
  )
})

mortality_grid <- pmap_dfr(mortality_definitions,function(Mortality_definition,Definition_label,Publication_eligible,Definition_note) {
  imap_dfr(adjustment_specs,function(covs,adj_lab) {
    map_dfr(sample_specs,function(sample_lab) {
      dd0 <- dat
      if(identical(sample_lab,"four_score_common")) {
        dd0 <- dd0 %>% filter(if_all(all_of(unname(grid_scores)),~is.finite(num(.x))))
      }
      imap_dfr(grid_scores,function(v,diet_lab) {
        fit_primary_cox(dd0,dd0[[v]],covs=covs,mortality_definition=Mortality_definition) %>%
          mutate(
            Baseline_wave=baseline_wave,Minimum_age=min_age,Diet=diet_lab,
            Mortality_definition=Mortality_definition,Definition_label=Definition_label,
            Publication_eligible=Publication_eligible,Definition_note=Definition_note,
            Adjustment=adj_lab,Sample_strategy=sample_lab,
            Event_target=mortality_event_target,Events_ge_target=Events>=mortality_event_target,
            .before=1
          )
      })
    })
  })
})

mortality_pattern <- mortality_grid %>%
  select(Baseline_wave,Minimum_age,Mortality_definition,Definition_label,Publication_eligible,
         Adjustment,Sample_strategy,Diet,HR,P,N,Events,Person_years) %>%
  pivot_wider(names_from=Diet,values_from=c(HR,P,N,Events,Person_years)) %>%
  mutate(
    Established_all_protective=is.finite(HR_DASH)&is.finite(HR_MIND)&is.finite(HR_MEDI)&
      HR_DASH<1&HR_MIND<1&HR_MEDI<1,
    Established_all_significant_protective=Established_all_protective&
      P_DASH<.05&P_MIND<.05&P_MEDI<.05,
    MAHA_attenuated_or_null_vs_all=is.finite(HR_MAHA)&is.finite(HR_DASH)&is.finite(HR_MIND)&is.finite(HR_MEDI)&
      HR_MAHA>=pmax(HR_DASH,HR_MIND,HR_MEDI),
    MAHA_null_or_risky=is.finite(HR_MAHA)&is.finite(P_MAHA)&(HR_MAHA>=1|P_MAHA>=.05),
    Descriptive_expected_pattern=Established_all_significant_protective&MAHA_null_or_risky&MAHA_attenuated_or_null_vs_all,
    Events_ge_target=Events_MAHA>=mortality_event_target,
    Selection_rule="Diagnostic only: HR direction is never used to choose a primary model"
  )

p_grid <- mortality_grid %>%
  filter(Adjustment=="Model5_bmi",Sample_strategy=="four_score_common",is.finite(HR)) %>%
  mutate(Definition_label=factor(Definition_label,levels=rev(mortality_definitions$Definition_label)),
         Diet=factor(Diet,levels=c("MAHA","DASH","MIND","MEDI"))) %>%
  ggplot(aes(HR,Definition_label,color=Diet))+
  geom_vline(xintercept=1,linetype=2,color="grey55")+
  geom_segment(aes(x=CI_low,xend=CI_high,yend=Definition_label),position=position_dodge(width=.62),linewidth=.65)+
  geom_point(position=position_dodge(width=.62),size=2.5)+
  labs(title=paste0("CHNS mortality-definition grid: baseline ",baseline_wave,", age >=",min_age),
       subtitle="Model 5 BMI sensitivity (no energy covariate); common sample across all four scores",
       x="HR per 1 SD healthier score",y=NULL,color=NULL)+
  theme_classic(11)+theme(legend.position="bottom")
save_plot(p_grid,"mortality_model_grid.png",11,7.5,420)
write_xlsx(list(flow=mortality_flow,results=mortality_grid,pattern_diagnostic=mortality_pattern),
           "mortality_model_grid.out.xlsx")

analysis_metadata <- tibble(
  schooling_definition="Estimated schooling, y: A11 0=0; 11-16=1-6; 21-26=7-12; 27-29=10-12; 31-36=13-18; other codes missing. EDUYRS retained as recorded. A11-derived years are estimated, not individually measured.",
  primary_model=analysis_spec,
  baseline_wave=baseline_wave,
  followup_end=followup_end,
  min_age=min_age,
  min_diet_days=min_diet_days,
  association_method="Cox",
  association_measure="HR",
  primary_covariates=paste(primary_covariates,collapse=" + "),
  primary_complete_case_N=mort$N[mort$Diet=="MAHA"][1],
  primary_deaths=mort$Events[mort$Diet=="MAHA"][1]
)

benchmark <- tibble(
  N=analysis_metadata$primary_complete_case_N,
  Deaths=analysis_metadata$primary_deaths,
  Baseline=baseline_wave,
  Followup_end=followup_end
)

write_xlsx(
  list(
    mortality=mort,
    quartiles=quartiles,
    precision=pub_precision,
    head2head=pub_head2head,
    sensitivity=mort_sens,
    benchmark=benchmark,
    reproducibility_audit=primary_reproducibility_audit,
    analysis_metadata=analysis_metadata
  ),
  "publication_inputs.out.xlsx"
)


# 🚩 16) Follow-up QC
fu <- dat %>% summarise(
  Baseline_N=n(),
  Known_followup=sum(followup_known==1,na.rm=TRUE),
  Observed_in_end_wave_roster=sum(observed_at_followup==1,na.rm=TRUE),
  Observed_after_baseline=sum(observed_after_baseline==1,na.rm=TRUE),
  Deaths=sum(death_reported==1,na.rm=TRUE),
  Unknown_at_end_wave=sum(followup_known!=1|is.na(followup_known),na.rm=TRUE),
  Note="Roster presence is contact evidence, not a direct alive-status indicator"
)
fu_prov <- dat %>% group_by(province) %>% summarise(
  N=n(),known=sum(followup_known==1,na.rm=TRUE),known_pct=known/N,
  deaths=sum(death_reported==1,na.rm=TRUE),.groups="drop"
)
p9 <- ggplot(fu_prov,aes(known_pct,reorder(as.character(province),known_pct)))+
  geom_col()+scale_x_continuous(labels=percent)+
  labs(title="Follow-up ascertainment by province",x="Known vital status",y="Province code")+theme_classic(11)
save_plot(p9,"FigS9.followup_qc.png",8,6.5,420)
write_xlsx(list(overall=fu,by_province=fu_prov),"FigS9.followup_qc.out.xlsx")


# 🚩 17) Mortality sensitivity figure
p10 <- mort_sens %>%
  filter(Diet %in% c("MAHA","DASH","MIND","MEDI"),is.finite(HR)) %>%
  mutate(Analysis=factor(Analysis,levels=rev(sensitivity_specs$Analysis))) %>%
  ggplot(aes(HR,Analysis,color=Diet))+
  geom_vline(xintercept=1,linetype=2,color="grey55")+
  geom_segment(aes(x=CI_low,xend=CI_high,yend=Analysis),position=position_dodge(width=.55),linewidth=.7)+
  geom_point(position=position_dodge(width=.55),size=2.5)+
  labs(title="CHNS all-cause mortality sensitivity analyses",x="Adjusted HR per 1 SD healthier score",y=NULL,color=NULL)+
  theme_classic(11)+theme(legend.position="bottom")
save_plot(p10,"FigS10.mortality_sensitivity.png",10.5,7.5,420)
write_xlsx(list(results=mort_sens,metadata=analysis_metadata),"FigS10.mortality_sensitivity.out.xlsx")


# 🚩 18) Table 1, metadata, reusable dataset
active_3c <- paste0(sub("\\.sum$","",active_maha),".3c")
tab1 <- dat %>% mutate(Diet_group=.data[[active_3c]]) %>% group_by(Diet_group) %>% summarise(
  N=n(),age_mean=mean(age,na.rm=TRUE),female_pct=mean(female,na.rm=TRUE),BMI_mean=mean(bmi,na.rm=TRUE),
  kcal_mean=mean(kcal,na.rm=TRUE),income_mean=mean(income_pc,na.rm=TRUE),education_mean=mean(education_years,na.rm=TRUE),
  known_followup_pct=mean(followup_known==1,na.rm=TRUE),mortality_pct=mean(death,na.rm=TRUE),.groups="drop"
)
meta <- tibble(
  schooling_definition=analysis_metadata$schooling_definition,
  primary_model=analysis_spec,
  baseline_wave=baseline_wave,followup_end=followup_end,followup_gap=followup_end-baseline_wave,
  min_age=min_age,min_diet_days=min_diet_days,association_method="Cox",active_score=active_label,
  nutr3=src$NUTR3[["path"]],c12diet=src$C12DIET[["path"]],master=src$MASTER[["path"]],
  roster=src$ROSTER[["path"]],pexam=src$PEXAM[["path"]],pact=src$PACT[["path"]],
  food_map=ifelse(is.na(food_map_file),"",food_map_file),
  nutrient_map=ifelse(is.na(nutrient_map_file),"",nutrient_map_file),
  classification_mapping_coverage=class_cov,foodcode_mapping_coverage=full_map_coverage,
  code_map_available=code_map_available,enhanced_map_mode=enhanced_map_mode,full_map_mode=full_map_mode,
  foodcode_field_threshold=full_map_min_cov
)
write_xlsx(list(Table1=tab1,metadata=meta,input_sources=source_manifest,benchmark=benchmark),"Table1.out.xlsx")
saveRDS(dat,cohort_file("analysis_dataset.rds"))


# 🚩 19) Focused high-event cumulative-average sensitivity
# Only age >=40 Models 4 and 5 are retained.  The earlier exploratory age,
# first-entry, landmark, minimal-adjustment and all-at-end grids are not fitted.
if(run_mortality_sweep==1L && !mortality_sweep_child) {
  step_header("CHNS focused high-event mortality sensitivity")
  stale_focused_outputs <- cohort_file(c(
    "mortality_design_sweep.out.xlsx","mortality_design_sweep.png",
    "mortality_competing_models.out.xlsx","mortality_competing_all4.png",
    "mortality_competing_pairwise.png"
  ))
  unlink(stale_focused_outputs[file.exists(stale_focused_outputs)],force=TRUE)

  parse_int_list <- function(x) {
    z <- suppressWarnings(as.integer(trimws(strsplit(x,",",fixed=TRUE)[[1]])))
    sort(unique(z[is.finite(z)]))
  }
  eligible_sweep_waves <- common_baseline[common_baseline<followup_end]
  if(!identical(tolower(mortality_sweep_wave_req),"auto")) {
    eligible_sweep_waves <- intersect(eligible_sweep_waves,parse_int_list(mortality_sweep_wave_req))
  }
  sweep_ages <- mortality_model_age
  if(!length(eligible_sweep_waves)) {
    warning("No eligible common waves for CHNS mortality sweep.")
  } else {
    sweep_root <- file.path(maha_auxdir,"mortality_sweep_runs")
    dir.create(sweep_root,recursive=TRUE,showWarnings=FALSE)
    rscript_bin <- Sys.getenv("RSCRIPT_BIN",unset=Sys.which("Rscript"))
    if(!nzchar(rscript_bin)) stop("Cannot locate Rscript for mortality sweep.",call.=FALSE)

    sweep_runs <- list()
    run_manifest <- list()
    md5_or_blank <- function(x) if(length(x)==1L&&!is.na(x)&&file.exists(x))unname(tools::md5sum(x))else""
    kk <- 0L
    for(aa in sweep_ages) for(ww in eligible_sweep_waves) {
      kk <- kk+1L
      run_id <- paste0("wave",ww,"_age",aa)
      child_root <- file.path(sweep_root,run_id)
      child_aux <- file.path(child_root,"chns")
      child_grid_file <- file.path(child_aux,"mortality_model_grid.out.xlsx")
      child_data_file <- file.path(child_aux,"analysis_dataset.rds")
      child_log <- file.path(sweep_root,paste0(run_id,".log"))
      signature_file <- file.path(child_root,"run_signature.txt")
      run_signature <- paste(md5_or_blank(.this_file),ww,aa,followup_end,min_diet_days,full_map_min_cov,
                             require_nutrient_map,require_enhanced_map,
                             md5_or_blank(food_map_file),md5_or_blank(nutrient_map_file),sep="|")
      cache_ok <- file.exists(child_grid_file)&&file.exists(child_data_file)&&file.exists(signature_file)&&
        identical(readLines(signature_file,warn=FALSE,n=1L),run_signature)

      use_parent <- identical(as.integer(ww),as.integer(baseline_wave)) && identical(as.integer(aa),as.integer(min_age))
      status <- 0L
      if(use_parent) {
        child_grid_file <- cohort_file("mortality_model_grid.out.xlsx")
        child_data_file <- cohort_file("analysis_dataset.rds")
      } else if(!cache_ok) {
        dir.create(child_root,recursive=TRUE,showWarnings=FALSE)
        env_child <- c(
          paste0("MAHA_OUTDIR=",shQuote(child_root)),
          paste0("MAHA_AUXDIR=",shQuote(child_aux)),
          paste0("CHNS_BASELINE_WAVE=",ww),
          paste0("CHNS_FOLLOWUP_END=",followup_end),
          paste0("CHNS_MIN_AGE=",aa),
          paste0("CHNS_MIN_DIET_DAYS=",min_diet_days),
          paste0("CHNS_MIN_MAP_COVERAGE=",full_map_min_cov),
          paste0("CHNS_REQUIRE_NUTRIENT_MAP=",require_nutrient_map),
          paste0("CHNS_REQUIRE_ENHANCED_MAP=",require_enhanced_map),
          paste0("CHNS_PROJECT_DIR=",shQuote(project_dir)),
          paste0("CHNS_FOOD_MAP=",shQuote(food_map_file)),
          paste0("CHNS_NUTRIENT_MAP=",shQuote(nutrient_map_file)),
          paste0("RSCRIPT_BIN=",shQuote(rscript_bin)),
          "CHNS_MORTALITY_SWEEP_CHILD=1",
          "CHNS_RUN_MORTALITY_SWEEP=0",
          "MAHA_NULL_MAX=1"
        )
        status <- tryCatch(
          system2(rscript_bin,args=shQuote(.this_file),env=env_child,stdout=child_log,stderr=child_log),
          error=function(e) 999L
        )
        if(identical(as.integer(status),0L)&&file.exists(child_grid_file)&&file.exists(child_data_file))
          writeLines(run_signature,signature_file,useBytes=TRUE)
      }
      ok <- identical(as.integer(status),0L) && file.exists(child_grid_file) && file.exists(child_data_file)
      failure_tail <- if(!ok && file.exists(child_log)) {
        paste(utils::tail(readLines(child_log,warn=FALSE),20L),collapse=" | ")
      } else ""
      run_manifest[[kk]] <- tibble(
        Run=run_id,Baseline_wave=ww,Minimum_age=aa,Followup_end=followup_end,
        Status=ifelse(ok,"complete",paste0("failed_status_",status)),
        Grid_file=child_grid_file,Dataset_file=child_data_file,Log_file=child_log,
        Failure_tail=failure_tail
      )
      if(ok) sweep_runs[[run_id]] <- list(wave=ww,age=aa,grid=child_grid_file,data=child_data_file)
      cat("[MORTALITY SWEEP] ",run_id," : ",ifelse(ok,"complete","failed"),"\n",sep="")
      if(!ok && nzchar(failure_tail)) cat("[MORTALITY SWEEP FAILURE] ",run_id," : ",failure_tail,"\n",sep="")
    }
    run_manifest <- bind_rows(run_manifest)

    # Build the single prespecified age >=40 repeated-exposure panel.
    panel_by_age <- set_names(sweep_ages,as.character(sweep_ages)) %>% map(function(panel_age_i) {
      panel_runs_i <- sweep_runs[vapply(sweep_runs,function(x)x$age==panel_age_i,logical(1))]
      imap_dfr(panel_runs_i,function(rr,run_id) {
        readRDS(rr$data) %>% mutate(exposure_wave=as.integer(rr$wave),sweep_run=run_id)
      }) %>% arrange(IDIND,exposure_wave) %>% distinct(IDIND,exposure_wave,.keep_all=TRUE)
    })
    panel_by_age <- panel_by_age[vapply(panel_by_age,nrow,integer(1))>0L]
    if(!length(panel_by_age)) stop("No successful datasets were available for a repeated-wave panel.",call.=FALSE)
    available_panel_ages <- as.integer(names(panel_by_age))
    if(!mortality_model_age%in%available_panel_ages) {
      replacement_age <- available_panel_ages[which.min(abs(available_panel_ages-mortality_model_age))]
      warning("Requested named-model age ",mortality_model_age," unavailable; using age ",replacement_age,".")
      panel_age <- replacement_age
    } else panel_age <- mortality_model_age
    panel <- panel_by_age[[as.character(panel_age)]]

    cummean_known <- function(x) {
      x <- num(x); n <- cumsum(is.finite(x)); s <- cumsum(ifelse(is.finite(x),x,0));
      out <- s/pmax(n,1); out[n==0] <- NA_real_; out
    }
    focused_algorithm_score_set <- algorithm_score_set[c(
      "MAHA","MAHA-density","MAHA-residual","DASH-reference","MIND","MEDI"
    )]
    panel_cum_by_age <- map(panel_by_age,function(panel_i) {
      sweep_score_vars_i <- unique(c(unname(grid_scores),unname(focused_algorithm_score_set)))
      sweep_score_vars_i <- sweep_score_vars_i[sweep_score_vars_i %in% names(panel_i)]
      panel_i %>% group_by(IDIND) %>% arrange(exposure_wave,.by_group=TRUE) %>%
        mutate(across(all_of(sweep_score_vars_i),cummean_known,.names="{.col}_cum")) %>% ungroup()
    })
    panel_cum <- panel_cum_by_age[[as.character(panel_age)]]

    build_time_updated <- cache_identical_inputs(function(d,mortality_definition) {
      d %>% group_by(IDIND) %>% arrange(exposure_wave,.by_group=TRUE) %>% group_modify(~{
        g <- .x
        dd <- c(num(g$DOD_Y),num(g$DOD_RPT)); dd <- dd[is.finite(dd)&dd>min(g$exposure_wave)&dd<=followup_end]
        dy <- if(length(dd))min(dd)else NA_real_
        lc0 <- num(g$last_contact_year); lc <- if(any(is.finite(lc0)))max(lc0,na.rm=TRUE)else NA_real_
        if(identical(mortality_definition,"end_wave")) {
          eligible <- any(g$observed_at_followup==1,na.rm=TRUE)|is.finite(dy); analysis_end <- followup_end
        } else if(identical(mortality_definition,"last_contact")) {
          eligible <- is.finite(dy)|(is.finite(lc)&lc>min(g$exposure_wave)); analysis_end <- if(is.finite(dy))dy else lc
        } else stop("Unknown mortality definition: ",mortality_definition,call.=FALSE)
        if(!eligible||!is.finite(analysis_end)) return(g[0,,drop=FALSE])
        if(is.finite(dy)) analysis_end <- min(analysis_end,dy)
        g <- g %>% filter(exposure_wave<analysis_end)
        if(!nrow(g)) return(g)
        nxt <- dplyr::lead(g$exposure_wave,default=analysis_end)
        g$tstart <- g$exposure_wave
        g$tstop <- pmin(nxt,analysis_end)
        g$surv_event <- as.integer(is.finite(dy)&g$tstart<dy&g$tstop==dy)
        g %>% filter(tstop>tstart)
      }) %>% ungroup()
    })

    fit_time_updated_score <- function(d,score_var,covs,mortality_definition) {
      dd <- build_time_updated(d,mortality_definition)
      dd$score_z <- zstd(dd[[paste0(score_var,"_cum")]])
      covs_use <- good_covariates(dd,intersect(covs,names(dd)))
      dd <- dd %>% select(all_of(c("IDIND","tstart","tstop","surv_event","score_z",covs_use))) %>% drop_na()
      empty <- function() tibble(N=n_distinct(dd$IDIND),Rows=nrow(dd),Events=sum(dd$surv_event==1),
        Person_years=sum(dd$tstop-dd$tstart),beta=NA_real_,se=NA_real_,HR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,
        Covariates=paste(covs_use,collapse=" + "))
      if(nrow(dd)<100||n_distinct(dd$IDIND)<100||sum(dd$surv_event==1)<10||length(unique(dd$score_z))<2)return(empty())
      fm <- as.formula(paste0("Surv(tstart,tstop,surv_event) ~ score_z",
        if(length(covs_use))paste0(" + ",paste(covs_use,collapse=" + "))else""," + cluster(IDIND)"))
      fit <- tryCatch(coxph(fm,data=dd,ties="efron"),error=function(e)NULL)
      if(is.null(fit))return(empty())
      tt <- broom::tidy(fit)%>%filter(term=="score_z")
      if(!nrow(tt))return(empty())
      se_beta <- if ("robust.se" %in% names(tt)) {
        tt$robust.se[1]
      } else {
        tt$std.error[1]
      }
      tibble(N=n_distinct(dd$IDIND),Rows=nrow(dd),Events=sum(dd$surv_event==1),Person_years=sum(dd$tstop-dd$tstart),
        beta=tt$estimate[1],se=se_beta,HR=exp(tt$estimate[1]),
        CI_low=exp(tt$estimate[1]-1.96*se_beta),CI_high=exp(tt$estimate[1]+1.96*se_beta),
        P=tt$p.value[1],Covariates=paste(covs_use,collapse=" + "))
    }

    fit_time_updated_head2head <- function(d,maha_var,trad_var,trad_lab,covs,mortality_definition="end_wave") {
      dd <- build_time_updated(d,mortality_definition)
      maha_cum <- paste0(maha_var,"_cum"); trad_cum <- paste0(trad_var,"_cum")
      if(!all(c(maha_cum,trad_cum)%in%names(dd))) return(tibble())
      dd$maha_z <- zstd(dd[[maha_cum]]); dd$trad_z <- zstd(dd[[trad_cum]])
      covs_use <- good_covariates(dd,intersect(covs,names(dd)))
      dd <- dd %>% select(all_of(c("IDIND","tstart","tstop","surv_event","maha_z","trad_z",covs_use))) %>% drop_na()
      if(nrow(dd)<100||n_distinct(dd$IDIND)<100||sum(dd$surv_event==1)<10)return(tibble())
      fm <- as.formula(paste0("Surv(tstart,tstop,surv_event) ~ maha_z + trad_z",
        if(length(covs_use))paste0(" + ",paste(covs_use,collapse=" + "))else""," + cluster(IDIND)"))
      fit <- tryCatch(coxph(fm,data=dd,ties="efron"),error=function(e)NULL)
      if(is.null(fit))return(tibble())
      cf <- coef(fit); vc <- vcov(fit)
      if(!all(c("maha_z","trad_z")%in%names(cf)))return(tibble())
      b1 <- unname(cf[["maha_z"]]); b2 <- unname(cf[["trad_z"]]); dff <- b1-b2
      se_diff <- sqrt(vc["maha_z","maha_z"]+vc["trad_z","trad_z"]-2*vc["maha_z","trad_z"])
      tibble(Comparison=paste0("MAHA vs ",trad_lab),Comparator=trad_lab,
        N=n_distinct(dd$IDIND),Events=sum(dd$surv_event==1),HR_MAHA_mutual=exp(b1),HR_comparator_mutual=exp(b2),
        HR_ratio_MAHA_vs_comparator=exp(dff),CI_low_ratio=exp(dff-1.96*se_diff),CI_high_ratio=exp(dff+1.96*se_diff),
        P_equal_coefficients=2*pnorm(-abs(dff/se_diff)),Covariates=paste(covs_use,collapse=" + "))
    }

    time_updated_results <- imap_dfr(panel_cum_by_age,function(panel_cum_i,age_lab) {
      pmap_dfr(mortality_definitions,function(Mortality_definition,Definition_label,Publication_eligible,Definition_note) {
        imap_dfr(adjustment_specs,function(covs,adj_lab) {
          map_dfr(sample_specs,function(sample_lab) {
            dd0 <- panel_cum_i
            if(identical(sample_lab,"four_score_common")) dd0 <- dd0 %>%
              filter(if_all(all_of(paste0(unname(grid_scores),"_cum")),~is.finite(num(.x))))
            imap_dfr(grid_scores,function(v,diet_lab) {
              fit_time_updated_score(dd0,v,covs,Mortality_definition) %>% mutate(
                Model_id=paste0("time_updated_age",age_lab,"_",Mortality_definition,"_",adj_lab,"_",sample_lab),
                Design="Chen-like time-updated cumulative average",Exposure_waves=paste(sort(unique(panel_cum_i$exposure_wave)),collapse=","),
                Baseline_wave=NA_integer_,Minimum_age=as.integer(age_lab),Diet=diet_lab,
                Mortality_definition=Mortality_definition,Definition_label=Definition_label,
                Publication_eligible=Publication_eligible,Definition_note=Definition_note,
                Adjustment=adj_lab,Sample_strategy=sample_lab,Event_target=mortality_event_target,
                Events_ge_target=Events>=mortality_event_target,
                Method_note="Start-stop Cox; cumulative mean updated only with current/prior waves; robust SE clustered by participant",.before=1)
            })
          })
        })
      })
    })

    # Model 1 is the submitted 2011 analysis; Models 4-5 are the only retained
    # high-event cumulative-average sensitivities.
    model0_name <- "Model1_default_reference"
    model_definitions <- bind_rows(
      tibble(Model=model0_name,Design="Fixed 2011 baseline",Adjustment_order=1L,
        Minimum_age=min_age,Exposure_waves=as.character(baseline_wave),
        Covariates=paste(primary_covariates,collapse=" + "),Energy_covariate=TRUE,
        Role="Default / submitted primary model",
        Note=paste0("Exact reproduction of Fig6; reference score algorithms; ",reference_event_n,"-death benchmark")),
      imap_dfr(adjustment_specs,function(covs,lab) tibble(
        Model=lab,Design="Chen-like time-updated cumulative average",Adjustment_order=match(lab,names(adjustment_specs))+1L,
        Minimum_age=panel_age,Exposure_waves=paste(sort(unique(panel_cum$exposure_wave)),collapse=","),
        Covariates=paste(covs,collapse=" + "),Energy_covariate=FALSE,
        Role=unname(final_model_roles[[lab]]),Note=unname(final_model_notes[[lab]])
      ))
    ) %>% mutate(
      Selection_rule="Retained before final reporting for design coverage and interpretability; never selected from HR direction or P value")
    if(nrow(model_definitions)>6L)
      stop("Focused CHNS analysis may retain no more than six named models.",call.=FALSE)

    model0_results <- mort %>% filter(Diet%in%c("MAHA","DASH","MIND","MEDI")) %>%
      transmute(Model=model0_name,Design="Fixed 2011 baseline",
        Score_algorithm=ifelse(Diet=="DASH","DASH-reference",Diet),Diet,N,Events,Person_years,
        beta,se,HR,CI_low,CI_high,P,Covariates,Rows=NA_integer_)

    # Anchor alternatives to one reference-comparable starting cohort.
    chen_core <- c(MAHA="diet.maha.sum",DASH="diet.dash.sum",MIND="diet.mind.sum",MEDI="diet.medi.sum")
    chen_core_cum <- paste0(unname(chen_core),"_cum")
    chen_common <- panel_cum %>% filter(if_all(all_of(chen_core_cum),~is.finite(num(.x))))
    chen_model_results <- imap_dfr(adjustment_specs,function(covs,model_lab) {
      imap_dfr(focused_algorithm_score_set,function(v,score_lab) {
        fit_time_updated_score(chen_common,v,covs,"end_wave") %>% mutate(
          Model=model_lab,Design="Chen-like time-updated cumulative average",Score_algorithm=score_lab,
          Diet=case_when(startsWith(score_lab,"MAHA")~"MAHA",startsWith(score_lab,"DASH")~"DASH",TRUE~score_lab),
          .before=1)
      })
    })
    model_results_0N <- bind_rows(model0_results,chen_model_results)

    model_primary_results <- model_results_0N %>% filter(
      (Model==model0_name & Score_algorithm%in%c("MAHA","DASH-reference","MIND","MEDI")) |
      (Model!=model0_name & Score_algorithm%in%c("MAHA","DASH-reference","MIND","MEDI"))
    )
    maha_algorithm_results <- chen_model_results %>%
      filter(Score_algorithm%in%c("MAHA","MAHA-density","MAHA-residual"))
    last_contact_results <- time_updated_results %>%
      filter(Mortality_definition=="last_contact",Sample_strategy=="four_score_common",
             Adjustment%in%names(adjustment_specs))

    model0_head2head <- pub_head2head %>% transmute(
      Model=model0_name,Comparison,Comparator=Comp2,N=N_total,Events=N_event,
      HR_MAHA_mutual=HR_maha_mutual,HR_comparator_mutual=HR_trad_mutual,
      HR_ratio_MAHA_vs_comparator=HR_ratio_maha_vs_trad,CI_low_ratio=LCI_ratio,CI_high_ratio=UCI_ratio,
      P_equal_coefficients,Covariates=paste(primary_covariates,collapse=" + "))
    chen_head2head <- imap_dfr(adjustment_specs,function(covs,model_lab) {
      bind_rows(
        fit_time_updated_head2head(chen_common,"diet.maha.sum","diet.dash.sum","DASH-reference",covs),
        fit_time_updated_head2head(chen_common,"diet.maha.sum","diet.mind.sum","MIND",covs),
        fit_time_updated_head2head(chen_common,"diet.maha.sum","diet.medi.sum","MEDI",covs)
      ) %>% mutate(Model=model_lab,.before=1)
    })
    model_head2head <- bind_rows(model0_head2head,chen_head2head)

    model_pattern <- model_primary_results %>%
      select(Model,Diet,HR,P,N,Events) %>%
      pivot_wider(names_from=Diet,values_from=c(HR,P,N,Events)) %>%
      mutate(
        Traditional_all_protective=is.finite(HR_DASH)&is.finite(HR_MIND)&is.finite(HR_MEDI)&
          HR_DASH<1&HR_MIND<1&HR_MEDI<1,
        Traditional_all_significant_protective=Traditional_all_protective&
          P_DASH<.05&P_MIND<.05&P_MEDI<.05,
        MAHA_null_or_risky=is.finite(HR_MAHA)&is.finite(P_MAHA)&(HR_MAHA>=1|P_MAHA>=.05),
        Descriptive_target_pattern=Traditional_all_significant_protective&MAHA_null_or_risky
      ) %>% left_join(
        model_head2head %>% filter(Comparator%in%c("DASH","DASH-reference","MIND","MEDI")) %>% group_by(Model) %>% summarise(
          MAHA_formally_diff_from_all=all(is.finite(P_equal_coefficients)&P_equal_coefficients<.05),
          .groups="drop"),by="Model") %>%
      mutate(Selection_rule="Diagnostic only; significance and direction never select a reported model")

    model_0N_labels <- model_primary_results %>% group_by(Model) %>% summarise(
      Events=first(Events),Event_counts_consistent=n_distinct(Events)==1L,.groups="drop") %>%
      right_join(model_definitions%>%select(Model,Adjustment_order),by="Model") %>%
      arrange(Adjustment_order) %>% mutate(Model_label=case_when(
        Model==model0_name~paste0("Primary: fixed 2011 (deaths=",Events,")"),
        Model=="Model4_socioeconomic"~paste0("High-event: socioeconomic (deaths=",Events,")"),
        Model=="Model5_bmi"~paste0("High-event: + BMI (deaths=",Events,")"),
        TRUE~paste0(Model," (deaths=",Events,")")
      ))
    model_0N_label_levels <- rev(model_0N_labels$Model_label)

    p_models_0N <- model_primary_results %>%
      left_join(model_0N_labels%>%select(Model,Model_label),by="Model") %>% mutate(
      Model_label=factor(Model_label,levels=model_0N_label_levels),
      Diet=factor(Diet,levels=c("DASH","MEDI","MIND","MAHA"))) %>%
      ggplot(aes(HR,Model_label,color=Diet))+
      geom_vline(xintercept=1,linetype=2,color="grey55")+
      geom_segment(aes(x=CI_low,xend=CI_high,yend=Model_label),position=position_dodge(width=.68),linewidth=.6)+
      geom_point(aes(size=Events),position=position_dodge(width=.68))+
      scale_color_manual(values=c(MAHA="#F26D60",DASH="#16B9C0",MIND="#B97AF7",MEDI="#7CAE00"))+
      labs(title="CHNS mortality: selected primary and high-event models",
        subtitle="Fixed 2011 primary analysis and age >=40 cumulative-average sensitivities",
        x="HR per 1 SD healthier score",y=NULL,color=NULL,size="Deaths")+
      theme_classic(10)+theme(legend.position="bottom")
    save_plot(p_models_0N,"mortality_models_final.png",12,6.5,420)

    all_sweep_results <- time_updated_results %>%
      filter(Adjustment%in%names(adjustment_specs),Sample_strategy=="four_score_common")
    if(n_distinct(all_sweep_results$Model_id)>6L)
      stop("Focused CHNS high-event output may contain no more than six model combinations.",call.=FALSE)
    sweep_pattern <- all_sweep_results %>%
      select(Model_id,Design,Exposure_waves,Baseline_wave,Minimum_age,Mortality_definition,Definition_label,
             Publication_eligible,Adjustment,Sample_strategy,Diet,HR,P,N,Events,Person_years) %>%
      pivot_wider(names_from=Diet,values_from=c(HR,P,N,Events,Person_years)) %>%
      mutate(
        Established_all_protective=is.finite(HR_DASH)&is.finite(HR_MIND)&is.finite(HR_MEDI)&HR_DASH<1&HR_MIND<1&HR_MEDI<1,
        Established_all_significant_protective=Established_all_protective&P_DASH<.05&P_MIND<.05&P_MEDI<.05,
        MAHA_attenuated_or_null_vs_all=is.finite(HR_MAHA)&is.finite(HR_DASH)&is.finite(HR_MIND)&is.finite(HR_MEDI)&
          HR_MAHA>=pmax(HR_DASH,HR_MIND,HR_MEDI),
        MAHA_null_or_risky=is.finite(HR_MAHA)&is.finite(P_MAHA)&(HR_MAHA>=1|P_MAHA>=.05),
        Descriptive_expected_pattern=Established_all_significant_protective&MAHA_null_or_risky&MAHA_attenuated_or_null_vs_all,
        Events_ge_target=Events_MAHA>=mortality_event_target,
        Selection_rule="Diagnostic only: direction is not a primary-model selection criterion"
      )

    sweep_sources <- tibble(
      Design="Chen-like time-updated cumulative average",
      Reference="Chen et al., Nutrition 2021;81:110902",
      Reported_sample="14,117 participants; 1,007 deaths",
      URL="https://doi.org/10.1016/j.nut.2020.110902",
      Replication_scope="Time-updated cumulative exposure structure with age >=40 and the present score definitions"
    )
    sweep_configuration <- tibble(
      Setting=c("Exposure waves","Minimum age","Retained models","Retained named-model count",
                "High-event model-combination count","Event-count flag","Selection policy","Schooling definition"),
      Value=c(paste(eligible_sweep_waves,collapse=","),as.character(panel_age),
              "Model4_socioeconomic,Model5_bmi",as.character(nrow(model_definitions)),
              as.character(n_distinct(all_sweep_results$Model_id)),as.character(mortality_event_target),
              "No model is selected from coefficient direction or P value",analysis_metadata$schooling_definition)
    )
    final_workbook <- list(
      model_definitions=model_definitions,
      primary_results=model_primary_results,
      last_contact_results=last_contact_results,
      maha_algorithm_results=maha_algorithm_results,
      head2head=model_head2head,
      pattern_diagnostic=model_pattern,
      configuration=sweep_configuration,
      run_manifest=run_manifest
    )
    write_xlsx(final_workbook,"mortality_models_final.out.xlsx")
    write_xlsx(list(configuration=sweep_configuration,sources=sweep_sources,
                    end_and_last_contact=all_sweep_results,pattern_diagnostic=sweep_pattern,
                    run_manifest=run_manifest),
               "mortality_high_event_sensitivity.out.xlsx")

    p_sweep <- all_sweep_results %>%
      filter(is.finite(HR)) %>%
      mutate(Model_label=case_when(
        Adjustment=="Model4_socioeconomic"~"Model 4: socioeconomic",
        Adjustment=="Model5_bmi"~"Model 5: + BMI",
        TRUE~Adjustment
      ),Diet=factor(Diet,levels=c("DASH","MEDI","MIND","MAHA"))) %>%
      ggplot(aes(HR,Model_label,color=Diet))+
      geom_vline(xintercept=1,linetype=2,color="grey55")+
      geom_segment(aes(x=CI_low,xend=CI_high,yend=Model_label),
                   position=position_dodge(width=.65),linewidth=.55)+
      geom_point(aes(size=Events),position=position_dodge(width=.65))+
      facet_wrap(~Definition_label,nrow=1)+
      scale_color_manual(values=c(MAHA="#F26D60",DASH="#16B9C0",MIND="#B97AF7",MEDI="#7CAE00"))+
      labs(title="CHNS high-event cumulative-average sensitivity",
           subtitle=paste0("Age >=",panel_age,"; exposure waves ",paste(eligible_sweep_waves,collapse=","),"; common-score samples"),
           x="HR per 1 SD healthier score",y=NULL,color=NULL,size="Deaths")+
      theme_classic(10)+theme(legend.position="bottom")
    save_plot(p_sweep,"mortality_high_event_sensitivity.png",12,6.5,420)
  }
}

cat("\n[CHNS] Completed.\n")
cat("[CHNS] Primary fixed-wave analysis uses baseline ",baseline_wave," and Cox regression.\n",sep="")
cat("[CHNS] Publication inputs: ",cohort_file("publication_inputs.out.xlsx"),"\n",sep="")
