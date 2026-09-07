pacman::p_load(tidyverse, data.table, survival)

dir0 <- Sys.getenv("DIR0", unset = ifelse(Sys.info()[["sysname"]] == "Windows", "D:/", "/mnt/d"))
indir <- Sys.getenv("PHEDIR", unset = paste0(dir0, "/data/ukb/phe"))
outdir <- Sys.getenv("UKB_OUT", unset = paste0(dir0, "/analysis/ukb/phe"))
script_dir <- Sys.getenv("SCRIPT_DIR", unset = file.path(dir0, "scripts"))
gen_dir <- Sys.getenv("GEN_DIR", unset = file.path(dir0, "data/ukb/gen"))
map_tile_cache <- normalizePath(
	file.path(dir0, "data", "ukb", "map", "maha", "rosm.cache"),
	winslash = "/", mustWork = TRUE
)
invisible(lapply(c("0phe.f.R", "mr.f.R", "assoc.f.R", "plot.f.R", "pred.f.R"), \(f) source(file.path(script_dir, "0f", f))))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd2(outdir)
analysis_file <- function(...) file.path(outdir, ...)

step_order <- c(
	"data_qc", "vip_phe", "ses", "move", "death", "srd", "icd", "audit_i67", "fod_ref", "gp_prep", "gp",
	"biom", "hla_pca", "pgs", "merge", "audit_dates",
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
get_first_existing_col <- function(dat, trait, prefixes = c("", "fod_", "date_")) {
	hit <- get_existing_date_cols(dat, trait, prefixes)
	if (length(hit)) hit[1] else NA_character_
}
get_date_by_eid <- function(dat, trait, eids, prefixes = c("", "fod_", "date_")) {
	col <- get_first_existing_col(dat, trait, prefixes)
	if (is.na(col)) return(rep(as.Date(NA), length(eids)))
	asDate2(dat[[col]])[match(as.character(eids), as.character(dat$eid))]
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


# CKM component helpers for phe.R
# Baseline UACR, baseline/pre-baseline GP-derived CAC/ABI/pre-HF components,
# and imaging-visit CMR/carotid auxiliary phenotypes.

ckm_num <- function(x) {
  z <- suppressWarnings(as.numeric(as.character(x)))
  z[!is.finite(z)] <- NA_real_
  z
}

ckm_date <- function(x) {
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

  alb <- ckm_num(get_col(alb_col))
  cr <- ckm_num(get_col(cr_col))
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

read_gp_ckm_spec <- function(file, code_system) {
  if (!file.exists(file)) return(data.table::data.table())
  z <- data.table::fread(file, sep = "\t", header = TRUE, fill = TRUE,
                         colClasses = "character", showProgress = FALSE)
  req <- c("codes_csv", "component", "kind")
  miss <- setdiff(req, names(z))
  if (length(miss)) stop("Missing GP CKM list columns in ", file, ": ", paste(miss, collapse = ", "))
  for (nm in c("plausible_min", "plausible_max", "positive_threshold")) {
    if (!nm %in% names(z)) z[[nm]] <- NA_character_
  }
  if (!"positive_op" %in% names(z)) z[, positive_op := NA_character_]
  if (!"primary_use" %in% names(z)) z[, primary_use := NA_character_]
  z[, `:=`(
    codes_csv = trimws(as.character(codes_csv)),
    component = trimws(as.character(component)),
    kind = trimws(as.character(kind)),
    plausible_min = suppressWarnings(as.numeric(plausible_min)),
    plausible_max = suppressWarnings(as.numeric(plausible_max)),
    positive_threshold = suppressWarnings(as.numeric(positive_threshold)),
    positive_op = trimws(as.character(positive_op)),
    primary_use = trimws(as.character(primary_use)),
    code_system = code_system
  )]
  z <- z[!is.na(codes_csv) & nzchar(codes_csv) & !is.na(component) & nzchar(component)]
  out <- z[, .(
    code = trimws(unlist(strsplit(codes_csv, ",", fixed = TRUE)))
  ), by = .(component, kind, plausible_min, plausible_max,
            positive_op, positive_threshold, primary_use, code_system)]
  unique(out[!is.na(code) & nzchar(code)])
}

# Choose the first value1/value2/value3 entry compatible with the component's
# plausible range. The GP release does not provide a universal unit column;
# component-specific distributions must still be audited by data_provider.
first_plausible_gp_value <- function(v1, v2, v3, lo, hi) {
  v <- c(ckm_num(v1), ckm_num(v2), ckm_num(v3))
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  if (is.finite(lo)) v <- v[v >= lo]
  if (is.finite(hi)) v <- v[v <= hi]
  if (!length(v)) NA_real_ else v[1]
}

# Build baseline GP components from the active-list specification.
make_gp_ckm_from_fst <- function(fst_file, phe_file, gp2_list, gp3_list) {
  if (!requireNamespace("fst", quietly = TRUE)) stop("Package 'fst' is required")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package 'data.table' is required")
  DT <- data.table::as.data.table

  base <- readRDS(phe_file) |> as.data.frame()
  if (!all(c("eid", "date_attend") %in% names(base))) {
    stop("phe.rds must contain eid and date_attend")
  }
  base <- DT(base[, c("eid", "date_attend")])
  base[, `:=`(eid = as.character(eid), date_attend = ckm_date(date_attend))]
  base <- unique(base, by = "eid")

  need <- c("eid", "data_provider", "event_dt", "read_2", "read_3",
            "value1", "value2", "value3")
  dt <- DT(fst::read_fst(fst_file, columns = need))
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("gp_clinical.full.fst is missing: ", paste(miss, collapse = ", "))
  dt[, `:=`(
    eid = as.character(eid),
    data_provider = as.integer(data_provider),
    event_dt = ckm_date(event_dt),
    read_2 = trimws(as.character(read_2)),
    read_3 = trimws(as.character(read_3))
  )]
  dt <- merge(dt, base, by = "eid", all.x = FALSE, all.y = FALSE, sort = FALSE)
  bad_dates <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02",
                         "1903-03-03", "1909-09-09", "2037-07-07"))
  dt <- dt[!is.na(event_dt) & !is.na(date_attend) &
             event_dt <= date_attend & !event_dt %in% bad_dates]

  s2 <- read_gp_ckm_spec(gp2_list, "read2")
  s3 <- read_gp_ckm_spec(gp3_list, "ctv3")
  hit <- list()
  if (nrow(s2)) {
    d2 <- dt[!is.na(read_2) & nzchar(read_2) & read_2 %chin% s2$code]
    if (nrow(d2)) {
      h2 <- merge(d2, s2, by.x = "read_2", by.y = "code",
                  allow.cartesian = TRUE, sort = FALSE)
      h2[, code := read_2]
      hit[["read2"]] <- h2
    }
  }
  if (nrow(s3)) {
    d3 <- dt[!is.na(read_3) & nzchar(read_3) & read_3 %chin% s3$code]
    if (nrow(d3)) {
      h3 <- merge(d3, s3, by.x = "read_3", by.y = "code",
                  allow.cartesian = TRUE, sort = FALSE)
      h3[, code := read_3]
      hit[["ctv3"]] <- h3
    }
  }
  h <- data.table::rbindlist(hit, use.names = TRUE, fill = TRUE)

  all_out <- data.table::data.table(eid = base$eid)
  if (!nrow(h)) {
    warning("No GP CKM codes matched active gp2.ckm.lst/gp3.ckm.lst")
    all_out[, `:=`(
      subclinical_ascvd = NA_integer_,
      subclinical_ascvd_observed = FALSE,
      date_subclinical_ascvd = as.Date(NA),
      date_subclinical_ascvd_assessed = as.Date(NA),
      subclinical_ascvd_source = NA_character_,
      pre_hf = NA_integer_,
      pre_hf_observed = FALSE,
      date_pre_hf = as.Date(NA),
      date_pre_hf_assessed = as.Date(NA),
      pre_hf_source = NA_character_
    )]
    return(as.data.frame(all_out))
  }

  h[, value := mapply(
    first_plausible_gp_value,
    value1, value2, value3, plausible_min, plausible_max,
    SIMPLIFY = TRUE, USE.NAMES = FALSE
  )]

  # Latest valid pre-baseline numeric result per component.
  num <- h[kind == "numeric" & is.finite(value)]
  if (nrow(num)) {
    data.table::setorder(num, eid, component, event_dt)
    num_latest <- num[, .SD[.N], by = .(eid, component)]
    num_val <- data.table::dcast(num_latest, eid ~ component, value.var = "value")
    num_date <- data.table::dcast(num_latest, eid ~ component, value.var = "event_dt")
    data.table::setnames(num_val, setdiff(names(num_val), "eid"),
                         paste0("gp_", setdiff(names(num_val), "eid")))
    data.table::setnames(num_date, setdiff(names(num_date), "eid"),
                         paste0("date_gp_", setdiff(names(num_date), "eid")))
    all_out <- merge(all_out, num_val, by = "eid", all.x = TRUE, sort = FALSE)
    all_out <- merge(all_out, num_date, by = "eid", all.x = TRUE, sort = FALSE)
  }

  # Positive structural/abnormality codes.
  pos <- h[kind == "positive_code", .(
    positive = TRUE,
    first_date = min(event_dt),
    last_date = max(event_dt)
  ), by = .(eid, component)]
  if (nrow(pos)) {
    pos_w <- data.table::dcast(pos, eid ~ component, value.var = "positive", fill = FALSE)
    pos_d <- data.table::dcast(pos, eid ~ component, value.var = "first_date")
    data.table::setnames(pos_w, setdiff(names(pos_w), "eid"),
                         paste0("gp_code_", setdiff(names(pos_w), "eid")))
    data.table::setnames(pos_d, setdiff(names(pos_d), "eid"),
                         paste0("date_gp_code_", setdiff(names(pos_d), "eid")))
    all_out <- merge(all_out, pos_w, by = "eid", all.x = TRUE, sort = FALSE)
    all_out <- merge(all_out, pos_d, by = "eid", all.x = TRUE, sort = FALSE)
  }

  assessed <- unique(h[, .(eid, component)])
  if (nrow(assessed)) {
    assessed[, assessed := TRUE]
    ass_w <- data.table::dcast(assessed, eid ~ component,
                               value.var = "assessed", fill = FALSE)
    data.table::setnames(ass_w, setdiff(names(ass_w), "eid"),
                         paste0("gp_assessed_", setdiff(names(ass_w), "eid")))
    all_out <- merge(all_out, ass_w, by = "eid", all.x = TRUE, sort = FALSE)
  }

  get_num <- function(nm) {
    if (nm %in% names(all_out)) ckm_num(all_out[[nm]]) else rep(NA_real_, nrow(all_out))
  }
  get_log <- function(nm) {
    if (nm %in% names(all_out)) all_out[[nm]] %in% TRUE else rep(FALSE, nrow(all_out))
  }
  get_date <- function(nm) {
    if (nm %in% names(all_out)) ckm_date(all_out[[nm]]) else as.Date(rep(NA, nrow(all_out)))
  }
  keep_date <- function(x, keep) {
    z <- ckm_date(x)
    z[!(keep %in% TRUE)] <- as.Date(NA)
    z
  }
  min_date_vec <- function(...) {
    z <- list(...)
    z <- lapply(z, ckm_date)
    if (!length(z)) return(as.Date(rep(NA, nrow(all_out))))
    ans <- do.call(pmin, c(lapply(z, as.numeric), list(na.rm = TRUE)))
    ans[!is.finite(ans)] <- NA_real_
    as.Date(ans, origin = "1970-01-01")
  }

  abi <- get_num("gp_abi")
  cac <- get_num("gp_cac")
  lvef <- get_num("gp_lvef")
  date_abi <- get_date("date_gp_abi")
  date_cac <- get_date("date_gp_cac")
  date_lvef <- get_date("date_gp_lvef")

  abi_pos <- is.finite(abi) & abi <= 0.90
  cac_pos <- is.finite(cac) & cac >= 100
  sub_pos <- abi_pos | cac_pos
  sub_obs <- is.finite(abi) | is.finite(cac)
  date_sub_assessed <- min_date_vec(date_abi, date_cac)
  date_sub_positive <- min_date_vec(
    keep_date(date_abi, abi_pos),
    keep_date(date_cac, cac_pos)
  )
  all_out[, `:=`(
    subclinical_ascvd = ifelse(sub_pos, 1L, ifelse(sub_obs, 0L, NA_integer_)),
    subclinical_ascvd_observed = sub_obs,
    date_subclinical_ascvd = date_sub_positive,
    date_subclinical_ascvd_assessed = date_sub_assessed,
    subclinical_ascvd_source = fifelse(
      cac_pos & abi_pos, "GP CAC+ABI",
      fifelse(cac_pos, "GP CAC",
              fifelse(abi_pos, "GP ABI",
                      fifelse(sub_obs, "GP measured negative", NA_character_)))
    )
  )]

  structural_components <- c(
    "lv_systolic_dysfunction", "lv_hypertrophy",
    "chamber_enlargement", "valvular_disease"
  )
  structural_cols <- paste0("gp_code_", structural_components)
  structural_pos <- Reduce(`|`, lapply(structural_cols, get_log),
                           init = rep(FALSE, nrow(all_out)))
  structural_dates <- lapply(
    paste0("date_gp_code_", structural_components),
    get_date
  )
  structural_date <- do.call(min_date_vec, structural_dates)
  lvef_pos <- is.finite(lvef) & lvef < 40
  pre_pos <- lvef_pos | structural_pos
  # A normal LVEF is a valid negative for the reduced-EF route. Absence of a
  # structural diagnosis code is not interpreted as a negative imaging test.
  pre_obs <- is.finite(lvef) | structural_pos
  date_pre_positive <- min_date_vec(
    keep_date(date_lvef, lvef_pos),
    structural_date
  )
  date_pre_assessed <- min_date_vec(date_lvef, structural_date)
  all_out[, `:=`(
    pre_hf = ifelse(pre_pos, 1L, ifelse(is.finite(lvef), 0L, NA_integer_)),
    pre_hf_observed = pre_obs,
    date_pre_hf = date_pre_positive,
    date_pre_hf_assessed = date_pre_assessed,
    pre_hf_source = fifelse(
      lvef_pos & structural_pos, "GP LVEF+structural code",
      fifelse(lvef_pos, "GP LVEF<40%",
              fifelse(structural_pos, "GP structural code",
                      fifelse(is.finite(lvef), "GP LVEF measured negative", NA_character_)))
    ),
    pre_hf_biomarker_available = is.finite(get_num("gp_bnp")) |
      is.finite(get_num("gp_ntprobnp")) | is.finite(get_num("gp_troponin")) |
      get_log("gp_code_elevated_troponin")
  )]

  # Audit GP values by provider before using BNP/NT-proBNP/troponin thresholds.
  attr(all_out, "gp_ckm_long") <- h[, .(
    eid, data_provider, event_dt, code_system, code, component, kind,
    value1, value2, value3, value, primary_use
  )]
  as.data.frame(all_out)
}


make_img_ckm <- function(img_file, img_list) {
  if (!file.exists(img_file)) stop("Cannot find imaging table: ", img_file)
  if (!file.exists(img_list)) stop("Cannot find img.ckm.lst: ", img_list)
  # Follow the existing UKB extraction-list convention:
  # field, renamed variable, i0/all. Imaging rows use "all" so phe.sh cuts all
  # available instances; CKM derivation then explicitly retains visits i2/i3.
  mp0 <- data.table::fread(
    img_list, sep = "\t", header = FALSE, select = 1:3, fill = TRUE,
    colClasses = "character", showProgress = FALSE
  )
  data.table::setnames(mp0, c("field", "variable", "scope"))
  mp0 <- mp0[grepl("^p[0-9]+$", field) & nzchar(variable)]
  if (!nrow(mp0)) stop("No valid rows in img.ckm.lst")
  mp <- data.table::rbindlist(lapply(c(2L, 3L), function(inst) {
    mp0[, .(
      field,
      variable = paste0(variable, "_i", inst),
      instance = as.character(inst)
    )]
  }))
  mp[, raw := paste0(field, "_i", instance)]
  if (anyDuplicated(mp$raw)) stop("Duplicated raw imaging fields in img.ckm.lst")
  if (anyDuplicated(mp$variable)) stop("Duplicated variable names in img.ckm.lst")

  header_names <- names(data.table::fread(img_file, nrows = 0, showProgress = FALSE))
  selected <- intersect(mp$raw, header_names)
  if (!length(selected)) stop("No img.ckm.lst fields were found in imaging table")
  x <- data.table::fread(
    img_file, sep = "\t", header = TRUE,
    select = c("eid", selected), showProgress = TRUE
  )
  data.table::setnames(x, selected, mp$variable[match(selected, mp$raw)])
  x[, eid := as.character(eid)]

  date_cols <- grep("^date_img", names(x), value = TRUE)
  for (cc in date_cols) data.table::set(x, j = cc, value = ckm_date(x[[cc]]))
  num_cols <- setdiff(names(x), c("eid", date_cols))
  for (cc in num_cols) data.table::set(x, j = cc, value = ckm_num(x[[cc]]))

  getv <- function(nm) {
    if (nm %in% names(x)) ckm_num(x[[nm]]) else rep(NA_real_, nrow(x))
  }
  getd <- function(nm) {
    if (nm %in% names(x)) ckm_date(x[[nm]]) else as.Date(rep(NA, nrow(x)))
  }

  derive_one_instance <- function(inst) {
    suf <- paste0("_i", inst)
    bsa <- getv(paste0("cmr_bsa_m2", suf))
    index_map <- c(
      cmr_lvedv_ml = "cmr_lvedv_index",
      cmr_lvesv_ml = "cmr_lvesv_index",
      cmr_lv_mass_g = "cmr_lv_mass_index",
      cmr_rvedv_ml = "cmr_rvedv_index",
      cmr_rvesv_ml = "cmr_rvesv_index",
      cmr_la_max_ml = "cmr_la_max_index",
      cmr_ra_max_ml = "cmr_ra_max_index"
    )
    for (src in names(index_map)) {
      src_nm <- paste0(src, suf)
      if (src_nm %in% names(x)) {
        out_nm <- paste0(index_map[[src]], suf)
        x[[out_nm]] <<- ifelse(
          is.finite(bsa) & bsa > 0,
          getv(src_nm) / bsa,
          NA_real_
        )
      }
    }

    lvef <- getv(paste0("cmr_lvef_pct", suf))
    x[[paste0("pre_hf_img_hfref", suf)]] <<-
      ifelse(is.finite(lvef), as.integer(lvef < 40), NA_integer_)
    x[[paste0("pre_hf_img_broad_sensitivity", suf)]] <<-
      ifelse(is.finite(lvef), as.integer(lvef < 50), NA_integer_)

    imt_cols <- intersect(
      paste0(c(
        "carotid_imt_mean_120_um", "carotid_imt_mean_150_um",
        "carotid_imt_mean_210_um", "carotid_imt_mean_240_um"
      ), suf),
      names(x)
    )
    imt_name <- paste0("carotid_imt_mean_um", suf)
    if (length(imt_cols)) {
      # UKB carotid IMT uses 0 for an invalid measurement.
      imt_dat <- as.data.frame(x[, ..imt_cols])
      imt_dat[] <- lapply(imt_dat, function(z) {
        z <- ckm_num(z)
        z[z <= 0] <- NA_real_
        z
      })
      x[[imt_name]] <<- rowMeans(imt_dat, na.rm = TRUE)
      x[[imt_name]][!is.finite(x[[imt_name]])] <<- NA_real_
    } else {
      x[[imt_name]] <<- NA_real_
    }
    x[[paste0("carotid_imt_assessed", suf)]] <<- is.finite(x[[imt_name]])

    date_img <- getd(paste0("date_img", suf))
    hfref <- x[[paste0("pre_hf_img_hfref", suf)]]
    broad <- x[[paste0("pre_hf_img_broad_sensitivity", suf)]]
    x[[paste0("date_pre_hf_img_hfref", suf)]] <<- date_img
    x[[paste0("date_pre_hf_img_hfref", suf)]][!(hfref == 1L)] <<- as.Date(NA)
    x[[paste0("date_pre_hf_img_broad", suf)]] <<- date_img
    x[[paste0("date_pre_hf_img_broad", suf)]][!(broad == 1L)] <<- as.Date(NA)
  }

  for (inst in intersect(c("2", "3"), unique(mp$instance))) derive_one_instance(inst)

  # Generic aliases refer to the first imaging visit (Instance 2).
  alias_map <- c(
    date_img = "date_img_i2",
    pre_hf_img_hfref = "pre_hf_img_hfref_i2",
    pre_hf_img_broad_sensitivity = "pre_hf_img_broad_sensitivity_i2",
    date_pre_hf_img_hfref = "date_pre_hf_img_hfref_i2",
    date_pre_hf_img_broad = "date_pre_hf_img_broad_i2",
    carotid_imt_mean_um = "carotid_imt_mean_um_i2",
    carotid_imt_assessed = "carotid_imt_assessed_i2"
  )
  for (to in names(alias_map)) {
    from <- alias_map[[to]]
    if (from %in% names(x)) x[[to]] <- x[[from]]
  }

  # Repeat-imaging transition endpoints. These are future/landmark outcomes,
  # never baseline-stage definitions.
  if (all(c("pre_hf_img_hfref_i2", "pre_hf_img_hfref_i3") %in% names(x))) {
    x[, pre_hf_img_incident_i3 := as.integer(
      pre_hf_img_hfref_i2 == 0L & pre_hf_img_hfref_i3 == 1L
    )]
    x[is.na(pre_hf_img_hfref_i2) | is.na(pre_hf_img_hfref_i3),
      pre_hf_img_incident_i3 := NA_integer_]
    x[, date_pre_hf_img_incident_i3 := ckm_date(date_img_i3)]
    x[pre_hf_img_incident_i3 != 1L | is.na(pre_hf_img_incident_i3),
      date_pre_hf_img_incident_i3 := as.Date(NA)]
  }

  # cIMT is retained as an atherosclerosis proxy. It is not renamed to strict
  # subclinical ASCVD because the AHA Stage-3 criterion is CAC/ABI based.
  x[, subclinical_atherosclerosis_proxy_img := NA_integer_]
  as.data.frame(x)
}


run_if("data_qc", {
# 🚩 data QC
fn	<- fread(paste0(indir, "/common/vip.lst"), sep = "\t", select = 1:3, fill = Inf, header = FALSE) %>% as.data.frame()
fn$V1[duplicated(fn$V1)]; fn$V2[duplicated(fn$V2)]
phe0 <- fread(paste0(indir, "/rap/vip.tab.gz"), sep = "\t", header = TRUE) %>% as.data.frame()
names(phe0)	<- names(phe0) %>% gsub("_i0$|_i0_a0$", "", .)
setdiff(fn$V1, names(phe0)); setdiff(names(phe0), fn$V1)
names(phe0)	<- sapply(names(phe0), replace_code, sep = "_", mapping_df = fn) %>% as.character()

phe0 <- phe0 %>% mutate(across(grep("date", names(phe0), value = TRUE), as.Date), edu = ifelse(edu == -7, 7, edu))
#special negatvie values
phe0 <- neg2NA(phe0, c(-1,-2,-3,-7,-10,-121,-818)); table(phe0$social) # 
#special dates
date_4na <- as.Date(c("1900-01-01","1901-01-01","1902-02-02","1903-03-03","1909-09-09","2037-07-07"))
date_cols <- grep("date", names(phe0), ignore.case = TRUE, value = TRUE)
purrr::walk(date_cols, \(cc) print(summary(phe0[[cc]])))
phe0 <- phe0 %>% mutate(across(all_of(date_cols), \(x) { x[x %in% date_4na] <- NA; x }))
#combo variables
phe0 %>% dplyr::select(where(~ any(stringr::str_detect(as.character(.x), stringr::fixed("|")), na.rm = TRUE))) %>% names()
saveRDS(phe0, paste0(indir, "/Rdata/phe0.rds"))
})


run_if("vip_phe", {
# 🚩 VIP Phe
if (!exists("phe0")) phe0 <- read_rds("phe0")
phe	<- phe0 %>% mutate( # age_reception [p21003], not age_recruit [p21022]
	center = paste0("c", center),
	abo = recode_factor(blood_type, "AA"="A", "BB"="B", "OO"="O", "AO"="A", "BO"="B"),
	ethnic.c = ifelse(ethnicity==2002, 12, ifelse(ethnicity==2003, 13, ifelse(ethnicity %in% c(2004, 3004), 6, ifelse(grepl("^1", ethnicity), 1, ifelse(grepl("^2", ethnicity), 2, ifelse(grepl("^3", ethnicity), 3, ifelse(grepl("^4", ethnicity), 4, ethnicity))))))),
	ethnic.c = factor(ethnic.c, levels = c(1:6, 12:13), labels = c("White", "Mixed", "Asian", "Black", "Chinese", "Other", "White-Black", "White-Asian")),
	tdi.i = scale(inormal(-tdi)), 
	age.c = cut(age, breaks = seq(35, 75, 5)), age2 = age * age,
	birth_date = as.Date(paste0(birth_year, "-", birth_month, "-15")),
	sex.b = ifelse(sex==0, "female", ifelse(sex==1, "male", NA)),
	sex.f = factor(sex, c(0, 1), c("Female", "Male")),
	sex_first.c = factor(ifelse(sex_first<12, NA, ifelse(sex_first<16, "early", ifelse(sex_first<22, "average", "late"))), levels = c("early", "average", "late")),
	sex_par.c = factor(ifelse(sex_par<=1, "low", ifelse(sex_par<=10, "average", "high")), levels = c("low", "average", "high")),
	whr = waist / hip,  
	leg = height - height_sit, lhr = leg / height, lhr = residuals(lm(lhr ~ height, na.action = na.exclude)),
	bald.c = factor(bald, levels = 1:4, labels = c("ctrl", "bald2", "bald3", "bald4")),
	bald12 = ifelse(bald==1, 0, ifelse(bald==2, 1, NA)), bald13 = ifelse(bald==1, 0, ifelse(bald==3, 1, NA)), bald14 = ifelse(bald==1, 0, ifelse(bald==4, 1, NA)),
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
	) %>% dplyr::select(-center0) %>%
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
p3 <- plot_map(dat1, centers = c("Croydon", "Hounslow", "Bristol", "Edinburgh"), separate = TRUE, width = 50, nrow = 2, transparency = "tdi",  center.lst = center.lst, basemap = TRUE, basemap_type = "osm", tile_zoom = 11, tile_cachedir = map_tile_cache)
	ggsave("Fig3_4centers_tdi.png", p3, width = 10, height = 10, dpi = 300, bg = "white")

dat1 <- dat1 %>% dplyr::select(eid, lon.birth, lat.birth, lon.home, lat.home, move, tdi.birth, tdi.birth.nn)
saveRDS(dat1, paste0(indir, "/Rdata/move.rds"))
})


run_if("death", {
# 🚩 Death
fn <- fread(paste0(indir,"/common/death.lst"), sep = "\t", header = FALSE, select = 1:2) %>% as.data.frame()
dat0 <- fread(paste0(indir,"/rap/death.tab.gz"), sep = "\t", header = TRUE) %>% as.data.frame()
concat <- function(x) {
	x <- as.character(x); x <- x[!is.na(x) & x != ""]
	if (!length(x)) NA_character_ else paste(x, collapse = " ") 
}

dat <- dat0 %>% mutate(
	p40000 = as.Date(coalesce(p40000_i0, p40000_i1)),
	p40001 = apply(as.data.frame(dplyr::select(., matches("^p40001_i"))), 1, concat),
	p40002 = apply(as.data.frame(dplyr::select(., matches("^p40002_i[01]_a"))), 1, concat),
	p40007 = coalesce(p40007_i0, p40007_i1),
	p40010 = as.character(p40010_i0)
) %>% dplyr::select(eid, any_of(fn$V1)) %>% rename_with(~ fn$V2[match(.x, fn$V1)], any_of(fn$V1)) %>% dplyr::select(-death_cause) %>% mutate(
	date_death = if_else(date_death > date_follow_end, date_follow_end, date_death),
	date_cvd_death = if_else(!is.na(date_death) & grepl("(^|\\s)I", death_icd10), date_death, as.Date(NA))
)
saveRDS(dat, paste0(indir, "/Rdata/death.rds"))
})


run_if("srd", {
# 🚩 SRD (self-report disease) (20002 & 87)
date0 <- fread(paste0(indir, "/rap/srd.time.tab.gz"), header = TRUE) %>% as.data.frame() %>% mutate(across(where(is.logical), as.integer))
birth <- readRDS(paste0(indir, "/Rdata/phe.rds")) %>% dplyr::select(eid, birth_year)
tmp <- merge_check(date0 = date0, birth = birth, by = "eid")
date0 <- merge(tmp$date0, tmp$birth, by = "eid", all.x = TRUE, sort = FALSE)

year_cap <- 2023L
bad_tbl <- date0 %>% dplyr::select(-birth_year) %>% pivot_longer(-eid, names_to = "var", values_to = "x") %>%
	left_join(birth, by = "eid") %>% mutate(year_from_age = birth_year + x) %>% filter(!is.na(x), x >= 0, x <= 100, year_from_age > year_cap) %>% arrange(year_from_age)
bad_tbl; bad_map <- bad_tbl %>% dplyr::select(eid, var) %>% dplyr::mutate(flag = TRUE)

date0 <- date0 %>% pivot_longer(-c(eid, birth_year), names_to = "var", values_to = "x") %>%
	left_join(bad_map, by = c("eid", "var")) %>% mutate(x = if_else(!is.na(flag) & flag, NA_real_, as.numeric(x))) %>%
	dplyr::select(-flag) %>% pivot_wider(names_from = "var", values_from = "x")

date0 <- date0 %>% mutate(
	across(-c(eid, birth_year), function(x) {
		res_year <- case_when(!is.na(x) & x >= 1930 & x <= year_cap ~ as.integer(x), !is.na(x) & x >= 0 & x <= 100 ~ as.integer(birth_year + x), TRUE ~ NA_integer_)
		if_else(is.na(res_year), NA_character_, paste0(res_year, "-07-01"))
	})
) %>% dplyr::select(-birth_year)
fwrite(date0, paste0(indir, "/rap/srd.date.tab.gz"), sep = "\t", na = "NA", quote = FALSE)
max(as.Date(unlist(date0[,-1])), na.rm = TRUE)
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
if (!exists("get_icd_fast", mode = "function")) stop("get_icd_fast() is missing. Update scripts/0f/0phe.f.R first.", call. = FALSE)

for (src in icd_sources) {
	date0 <- fread(paste0(indir, "/rap/", src, ".date.tab.gz"), header = TRUE) %>% as.data.frame()
	cols <- ncol(date0); colnames(date0) <- c("eid", paste0("d", seq_len(cols - 1)))
	date0$eid <- as.character(date0$eid)
	print(max(as.matrix(date0[ , -1]), na.rm = TRUE)) # 
	saveRDS(date0, paste0(indir, "/Rdata/", src, ".date0.rds"))
	code0 <- readLines(paste0(indir, "/rap/", src, ".code.tab.gz"))[-1] %>% strsplit(., split = "\\||\t")
	code0 <- do.call(rbind, lapply(code0, function(row) { length(row) <- cols; return(row) })) %>% as.data.frame()
	colnames(code0) <- c("eid", paste0("i", seq_len(cols - 1)))
	code0$eid <- as.character(code0$eid)
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
		warning("MAKE_PHEWAS=TRUE but icd10.code0.rds/icd10.date0.rds is missing; skipped PheWAS.rds.", call. = FALSE)
	} else if (!requireNamespace("PheWAS", quietly = TRUE)) {
		warning(
			"MAKE_PHEWAS=TRUE but package PheWAS is not installed; skipped PheWAS.rds. ",
			"This does not affect fod/lod/cnt/aodtxt or all.rds.",
			call. = FALSE
		)
	} else {
		code0 <- readRDS(icd10_code_rds) %>% as.data.frame()
		date0 <- readRDS(icd10_date_rds) %>% as.data.frame()
		stopifnot(identical(code0$eid, date0$eid))

		dat1 <- left_join(
			code0 %>%
				pivot_longer(-eid, names_pattern = "i(\\d+)", names_to = "pos", values_to = "code") %>%
				mutate(code = trimws(as.character(code)), pos = as.integer(pos)) %>%
				filter(!is.na(code), nzchar(code)),
			date0 %>%
				pivot_longer(-eid, names_pattern = "d(\\d+)", names_to = "pos", values_to = "date") %>%
				mutate(date = asDate2(date), pos = as.integer(pos)) %>%
				filter(!is.na(date)),
			by = c("eid", "pos")
		) %>%
			transmute(
				id = eid,
				vocabulary_id = "ICD10CM",
				code = sub("^([A-Z0-9]{3})([A-Z0-9]+)$", "\\1.\\2", toupper(gsub("\\.", "", code))),
				date
			) %>%
			PheWAS::mapCodesToPhecodes() %>%
			filter(!is.na(phecode)) %>%
			dplyr::count(id, phecode) %>%
			mutate(value = 1L) %>%
			dplyr::select(id, phecode, value) %>%
			tidyr::pivot_wider(names_from = phecode, values_from = value, values_fill = 0) %>%
			as.data.frame()

		saveRDS(dat1, paste0(indir, "/Rdata/PheWAS.rds"))
		cat("Wrote: ", paste0(indir, "/Rdata/PheWAS.rds"), "\n", sep = "")
	}
} else {
	message("MAKE_PHEWAS=FALSE; skipped PheWAS.rds. This is not needed for icd/merge/all.rds.")
}
})


run_if("audit_i67", {
# 🚩 Audit: I67 trajectory audit
pacman::p_load(dplyr, tidyr, data.table)
code0 <- read_rds("icd10.code0") %>% as.data.frame()
date0 <- read_rds("icd10.date0") %>% as.data.frame()
stopifnot(identical(code0$eid, date0$eid))

long_icd10 <- left_join(
	code0 %>% pivot_longer(-eid, names_pattern = "i(\\d+)", names_to = "pos", values_to = "code") %>% filter(!is.na(code), nzchar(code)) %>% mutate(pos = as.integer(pos)),
	date0 %>% pivot_longer(-eid, names_pattern = "d(\\d+)", names_to = "pos", values_to = "date") %>% filter(!is.na(date)) %>% mutate(pos = as.integer(pos)),
	by = c("eid", "pos")
) %>% mutate(date = asDate2(date))

traj <- long_icd10 %>%
	filter(grepl("^I67|^I60|^I61|^I63|^I64", code)) %>%
	transmute(
		eid, date,
		grp = case_when(
			grepl("^I678", code) ~ "I67.8",
			grepl("^I679", code) ~ "I67.9",
			grepl("^I67",  code) ~ "I67_other",
			grepl("^I63",  code) ~ "I63",
			grepl("^I61",  code) ~ "I61",
			grepl("^I60",  code) ~ "I60",
			grepl("^I64",  code) ~ "I64",
			TRUE ~ NA_character_
		)
	) %>% filter(!is.na(grp))

traj_wide <- traj %>% group_by(eid, grp) %>% summarise(date = min(date), .groups = "drop") %>%
	tidyr::pivot_wider(names_from = grp, values_from = date) %>%
	mutate(
		I67_any = do.call(pminDate2, across(any_of(c("I67.8", "I67.9", "I67_other")))),
		stroke_specific = do.call(pminDate2, across(any_of(c("I60", "I61", "I63", "I64")))),
		traj_class = case_when(
			!is.na(I67_any) &  is.na(stroke_specific) ~ "I67_only",
			 is.na(I67_any) & !is.na(stroke_specific) ~ "specific_only",
			!is.na(I67_any) & !is.na(stroke_specific) & I67_any < stroke_specific ~ "I67_before_specific",
			!is.na(I67_any) & !is.na(stroke_specific) & I67_any > stroke_specific ~ "specific_before_I67",
			!is.na(I67_any) & !is.na(stroke_specific) & I67_any == stroke_specific ~ "same_day",
			TRUE ~ NA_character_
		)
	)

fwrite(traj_wide, analysis_file("I67_trajectory_audit.csv"))
traj_wide %>% count(traj_class, sort = TRUE) %>% print(n = Inf)

})


run_if("fod_ref", {
# 🚩 FOD coming with UKB
date0 <- fread(paste0(indir, "/rap/fod.date.tab.gz"), header = TRUE) %>% as.data.frame()
saveRDS(date0, paste0(indir, "/Rdata/fod.date0.rds"))
fn <- fread(paste0(indir, "/common/fod.lst"), sep = "\t", header = FALSE, select = 1:2) %>% as.data.frame()
dat <- data.frame(eid = date0$eid)
for (i in 1:nrow(fn)) { 
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
  r2[, .(codes=paste(sort(unique(code)), collapse=",")), by=phecode],
  file.path(lstdir, "gp2.phecode.lst"),
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE
)
write.table(
  r3[, .(codes=paste(sort(unique(code)), collapse=",")), by=phecode],
  file.path(lstdir, "gp3.phecode.lst"),
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE
)

})


# 🚩 GP prescription CKM helper
make_gp_med_ckm <- function(gp_script_file) {
  library(data.table)
  min_hit_date <- function(date, hit) {
    keep <- !is.na(hit) & hit & !is.na(date)
    if (any(keep)) min(as.integer(date[keep])) else NA_integer_
  }
  x <- fread(
    gp_script_file,
    select=c("eid","issue_date","drug_name"),
    showProgress=TRUE
  )
  x[, eid := as.character(eid)]

  x[, gp_htn_med := grepl(
    "amlodipine|bisoprolol|atenolol|ramipril|lisinopril|losartan|valsartan|candesartan|indapamide|bendroflumethiazide|spironolactone",
    drug_name, ignore.case=TRUE)]

  x[, gp_dm_med := grepl(
    "metformin|insulin|gliclazide|glimepiride|pioglitazone|sitagliptin|dapagliflozin|empagliflozin",
    drug_name, ignore.case=TRUE)]

  x[, gp_lipid_med := grepl(
    "atorvastatin|simvastatin|rosuvastatin|pravastatin|ezetimibe|fenofibrate",
    drug_name, ignore.case=TRUE)]

  x[, issue_date := as.Date(issue_date)]

  ans <- x[, .(
    gp_htn_med=any(gp_htn_med),
    gp_dm_med=any(gp_dm_med),
    gp_lipid_med=any(gp_lipid_med),
    date_gp_htn_med=min_hit_date(issue_date, gp_htn_med),
    date_gp_dm_med=min_hit_date(issue_date, gp_dm_med),
    date_gp_lipid_med=min_hit_date(issue_date, gp_lipid_med)
  ), by=eid]
  date_cols <- grep("^date_", names(ans), value=TRUE)
  ans[, (date_cols) := lapply(.SD, as.Date, origin="1970-01-01"), .SDcols=date_cols]
  ans
}


run_if("gp", {
# 🚩 GP prescription-supported CKM evidence
gp.med.ckm <- make_gp_med_ckm(
  paste0(indir, "/rap/raw/gp_scripts.tsv.gz")
)
saveRDS(
  gp.med.ckm,
  paste0(indir, "/Rdata/gp.ckm.med.rds")
)


# 🚩 GP: disease FODs + CKM numeric/structural components
pacman::p_load(data.table, fst)


make_gp_med_ckm_v19 <- function(gp_script_file =
  paste0(indir, "/rap/raw/gp_scripts.tsv.gz")) {

  suppressPackageStartupMessages(library(data.table))
  min_hit_date <- function(date, hit) {
    keep <- !is.na(hit) & hit & !is.na(date)
    if (any(keep)) min(as.integer(date[keep])) else NA_integer_
  }

  x <- fread(gp_script_file,
             select=c("eid","issue_date","drug_name"),
             showProgress=TRUE)

  x[, eid := as.character(eid)]

  x[, gp_htn_med := grepl(
    "amlodipine|bisoprolol|atenolol|ramipril|lisinopril|losartan|valsartan|candesartan|indapamide|bendroflumethiazide|spironolactone",
    drug_name, ignore.case=TRUE)]

  x[, gp_dm_med := grepl(
    "metformin|insulin|gliclazide|glimepiride|pioglitazone|sitagliptin|dapagliflozin|empagliflozin",
    drug_name, ignore.case=TRUE)]

  x[, gp_lipid_med := grepl(
    "atorvastatin|simvastatin|rosuvastatin|pravastatin|ezetimibe|fenofibrate",
    drug_name, ignore.case=TRUE)]

  x[, issue_date := as.Date(issue_date)]

  ans <- x[, .(
    gp_htn_med = any(gp_htn_med),
    gp_dm_med = any(gp_dm_med),
    gp_lipid_med = any(gp_lipid_med),
    date_gp_htn_med=min_hit_date(issue_date, gp_htn_med),
    date_gp_dm_med=min_hit_date(issue_date, gp_dm_med),
    date_gp_lipid_med=min_hit_date(issue_date, gp_lipid_med)
  ), by=eid]
  date_cols <- grep("^date_", names(ans), value=TRUE)
  ans[, (date_cols) := lapply(.SD, as.Date, origin="1970-01-01"), .SDcols=date_cols]
  ans
}

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
  dt[, event_dt := as.integer(as.Date(event_dt))]
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
  dt[, event_dt := ckm_date(event_dt)]
  setnames(dt, code_col, "code")
  dt <- dt[!is.na(code) & nzchar(code) & !is.na(event_dt)]
  x <- long[dt, on = "code", nomatch = 0L, allow.cartesian = TRUE]
  ans <- x[, .(fod = min(event_dt)), by = .(eid, trait)]
  wide <- dcast(ans, eid ~ trait, value.var = "fod")
  old <- setdiff(names(wide), "eid")
  setnames(wide, old, paste0("fod_", src, "_", old))
  for (cc in grep("^fod_", names(wide), value = TRUE)) wide[, (cc) := ckm_date(get(cc))]
  wide
}

fod_gp2 <- make_fod_from_fst("gp2", fst_file)
fod_gp3 <- make_fod_from_fst("gp3", fst_file)
saveRDS(fod_gp2, paste0(out_dir, "/fod.gp2.rds"))
saveRDS(fod_gp3, paste0(out_dir, "/fod.gp3.rds"))

tmp <- merge_check(fod_gp2 = as.data.frame(fod_gp2), fod_gp3 = as.data.frame(fod_gp3), by = "eid")
fod_gp <- as.data.table(merge(tmp$fod_gp2, tmp$fod_gp3, by = "eid", all = TRUE))
traits <- intersect(
  sub("^fod_gp2_", "", grep("^fod_gp2_", names(fod_gp), value = TRUE)),
  sub("^fod_gp3_", "", grep("^fod_gp3_", names(fod_gp), value = TRUE))
)
for (t in traits) {
  a <- ckm_date(fod_gp[[paste0("fod_gp2_", t)]])
  b <- ckm_date(fod_gp[[paste0("fod_gp3_", t)]])
  z <- pmin(as.numeric(a), as.numeric(b), na.rm = TRUE)
  z[!is.finite(z)] <- NA_real_
  fod_gp[[paste0("fod_gp_", t)]] <- as.Date(z, origin = "1970-01-01")
}
for (cc in grep("^fod_", names(fod_gp), value = TRUE)) fod_gp[, (cc) := ckm_date(get(cc))]
saveRDS(fod_gp, paste0(out_dir, "/fod.gp.rds"))

})


run_if("biom", {
# 🚩 BBC, Prot, Met, Img + CKM components
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
  bbc_raw[[cc]] <- ckm_date(bbc_raw[[cc]])
}

# clean_biom() coerces non-numeric results such as "<6.7" to NA. Dates are
# kept outside numeric biomarker QC and merged back unchanged.
bbc_date_cols <- grep("^date_", names(bbc_raw), value = TRUE)
bbc_dates <- bbc_raw[, c("eid", bbc_date_cols), drop = FALSE]
uacr <- derive_ukb_uacr(bbc_raw)
bbc_main <- bbc_raw[, setdiff(names(bbc_raw), bbc_date_cols), drop = FALSE]
bbc_main <- clean_biom(bbc_main, miss_col = 0.2)
bbc <- Reduce(function(x, y) merge(x, y, by = "eid", all = TRUE, sort = FALSE),
              list(bbc_main, bbc_dates, uacr))
saveRDS(bbc, paste0(indir, "/Rdata/bbc.rds"))
cat("UACR assessed from numeric values:",
    sum(bbc$bb_UACR_assessed %in% TRUE, na.rm = TRUE), "\n")

prot <- fread(paste0(indir, "/rap/raw/prot.tab.gz"), sep = "\t",
              header = TRUE) %>% as.data.frame()
names(prot)[-1] <- toupper(names(prot)[-1])
prot <- clean_biom(prot, miss_col = 0.2)
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
met <- clean_biom(as.data.frame(met))
saveRDS(met, paste0(indir, "/Rdata/met.rds"))

# Preserve the original full numeric imaging object for other analyses.
img_file <- paste0(indir, "/rap/img.tab.gz")
img0 <- fread(img_file, sep = "\t", header = TRUE) %>% as.data.frame()
img_cols <- grep("_i2$", names(img0), value = TRUE)
img_num_cols <- img_cols[vapply(img0[img_cols], function(x) is.numeric(x) || is.integer(x), logical(1))]
img <- img0[, c("eid", img_num_cols), drop = FALSE]
names(img)[-1] <- paste0("img_", sub("_i2$", "", names(img)[-1]))
img <- img[rowSums(!is.na(img[, -1, drop = FALSE])) > 0, , drop = FALSE]
img <- clean_biom(img)
saveRDS(img, paste0(indir, "/Rdata/img.rds"))
rm(img0); invisible(gc())

})


run_if("hla_pca", {
# 🚩 HLA & PCA
hla <- fread(paste0(indir, "/rap/hla.tab.gz"), sep = "\t", header = TRUE) %>% as.data.frame() # 362 commas
	alleles <- fread(paste0(indir, "/common/ukb_hla_v2.txt"), header = TRUE) %>% as.data.frame() # 1 row, 362 columns
	hla <- hla %>% separate(p22182, into = names(alleles), remove = TRUE, convert = TRUE, sep = ",") %>% drop_na()
	names(hla)[-1] = paste0("hla_", names(hla)[-1])
	pc <- prcomp(hla[, -1], center = TRUE, scale. = TRUE)
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
	apoe = factor(apoe, levels=c("e3e3", "e2e2", "e2e3", "e2e4", "e3e4", "e4e4")),
	apoe.e4 = stringr::str_count(apoe, "4"),
	sp1 = ifelse(sp1.Z.rs28929474_C == 0, "ZZ", ifelse (sp1.S.rs17580_T == 0, "SS", 
		ifelse( (sp1.Z.rs28929474_C == 1 & sp1.S.rs17580_T == 1), "SZ", ifelse( (sp1.Z.rs28929474_C == 1 & sp1.S.rs17580_T == 2), "MZ", 
		ifelse( (sp1.Z.rs28929474_C == 2 & sp1.S.rs17580_T == 1), "MS", ifelse( (sp1.Z.rs28929474_C == 2 & sp1.S.rs17580_T == 2), "MM", 
		NA)))))),
	sp1.M = stringr::str_count(sp1, "M"), sp1.S = stringr::str_count(sp1, "S"), sp1.Z = stringr::str_count(sp1, "Z")		
)
saveRDS(gen, paste0(indir, "/Rdata/gen.rds"))

pgs <- fread(file.path(dir0, "data/ukb/pgs/all.pgs.gz"), header = TRUE) %>% as.data.frame()
names(pgs) <- sub("\\.score_sum$", ".pgs", names(pgs))
saveRDS(pgs, paste0(indir, "/Rdata/pgs.rds"))

#rel <- fread(file.path(gen_dir, "typ/ukb_rel.dat"), header = TRUE) %>% as.data.frame()
#rel.id <- ukbtools::ukb_gen_samples_to_remove(data = rel, ukb_with_data = dat$eid, cutoff = 0.0882) 
#rel.res <- rel_seek(rel = rel, eid = unique(c(rel$ID1, rel$ID2)), cutoff = 0.0882)
})


run_if("merge", {
# 🚩 merge (diet.rds must exist)
pacman::p_load(matrixStats)
lst <- c("phe", "bbc", "diet", "move", "ses", "death", "fod.icd10", "fod.icd9", "fod.opcs4", "fod.srd", "fod.gp", "gp.ckm.med", "img.ckm", "fod.ref", "pca", "gen", "pgs")
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
	date_mi = do.call(pminDate2, across(any_of(c("fod_icd10_cvd_cad", "fod_icd9_cvd_cad", "fod_srd_cvd_cad")))),
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
	htn.yes = if_else(!is.na(fod_srd_cvd_htn) | !is.na(fod_icd10_cvd_htn)  | !is.na(fod_ref_cvd_htn) | drug.htn == 1, 1L, 0L),
	bp_base = case_when(sbp >= 160 | dbp >= 100 ~ 0, (sbp >= 140 & sbp < 160) | (dbp >= 90 & dbp < 100) ~ 25, (sbp >= 130 & sbp < 140) | (dbp >= 80 & dbp < 90) ~ 50, (sbp >= 120 & sbp < 130) & (dbp > 0 & dbp < 80) ~ 75, (sbp > 0 & sbp < 120) & (dbp > 0 & dbp < 80) ~ 100, TRUE ~ NA_real_),
	bp.pts = if_else(coalesce(drug.htn, 0) == 1 & bp_base > 0, pmax(bp_base - 20, 0), bp_base),
	# bp1.pts = case_when(is.na(sbp) | is.na(dbp) ~ NA_real_, TRUE ~ pmax(pmin(100 - 2 * pmax(sbp - 115, 0) - 2 * pmax(dbp - 75, 0), 100), 0)),
	# bp2.pts = case_when(is.na(htn.yes) & is.na(drug.htn) ~ NA_real_, TRUE ~ pmax(100 - 30 * coalesce(htn.yes, 0) - 10 * coalesce(drug.htn, 0), 0))	
	# 
	alcohol.c = factor(alcohol_status, levels = 0:2, labels = c("never", "previous", "current")),
	alcohol_days = case_match(alcohol_freq, 1~7, 2~3.5, 3~1.5, 4~0.5, 5~0.1, 6~0, .default=NA),
	alcohol_units = case_match(alcohol_amount, 1~1.5, 2~3.5, 3~5.5, 4~8, 5~12, .default=NA),
	# alcohol.pts = case_when(alcohol_status %in% c(0, 1) ~ 100, sex == 0 ~ as.numeric(as.character(cut(alcohol_days * alcohol_units * 8 / 7, c(0, 7, 14, 21, 28, Inf), c(100, 80, 50, 20, 0), right = FALSE))), sex == 1 ~ as.numeric(as.character(cut(alcohol_days * alcohol_units * 8 / 7, c(0, 14, 28, 42, 56, Inf), c(100, 80, 50, 20, 0), right = FALSE))), TRUE ~ NA_real_),
	alcohol.pts = case_when(alcohol_status %in% c(0, 1) ~ 100, is.na(alcohol_status) ~ NA_real_, sex == 0 ~ as.numeric(as.character(cut(coalesce(alcohol_days * alcohol_units * 8 / 7, 14.1), c(0, 7, 14, 21, 28, Inf), c(100, 80, 50, 20, 0), right = FALSE))), sex == 1 ~ as.numeric(as.character(cut(coalesce(alcohol_days * alcohol_units * 8 / 7, 28.1), c(0, 14, 28, 42, 56, Inf), c(100, 80, 50, 20, 0), right = FALSE))), TRUE ~ NA_real_)
) %>% mutate(
	le8.sco = rowMeans1(across(all_of(vars.le8)), 0.7), # 
	le8r.sco = -le8.sco, le8r.3c = f3c(le8r.sco),	
	cad.pgs = cad.6m.pgs, cad.pgs.3c = f3c(cad.pgs), # 
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


run_if("audit_dates", {
# 🚩 Audit: dates across ref / ICD / GP / SRD / me
pacman::p_load(data.table, dplyr, purrr)

audit_traits <- c("cvd_cad", "cvd_stroke_i", "cvd_stroke_ih", "cvd_stroke_sh", "cvd_htn", "t2dm")
ref  <- read_rds("fod.ref")
ic10 <- read_rds("fod.icd10")
ic9  <- read_rds("fod.icd9")
op4  <- read_rds("fod.opcs4")
srd0 <- read_rds("fod.srd")
gp0  <- read_rds("fod.gp")
eids <- as.character(ref$eid)

dates_audit <- lapply(audit_traits, function(tr) {
	dat <- data.frame(eid = eids)
	dat$ref   <- get_date_by_eid(ref,  tr, eids, c("", "fod_ref_", "date_"))
	dat$icd10 <- get_date_by_eid(ic10, tr, eids, c("", "fod_icd10_", "date_"))
	dat$icd9  <- get_date_by_eid(ic9,  tr, eids, c("", "fod_icd9_", "date_"))
	dat$opcs4 <- get_date_by_eid(op4,  tr, eids, c("", "fod_opcs4_", "date_"))
	dat$srd   <- get_date_by_eid(srd0, tr, eids, c("", "fod_srd_", "date_"))
	dat$gp    <- get_date_by_eid(gp0,  tr, eids, c("", "fod_gp_", "date_"))
	dat$me    <- do.call(pminDate2, dat[c("icd10", "icd9", "opcs4", "srd", "gp")])
	dat$trait <- tr
	dat$diff_me_vs_ref <- as.numeric(dat$me  - dat$ref)
	dat$diff_gp_vs_ref <- as.numeric(dat$gp  - dat$ref)
	dat$diff_icd10_vs_ref <- as.numeric(dat$icd10 - dat$ref)
	dat
}) %>% bind_rows()

saveRDS(dates_audit, analysis_file("dates_audit.rds"))
fwrite(dates_audit, analysis_file("dates_audit.csv"))

safe_stat <- function(x, fun) {
	x <- x[is.finite(x)]
	if (!length(x)) return(NA_real_)
	as.numeric(fun(x))
}
dates_audit %>% mutate(diff_me_vs_ref_overlap = if_else(!is.na(ref) & !is.na(me), diff_me_vs_ref, NA_real_)) %>%
	group_by(trait) %>% summarise(
		n_ref = sum(!is.na(ref)),
		n_me = sum(!is.na(me)),
		n_overlap = sum(!is.na(ref) & !is.na(me)),
		median_diff_me_vs_ref = safe_stat(diff_me_vs_ref_overlap, median),
		p25 = safe_stat(diff_me_vs_ref_overlap, \(x) quantile(x, 0.25, names = FALSE)),
		p75 = safe_stat(diff_me_vs_ref_overlap, \(x) quantile(x, 0.75, names = FALSE)),
		.groups = "drop"
	) %>% as.data.frame() %>% print(row.names = FALSE)

dat_all <- read_rds("all")
#  FOD  
Y <- "cvd_afib" 
cols <- c("ref", "icd10", "icd9", "srd", "gp")
dat1 <- dat_all %>% dplyr::select(eid, any_of(paste0("fod_", cols, "_", Y))) %>%
	rename_with(~ gsub(paste0("^fod_|_", Y, "$"), "", .x), -eid) %>%
	mutate(across(any_of(cols), asDate2)) %>% drop_na(ref)
date_cols_present <- intersect(cols, names(dat1))
dat1$me <- if (length(date_cols_present)) do.call(pminDate2, dat1[date_cols_present]) else as.Date(NA)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1)); on.exit(par(op), add = TRUE)
spec <- list(
	list("ref", "icd10",c(38, 139, 210),  "REF", "ICD10"),
	list("ref", "srd",	c(220, 50, 47),   "REF", "SRD"),
	list("ref", "gp",	c(133, 153, 0),   "REF", "GP"),
	list("ref", "me",	c(108, 113, 196), "REF", "Me")
)

plot_date_diff <- function(dat, date1, date2, col = c(38, 139, 210), alpha = 80, cex = 0.3, pch = 16, date1_lab = date1, date2_lab = date2) {
	d1 <- as.Date(dat[[date1]])
	d2 <- as.Date(dat[[date2]])
	ok <- !is.na(d1) & !is.na(d2)
	d1 <- d1[ok]; d2 <- d2[ok]
	plot(d1, as.numeric(d1 - d2),
		 pch = pch, cex = cex, col = rgb(col[1], col[2], col[3], alpha, maxColorValue = 255),
		 xlab = paste0(date1_lab, " (reference)"),
		 ylab = paste0(date1_lab, " - ", date2_lab, " (days)"),
		 main = paste0(date1_lab, " vs ", date2_lab))
	abline(h = 0, lty = 2)
	invisible(NULL)
}
for (p in spec) { plot_date_diff(dat1, p[[1]], p[[2]], col = p[[3]], date1_lab = p[[4]], date2_lab = p[[5]]) }

#  BP  
lt_date <- function(x, y) !is.na(x) & !is.na(y) & x < y
gt_date <- function(x, y) !is.na(x) & !is.na(y) & x > y
dat1 <- dat_all %>% mutate(
    htn_pre_1  = lt_date(fod_icd10_cvd_htn, date_attend) | lt_date(fod_srd_cvd_htn, date_attend) | (drug.htn == 1), htn_post_1 = !htn_pre_1 & (gt_date(fod_icd10_cvd_htn, date_attend) | gt_date(fod_srd_cvd_htn, date_attend)),
    htn_pre_2  = lt_date(fod_icd10_cvd_htn, date_attend) | (drug.htn == 1), htn_post_2 = !htn_pre_2 & gt_date(fod_icd10_cvd_htn, date_attend),
    htn_pre_3  = lt_date(fod_icd10_cvd_htn, date_attend), htn_post_3 = !htn_pre_3 & gt_date(fod_icd10_cvd_htn, date_attend)
)
table(dat1$htn_pre_1, dat1$htn_post_1, useNA = "always") # ICD10 + SRD + drug.htn
table(dat1$htn_pre_2, dat1$htn_post_2, useNA = "always") # ICD10 + drug.htn
table(dat1$htn_pre_3, dat1$htn_post_3, useNA = "always") # ICD10 only

})


run_if("drug", {
# 🚩 drug⭕
prot <- fread(paste0(dir0, '/data.BIG/gwas/prot/map_3k_v1.tsv'), sep = '\t', select = c(1,3,5), header = TRUE, fill = TRUE) %>% as.data.frame()
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
dat1 <- dat %>% filter(ethnic.c == "White") %>% rename(IID = eid) %>% mutate(FID = IID)
for (Y in c("cvd_cad", "stroke", "cvd_stroke_i", "stroke_o")) {
	dat1 <- t2e(dat1, "cvd", paste0("fod_icd10_", Y), "birth_date", "date_attend", "date_lost", "date_death", date_follow_end, Y, "year")
}
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
	dplyr::select(FID, IID, ethnic.c, center, tdi, edu.sco, age, sex, bmi, height, bald, matches("^bald1|^cvd_cad|^stroke_|^hap|t2e$|^PC[1-4]$"), ends_with(".adu"), -starts_with("happy_"))

# mm <- model.matrix(~ center - 1, data = dat1)
# colnames(mm) <- paste0("CR", seq_len(ncol(mm)))
# dat1 <- cbind(dat1, as.data.table(mm))
dir.create(file.path(indir, "common"), recursive = TRUE, showWarnings = FALSE)
# Stage the completed table beside the destination; errors cannot truncate the
# existing phenotype file. Rename only after ADuLT and the entire write succeed.
phe_destination <- file.path(indir, "common", "ukb.phe")
phe_staged <- tempfile("ukb.phe.", tmpdir = dirname(phe_destination))
write.table(dat1, phe_staged, na = "NA", append = FALSE, quote = FALSE, row.names = FALSE)
if (!file.rename(phe_staged, phe_destination)) {
	unlink(phe_staged)
	stop("Could not replace ", phe_destination)
}
})
