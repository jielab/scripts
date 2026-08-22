#!/usr/bin/env Rscript


get_arg <- function(key, default = NA_character_) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- match(key, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}

parse_vec <- function(x, default = character()) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  z <- trimws(unlist(strsplit(x, "[,;[:space:]]+")))
  z[nzchar(z)]
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(as.character(x))) return(default)
  toupper(as.character(x)[1]) %in% c("TRUE", "T", "1", "YES", "Y")
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b

fmt_int <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)

elapsed_since <- function(t0) {
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.finite(sec)) return("unknown")
  if (sec < 60) return(sprintf("%.1fs", sec))
  if (sec < 3600) return(sprintf("%.1fmin", sec / 60))
  sprintf("%.2fh", sec / 3600)
}

progress_log <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

timed_eval <- function(label, expr) {
  t0 <- Sys.time()
  progress_log(label, " ...")
  out <- eval(substitute(expr), envir = parent.frame())
  progress_log(label, " done in ", elapsed_since(t0))
  out
}

ckm_make_gp_ckm_lists <- function(args) {
  suppressPackageStartupMessages(library(data.table))
  # Generate candidate Read v2 / CTV3 code lists for CKM components from
  # UK Biobank Resource 592 (primary-care coding lookups).
  #
  # Usage:
  #   Rscript ckm.R --make-gp-ckm-lists \
  #     --lookup-dir D:/data/ukb/phe/common/map/primarycare_codings \
  #     --terms D:/data/ukb/phe/common/gp_ckm_terms.tsv \
  #     --out-dir D:/data/ukb/phe/common/gp_ckm_review
  #
  # After manual review, set review_keep=1 in gp_ckm.candidates.tsv:
  #   Rscript ckm.R --make-gp-ckm-lists \
  #     --approved D:/data/ukb/phe/common/gp_ckm_review/gp_ckm.candidates.tsv \
  #     --out-dir D:/data/ukb/phe/common
  #
  # IMPORTANT: The generated *.candidate.lst files are discovery outputs.
  # Review descriptions and remove false-positive codes before copying them to:
  #   common/gp2.ckm.lst
  #   common/gp3.ckm.lst
  
  get_arg <- function(key, default = NA_character_) {
    a <- args
    i <- match(key, a)
    if (!is.na(i) && i < length(a)) a[i + 1] else default
  }
  
  lookup_dir <- get_arg("--lookup-dir", Sys.getenv("GP_CODING_DIR", unset = ""))
  terms_file <- get_arg("--terms", Sys.getenv("GP_CKM_TERMS", unset = "gp_ckm_terms.tsv"))
  approved_file <- get_arg("--approved", Sys.getenv("GP_CKM_APPROVED", unset = ""))
  out_dir <- get_arg("--out-dir", Sys.getenv("GP_CKM_REVIEW_DIR", unset = "gp_ckm_review"))
  auto_activate <- "--auto" %in% args
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_active_lists <- function(review_file, out_dir) {
    z <- fread(review_file, sep = "\t", header = TRUE, fill = TRUE,
               colClasses = "character")
    req <- c("code_system", "code", "component", "kind", "plausible_min",
             "plausible_max", "positive_op", "positive_threshold",
             "primary_use", "review_keep")
    miss <- setdiff(req, names(z))
    if (length(miss)) stop("Reviewed candidate file missing: ", paste(miss, collapse = ", "))
    keep <- tolower(trimws(z$review_keep)) %in% c("1", "true", "t", "yes", "y", "keep")
    z <- z[keep & code_system %in% c("read2", "ctv3") &
             !is.na(code) & nzchar(code)]
    if (!nrow(z)) stop("No rows have review_keep=1/yes/true in: ", review_file)
  
    for (sys in c("read2", "ctv3")) {
      zz <- z[code_system == sys]
      if (!nrow(zz)) next
      out <- zz[, .(
        codes_csv = paste(sort(unique(code)), collapse = ","),
        kind = unique(kind)[1],
        plausible_min = unique(plausible_min)[1],
        plausible_max = unique(plausible_max)[1],
        positive_op = unique(positive_op)[1],
        positive_threshold = unique(positive_threshold)[1],
        primary_use = unique(primary_use)[1]
      ), by = component]
      setcolorder(out, c("codes_csv", "component", "kind", "plausible_min",
                         "plausible_max", "positive_op", "positive_threshold",
                         "primary_use"))
      out_name <- if (sys == "read2") "gp2.ckm.lst" else "gp3.ckm.lst"
      fwrite(out, file.path(out_dir, out_name), sep = "\t",
             quote = FALSE, na = "NA")
    }
    keep_cols <- intersect(c("code_system", "code", "description", "component",
                             "kind", "primary_use", "review_keep", "review_note"),
                           names(z))
    fwrite(z[, ..keep_cols], file.path(out_dir, "gp_ckm.approved_codes.tsv"),
           sep = "\t", quote = FALSE)
    cat("Wrote manually approved active GP CKM lists to: ",
        normalizePath(out_dir, winslash = "/"), "\n", sep = "")
  }
  
  if (nzchar(approved_file)) {
    if (!file.exists(approved_file)) stop("Cannot find --approved file: ", approved_file)
    write_active_lists(approved_file, out_dir)
    return(invisible(0L))
  }
  
  if (!nzchar(lookup_dir) || !dir.exists(lookup_dir)) stop("Missing/invalid --lookup-dir")
  if (!file.exists(terms_file)) stop("Cannot find terms file: ", terms_file)
  
  terms <- fread(terms_file, sep = "\t", header = TRUE, fill = TRUE,
                 colClasses = "character")
  need_terms <- c("component", "kind", "include_regex", "exclude_regex",
                  "plausible_min", "plausible_max", "positive_op",
                  "positive_threshold", "primary_use")
  miss <- setdiff(need_terms, names(terms))
  if (length(miss)) stop("Terms file missing columns: ", paste(miss, collapse = ", "))
  
  read_any <- function(f) {
    ext <- tolower(tools::file_ext(f))
    tryCatch({
      if (ext %in% c("csv", "tsv", "txt", "tab")) {
        sep <- if (ext == "csv") "," else "\t"
        fread(f, sep = sep, header = TRUE, fill = TRUE,
              colClasses = "character", showProgress = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("Package readxl is required for Excel lookup files")
        }
        sh <- readxl::excel_sheets(f)
        rbindlist(lapply(sh, function(s) {
          z <- as.data.table(readxl::read_excel(f, sheet = s, col_types = "text"))
          z[, source_sheet := s]
          z
        }), fill = TRUE)
      } else NULL
    }, error = function(e) {
      message("Skip unreadable lookup file: ", f, " | ", conditionMessage(e))
      NULL
    })
  }
  
  files <- list.files(lookup_dir, recursive = TRUE, full.names = TRUE,
                      pattern = "\\.(csv|tsv|txt|tab|xlsx|xls)$", ignore.case = TRUE)
  if (!length(files)) stop("No lookup files found under: ", lookup_dir)
  
  # Resource 592 v4 stores the authoritative descriptions in two specific
  # worksheets. Read these directly; combining every heterogeneous worksheet
  # can misclassify the coding system and discard one of the code columns.
  resource592 <- files[tolower(basename(files)) == "all_lkps_maps_v4.xlsx"]
  read_resource592 <- function(f) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Package readxl is required")
    read_sheet <- function(sheet, system) {
      z <- as.data.table(readxl::read_excel(f, sheet = sheet, col_types = "text"))
      req <- c("read_code", "term_description")
      miss <- setdiff(req, names(z))
      if (length(miss)) stop(sheet, " missing: ", paste(miss, collapse = ", "))
      unique(z[
        !is.na(read_code) & nzchar(trimws(read_code)) &
          !is.na(term_description) & nzchar(trimws(term_description)),
        .(
          code = trimws(read_code),
          description = trimws(term_description),
          code_system = system,
          source_file = f,
          source_sheet = sheet
        )
      ])
    }
    rbindlist(list(
      read_sheet("read_v2_lkp", "read2"),
      read_sheet("read_ctv3_lkp", "ctv3")
    ), use.names = TRUE, fill = TRUE)
  }
  
  find_col <- function(nms, patterns) {
    n0 <- tolower(gsub("[^a-z0-9]+", "_", nms))
    for (p in patterns) {
      i <- grep(p, n0, perl = TRUE)
      if (length(i)) return(nms[i[1]])
    }
    NA_character_
  }
  
  infer_system <- function(f, z) {
    s <- tolower(paste(basename(f), paste(names(z), collapse = " "),
                       if ("source_sheet" %in% names(z)) paste(unique(z$source_sheet), collapse = " ") else ""))
    if (grepl("ctv3|read[ _-]?3|readv3", s)) return("ctv3")
    if (grepl("read[ _-]?2|readv2", s)) return("read2")
    NA_character_
  }
  
  normalize_one <- function(f) {
    z <- read_any(f)
    if (is.null(z) || !nrow(z)) return(NULL)
    sys <- infer_system(f, z)
    code_col <- find_col(names(z), c(
      "^read_?code$", "^read_?2_?code$", "^ctv3_?code$",
      "^concept_?id$", "^code$", "^term_?code$"
    ))
    desc_col <- find_col(names(z), c(
      "^description$", "^term$", "^term_?description$", "^preferred_?term$",
      "^meaning$", "^rubric$", "^concept_?description$", "^text$"
    ))
    if (is.na(code_col) || is.na(desc_col)) return(NULL)
    out <- z[, .(
      code = trimws(as.character(get(code_col))),
      description = trimws(as.character(get(desc_col)))
    )]
    out[, `:=`(code_system = sys, source_file = f)]
    out[!is.na(code) & nzchar(code) & !is.na(description) & nzchar(description)]
  }
  
  look <- if (length(resource592)) {
    rbindlist(lapply(resource592, read_resource592), use.names = TRUE, fill = TRUE)
  } else {
    rbindlist(lapply(files, normalize_one), fill = TRUE)
  }
  if (!nrow(look)) {
    stop("No readable code/description lookup table was detected. Inspect Resource 592 column names and extend find_col().")
  }
  look <- unique(look, by = c("code_system", "code", "description"))
  
  # If file-name inference was ambiguous, infer Read v2 vs CTV3 from code shape.
  look[is.na(code_system), code_system := fifelse(
    grepl("^[A-Za-z0-9.]{5}$", code), "read2",
    fifelse(grepl("^X[a-zA-Z0-9]{4}$|^Xa[a-zA-Z0-9]{3}$", code), "ctv3", NA_character_)
  )]
  look <- look[code_system %in% c("read2", "ctv3")]
  
  cand <- rbindlist(lapply(seq_len(nrow(terms)), function(i) {
    inc <- terms$include_regex[i]
    exc <- terms$exclude_regex[i]
    hit <- grepl(inc, look$description, ignore.case = TRUE, perl = TRUE)
    if (!is.na(exc) && nzchar(exc)) {
      hit <- hit & !grepl(exc, look$description, ignore.case = TRUE, perl = TRUE)
    }
    if (!any(hit)) return(NULL)
    x <- copy(look[hit])
    x[, `:=`(
      component = terms$component[i],
      kind = terms$kind[i],
      plausible_min = terms$plausible_min[i],
      plausible_max = terms$plausible_max[i],
      positive_op = terms$positive_op[i],
      positive_threshold = terms$positive_threshold[i],
      primary_use = terms$primary_use[i],
      include_regex = inc,
      exclude_regex = exc,
      review_keep = NA_character_,
      review_note = NA_character_
    )]
    x
  }), fill = TRUE)
  
  if (!nrow(cand)) stop("No code descriptions matched gp_ckm_terms.tsv")
  setorder(cand, code_system, component, description, code)
  fwrite(cand, file.path(out_dir, "gp_ckm.candidates.tsv"), sep = "\t", quote = FALSE)
  
  make_candidate_lst <- function(system, filename) {
    z <- cand[code_system == system]
    if (!nrow(z)) return(invisible(NULL))
    # This intentionally keeps all discovery candidates. Do not activate it;
    # edit gp_ckm.candidates.tsv and use --approved to write active lists.
    out <- z[, .(
      codes_csv = paste(sort(unique(code)), collapse = ","),
      kind = unique(kind)[1],
      plausible_min = unique(plausible_min)[1],
      plausible_max = unique(plausible_max)[1],
      positive_op = unique(positive_op)[1],
      positive_threshold = unique(positive_threshold)[1],
      primary_use = unique(primary_use)[1]
    ), by = component]
    setcolorder(out, c("codes_csv", "component", "kind", "plausible_min",
                       "plausible_max", "positive_op", "positive_threshold", "primary_use"))
    fwrite(out, file.path(out_dir, filename), sep = "\t", quote = FALSE, na = "NA")
  }
  
  make_candidate_lst("read2", "gp2.ckm.candidate.lst")
  make_candidate_lst("ctv3", "gp3.ckm.candidate.lst")
  
  summary <- cand[, .(
    candidate_codes = uniqueN(code),
    candidate_descriptions = uniqueN(description)
  ), by = .(code_system, component, kind, primary_use)]
  fwrite(summary, file.path(out_dir, "gp_ckm.candidate_summary.tsv"),
         sep = "\t", quote = FALSE)
  
  cat("Wrote candidate review files to: ", normalizePath(out_dir, winslash = "/"), "\n", sep = "")
  if (auto_activate) {
    cand[, review_keep := "1"]
    auto_file <- file.path(out_dir, "gp_ckm.auto_codes.tsv")
    fwrite(cand, auto_file, sep = "\t", quote = FALSE)
    write_active_lists(auto_file, out_dir)
    cat("Automatic activation requested: all include-regex matches surviving exclude-regex were activated.\n")
  } else {
    cat("Review gp_ckm.candidates.tsv, set review_keep=1 only for valid codes, then rerun with --approved.\n")
  }
  
  invisible(0L)
}

ckm_check_stage3_components <- function(args) {
  suppressPackageStartupMessages(library(data.table))
  get_arg <- function(key, default = NA_character_) {
    a <- args
    i <- match(key, a)
    if (!is.na(i) && i < length(a)) a[i + 1] else default
  }
  
  dir0 <- get_arg("--dir0", Sys.getenv("DIR0", unset = ifelse(.Platform$OS.type == "windows", "D:", "/mnt/d")))
  phedir <- get_arg("--phedir", Sys.getenv("PHEDIR", unset = file.path(dir0, "data", "ukb", "phe")))
  all_file <- get_arg("--all", file.path(phedir, "Rdata", "all.rds"))
  ckm_file <- get_arg("--ckm", file.path(phedir, "Rdata", "ckm.rds"))
  img_file <- get_arg("--ckm-img", file.path(phedir, "Rdata", "img.ckm.rds"))
  
  if (!file.exists(all_file)) stop("Missing: ", all_file)
  dat <- as.data.table(readRDS(all_file))
  cat("all.rds: ", format(nrow(dat), big.mark = ","), " rows; ", ncol(dat), " columns\n", sep = "")
  
  need <- c(
    "bb_UALB", "bb_UALB_flag", "bb_UCRE", "bb_UCRE_flag",
    "bb_UACR", "bb_UACR_category", "bb_UACR_assessed",
    "subclinical_ascvd", "subclinical_ascvd_observed",
    "pre_hf", "pre_hf_observed"
  )
  cat("\nColumns present\n")
  print(data.table(variable = need, present = need %in% names(dat)))
  
  safe_n <- function(expr) sum(expr %in% TRUE, na.rm = TRUE)
  if ("bb_UACR_assessed" %in% names(dat)) {
    cat("\nUACR\n")
    print(data.table(
      item = c("assessed", "exact numeric", "albumin below LOD", "interval uncertain", "A2/A3"),
      N = c(
        safe_n(dat$bb_UACR_assessed),
        if ("bb_UACR_exact" %in% names(dat)) sum(is.finite(dat$bb_UACR_exact)) else NA_integer_,
        if ("bb_UALB_below_lod" %in% names(dat)) safe_n(dat$bb_UALB_below_lod) else NA_integer_,
        if ("bb_UACR_interval_uncertain" %in% names(dat)) safe_n(dat$bb_UACR_interval_uncertain) else NA_integer_,
        if ("bb_UACR" %in% names(dat)) sum(is.finite(dat$bb_UACR) & dat$bb_UACR >= 30) else NA_integer_
      )
    ))
    if ("bb_UACR_category" %in% names(dat)) print(table(dat$bb_UACR_category, useNA = "always"))
  }
  
  cat("\nGP Stage-3 components\n")
  for (v in c("subclinical_ascvd", "pre_hf")) {
    obs <- paste0(v, "_observed")
    if (v %in% names(dat)) {
      cat(v, ":\n", sep = "")
      print(table(dat[[v]], useNA = "always"))
      if (obs %in% names(dat)) cat("observed N=", safe_n(dat[[obs]]), "\n", sep = "")
    }
  }
  
  if (file.exists(ckm_file)) {
    ckm <- as.data.table(readRDS(ckm_file))
    cat("\nckm.rds: ", format(nrow(ckm), big.mark = ","), " rows; ", ncol(ckm), " columns\n", sep = "")
    if ("ckm.stage" %in% names(ckm)) print(table(ckm$ckm.stage, useNA = "always"))
    if ("ckm_stage3_subtype" %in% names(ckm)) print(table(ckm$ckm_stage3_subtype, useNA = "always"))
    component_cols <- c(
      "uacr_mg_g", "uacr_assessed", "ckd_risk",
      "subclinical_ascvd", "subclinical_ascvd_observed",
      "pre_hf", "pre_hf_observed",
      "stage3_prevent_component", "stage3_ckd_component",
      "stage3_subclinical_ascvd_component", "stage3_prehf_component"
    )
    cat("Expanded CKM columns present: ", sum(component_cols %in% names(ckm)), "/", length(component_cols), "\n", sep = "")
  }
  
  if (file.exists(img_file)) {
    im <- as.data.table(readRDS(img_file))
    cat("\nimg.ckm.rds: ", format(nrow(im), big.mark = ","), " rows; ", ncol(im), " columns\n", sep = "")
    for (v in intersect(c("pre_hf_img_hfref_i2", "pre_hf_img_hfref_i3", "pre_hf_img_incident_i3", "stage3_prehf_interval_i2"), names(im))) {
      cat(v, ":\n", sep = "")
      print(table(im[[v]], useNA = "always"))
    }
  }
  
  invisible(0L)
}

ckm_cli_args <- commandArgs(trailingOnly = TRUE)
if ("--make-gp-ckm-lists" %in% ckm_cli_args) {
  ckm_make_gp_ckm_lists(setdiff(ckm_cli_args, "--make-gp-ckm-lists"))
  quit(save = "no", status = 0)
}
if ("--check-stage3-components" %in% ckm_cli_args) {
  ckm_check_stage3_components(setdiff(ckm_cli_args, "--check-stage3-components"))
  quit(save = "no", status = 0)
}


step_order <- c(
  "build_ckm", "summarize_ckm", "stage_clock", "data_prep", "marker_map", "biom_refine", "yy_timeline",
  "stage_specific", "causal_triage", "prediction", "consolidate"
)
start_step <- Sys.getenv("START_STEP", unset = get_arg("--start-step", get_arg("--start_step", "stage_clock")))
end_step <- Sys.getenv("END_STEP", unset = get_arg("--stop-step", get_arg("--end-step", tail(step_order, 1))))
if (is.na(start_step) || !nzchar(start_step) || tolower(start_step) == "all") start_step <- step_order[1]
if (is.na(end_step) || !nzchar(end_step) || tolower(end_step) == "all") end_step <- tail(step_order, 1)
if (!start_step %in% step_order) stop("Unknown start step: ", start_step, ". Available: ", paste(step_order, collapse = ", "), call. = FALSE)
if (!end_step %in% step_order) stop("Unknown stop step: ", end_step, ". Available: ", paste(step_order, collapse = ", "), call. = FALSE)
if (match(end_step, step_order) < match(start_step, step_order)) {
  stop("stop-step must be the same as or later than start-step", call. = FALSE)
}
run_step <- function(step) {
  i <- match(step, step_order)
  i >= match(start_step, step_order) && i <= match(end_step, step_order)
}
run_if <- function(step, expr) {
  if (run_step(step)) {
    t0 <- Sys.time()
    progress_log("========== RUN STEP: ", step, " ==========")
    eval(substitute(expr), envir = parent.frame())
    invisible(gc())
    progress_log("========== DONE STEP: ", step, " (", elapsed_since(t0), ") ==========")
  } else {
    progress_log("Skip step: ", step)
  }
}

omic_steps <- c(
  "data_prep", "marker_map", "biom_refine", "yy_timeline",
  "stage_specific", "causal_triage", "prediction"
)
analysis_requires_biom <- any(vapply(omic_steps, run_step, logical(1)))
biom_outputs_active <- analysis_requires_biom || run_step("consolidate")

options(stringsAsFactors = FALSE, timeout = max(1200, getOption("timeout")))
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(pacman::p_load(
  data.table, dplyr, tidyr, purrr, ggplot2, patchwork, ggrepel,
  survival, splines, MASS, writexl, scales, forcats
))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Global settings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir0 <- Sys.getenv("DIR0", unset = ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d"))
indir <- Sys.getenv("PHEDIR", unset = file.path(dir0, "data", "ukb", "phe"))
script_dir <- Sys.getenv("SCRIPT_DIR", unset = file.path(dir0, "scripts", "0f"))
script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file_dir <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), mustWork = FALSE))
} else {
  getwd()
}
ckm_module_dir <- Sys.getenv("CKM_SCRIPT_DIR_THIS", unset = script_file_dir)
helper_file <- file.path(script_dir, "0phe.f.R")
if (!file.exists(helper_file)) stop("Cannot find 0phe.f.R: ", helper_file, call. = FALSE)
source(helper_file)
if (!exists("vars.basic", inherits = TRUE)) stop("0phe.f.R must define vars.basic", call. = FALSE)
if (!exists("vars.le8", inherits = TRUE)) stop("0phe.f.R must define vars.le8", call. = FALSE)
ckm_vars_basic <- unique(as.character(vars.basic))
ckm_vars_le8 <- unique(as.character(vars.le8))

# The CKM pipeline has one fixed target: incident clinical CKM stage 4
# (first CHD/HF/stroke/PAD/AF) among baseline stages 0-3. Baseline stage-4
# participants are retained as Yang evidence. There is deliberately no --trait option.
trait <- "ckm4"
trait_label <- "clinical CKM stage 4"
biom_request <- tolower(gsub("[[:space:]]+", "", Sys.getenv("CKM_BIOM", unset = get_arg("--biom", "prot"))))
if (biom_request %in% c("met+prot", "met.prot", "prot.met")) biom_request <- "prot+met"
if (!biom_request %in% c("prot", "met", "prot+met")) stop("ckm.R accepts one biomarker mode per process: prot, met, or prot+met", call. = FALSE)
biom <- biom_request
biom_layers <- if (biom == "prot+met") c("prot", "met") else biom
omic_label <- switch(biom, prot = "proteomic", met = "metabolomic", `prot+met` = "multi-omic", "omic")
scan_all_biom <- as_bool(Sys.getenv("CKM_SCAN_ALL_BIOM", unset = "FALSE"))
max_biom <- as.integer(Sys.getenv("CKM_MAX_BIOM", unset = "0"))
if (!is.finite(max_biom) || max_biom < 0) stop("CKM_MAX_BIOM must be a nonnegative integer", call. = FALSE)
requested_biom <- parse_vec(Sys.getenv("CKM_BIOMARKERS", unset = ""))
ethnic_keep <- Sys.getenv("CKM_ETHNIC", unset = "all")
cli_args <- commandArgs(trailingOnly = TRUE)
# The explicit CLI flag is authoritative. The environment variable remains
# supported for direct Rscript invocations and backward compatibility.
rebuild_ckm <- "--rebuild-ckm" %in% cli_args ||
  as_bool(Sys.getenv("CKM_REBUILD", unset = "FALSE"))
merge_ckm_to_all <- as_bool(Sys.getenv("CKM_MERGE_ALL", unset = "FALSE"))
stage4_require_risk <- as_bool(Sys.getenv("CKM_STAGE4_REQUIRE_RISK", unset = "TRUE"), TRUE)
stage4_source_primary <- tolower(Sys.getenv("CKM_STAGE4_SOURCE", unset = "all"))
if (!stage4_source_primary %in% c("all", "icd10")) {
  stop("CKM_STAGE4_SOURCE must be all or icd10", call. = FALSE)
}
stage3_method <- toupper(Sys.getenv("CKM_STAGE3_METHOD", unset = "SCORE2"))
if (!stage3_method %in% c("SCORE2", "PREVENT")) stop("CKM_STAGE3_METHOD must be SCORE2 or PREVENT", call. = FALSE)
prevent_cutoff <- as.numeric(Sys.getenv("CKM_PREVENT_CUTOFF", unset = "0.20"))
if (!is.finite(prevent_cutoff) || prevent_cutoff <= 0 || prevent_cutoff >= 1) {
  stop("CKM_PREVENT_CUTOFF must be a probability strictly between 0 and 1", call. = FALSE)
}
allow_unknown_fasting_glucose <- as_bool(Sys.getenv("CKM_ALLOW_UNKNOWN_FASTING_GLUCOSE", unset = "FALSE"))
allow_unknown_fasting_tg <- as_bool(Sys.getenv("CKM_ALLOW_UNKNOWN_FASTING_TG", unset = "FALSE"))

tg_fasting_cutoff <- as.numeric(Sys.getenv("CKM_TG_FASTING_CUTOFF", unset = "1.69"))
tg_nonfasting_cutoff <- as.numeric(Sys.getenv("CKM_TG_NONFASTING_CUTOFF", unset = "1.98"))
if (!is.finite(tg_fasting_cutoff) || tg_fasting_cutoff <= 0 ||
    !is.finite(tg_nonfasting_cutoff) || tg_nonfasting_cutoff <= 0) {
  stop("TG cutoffs must be positive mmol/L values", call. = FALSE)
}
if (tg_nonfasting_cutoff < tg_fasting_cutoff) {
  stop("The nonfasting TG cutoff cannot be lower than the fasting cutoff", call. = FALSE)
}
prevent_model_primary <- tolower(Sys.getenv("CKM_PREVENT_MODEL", unset = "base"))
if (!prevent_model_primary %in% c("base", "best_available")) {
  stop("CKM_PREVENT_MODEL must be base or best_available", call. = FALSE)
}
risk_age_reference_age <- as.numeric(Sys.getenv("CKM_RISK_AGE_REFERENCE", unset = "60"))
if (!is.finite(risk_age_reference_age) || risk_age_reference_age < 30 || risk_age_reference_age > 79) risk_age_reference_age <- 60
marker_min_stage_beta <- as.numeric(Sys.getenv("CKM_MARKER_MIN_STAGE_BETA", unset = "0.10"))
marker_min_cox_loghr <- as.numeric(Sys.getenv("CKM_MARKER_MIN_COX_LOGHR", unset = as.character(log(1.05))))
marker_min_prevalence_beta <- as.numeric(Sys.getenv("CKM_MARKER_MIN_PREV_BETA", unset = "0.10"))
marker_min_proximity_beta <- as.numeric(Sys.getenv("CKM_MARKER_MIN_PROX_BETA", unset = "0.05"))
for (nm in c("marker_min_stage_beta", "marker_min_cox_loghr", "marker_min_prevalence_beta", "marker_min_proximity_beta")) {
  vv <- get(nm)
  if (!is.finite(vv) || vv < 0) stop(nm, " must be finite and nonnegative", call. = FALSE)
}
lag_min_retention <- as.numeric(Sys.getenv("CKM_LAG_MIN_RETENTION", unset = "0.50"))
if (!is.finite(lag_min_retention) || lag_min_retention < 0 || lag_min_retention > 2) lag_min_retention <- 0.50
clock_reference <- tolower(Sys.getenv("CKM_CLOCK_REFERENCE", unset = "healthy_le8"))
if (!clock_reference %in% c("healthy_le8", "population")) {
  stop("CKM_CLOCK_REFERENCE must be healthy_le8 or population", call. = FALSE)
}
clock_adjust_mode <- tolower(Sys.getenv("CKM_CLOCK_ADJUST_MODE", unset = "le8_behavior"))
if (!clock_adjust_mode %in% c("crude", "demog", "le8_behavior", "le8_full")) {
  stop("CKM_CLOCK_ADJUST_MODE must be crude, demog, le8_behavior, or le8_full", call. = FALSE)
}
clock_behavior_vars <- intersect(c("diet.pts", "pa.pts", "smoke.pts", "sleep.pts"), ckm_vars_le8)
clock_draws <- as.integer(Sys.getenv("CKM_CLOCK_DRAWS", unset = "500"))
clock_ref_n <- as.integer(Sys.getenv("CKM_CLOCK_REF_N", unset = "5000"))
clock_draw_ref_n <- as.integer(Sys.getenv("CKM_CLOCK_DRAW_REF_N", unset = "1000"))
if (!is.finite(clock_draws) || clock_draws < 0) stop("CKM_CLOCK_DRAWS must be a nonnegative integer", call. = FALSE)
if (!is.finite(clock_ref_n) || clock_ref_n < 100) stop("CKM_CLOCK_REF_N must be at least 100", call. = FALSE)
if (!is.finite(clock_draw_ref_n) || clock_draw_ref_n < 100) clock_draw_ref_n <- 1000L
clock_tau <- as.numeric(Sys.getenv("CKM_CLOCK_TAU", unset = "15"))
if (!is.finite(clock_tau) || clock_tau <= 1) stop("CKM_CLOCK_TAU must be greater than 1 year", call. = FALSE)
pred_base_mode <- tolower(Sys.getenv("CKM_PRED_BASE_MODE", unset = "behavior"))
if (!pred_base_mode %in% c("demog", "behavior", "le8", "full")) pred_base_mode <- "behavior"
K_outer <- as.integer(Sys.getenv("CKM_PRED_FOLDS", unset = "5"))
pred_screen_n <- as.integer(Sys.getenv("CKM_PRED_SCREEN_N", unset = "300"))
pred_top_n <- as.integer(Sys.getenv("CKM_PRED_TOP_N", unset = "50"))
pred_min_yin_n <- as.integer(Sys.getenv("CKM_PRED_MIN_YIN_N", unset = "1000"))
pred_min_events <- as.integer(Sys.getenv("CKM_PRED_MIN_EVENTS", unset = "50"))
if (!is.finite(K_outer) || K_outer < 2) stop("CKM_PRED_FOLDS must be at least 2", call. = FALSE)
if (!is.finite(pred_screen_n) || pred_screen_n < 5) stop("CKM_PRED_SCREEN_N must be at least 5", call. = FALSE)
if (!is.finite(pred_top_n) || pred_top_n < 1 || pred_top_n > pred_screen_n) {
  stop("CKM_PRED_TOP_N must be between 1 and CKM_PRED_SCREEN_N", call. = FALSE)
}
if (!is.finite(pred_min_yin_n) || pred_min_yin_n < 100) stop("CKM_PRED_MIN_YIN_N must be at least 100", call. = FALSE)
if (!is.finite(pred_min_events) || pred_min_events < 10) stop("CKM_PRED_MIN_EVENTS must be at least 10", call. = FALSE)
if (!is.finite(K_outer) || K_outer < 2 || K_outer > 20) K_outer <- 5L
if (!is.finite(pred_screen_n) || pred_screen_n < 20) pred_screen_n <- 300L
if (!is.finite(pred_top_n) || pred_top_n < 5) pred_top_n <- 50L
pred_top_n <- min(pred_top_n, pred_screen_n)
if (!is.finite(pred_min_yin_n) || pred_min_yin_n < 100) pred_min_yin_n <- 1000L
if (!is.finite(pred_min_events) || pred_min_events < 10) pred_min_events <- 50L
lag_years <- suppressWarnings(as.numeric(parse_vec(Sys.getenv("CKM_LAG_YEARS", unset = "0,2,5,10"))))
lag_years <- unique(lag_years[is.finite(lag_years) & lag_years >= 0])
if (!5 %in% lag_years) lag_years <- sort(unique(c(lag_years, 5)))
skip_prediction <- as_bool(Sys.getenv("CKM_SKIP_PREDICTION", unset = "FALSE"))
benchmark_glmnet <- as_bool(Sys.getenv("CKM_BENCHMARK_GLMNET", unset = "FALSE"))
seed <- as.integer(Sys.getenv("SEED", unset = "2026"))
if (!is.finite(seed)) stop("SEED must be an integer", call. = FALSE)
set.seed(seed)
met_strict_exclude_lipid_block <- as_bool(Sys.getenv("CKM_MET_STRICT_EXCLUDE_LIPID_BLOCK", unset = "TRUE"), TRUE)
yy_top_show <- as.integer(Sys.getenv("CKM_YY_TOP_SHOW", unset = "18"))
if (!is.finite(yy_top_show) || yy_top_show < 8) yy_top_show <- 18L

outbase <- Sys.getenv("CKM_OUTDIR", unset = file.path(dir0, "analysis", "ckm"))
global_outroot <- outbase
global_rawdir <- file.path(global_outroot, "raw")
outroot <- if (biom_outputs_active) file.path(outbase, biom) else global_outroot
rawdir <- file.path(outroot, "raw")
auditdir <- Sys.getenv("CKM_AUDIT_DIR", unset = file.path(outbase, "audit"))
run_shared_outputs <- as_bool(Sys.getenv("CKM_RUN_SHARED_OUTPUTS", unset = "TRUE"), TRUE)
ensure_r_directory <- function(path, label = "directory") {
  if (file.exists(path) && !dir.exists(path)) {
    stop(label, " path exists but is not a directory: ", path,
         ". Inspect it with `ls -ld` and move or rename the blocking file.",
         call. = FALSE)
  }
  dir.create(path, recursive = TRUE, showWarnings = TRUE)
  if (!dir.exists(path)) {
    stop("Cannot create ", label, ": ", path, call. = FALSE)
  }
  test_file <- file.path(path, paste0(".ckm_write_test_", Sys.getpid()))
  writable <- try(file.create(test_file), silent = TRUE)
  if (inherits(writable, "try-error") || !isTRUE(writable)) {
    stop(label, " is not writable: ", path, call. = FALSE)
  }
  unlink(test_file)
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

ensure_r_directory(outbase, "CKM output base directory")
if (!analysis_requires_biom || run_step("stage_clock")) {
  ensure_r_directory(global_rawdir, "CKM shared raw-output directory")
}
if (biom_outputs_active) {
  ensure_r_directory(outroot, "CKM biomarker-output directory")
  ensure_r_directory(rawdir, "CKM biomarker raw-output directory")
} else {
  ensure_r_directory(rawdir, "CKM shared raw-output directory")
}
setwd(outroot)

prot_meta_file <- Sys.getenv("CKM_PROT_META", unset = file.path(dir0, "data.BIG", "gwas", "ppp", "ppp_3k_b38.bed"))
met_meta_file <- Sys.getenv("CKM_MET_META", unset = file.path(indir, "common", "met.lst"))
icd10_list_file <- Sys.getenv("CKM_ICD10_LIST", unset = file.path(indir, "common", "icd10.lst"))
ckm_file <- file.path(indir, "Rdata", "ckm.rds")
all_file <- file.path(indir, "Rdata", "all.rds")
if (!file.exists(all_file)) stop("Cannot find all.rds: ", all_file, call. = FALSE)

# ckm.rds stores both SCORE2 and PREVENT stage assignments. stage3_method only
# selects which stored assignment is exposed as ckm.stage for this run.

date_follow_end_global <- if (exists("date_follow_end", inherits = TRUE)) {
  as.Date(get("date_follow_end", inherits = TRUE))
} else {
  as.Date(Sys.getenv("CKM_FOLLOW_END", unset = "2025-12-31"))
}

progress_log("CKM risk-age staging Yin-Yang pipeline")
progress_log("rebuild_ckm = ", rebuild_ckm,
             if ("--rebuild-ckm" %in% cli_args) " (explicit CLI request)" else "")
progress_log("fixed CKM target = ", trait, " (", trait_label, "); biom = ", biom,
             "; start_step = ", start_step, "; stop_step = ", end_step)
progress_log("shared_outroot = ", global_outroot)
if (biom_outputs_active) progress_log("biomarker_outroot = ", outroot)
progress_log("active_outroot = ", outroot)
progress_log("auditdir = ", auditdir)
progress_log("stage-4 primary date source = ", stage4_source_primary, "; exact ICD-10 is retained as a sensitivity definition")
progress_log("stage4_require_risk = ", stage4_require_risk, "; stage3_method = ", stage3_method, "; stage4_source = ", stage4_source_primary, "; PREVENT cutoff = ", prevent_cutoff,
             "; PREVENT primary model = ", prevent_model_primary,
             "; TG cutoffs fasting/nonfasting = ", tg_fasting_cutoff, "/", tg_nonfasting_cutoff,
             "; risk-age reference = ", risk_age_reference_age)
progress_log("clock_reference = ", clock_reference, "; clock_adjust_mode = ", clock_adjust_mode,
             "; clock_behavior_vars = ", paste(clock_behavior_vars, collapse = ", "),
             "; clock_draws = ", clock_draws,
             "; clock_ref_n = ", clock_ref_n, "; clock_draw_ref_n = ", clock_draw_ref_n,
             "; clock_tau = ", clock_tau)
progress_log("pred_base_mode = ", pred_base_mode)
progress_log("met strict lipoprotein-lipid exclusion sensitivity = ", met_strict_exclude_lipid_block, "; YY top-show = ", yy_top_show)
progress_log("vars.basic = ", paste(ckm_vars_basic, collapse = ", "))
progress_log("vars.le8 = ", paste(ckm_vars_le8, collapse = ", "))

# Shared helpers are maintained separately.
source(file.path(ckm_module_dir, "f", "ckm-helpers.R"), local = TRUE)

# CKM construction functions are maintained separately.
source(file.path(ckm_module_dir, "f", "ckm-construction.R"), local = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Biomarker data and metadata
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

flag_met_lipid_block <- function(group, subgroup, label) {
  norm_text <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    tolower(gsub("[[:space:]_]+", " ", x))
  }
  group0 <- norm_text(group)
  subgroup0 <- norm_text(subgroup)
  label0 <- tolower(as.character(label))
  label0[is.na(label0)] <- ""
  group_hit <- grepl(
    "lipoprotein|apolipoprotein|cholesterol|triglyceride|cholesteryl|free cholesterol|phospholipid|total lipids",
    group0
  )
  label_hit <- grepl(
    "(^|_)(apoa1|apob|apob2apoa1|hdl|ldl|vldl|idl)|clinical_ldl|non_hdl|remnant|total_(c|l|pl|tg)|_size$",
    label0
  )
  subgroup_hit <- grepl("hdl|ldl|vldl|idl|chylomicron|lipoprotein", subgroup0)
  group_hit | label_hit | subgroup_hit
}

read_layer_metadata <- function(layer, raw_features, rename_map) {
  if (layer == "prot") {
    out <- data.table(label = raw_features, category = NA_character_, subgroup = NA_character_)
    if (file.exists(prot_meta_file)) {
      mm <- try(fread(prot_meta_file, header = FALSE, select = c(4, 5),
                      col.names = c("label", "category"), showProgress = FALSE), silent = TRUE)
      if (!inherits(mm, "try-error") && nrow(mm)) {
        mm[, label := toupper(as.character(label))]
        mm <- mm[!is.na(label) & nzchar(label), .(
          category = {z <- as.character(category); z <- z[!is.na(z) & nzchar(z)]; if (length(z)) z[1] else NA_character_}
        ), by = label]
        out <- merge(data.table(label = raw_features), mm, by = "label", all.x = TRUE, sort = FALSE)
        out[, subgroup := NA_character_]
      }
    }
    out[is.na(category) | !nzchar(category), category := "Protein/unclassified"]
    out[, strict_exclude_lipoprotein_block := FALSE]
  } else {
    out <- data.table(
      label = raw_features,
      category = "Metabolite/unclassified",
      subgroup = NA_character_,
      strict_exclude_lipoprotein_block = FALSE
    )
    if (file.exists(met_meta_file)) {
      mm <- try(fread(met_meta_file, header = TRUE, fill = TRUE, showProgress = FALSE), silent = TRUE)
      if (!inherits(mm, "try-error") && nrow(mm)) {
        nn <- tolower(names(mm))
        named_layout <- "met_name" %in% nn
        if (!named_layout) {
          # Headerless/nonstandard files are interpreted using the documented
          # column positions: met name=2, group=5, subgroup=6. Re-reading with
          # header=FALSE prevents the first metabolite row from being lost.
          mm0 <- try(fread(met_meta_file, header = FALSE, fill = TRUE, showProgress = FALSE), silent = TRUE)
          if (!inherits(mm0, "try-error") && nrow(mm0)) mm <- mm0
          name_i <- if (ncol(mm) >= 2) 2L else NA_integer_
          group_i <- if (ncol(mm) >= 5) 5L else NA_integer_
          subgroup_i <- if (ncol(mm) >= 6) 6L else NA_integer_
        } else {
          nn <- tolower(names(mm))
          name_i <- match("met_name", nn)
          group_i <- match("group", nn)
          subgroup_i <- match("subgroup", nn)
          if (is.na(group_i) && ncol(mm) >= 5) group_i <- 5L
          if (is.na(subgroup_i) && ncol(mm) >= 6) subgroup_i <- 6L
        }
        if (!is.na(name_i)) {
          group_raw <- if (!is.na(group_i)) as.character(mm[[group_i]]) else rep(NA_character_, nrow(mm))
          subgroup_raw <- if (!is.na(subgroup_i)) as.character(mm[[subgroup_i]]) else rep(NA_character_, nrow(mm))
          mm2 <- data.table(
            label = as.character(mm[[name_i]]),
            category = group_raw,
            subgroup = subgroup_raw
          )
          mm2[, category := gsub("_", " ", category, fixed = TRUE)]
          mm2[, subgroup := gsub("_", " ", subgroup, fixed = TRUE)]
          mm2[, strict_exclude_lipoprotein_block := flag_met_lipid_block(category, subgroup, label)]
          out <- merge(data.table(label = raw_features), unique(mm2, by = "label"), by = "label", all.x = TRUE, sort = FALSE)
          out[is.na(category) | !nzchar(category), category := "Metabolite/unclassified"]
          out[is.na(subgroup) | !nzchar(subgroup), subgroup := NA_character_]
          out[is.na(strict_exclude_lipoprotein_block), strict_exclude_lipoprotein_block := FALSE]
        }
      }
    }
  }
  out[, layer := layer]
  out[, biomarker := unname(rename_map[label])]
  out[is.na(biomarker), biomarker := label]
  out[, .(biomarker, label, layer, category, subgroup, strict_exclude_lipoprotein_block)]
}

load_analysis_data <- function() {
  t_load <- Sys.time()
  progress_log("data load: starting for biom mode=", biom)

  all0 <- timed_eval("data load: reading all.rds", readRDS(all_file) %>% as.data.frame())
  progress_log("data load: all.rds rows=", fmt_int(nrow(all0)), "; cols=", fmt_int(ncol(all0)))
  idate_cols <- vapply(all0, inherits, logical(1), "IDate")
  for (cc in names(all0)[idate_cols]) {
    all0[[cc]] <- as.Date(as.numeric(unclass(all0[[cc]])), origin = "1970-01-01")
  }
  all0$eid <- as.character(all0$eid)

  ckm0 <- timed_eval("data load: reading compact ckm.rds", readRDS(ckm_file) %>% as.data.frame())
  required_ckm <- c(
    "eid", "ckm.stage.SCORE2", "ckm.stage.PREVENT", "ckm_definition",
    "ckm_stage_certainty_SCORE2", "ckm_stage_certainty_PREVENT",
    "ckm_stage3_subtype_SCORE2", "ckm_stage3_subtype_PREVENT",
    "date_first_clinical_cvd", "date_ckm4", "date_ckm_death"
  )
  missing_ckm <- setdiff(required_ckm, names(ckm0))
  if (length(missing_ckm)) {
    stop("Existing ckm.rds is missing required compact fields: ", paste(missing_ckm, collapse = ", "),
         ". Re-run ./ckm.sh --rebuild-ckm ...", call. = FALSE)
  }
  optional_ckm <- c(
    "date_ckm4_icd10", "date_ckm4_allsource",
    "ckm_stage_certainty_SCORE2", "ckm_stage_certainty_PREVENT",
    "ckm_stage2_subtype", "ckm_stage3_subtype_SCORE2", "ckm_stage3_subtype_PREVENT", "ckm_stage4_substage"
  )
  keep_ckm <- unique(c(required_ckm, intersect(optional_ckm, names(ckm0))))
  ckm0 <- dplyr::select(ckm0, dplyr::all_of(keep_ckm))
  ckm0$eid <- as.character(ckm0$eid)
  ckm0$ckm.stage <- as.integer(ckm0[[paste0("ckm.stage.", stage3_method)]])
  ckm0$ckm_stage_certainty <- as.character(ckm0[[paste0("ckm_stage_certainty_", stage3_method)]])
  ckm0$ckm_stage3_subtype <- as.character(ckm0[[paste0("ckm_stage3_subtype_", stage3_method)]])
  if (any(is.na(ckm0$ckm_definition)) || !all(ckm0$ckm_definition == "2026_AHA")) {
    stop("Existing ckm.rds lacks the current 2026 AHA definition stamp. Rebuild with ./ckm.sh --rebuild-ckm ...", call. = FALSE)
  }
  legacy_yang_missing <- sum(ckm0$ckm.stage == 4 & is.na(ckm0$date_ckm4), na.rm = TRUE)
  if (legacy_yang_missing > 0) {
    stop("Existing ckm.rds appears to be the legacy incident-date file: ", legacy_yang_missing,
         " baseline stage-4 rows have missing date_ckm4. Rebuild with ./ckm.sh --rebuild-ckm ...",
         call. = FALSE)
  }

  all0 <- timed_eval(
    "data load: joining CKM fields into all.rds frame",
    dplyr::left_join(
      dplyr::select(all0, -dplyr::any_of(c(
        "ckm.stage", "ckm.stage.SCORE2", "ckm.stage.PREVENT", "ckm_definition", "date_first_clinical_cvd", "date_ckm4", "date_ckm4_icd10", "date_ckm4_allsource", "date_ckm_death",
        "ckm_stage_certainty", "ckm_stage_certainty_SCORE2", "ckm_stage_certainty_PREVENT", "ckm_stage2_subtype", "ckm_stage3_subtype", "ckm_stage3_subtype_SCORE2", "ckm_stage3_subtype_PREVENT", "ckm_stage4_substage",
    "egfr", "uacr_mg_g", "uacr_category", "uacr_assessed",
    "uacr_lower_mg_g", "uacr_upper_mg_g", "uacr_interval_uncertain", "ckd_risk",
    "subclinical_ascvd", "subclinical_ascvd_observed", "date_subclinical_ascvd", "subclinical_ascvd_source",
    "pre_hf", "pre_hf_observed", "date_pre_hf", "pre_hf_source",
    "stage3_prevent_component", "stage3_score2_component", "stage3_ckd_component",
    "stage3_subclinical_ascvd_component", "stage3_prehf_component"
      ))),
      ckm0,
      by = "eid"
    )
  )
  rm(ckm0)
  invisible(gc())

  t0 <- Sys.time()
  progress_log("data load: converting analysis date columns")
  date_candidates <- unique(c(
    "date_attend", "date_lost", "date_death", "birth_date",
    "date_first_clinical_cvd", "date_ckm4", "date_ckm4_icd10",
    "date_ckm4_allsource", "date_ckm_death"
  ))
  for (v in intersect(date_candidates, names(all0))) all0[[v]] <- as_date2(all0[[v]])
  progress_log("data load: date conversion done in ", elapsed_since(t0))

  if (tolower(ethnic_keep) != "all" && "ethnic.c" %in% names(all0)) {
    all0 <- all0[as.character(all0$ethnic.c) == ethnic_keep, , drop = FALSE]
    progress_log("data load: ethnic filter ", ethnic_keep, " retained rows=", fmt_int(nrow(all0)))
  }

  # The quantitative CKM clock is a population construct and should not be
  # estimated only in the proteomic/metabolomic subset.  Retain a compact
  # full-UKB clock cohort before restricting to participants with biomarkers.
  clock_keep <- unique(c(
    "eid", "ckm.stage", "ckm.stage.SCORE2", "ckm.stage.PREVENT", "ckm_definition", "date_first_clinical_cvd", "date_ckm4", "date_ckm4_icd10", "date_ckm4_allsource", "date_ckm_death",
    "ckm_stage_certainty", "ckm_stage_certainty_SCORE2", "ckm_stage_certainty_PREVENT", "ckm_stage2_subtype", "ckm_stage3_subtype", "ckm_stage3_subtype_SCORE2", "ckm_stage3_subtype_PREVENT", "ckm_stage4_substage",
    "egfr", "uacr_mg_g", "uacr_category", "uacr_assessed",
    "uacr_lower_mg_g", "uacr_upper_mg_g", "uacr_interval_uncertain", "ckd_risk",
    "subclinical_ascvd", "subclinical_ascvd_observed", "date_subclinical_ascvd", "subclinical_ascvd_source",
    "pre_hf", "pre_hf_observed", "date_pre_hf", "pre_hf_source",
    "stage3_prevent_component", "stage3_score2_component", "stage3_ckd_component",
    "stage3_subclinical_ascvd_component", "stage3_prehf_component",
    "date_attend", "date_lost", "date_death", "birth_date",
    ckm_vars_basic, ckm_vars_le8
  ))
  clock_missing <- setdiff(unique(c("eid", "ckm.stage", "date_ckm4", "date_attend", ckm_vars_basic, ckm_vars_le8)), names(all0))
  if (length(clock_missing)) {
    stop("Variables required for the full-cohort stage clock are missing: ",
         paste(clock_missing, collapse = ", "), call. = FALSE)
  }
  clock_dat <- all0[, intersect(clock_keep, names(all0)), drop = FALSE]
  progress_log("data load: compact full-cohort clock rows=", fmt_int(nrow(clock_dat)),
               "; cols=", fmt_int(ncol(clock_dat)))

  if (!analysis_requires_biom) {
    rm(all0)
    invisible(gc())
    progress_log("data load: omic layers are not required for requested step range; clock-only load completed")
    return(list(
      dat = clock_dat,
      clock_dat = clock_dat,
      features = character(),
      metadata = data.table(),
      layer_counts = data.table(),
      vars_by_layer = list()
    ))
  }

  # Keep only variables needed downstream before adding thousands of omic
  # columns.  This avoids retaining all ~1,300 all.rds columns in the 41k
  # biomarker cohort and substantially lowers peak RAM.
  target_date_cols <- character()
  genetic_cols <- if (pred_base_mode == "full") {
    grep("prs|polygenic|apoe|rare|lof|plof", names(all0), value = TRUE, ignore.case = TRUE)
  } else character()
  analysis_keep <- unique(c(
    clock_keep, "ethnic.c", "drug.lipid", "drug.htn", "drug.dm", target_date_cols, genetic_cols
  ))
  analysis_keep <- intersect(analysis_keep, names(all0))
  all_small <- all0[, analysis_keep, drop = FALSE]
  rm(all0)
  invisible(gc())
  progress_log("data load: reduced phenotype frame cols=", fmt_int(ncol(all_small)))

  layer_data <- list()
  raw_features <- list()
  rename_maps <- list()
  layer_counts <- data.table()
  for (layer in biom_layers) {
    f <- file.path(indir, "Rdata", paste0(layer, ".rds"))
    if (!file.exists(f)) stop("Cannot find biomarker data: ", f, call. = FALSE)
    dd <- timed_eval(paste0("data load: reading ", layer, ".rds"), readRDS(f) %>% as.data.frame())
    dd$eid <- as.character(dd$eid)
    dd <- dd[!duplicated(dd$eid), , drop = FALSE]
    if (layer == "prot") names(dd)[-1] <- toupper(names(dd)[-1])
    layer_data[[layer]] <- dd
    raw_features[[layer]] <- setdiff(names(dd), "eid")
    rename_maps[[layer]] <- setNames(raw_features[[layer]], raw_features[[layer]])
    layer_counts <- rbind(
      layer_counts,
      data.table(layer, rds_N = nrow(dd), features = length(raw_features[[layer]]))
    )
    progress_log("data load: layer=", layer, "; rows=", fmt_int(nrow(dd)),
                 "; raw features=", fmt_int(length(raw_features[[layer]])))
  }

  if (length(biom_layers) > 1) {
    collisions <- Reduce(intersect, raw_features[biom_layers])
    if (length(collisions)) {
      progress_log("data load: resolving cross-layer feature-name collisions=", fmt_int(length(collisions)))
      for (layer in biom_layers) {
        old <- collisions
        new <- paste0(collisions, "__", layer)
        names(layer_data[[layer]])[match(old, names(layer_data[[layer]]))] <- new
        rename_maps[[layer]][old] <- new
      }
    }
  }

  for (layer in biom_layers) {
    current <- unname(rename_maps[[layer]])
    hits <- intersect(current, setdiff(names(all_small), "eid"))
    if (length(hits)) {
      for (h in hits) {
        original <- names(rename_maps[[layer]])[match(h, rename_maps[[layer]])]
        new <- paste0(h, "__", layer)
        names(layer_data[[layer]])[match(h, names(layer_data[[layer]]))] <- new
        rename_maps[[layer]][original] <- new
      }
    }
  }

  layer_eids <- lapply(layer_data, function(x) intersect(all_small$eid, x$eid))
  cohort_eids <- if (length(layer_eids) == 1) layer_eids[[1]] else Reduce(intersect, layer_eids)
  progress_log("data load: common biomarker cohort rows=", fmt_int(length(cohort_eids)))

  dat <- all_small[all_small$eid %in% cohort_eids, , drop = FALSE]
  rm(all_small)
  invisible(gc())

  for (layer in biom_layers) {
    dat <- timed_eval(
      paste0("data load: merging ", layer, " biomarkers"),
      merge(dat, layer_data[[layer]], by = "eid", all.x = TRUE, sort = FALSE)
    )
    layer_data[[layer]] <- NULL
    invisible(gc())
    progress_log("data load: after ", layer, " merge rows=", fmt_int(nrow(dat)),
                 "; cols=", fmt_int(ncol(dat)))
  }


  vars_by_layer <- lapply(
    biom_layers,
    function(layer) unname(rename_maps[[layer]][raw_features[[layer]]])
  )
  names(vars_by_layer) <- biom_layers
  features <- unique(unlist(vars_by_layer, use.names = FALSE))
  meta <- rbindlist(
    lapply(
      biom_layers,
      function(layer) read_layer_metadata(layer, raw_features[[layer]], rename_maps[[layer]])
    ),
    fill = TRUE
  )

  # Biomarker QC: retain features measured in at least 80% of the biomarker cohort.
  t0 <- Sys.time()
  progress_log("data load: biomarker QC starting; candidate features=", fmt_int(length(features)))
  nonmissing <- vapply(intersect(features, names(dat)), function(v) mean(!is.na(dat[[v]])), numeric(1))
  features <- names(nonmissing)[nonmissing >= 0.80]
  if (max_biom > 0 && length(features) > max_biom) features <- head(features, max_biom)
  meta <- meta[biomarker %in% features]
  progress_log("data load: biomarker QC retained features=", fmt_int(length(features)),
               " in ", elapsed_since(t0))

  invisible(gc())
  progress_log("data load: completed in ", elapsed_since(t_load))
  list(
    dat = dat,
    clock_dat = clock_dat,
    features = features,
    metadata = meta,
    layer_counts = layer_counts,
    vars_by_layer = vars_by_layer
  )
}

resolve_display_biomarkers <- function(features, meta, marker_scan = NULL, n = 6L) {
  ans <- character()
  if (length(requested_biom)) {
    for (x in requested_biom) {
      hit <- meta[tolower(biomarker) == tolower(x) | tolower(label) == tolower(x), biomarker]
      if (length(hit)) ans <- c(ans, hit[1])
    }
  }
  defaults <- if (biom == "prot") c("GDF15", "PCSK9", "APOE", "ANGPTL3", "FGA", "IL6") else
    if (biom == "met") c("GlycA", "ApoB", "ApoA1", "LDL_C", "HDL_C", "Glucose") else
      c("GDF15", "PCSK9", "APOE", "GlycA", "ApoB", "LDL_C")
  for (x in defaults) {
    hit <- meta[tolower(biomarker) == tolower(x) | tolower(label) == tolower(x), biomarker]
    if (length(hit)) ans <- c(ans, hit[1])
  }
  if (!is.null(marker_scan) && nrow(marker_scan)) {
    ms <- as.data.table(marker_scan)
    ms[, rank_p := pmin(stage_p, yin_cox_p, yang_proximity_p, prev_vs_yin_p, na.rm = TRUE)]
    ans <- c(ans, ms[is.finite(rank_p)][order(rank_p)]$biomarker)
  }
  ans <- unique(intersect(ans, features))
  if (length(ans) < n) ans <- unique(c(ans, head(setdiff(features, ans), n - length(ans))))
  head(ans, n)
}

prepare_target_data <- function(dat, y = "ckm4") {
  if (!identical(y, "ckm4")) {
    stop("This pipeline has one fixed target: incident clinical CKM stage 4 (ckm4).", call. = FALSE)
  }
  d <- dat
  event_date <- if ("date_first_clinical_cvd" %in% names(d)) {
    as_date2(d$date_first_clinical_cvd)
  } else {
    as_date2(d$date_ckm4)
  }
  baseline <- as_date2(d$date_attend)
  lost <- if ("date_lost" %in% names(d)) as_date2(d$date_lost) else as.Date(rep(NA, nrow(d)))
  death <- if ("date_death" %in% names(d)) as_date2(d$date_death) else as.Date(rep(NA, nrow(d)))
  follow <- rep(date_follow_end_global, nrow(d))
  censor <- pmin_date2(data.frame(lost = lost, death = death, follow = follow))
  censor[is.na(censor)] <- date_follow_end_global

  # Match 0phe.f.R::t2e(): an event recorded on the assessment date is not
  # assigned a positive or negative b2e. It remains a baseline stage-4 record
  # for categorical staging, but it is excluded from event-proximity YY curves.
  same_day_stage4 <- !is.na(event_date) & !is.na(baseline) & event_date == baseline
  baseline_stage4 <- !is.na(event_date) & !is.na(baseline) & event_date <= baseline
  yang <- !is.na(event_date) & !is.na(baseline) & event_date < baseline
  event_observed <- !is.na(event_date) & !is.na(baseline) & event_date > baseline & event_date <= censor

  stop_date <- censor
  stop_date[event_observed] <- event_date[event_observed]
  time <- as.numeric(stop_date - baseline) / 365.25
  b2e <- as.numeric(event_date - baseline) / 365.25
  b2e[same_day_stage4] <- NA_real_
  yin_eligible <- !baseline_stage4 & d$ckm.stage %in% 0:3 & is.finite(time) & time > 0

  d$target_date <- event_date
  d$censor_date <- censor
  d$b2e <- b2e
  d$baseline_stage4 <- baseline_stage4
  d$same_day_stage4 <- same_day_stage4
  d$yang <- yang
  d$yin <- yin_eligible
  d$event <- as.integer(event_observed & yin_eligible)
  d$time <- time
  d$yin_control <- yin_eligible & d$event == 0
  d
}


# Survival and CKM stage-clock functions are maintained separately.
source(file.path(ckm_module_dir, "f", "ckm-stage-clock.R"), local = TRUE)

# CKM summary figures and workbooks are maintained separately.
source(file.path(ckm_module_dir, "f", "ckm-summary.R"), local = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Marker map and Yin-Yang timeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ncore <- as.integer(Sys.getenv("CKM_NCORE", unset = as.character(min(2, max(1, parallel::detectCores() - 2)))))
if (!is.finite(ncore) || ncore < 1) ncore <- 1L
progress_log("parallel feature-scan cores = ", ncore)
feature_chunk_size <- suppressWarnings(as.integer(Sys.getenv("CKM_FEATURE_CHUNK_SIZE", unset = "0")))
map_features <- function(features, FUN, label = "feature scan") {
  total <- length(features)
  if (!total) return(list())
  chunk_size <- feature_chunk_size
  if (!is.finite(chunk_size) || chunk_size <= 0) {
    chunk_size <- min(200L, max(1L, ncore, ceiling(total / 40)))
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  chunks <- split(features, ceiling(seq_along(features) / chunk_size))
  progress_log(label, ": starting ", fmt_int(total), " feature(s) in ", length(chunks),
               " chunk(s); chunk_size=", fmt_int(chunk_size), "; ncore=", ncore)
  out <- vector("list", total)
  done <- 0L
  t_all <- Sys.time()
  for (i in seq_along(chunks)) {
    ff <- chunks[[i]]
    t0 <- Sys.time()
    part <- if (.Platform$OS.type != "windows" && ncore > 1 && length(ff) > 1) {
      parallel::mclapply(ff, FUN, mc.cores = min(ncore, length(ff)), mc.preschedule = TRUE)
    } else {
      lapply(ff, FUN)
    }
    idx <- done + seq_along(part)
    out[idx] <- part
    done <- done + length(ff)
    rm(part)
    invisible(gc())
    progress_log(label, ": chunk ", i, "/", length(chunks), " done in ", elapsed_since(t0),
                 "; completed=", fmt_int(done), "/", fmt_int(total))
  }
  progress_log(label, ": completed in ", elapsed_since(t_all))
  out
}

scan_marker_one <- function(dat, biomarker, position_map, covars, covars_med = covars) {
  cols <- unique(c("ckm.stage", "yang", "yin", "event", "time", "b2e", "yin_control", biomarker, covars, covars_med))
  d <- dat[, intersect(cols, names(dat)), drop = FALSE]
  d$.z <- standardize_ref(d[[biomarker]], d$ckm.stage == 0 & !d$yang)
  d$stage_position <- unname(position_map[as.character(d$ckm.stage)])
  stage_ok <- d$ckm.stage %in% 0:3 & is.finite(d$stage_position)
  d$stage_position_z <- NA_real_
  if (sum(stage_ok) > 100) {
    mu_pos <- mean(d$stage_position[stage_ok], na.rm = TRUE)
    sd_pos <- sd(d$stage_position[stage_ok], na.rm = TRUE)
    if (is.finite(sd_pos) && sd_pos > 0) d$stage_position_z[stage_ok] <- (d$stage_position[stage_ok] - mu_pos) / sd_pos
  }

  stage_dat <- d[stage_ok & is.finite(d$stage_position_z), , drop = FALSE]
  stage <- lm_one(stage_dat, ".z", "stage_position_z", covars)
  yin <- cox_one(d[d$yin, , drop = FALSE], ".z", covars, lag = 0)

  yd <- d[d$yang & is.finite(d$b2e) & d$b2e >= -20 & d$b2e <= 0, , drop = FALSE]
  yd$proximity10 <- yd$b2e / 10
  prox <- lm_one(yd, ".z", "proximity10", covars, min_n = 100)
  prox_med <- lm_one(yd, ".z", "proximity10", covars_med, min_n = 100)

  cd <- d[(d$yang | d$yin_control) & is.finite(d$.z), , drop = FALSE]
  cd$prevalent_yang <- as.integer(cd$yang)
  prev <- lm_one(cd, ".z", "prevalent_yang", covars, min_n = 300)
  prev_med <- lm_one(cd, ".z", "prevalent_yang", covars_med, min_n = 300)

  data.table(
    biomarker,
    stage_beta = stage$beta, stage_se = stage$se, stage_p = stage$p, stage_n = stage$n,
    yin_cox_beta = yin$beta, yin_cox_se = yin$se, yin_cox_p = yin$p, yin_cox_HR = yin$HR,
    yin_n = yin$n, yin_events = yin$events,
    yang_proximity_beta10 = prox$beta, yang_proximity_se = prox$se, yang_proximity_p = prox$p,
    yang_proximity_beta10_med = prox_med$beta, yang_proximity_p_med = prox_med$p, yang_n = prox$n,
    prev_vs_yin_beta = prev$beta, prev_vs_yin_se = prev$se, prev_vs_yin_p = prev$p,
    prev_vs_yin_beta_med = prev_med$beta, prev_vs_yin_p_med = prev_med$p, prev_n = prev$n
  )
}

make_marker_map <- function(dat, features, meta, clock_obj) {
  t_step <- Sys.time()
  position_map <- get_clock_positions(clock_obj)
  covars <- get_clock_covars(dat, "le8_behavior")
  med_covars <- safe_covars(dat, unique(c(covars, "drug.lipid", "drug.htn", "drug.dm")))
  scan_features <- if (scan_all_biom) features else resolve_display_biomarkers(features, meta, n = min(20, length(features)))
  progress_log("Marker scan features: ", fmt_int(length(scan_features)))
  res <- rbindlist(map_features(scan_features, function(p) scan_marker_one(dat, p, position_map, covars, med_covars),
                                label = "marker_map scan"), fill = TRUE)
  res[, biomarker := as.character(biomarker)]
  meta <- as.data.table(copy(meta))
  meta[, biomarker := as.character(biomarker)]
  res <- merge(res, meta, by = "biomarker", all.x = TRUE, sort = FALSE)
  res[, `:=`(
    stage_q = p.adjust(stage_p, method = "BH"),
    yin_cox_q = p.adjust(yin_cox_p, method = "BH"),
    yang_proximity_q = p.adjust(yang_proximity_p, method = "BH"),
    prev_vs_yin_q = p.adjust(prev_vs_yin_p, method = "BH"),
    yang_proximity_q_med = p.adjust(yang_proximity_p_med, method = "BH"),
    prev_vs_yin_q_med = p.adjust(prev_vs_yin_p_med, method = "BH")
  )]
  res[, medication_robust_yang :=
    is.finite(yang_proximity_q) & yang_proximity_q < .05 &
    is.finite(yang_proximity_q_med) & yang_proximity_q_med < .05 &
    sign(yang_proximity_beta10) == sign(yang_proximity_beta10_med) &
    is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 &
    is.finite(prev_vs_yin_q_med) & prev_vs_yin_q_med < .05 &
    sign(prev_vs_yin_beta) == sign(prev_vs_yin_beta_med)
  ]
  res[, marker_class := fcase(
    is.finite(yin_cox_q) & yin_cox_q < .05 & abs(yin_cox_beta) >= marker_min_cox_loghr &
      is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 & abs(prev_vs_yin_beta) >= marker_min_prevalence_beta &
      sign(yin_cox_beta) == sign(prev_vs_yin_beta), "prospective + disease-state",
    is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 & abs(prev_vs_yin_beta) >= marker_min_prevalence_beta &
      (!is.finite(yin_cox_q) | yin_cox_q >= .05 | abs(yin_cox_beta) < marker_min_cox_loghr), "prevalent marker",
    is.finite(yang_proximity_q) & yang_proximity_q < .05 & abs(yang_proximity_beta10) >= marker_min_proximity_beta &
      (!is.finite(yin_cox_q) | yin_cox_q >= .05 | abs(yin_cox_beta) < marker_min_cox_loghr), "reactive/proximity marker",
    is.finite(yin_cox_q) & yin_cox_q < .05 & abs(yin_cox_beta) >= marker_min_cox_loghr &
      (!is.finite(prev_vs_yin_q) | prev_vs_yin_q >= .05 | abs(prev_vs_yin_beta) < marker_min_prevalence_beta), "prospective-only",
    is.finite(stage_q) & stage_q < .05 & abs(stage_beta) >= marker_min_stage_beta, "stage marker",
    default = "weak/mixed"
  )]
  res[, evidence_p := pmin(stage_p, yin_cox_p, yang_proximity_p, prev_vs_yin_p, na.rm = TRUE)]
  res[!is.finite(evidence_p), evidence_p := NA_real_]
  setorder(res, evidence_p, na.last = TRUE)

  labels <- unique(c(requested_biom, res[seq_len(min(20, .N)), biomarker]))
  res[, plot_label := ifelse(biomarker %in% labels, biomarker, "")]
  class_cols <- c(
    `prospective + disease-state` = "#D95F02", `prevalent marker` = "#7570B3",
    `reactive/proximity marker` = "#E7298A", `prospective-only` = "#1B9E77",
    `stage marker` = "#E6AB02", `weak/mixed` = "grey70"
  )
  p1 <- ggplot(res, aes(stage_beta, yin_cox_beta, color = marker_class, label = plot_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(alpha = .8, size = 1.8) +
    ggrepel::geom_text_repel(size = 3.0, max.overlaps = 35, seed = seed) +
    scale_color_manual(values = class_cols) +
    labs(title = "a. Risk-age marker vs prospective risk",
         x = "Biomarker difference per 1-SD CKM risk-age position",
         y = "Incident clinical stage-4 Cox beta", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")
  p2 <- ggplot(res, aes(prev_vs_yin_beta, yang_proximity_beta10, color = marker_class, label = plot_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(alpha = .8, size = 1.8) +
    ggrepel::geom_text_repel(size = 3.0, max.overlaps = 35, seed = seed + 1) +
    scale_color_manual(values = class_cols) +
    labs(title = "b. Prevalence contrast vs Yang event proximity",
         x = "Yang vs disease-free Yin difference", y = "Yang proximity beta per 10 years", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")
  counts <- res[, .N, by = marker_class][order(-N)]
  p3 <- ggplot(counts, aes(N, fct_reorder(marker_class, N), fill = marker_class)) +
    geom_col(color = "white") + geom_text(aes(label = N), hjust = -.1, size = 3.5) +
    scale_fill_manual(values = class_cols, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, .18))) +
    labs(title = "c. Evidence classes", x = "Biomarkers", y = NULL) + theme_ckm(11)

  top <- res[seq_len(min(20, .N))]
  top[, biomarker_plot := factor(biomarker, levels = rev(biomarker))]
  top_long <- melt(top[, .(
    biomarker, biomarker_plot,
    Stage = sign(stage_beta) * pmin(-log10(pmax(stage_p, 1e-300)), 20),
    Prospective = sign(yin_cox_beta) * pmin(-log10(pmax(yin_cox_p, 1e-300)), 20),
    Prevalence = sign(prev_vs_yin_beta) * pmin(-log10(pmax(prev_vs_yin_p, 1e-300)), 20),
    Proximity = sign(yang_proximity_beta10) * pmin(-log10(pmax(yang_proximity_p, 1e-300)), 20)
  )], id.vars = c("biomarker", "biomarker_plot"), variable.name = "evidence_axis", value.name = "signed_logp")
  top_long[, biomarker_plot := factor(biomarker, levels = top$biomarker)]
  axis_cols <- c(Stage = "#2C7FB8", Prospective = "#D95F02", Prevalence = "#7570B3", Proximity = "#E7298A")
  p4 <- ggplot(top_long, aes(biomarker_plot, signed_logp, color = evidence_axis, group = evidence_axis)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_line(linewidth = .85) + geom_point(size = 2.0) +
    scale_color_manual(values = axis_cols) +
    labs(title = "d. Four-axis evidence profiles",
         subtitle = "Proteins are ordered by the combined evidence ranking; direction is encoded by the sign.",
         x = "Biomarker", y = "Signed -log10(P)", color = NULL) +
    theme_ckm(10) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 8.5), legend.position = "top")

  save_plot((p1 | p2) / (p3 | p4), biom_fig_filename("Fig3.marker_map.png"), 17, 12)
  write_tab(res, "marker_map.tsv.gz")
  marker_thresholds <- data.table(
    item = c("minimum absolute stage beta", "minimum absolute prospective log-HR", "minimum absolute Yang-vs-Yin beta", "minimum absolute Yang-proximity beta per 10 years"),
    value = c(marker_min_stage_beta, marker_min_cox_loghr, marker_min_prevalence_beta, marker_min_proximity_beta)
  )
  write_book(list(marker_thresholds = marker_thresholds, marker_map = res,
                  class_counts = counts, evidence_profiles = top_long),
             biom_fig_filename("Fig3.marker_map.out.xlsx"))
  saveRDS(res, file.path(rawdir, "marker_map.rds"), compress = "gzip")
  progress_log("marker_map: completed in ", elapsed_since(t_step))
  res
}


make_yy_timeline <- function(dat, features, meta, clock_obj, marker_scan = NULL) {
  display <- resolve_display_biomarkers(features, meta, marker_scan, n = 6)
  position_map <- get_clock_positions(clock_obj)
  ref <- dat$ckm.stage == 0 & !dat$yang
  cases <- dat[(dat$yang | (dat$yin & dat$event == 1)) & is.finite(dat$b2e) & dat$b2e >= -16 & dat$b2e <= 16, , drop = FALSE]
  if (!nrow(cases)) {
    warning("No target-event cases are available in the -16 to +16 year Yin-Yang window", call. = FALSE)
    return(list(display = display, bins = data.table(), bins_baseline = data.table(), stage_profile = data.table()))
  }
  cases$event_centered_year <- -cases$b2e
  cases$baseline_centered_year <- cases$b2e
  breaks <- seq(-16, 16, 2); mids <- head(breaks, -1) + 1

  build_bin_summary <- function(xvals, phase_fun) {
    rbindlist(lapply(display, function(p) {
      zall <- standardize_ref(dat[[p]], ref)
      z <- zall[match(cases$eid, dat$eid)]
      b <- cut(xvals, breaks = breaks, include.lowest = TRUE)
      x <- data.table(bin = b, z = z)[is.finite(z), .(mean_z = mean(z), se = sd(z) / sqrt(.N), n = .N), by = bin]
      full <- data.table(bin = factor(levels(b), levels = levels(b)), mid = mids)
      x <- merge(full, x, by = "bin", all.x = TRUE, sort = FALSE)
      x[n < 20, `:=`(mean_z = NA_real_, se = NA_real_)]
      x[, phase := phase_fun(mid)]
      x[, biomarker := p]
      x
    }), fill = TRUE)
  }

  bins <- build_bin_summary(cases$event_centered_year, function(mid) ifelse(mid < 0, "Yin (pre-event)", "Yang (post-event)"))
  bins_base <- build_bin_summary(cases$baseline_centered_year, function(mid) ifelse(mid < 0, "Yang before baseline", "Yin after baseline"))
  hist0 <- data.table(mid = mids, count = hist(cases$event_centered_year, breaks = breaks, plot = FALSE)$counts)
  histb <- data.table(mid = mids, count = hist(cases$baseline_centered_year, breaks = breaks, plot = FALSE)$counts)
  hist_max <- max(c(hist0$count, histb$count), na.rm = TRUE)
  if (!is.finite(hist_max) || hist_max <= 0) hist_max <- 1

  stage_profile <- rbindlist(lapply(display, function(p) {
    z <- standardize_ref(dat[[p]], ref)
    data.table(stage = dat$ckm.stage, z = z)[stage %in% 0:4 & is.finite(z), .(
      mean_z = mean(z), se = sd(z) / sqrt(.N), n = .N
    ), by = stage][, `:=`(position = unname(position_map[as.character(stage)]), biomarker = p)]
  }), fill = TRUE)

  legend_one_line <- guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE))
  panel_legend_theme <- theme(
    legend.position = "top", legend.justification = "left",
    legend.box.just = "left", legend.text = element_text(size = 8.3),
    legend.key.width = grid::unit(.45, "cm"), legend.spacing.x = grid::unit(.10, "cm")
  )
  draw_timeline_panel <- function(hist_dat, line_dat, xlab, title, subtitle, left_fill, right_fill) {
    ggplot() +
      geom_col(data = hist_dat[mid < 0], aes(mid, count), width = 1.9, fill = left_fill, color = "white") +
      geom_col(data = hist_dat[mid > 0], aes(mid, count), width = 1.9, fill = right_fill, color = "white") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      geom_line(data = line_dat, aes(mid, mean_z * hist_max / 2 + hist_max / 2, color = biomarker, group = biomarker), linewidth = .9, na.rm = TRUE) +
      geom_point(data = line_dat, aes(mid, mean_z * hist_max / 2 + hist_max / 2, color = biomarker), size = 1.9, na.rm = TRUE) +
      scale_y_continuous("Target-event cases", sec.axis = sec_axis(~ (. - hist_max / 2) * 2 / hist_max, name = "Biomarker z-score")) +
      labs(title = title, subtitle = subtitle, x = xlab, color = NULL) +
      theme_ckm(11) + legend_one_line + panel_legend_theme
  }

  p1 <- draw_timeline_panel(
    hist0, bins, "Years relative to target event",
    "a. Yin-Yang event-centered profile for CKM stage 4",
    "Bars: pre-event Yin (left) and post-event Yang (right). Connected lines compare cohorts, not repeated measures.",
    "#D9EAF7", "#F6DDD8"
  )
  p2 <- draw_timeline_panel(
    histb, bins_base, "Years relative to baseline",
    "b. Same proteins on the baseline-centered axis",
    "Negative: Yang years before baseline; positive: Yin years after baseline. Same-day records are excluded from b2e.",
    "#F6DDD8", "#D9EAF7"
  )
  p3 <- ggplot(stage_profile, aes(position, mean_z, color = biomarker, group = biomarker)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_line(linewidth = .9) + geom_point(size = 2.1) +
    geom_errorbar(aes(ymin = mean_z - 1.96 * se, ymax = mean_z + 1.96 * se), width = .08, alpha = .5) +
    scale_x_continuous(breaks = clock_obj$primary[stage %in% 0:3, position_from_stage0], labels = unname(stage_short_labels[as.character(0:3)])) +
    labs(title = "c. Proteins on the prospective CKM risk-age road",
         subtitle = "Risk-age is equivalent hazard age, not elapsed time to the next stage.",
         x = "CKM risk-age advancement vs s0 (years)", y = "Mean biomarker z-score", color = NULL) +
    theme_ckm(11) + legend_one_line + panel_legend_theme
  p4 <- ggplot(stage_profile, aes(factor(stage, levels = 0:4, labels = unname(stage_short_labels[as.character(0:4)])), mean_z, fill = biomarker)) +
    geom_col(position = "dodge") +
    labs(title = "d. Same proteins on the categorical CKM axis",
         subtitle = "s4 is the prevalent Yang disease state; it is not an observed s3→s4 transition.",
         x = "CKM stage", y = "Mean biomarker z-score", fill = NULL) +
    theme_ckm(11) + legend_one_line + panel_legend_theme

  save_plot((p1 | p2) / (p3 | p4), biom_fig_filename("Fig2.yy_timeline.png"), 17, 13.5)
  write_tab(bins, "yy_timeline.bins.tsv.gz")
  write_tab(bins_base, "yy_timeline.baseline_bins.tsv.gz")
  write_tab(stage_profile, "yy_timeline.stage_profile.tsv")
  write_book(list(timeline_bins = bins, baseline_bins = bins_base,
                  quantitative_stage_profile = stage_profile),
             biom_fig_filename("Fig2.yy_timeline.out.xlsx"))
  list(display = display, bins = bins, bins_baseline = bins_base, stage_profile = stage_profile)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stage-context-specific biomarker effects
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

interaction_p <- function(dat, zvar, covars) {
  dd <- dat[dat$yin & dat$ckm.stage %in% 0:3, unique(c("time", "event", "ckm.stage", zvar, covars)), drop = FALSE]
  dd$stage_f <- factor(dd$ckm.stage, levels = 0:3)
  dd <- dd[complete.cases(dd) & is.finite(dd$time) & dd$time > 0, , drop = FALSE]
  if (nrow(dd) < 500 || sum(dd$event, na.rm = TRUE) < 30) return(NA_real_)

  f0 <- as.formula(paste0(
    "Surv(time,event) ~ ", bt(zvar), " + stage_f",
    if (length(covars)) paste0(" + ", paste(bt(covars), collapse = " + ")) else ""
  ))
  f1 <- as.formula(paste0(
    "Surv(time,event) ~ ", bt(zvar), " * stage_f",
    if (length(covars)) paste0(" + ", paste(bt(covars), collapse = " + ")) else ""
  ))
  fit0 <- try(coxph(f0, data = dd, ties = "efron", model = FALSE, x = FALSE, y = FALSE), silent = TRUE)
  fit1 <- try(coxph(f1, data = dd, ties = "efron", model = FALSE, x = FALSE, y = FALSE), silent = TRUE)
  if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) return(NA_real_)

  ll0 <- logLik(fit0)
  ll1 <- logLik(fit1)
  df <- attr(ll1, "df") - attr(ll0, "df")
  stat <- 2 * (as.numeric(ll1) - as.numeric(ll0))
  if (length(stat) != 1L || !is.finite(stat) ||
      length(df) != 1L || !is.finite(df) || df <= 0) return(NA_real_)
  as.numeric(pchisq(max(stat, 0), df = df, lower.tail = FALSE))
}

stage_specific_one <- function(dat, biomarker, covars) {
  biomarker_id <- as.character(biomarker)[1]
  fail_result <- function(msg) {
    list(
      beta_rows = data.table(),
      summary_row = data.table(
        biomarker = biomarker_id,
        interaction_p = NA_real_,
        profile = NA_character_, preclinical_profile = NA_character_,
        yang_shift = NA_real_, yang_direction = NA_character_,
        delta01 = NA_real_, delta12 = NA_real_, delta23 = NA_real_, delta34 = NA_real_,
        analysis_ok = FALSE,
        error = as.character(msg)
      ),
      mean_rows = data.table()
    )
  }

  tryCatch({
    cols <- unique(c("ckm.stage", "yang", "yin", "event", "time", biomarker_id, covars))
    d <- dat[, intersect(cols, names(dat)), drop = FALSE]
    if (!biomarker_id %in% names(d)) return(fail_result("missing biomarker column"))

    z0 <- standardize_ref(d[[biomarker_id]], d$ckm.stage == 0 & !d$yang)
    if (sum(is.finite(z0)) < 200) return(fail_result("too few finite standardized values"))
    d$.z <- z0

    beta_rows <- rbindlist(lapply(0:3, function(s0) {
      rr0 <- cox_one(
        d[d$yin & d$ckm.stage == s0, , drop = FALSE],
        ".z",
        covars,
        min_n = 150,
        min_event = 8
      )
      rr0[, `:=`(
        biomarker = biomarker_id,
        baseline_stage = as.integer(s0)
      )]
      rr0
    }), fill = TRUE, use.names = TRUE)

    ip <- interaction_p(d, ".z", covars)
    mean_rows <- data.table(stage = d$ckm.stage, z = d$.z)[
      stage %in% 0:4 & is.finite(z),
      .(mean_z = mean(z), se = sd(z) / sqrt(.N), n = .N),
      by = stage
    ][order(stage)]
    mean_rows <- merge(
      data.table(stage = 0:4),
      mean_rows,
      by = "stage",
      all.x = TRUE,
      sort = FALSE
    )
    mean_rows[, delta := mean_z - shift(mean_z)]
    mean_rows[, biomarker := biomarker_id]

    deltas <- mean_rows[stage %in% 1:4, delta]
    if (length(deltas) < 4) deltas <- c(deltas, rep(NA_real_, 4 - length(deltas)))
    pre_deltas <- deltas[1:3]
    sig <- sign(pre_deltas[is.finite(pre_deltas) & abs(pre_deltas) > 0.05])
    profile <- if (!length(sig)) {
      "flat"
    } else if (all(sig >= 0)) {
      "monotone increasing"
    } else if (all(sig <= 0)) {
      "monotone decreasing"
    } else if (length(sig) >= 2 && sig[1] > 0 && tail(sig, 1) < 0) {
      "preclinical inverted-U"
    } else if (length(sig) >= 2 && sig[1] < 0 && tail(sig, 1) > 0) {
      "preclinical U-shaped"
    } else {
      "non-monotone"
    }
    yang_shift <- safe_num(deltas[4])
    yang_direction <- if (!is.finite(yang_shift) || abs(yang_shift) <= 0.05) "minimal" else if (yang_shift > 0) "higher in Yang" else "lower in Yang"

    summary_row <- data.table(
      biomarker = biomarker_id,
      interaction_p = safe_num(ip)[1],
      profile = profile,
      preclinical_profile = profile,
      yang_shift = yang_shift,
      yang_direction = yang_direction,
      delta01 = safe_num(deltas[1]),
      delta12 = safe_num(deltas[2]),
      delta23 = safe_num(deltas[3]),
      delta34 = safe_num(deltas[4]),
      analysis_ok = TRUE,
      error = NA_character_
    )

    list(
      beta_rows = beta_rows,
      summary_row = summary_row,
      mean_rows = mean_rows
    )
  }, error = function(e) fail_result(conditionMessage(e)))
}

normalize_stage_specific_result <- function(z, feature_id) {
  feature_id <- as.character(feature_id)[1]
  fallback <- list(
    beta_rows = data.table(),
    summary_row = data.table(
      biomarker = feature_id,
      interaction_p = NA_real_,
      profile = NA_character_, preclinical_profile = NA_character_,
      yang_shift = NA_real_, yang_direction = NA_character_,
      delta01 = NA_real_, delta12 = NA_real_, delta23 = NA_real_, delta34 = NA_real_,
      analysis_ok = FALSE,
      error = "malformed worker result"
    ),
    mean_rows = data.table()
  )
  if (is.null(z) || !is.list(z)) return(fallback)

  br <- z[["beta_rows"]]
  sr <- z[["summary_row"]]
  mr <- z[["mean_rows"]]

  if (is.null(sr) || !(is.data.frame(sr) || data.table::is.data.table(sr)) || nrow(sr) != 1) {
    return(fallback)
  }
  sr <- as.data.table(sr)
  sr[, biomarker := feature_id]
  if (is.null(br) || !(is.data.frame(br) || data.table::is.data.table(br))) br <- data.table()
  if (is.null(mr) || !(is.data.frame(mr) || data.table::is.data.table(mr))) mr <- data.table()
  br <- as.data.table(br)
  mr <- as.data.table(mr)
  if (nrow(br)) br[, biomarker := feature_id]
  if (nrow(mr)) mr[, biomarker := feature_id]

  list(beta_rows = br, summary_row = sr, mean_rows = mr)
}

make_stage_specific <- function(dat, features, meta, marker_scan) {
  t_step <- Sys.time()
  covars <- get_base_covars(dat)
  marker_scan <- as.data.table(copy(marker_scan)); marker_scan[, biomarker := as.character(biomarker)]
  meta <- as.data.table(copy(meta)); meta[, biomarker := as.character(biomarker)]
  scan_features <- if (scan_all_biom) as.character(features) else unique(as.character(c(
    resolve_display_biomarkers(features, meta, marker_scan, 12),
    marker_scan[seq_len(min(50, .N)), biomarker]
  )))
  scan_features <- intersect(scan_features, names(dat))
  if (!length(scan_features)) stop("stage_specific has no biomarkers to scan", call. = FALSE)
  progress_log("stage_specific: scan features=", fmt_int(length(scan_features)))

  probe <- normalize_stage_specific_result(stage_specific_one(dat, scan_features[1], covars), scan_features[1])
  if (!nrow(probe$summary_row)) stop("stage_specific probe returned no summary row", call. = FALSE)
  rr_raw <- map_features(scan_features, function(feature_id) stage_specific_one(dat, as.character(feature_id)[1], covars), label = "stage_specific scan")
  if (length(rr_raw) != length(scan_features)) {
    warning("stage_specific worker-result count mismatch: expected ", length(scan_features), ", obtained ", length(rr_raw), call. = FALSE)
    length(rr_raw) <- length(scan_features)
  }
  rr <- Map(function(z, feature_id) normalize_stage_specific_result(z, feature_id), rr_raw, scan_features)
  betas <- rbindlist(lapply(rr, `[[`, "beta_rows"), fill = TRUE, use.names = TRUE)
  sums <- rbindlist(lapply(rr, `[[`, "summary_row"), fill = TRUE, use.names = TRUE)
  means <- rbindlist(lapply(rr, `[[`, "mean_rows"), fill = TRUE, use.names = TRUE)
  if (nrow(sums) != length(scan_features)) stop("stage_specific internal collection failure: expected ", length(scan_features), " summary rows, obtained ", nrow(sums), call. = FALSE)

  sums[, biomarker := as.character(biomarker)]
  if (nrow(betas)) betas[, biomarker := as.character(biomarker)]
  if (nrow(means)) means[, biomarker := as.character(biomarker)]
  sums[, interaction_q := p.adjust(interaction_p, method = "BH")]
  sums[, stage_specific := analysis_ok %in% TRUE & is.finite(interaction_q) & interaction_q < .05]
  sums <- merge(sums, meta, by = "biomarker", all.x = TRUE, sort = FALSE)
  if (nrow(betas)) {
    betas[, stage_p_q := p.adjust(p, method = "BH"), by = baseline_stage]
    betas[, baseline_stage_label := factor(baseline_stage, levels = 0:3, labels = unname(stage_short_labels[as.character(0:3)]))]
  }
  audit <- sums[, .(scanned = .N, valid = sum(analysis_ok %in% TRUE), failed = sum(!(analysis_ok %in% TRUE)),
                     interaction_p_finite = sum(is.finite(interaction_p)), interaction_q_lt_005 = sum(is.finite(interaction_q) & interaction_q < .05))]
  error_counts <- sums[analysis_ok == FALSE & !is.na(error), .N, by = error][order(-N)]

  profile_counts <- sums[analysis_ok == TRUE & !is.na(profile), .N, by = profile][order(-N)]
  p1 <- if (nrow(profile_counts)) {
    ggplot(profile_counts, aes(N, fct_reorder(profile, N), fill = profile)) +
      geom_col(color = "white") + geom_text(aes(label = N), hjust = -.1, size = 3.6) +
      scale_x_continuous(expand = expansion(mult = c(0, .2))) +
      labs(title = "a. Preclinical s0-s3 protein-profile shapes", x = "Proteins", y = NULL, fill = NULL) +
      theme_ckm(12) + theme(legend.position = "none")
  } else ggplot() + annotate("text", x = 0, y = 0, label = "No valid profile summaries") + theme_void()

  ranked <- unique(c(
    sums[stage_specific == TRUE][order(interaction_q), biomarker],
    sums[analysis_ok == TRUE & is.finite(interaction_p)][order(interaction_p), biomarker],
    marker_scan[order(evidence_p), biomarker]
  ))
  ranked <- intersect(ranked, unique(c(betas$biomarker, means$biomarker)))
  top_interaction <- head(ranked, 8)
  beta_show <- betas[biomarker %in% top_interaction]
  if (nrow(beta_show)) beta_show[, biomarker := factor(biomarker, levels = top_interaction)]
  p2 <- if (nrow(beta_show)) {
    ggplot(beta_show, aes(baseline_stage_label, beta, group = 1)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_line(linewidth = .85, color = "#2C7FB8") + geom_point(size = 2.0, color = "#2C7FB8") +
      geom_errorbar(aes(ymin = beta - 1.96 * se, ymax = beta + 1.96 * se), width = .08, alpha = .5, color = "#2C7FB8") +
      facet_wrap(~ biomarker, ncol = 4, scales = "free_y") +
      labs(title = "b. Top stage-context-specific associations with future CKM stage 4",
           subtitle = "Two-row display. FDR hits are prioritized; otherwise the lowest interaction P values are shown.",
           x = "Baseline CKM stage", y = "Cox beta per SD") + theme_ckm(11)
  } else ggplot() + annotate("text", x = 0, y = 0, label = "No valid stage-specific Cox estimates") + theme_void()

  means2 <- merge(means, sums[, .(biomarker, profile, interaction_q)], by = "biomarker", all.x = TRUE, sort = FALSE)
  profile_order <- c("monotone increasing", "monotone decreasing", "preclinical inverted-U", "preclinical U-shaped", "non-monotone")
  representatives <- rbindlist(lapply(profile_order, function(pr) {
    z <- unique(means2[profile == pr & is.finite(interaction_q), .(biomarker, interaction_q)])[order(interaction_q)]
    head(z, 2)
  }), fill = TRUE)
  rep_ids <- unique(representatives$biomarker)
  mean_plot <- means2[biomarker %in% rep_ids]
  if (nrow(mean_plot)) {
    mean_plot[, stage_label := factor(stage, levels = 0:4, labels = unname(stage_short_labels[as.character(0:4)]))]
    mean_plot[, facet_label := paste0(biomarker, "\n", profile)]
  }
  p3 <- if (nrow(mean_plot)) {
    ggplot(mean_plot, aes(stage_label, mean_z, group = 1)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_line(linewidth = .85, color = "#D95F02") + geom_point(size = 2.0, color = "#D95F02") +
      geom_errorbar(aes(ymin = mean_z - 1.96 * se, ymax = mean_z + 1.96 * se), width = .08, alpha = .45, color = "#D95F02") +
      facet_wrap(~ facet_label, ncol = 5, scales = "free_y") +
      labs(title = "c. Representative preclinical protein trajectories",
           subtitle = "Two-row display. The s3→s4 point is a separate cross-sectional Yang disease-state contrast.",
           x = "CKM stage", y = "Mean protein z-score") + theme_ckm(10)
  } else ggplot() + annotate("text", x = 0, y = 0, label = "No representative stage trajectories") + theme_void()

  save_plot(p1 / p2 / p3 + plot_layout(heights = c(.55, 1.25, 1.25)), biom_fig_filename("Fig4.stage_specific.png"), 17, 17)
  write_tab(audit, "stage_specific.scan_audit.tsv")
  write_tab(error_counts, "stage_specific.error_counts.tsv")
  write_tab(betas, "stage_specific.cox.tsv.gz")
  write_tab(sums, "stage_specific.summary.tsv.gz")
  write_tab(means, "stage_specific.means.tsv.gz")
  write_book(list(scan_audit = audit, error_counts = error_counts, stage_specific_summary = sums,
                  stage_specific_cox = betas, stage_means = means),
             biom_fig_filename("Fig4.stage_specific.out.xlsx"))
  out <- list(betas = betas, summary = sums, means = means, audit = audit, error_counts = error_counts)
  saveRDS(out, file.path(rawdir, "stage_specific.rds"), compress = "gzip")
  progress_log("stage_specific: completed in ", elapsed_since(t_step))
  out
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Causal-compatible triage (not causal proof)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_causal_triage <- function(dat, features, meta, marker_scan, stage_specific) {
  t_step <- Sys.time()
  covars <- get_base_covars(dat)
  marker_scan <- as.data.table(copy(marker_scan)); marker_scan[, biomarker := as.character(biomarker)]
  if (!"medication_robust_yang" %in% names(marker_scan)) marker_scan[, medication_robust_yang := FALSE]
  meta <- as.data.table(copy(meta)); meta[, biomarker := as.character(biomarker)]
  stage_summary <- if (!is.null(stage_specific$summary)) as.data.table(copy(stage_specific$summary)) else data.table()
  if (nrow(stage_summary)) stage_summary[, biomarker := as.character(biomarker)]
  scan_features <- if (scan_all_biom) as.character(features) else {
    stage_hits <- if (nrow(stage_summary) && "stage_specific" %in% names(stage_summary)) stage_summary[stage_specific == TRUE, biomarker] else character()
    unique(as.character(c(marker_scan[seq_len(min(100, .N)), biomarker], stage_hits)))
  }
  scan_features <- intersect(scan_features, names(dat))
  progress_log("causal_triage: scan features=", fmt_int(length(scan_features)), "; lag years=", paste(lag_years, collapse = ","))

  lag_res <- rbindlist(map_features(scan_features, function(feature_id) {
    feature_id <- as.character(feature_id)[1]
    cols <- unique(c("ckm.stage", "yang", "yin", "event", "time", feature_id, covars))
    d <- dat[, intersect(cols, names(dat)), drop = FALSE]
    if (!feature_id %in% names(d)) return(data.table())
    d$.z <- standardize_ref(d[[feature_id]], d$ckm.stage == 0 & !d$yang)
    rbindlist(lapply(lag_years, function(lg) {
      rr0 <- cox_one(d[d$yin, , drop = FALSE], ".z", covars, lag = lg)
      rr0[, `:=`(biomarker = feature_id, lag_year = as.numeric(lg))]
      rr0
    }), fill = TRUE)
  }, label = "causal_triage lag scan"), fill = TRUE, use.names = TRUE)
  if (!nrow(lag_res)) stop("causal_triage produced no lagged Cox results", call. = FALSE)
  lag_res[, `:=`(biomarker = as.character(biomarker), lag_year = safe_num(lag_year))]
  lag_res[, q := p.adjust(p, method = "BH"), by = lag_year]
  wide <- data.table::dcast(lag_res, biomarker ~ lag_year, value.var = c("beta", "p", "q"))
  wide[, biomarker := as.character(biomarker)]
  setnames(wide, gsub("^beta_", "lag_beta_", names(wide)))
  setnames(wide, gsub("^p_", "lag_p_", names(wide)))
  setnames(wide, gsub("^q_", "lag_q_", names(wide)))
  key_audit <- data.table(
    object = c("marker_scan", "lag_wide", "stage_summary", "metadata"),
    class = c(class(marker_scan$biomarker)[1], class(wide$biomarker)[1], if (nrow(stage_summary)) class(stage_summary$biomarker)[1] else "empty", class(meta$biomarker)[1]),
    N = c(nrow(marker_scan), nrow(wide), nrow(stage_summary), nrow(meta)),
    unique_keys = c(uniqueN(marker_scan$biomarker), uniqueN(wide$biomarker), if (nrow(stage_summary)) uniqueN(stage_summary$biomarker) else 0L, uniqueN(meta$biomarker))
  )
  write_tab(key_audit, "causal_triage.key_audit.tsv")
  triage <- merge(marker_scan, wide, by = "biomarker", all.x = TRUE, sort = FALSE)
  stage_cols <- c("biomarker", "interaction_p", "interaction_q", "profile", "preclinical_profile", "yang_shift", "yang_direction", "delta01", "delta12", "delta23", "delta34", "analysis_ok")
  if (nrow(stage_summary)) {
    stage_add <- stage_summary[, intersect(stage_cols, names(stage_summary)), with = FALSE]
    triage <- merge(triage, stage_add, by = "biomarker", all.x = TRUE, sort = FALSE)
  } else {
    triage[, `:=`(interaction_p = NA_real_, interaction_q = NA_real_, profile = NA_character_, delta01 = NA_real_, delta12 = NA_real_, delta23 = NA_real_, delta34 = NA_real_, analysis_ok = NA)]
    warning("stage_specific summary is empty; causal triage will proceed without interaction evidence", call. = FALSE)
  }
  meta_cols <- setdiff(names(meta), names(triage))
  if (length(meta_cols)) triage <- merge(triage, meta[, c("biomarker", meta_cols), with = FALSE], by = "biomarker", all.x = TRUE, sort = FALSE)
  lag5b <- if ("lag_beta_5" %in% names(triage)) triage$lag_beta_5 else rep(NA_real_, nrow(triage))
  lag5p <- if ("lag_p_5" %in% names(triage)) triage$lag_p_5 else rep(NA_real_, nrow(triage))
  lag5q <- if ("lag_q_5" %in% names(triage)) triage$lag_q_5 else rep(NA_real_, nrow(triage))
  triage[, `:=`(lag5_beta = safe_num(lag5b), lag5_p = safe_num(lag5p), lag5_q = safe_num(lag5q))]
  triage[, retention5 := abs(lag5_beta) / pmax(abs(yin_cox_beta), 1e-8)]
  triage[, attenuation5 := 1 - retention5]
  triage[, lag_persistent :=
    is.finite(yin_cox_q) & yin_cox_q < .05 & is.finite(lag5_q) & lag5_q < .05 &
    is.finite(yin_cox_beta) & abs(yin_cox_beta) >= marker_min_cox_loghr &
    is.finite(lag5_beta) & abs(lag5_beta) >= marker_min_cox_loghr &
    sign(lag5_beta) == sign(yin_cox_beta) & is.finite(retention5) & retention5 >= lag_min_retention]
  triage[, evidence_class := fcase(
    lag_persistent %in% TRUE & is.finite(interaction_q) & interaction_q < .05, "stage-specific lag-persistent association",
    lag_persistent %in% TRUE, "lag-persistent prospective association",
    is.finite(yin_cox_beta) & is.finite(prev_vs_yin_beta) & sign(yin_cox_beta) != sign(prev_vs_yin_beta) &
      is.finite(yin_cox_q) & yin_cox_q < .05 & is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05, "mirror/opposing",
    medication_robust_yang %in% TRUE & is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 &
      is.finite(yang_proximity_q) & yang_proximity_q < .05 & !(lag_persistent %in% TRUE), "medication-robust reactive marker-like",
    !(medication_robust_yang %in% TRUE) & is.finite(prev_vs_yin_q) & prev_vs_yin_q < .05 &
      is.finite(yang_proximity_q) & yang_proximity_q < .05 & !(lag_persistent %in% TRUE), "treatment-sensitive Yang pattern",
    is.finite(interaction_q) & interaction_q < .05, "stage-context-specific",
    default = "mixed/uncertain")]
  triage[, triage_score := pmin(-log10(pmax(yin_cox_p, 1e-300)), 30) +
    .7 * pmin(-log10(pmax(lag5_p, 1e-300)), 20) +
    .5 * pmin(-log10(pmax(interaction_p, 1e-300)), 20) -
    .3 * pmin(-log10(pmax(yang_proximity_p, 1e-300)), 20) *
      (evidence_class %in% c("medication-robust reactive marker-like", "treatment-sensitive Yang pattern"))]
  triage[!is.finite(triage_score), triage_score := NA_real_]
  setorder(triage, -triage_score, na.last = TRUE)
  mr_candidates <- triage[evidence_class %in% c("stage-specific lag-persistent association", "lag-persistent prospective association", "stage-context-specific")]
  mr_candidates[, note := "Observational triage only; external cis-pQTL/mQTL MR, colocalization, and biological validation are required"]
  class_cols <- c(
    `stage-specific lag-persistent association` = "#D95F02", `lag-persistent prospective association` = "#1B9E77",
    `medication-robust reactive marker-like` = "#E7298A", `treatment-sensitive Yang pattern` = "#CC79A7",
    `mirror/opposing` = "#7570B3", `stage-context-specific` = "#E6AB02", `mixed/uncertain` = "grey70")

  show0 <- triage[is.finite(yin_cox_beta) & is.finite(lag5_beta)]
  if (nrow(show0) >= 10) {
    qx <- quantile(show0$yin_cox_beta, c(.02, .98), na.rm = TRUE, names = FALSE)
    qy <- quantile(show0$lag5_beta, c(.02, .98), na.rm = TRUE, names = FALSE)
    show <- show0[yin_cox_beta >= qx[1] & yin_cox_beta <= qx[2] & lag5_beta >= qy[1] & lag5_beta <= qy[2]]
  } else show <- copy(show0)
  n_trim <- nrow(show0) - nrow(show)
  show <- show[order(-abs(yin_cox_beta - lag5_beta))][seq_len(min(100, .N))]
  show[, plot_label := ifelse(rank(-abs(yin_cox_beta) - abs(lag5_beta), ties.method = "first") <= 24, biomarker, "")]
  p1 <- ggplot(show, aes(yin_cox_beta, lag5_beta, color = evidence_class, label = plot_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
    geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
    geom_point(size = 1.8) + ggrepel::geom_text_repel(size = 3.0, max.overlaps = 32, seed = seed) +
    scale_color_manual(values = class_cols) +
    labs(title = "a. Persistence after a 5-year lag",
         subtitle = if (n_trim > 0) paste0("Robust plotting window; ", n_trim, " extreme outlier(s) omitted from the panel only.") else NULL,
         x = "Main prospective Cox beta", y = "Lag-5-year Cox beta", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")

  context <- triage[is.finite(retention5) & is.finite(interaction_q)]
  context[, retention_plot := pmin(retention5, 2.5)]
  context[, interaction_evidence := -log10(pmax(interaction_q, 1e-300))]
  context[, plot_label := ifelse(rank(-triage_score, ties.method = "first") <= 22, biomarker, "")]
  p2 <- ggplot(context, aes(retention_plot, interaction_evidence, color = evidence_class, label = plot_label)) +
    geom_vline(xintercept = lag_min_retention, linetype = "dashed", color = "grey55") +
    geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "grey55") +
    geom_point(size = 1.8, alpha = .85) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 28, seed = seed + 2) +
    scale_color_manual(values = class_cols) +
    labs(title = "b. Persistence and stage-context evidence",
         subtitle = "Upper-right signals combine lag retention with biomarker-by-stage interaction support.",
         x = "Lag-5/main absolute effect ratio (capped at 2.5)", y = "-log10(stage-interaction FDR)", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")

  lag_rank <- lag_res[is.finite(beta), .(lag_range = diff(range(beta, na.rm = TRUE)), absmax = max(abs(beta), na.rm = TRUE)), by = biomarker][order(-lag_range, -absmax)]
  top_lag <- lag_rank[seq_len(min(16, .N)), biomarker]
  lag_show <- lag_res[biomarker %in% top_lag]
  p3 <- ggplot(lag_show, aes(lag_year, beta, color = biomarker, group = biomarker)) +
    geom_hline(yintercept = 0, linetype = "dashed") + geom_line(linewidth = .9) + geom_point(size = 1.9) +
    labs(title = "c. Lag-response trajectories",
         subtitle = "Signals with the largest change after excluding early events; shrinking curves are more compatible with reverse causation.",
         x = "Excluded early follow-up (years)", y = "Cox beta", color = NULL) +
    theme_ckm(10) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8.5))

  yy_contrast <- triage[is.finite(yin_cox_beta) & is.finite(yang_proximity_beta10)]
  yy_contrast[, plot_label := ifelse(rank(-triage_score, ties.method = "first") <= 24, biomarker, "")]
  p4 <- ggplot(yy_contrast, aes(yin_cox_beta, yang_proximity_beta10, color = evidence_class, label = plot_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(size = 1.8, alpha = .85) + ggrepel::geom_text_repel(size = 2.8, max.overlaps = 30, seed = seed + 3) +
    scale_color_manual(values = class_cols) +
    labs(title = "d. Prospective effect vs Yang event proximity",
         subtitle = "Large proximity-only effects suggest reactive disease-state biology rather than durable prospective information.",
         x = "Prospective Cox beta", y = "Yang proximity beta per 10 years", color = NULL) +
    theme_ckm(11) + guides(color = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")

  save_plot((p1 | p2) / (p3 | p4) + plot_layout(heights = c(.92, 1.25)), biom_fig_filename("Fig5.causal_triage.png"), 17, 14.5)
  triage_thresholds <- data.table(
    criterion = c("main FDR", "lag-5 FDR", "minimum absolute log-HR", "minimum lag/main effect retention", "direction"),
    threshold = c("<0.05", "<0.05", sprintf(">= %.4f [HR >= %.2f or <= %.2f]", marker_min_cox_loghr, exp(marker_min_cox_loghr), exp(-marker_min_cox_loghr)), sprintf(">= %.2f", lag_min_retention), "same")
  )
  write_tab(triage_thresholds, "causal_triage.thresholds.tsv")
  write_tab(lag_res, "causal_triage.lagged_cox.tsv.gz")
  write_tab(triage, "causal_triage.summary.tsv.gz")
  write_tab(mr_candidates, "causal_triage.mr_coloc_candidates.tsv")
  write_book(list(key_audit = key_audit, thresholds = triage_thresholds, triage = triage,
                  lagged_cox = lag_res, mr_coloc_candidates = mr_candidates),
             biom_fig_filename("Fig5.causal_triage.out.xlsx"))
  saveRDS(list(triage = triage, lagged_cox = lag_res, mr_candidates = mr_candidates), file.path(rawdir, "causal_triage.rds"), compress = "gzip")
  progress_log("causal_triage: completed in ", elapsed_since(t_step))
  list(triage = triage, lagged_cox = lag_res, mr_candidates = mr_candidates)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Explainable Yin-Yang prediction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_cov_design <- function(dat, train_rows, test_rows, covars) {
  covars <- safe_covars(dat, covars)
  if (!length(covars)) return(list(train = matrix(nrow = length(train_rows), ncol = 0), test = matrix(nrow = length(test_rows), ncol = 0), terms = character()))
  idx <- c(train_rows, test_rows)
  z <- dat[idx, covars, drop = FALSE]
  ntr <- length(train_rows)
  for (v in covars) {
    if (is.numeric(z[[v]]) || is.integer(z[[v]])) {
      x <- safe_num(z[[v]])
      med <- median(x[seq_len(ntr)], na.rm = TRUE); if (!is.finite(med)) med <- 0
      x[!is.finite(x)] <- med
      mu <- mean(x[seq_len(ntr)]); s <- sd(x[seq_len(ntr)]); if (!is.finite(s) || s == 0) s <- 1
      z[[v]] <- (x - mu) / s
    } else {
      tr <- as.character(z[[v]][seq_len(ntr)]); tr[is.na(tr) | !nzchar(tr)] <- "Missing"
      lev <- names(sort(table(tr), decreasing = TRUE)); if (!length(lev)) lev <- "Missing"
      x <- as.character(z[[v]]); x[is.na(x) | !nzchar(x)] <- "Missing"; x[!x %in% lev] <- "Other"
      z[[v]] <- factor(x, levels = unique(c(lev, "Other")))
    }
  }
  mm <- model.matrix(~ ., data = z)
  if (ncol(mm) > 1) mm <- mm[, -1, drop = FALSE] else mm <- matrix(nrow = nrow(z), ncol = 0)
  list(train = mm[seq_len(ntr), , drop = FALSE], test = mm[ntr + seq_along(test_rows), , drop = FALSE], terms = colnames(mm))
}

make_stage_design <- function(stage, position, train_n, categorical = FALSE) {
  if (categorical) {
    f <- factor(stage, levels = 0:3)
    mm <- model.matrix(~ f)
    mm <- mm[, -1, drop = FALSE]
    colnames(mm) <- paste0("stage_cat_", 1:3)
  } else {
    x <- safe_num(position)
    mu <- mean(x[seq_len(train_n)], na.rm = TRUE); s <- sd(x[seq_len(train_n)], na.rm = TRUE); if (!is.finite(s) || s == 0) s <- 1
    mm <- matrix((x - mu) / s, ncol = 1, dimnames = list(NULL, "stage_position"))
  }
  list(train = mm[seq_len(train_n), , drop = FALSE], test = mm[train_n + seq_len(nrow(mm) - train_n), , drop = FALSE])
}

fit_cox_matrix_predict <- function(time_tr, event_tr, xtr, xte) {
  xtr <- as.matrix(xtr); xte <- as.matrix(xte)
  if (!ncol(xtr)) return(rep(0, nrow(xte)))
  colnames(xtr) <- make.names(colnames(xtr), unique = TRUE)
  colnames(xte) <- colnames(xtr)
  dd <- data.frame(time = time_tr, event = event_tr, xtr, check.names = FALSE)
  ok <- complete.cases(dd) & is.finite(dd$time) & dd$time > 0
  if (sum(ok) < 300 || sum(dd$event[ok]) < 15) {
    return(rep(NA_real_, nrow(xte)))
  }
  fit <- try(coxph(Surv(time, event) ~ ., data = dd[ok, , drop = FALSE]), silent = TRUE)
  if (inherits(fit, "try-error")) {
    warning("Cox prediction model failed: ", as.character(fit), call. = FALSE)
    return(rep(NA_real_, nrow(xte)))
  }
  pr <- try(predict(fit, newdata = data.frame(xte, check.names = FALSE), type = "lp"), silent = TRUE)
  if (inherits(pr, "try-error")) {
    warning("Cox held-out prediction failed: ", as.character(pr), call. = FALSE)
    return(rep(NA_real_, nrow(xte)))
  }
  as.numeric(pr)
}


fit_lm_matrix_predict <- function(ytr, xtr, xte) {
  xtr <- as.matrix(xtr); xte <- as.matrix(xte)
  ok_y <- is.finite(ytr)
  if (!ncol(xtr)) return(rep(mean(ytr[ok_y], na.rm = TRUE), nrow(xte)))
  colnames(xtr) <- make.names(colnames(xtr), unique = TRUE)
  colnames(xte) <- colnames(xtr)
  dd <- data.frame(y = ytr, xtr, check.names = FALSE)
  ok <- complete.cases(dd) & is.finite(dd$y)
  if (sum(ok) < 100) return(rep(NA_real_, nrow(xte)))
  fit <- try(lm(y ~ ., data = dd[ok, , drop = FALSE]), silent = TRUE)
  if (inherits(fit, "try-error")) {
    warning("Linear held-out prediction model failed: ", as.character(fit), call. = FALSE)
    return(rep(NA_real_, nrow(xte)))
  }
  pr <- try(predict(fit, newdata = data.frame(xte, check.names = FALSE)), silent = TRUE)
  if (inherits(pr, "try-error")) {
    warning("Linear held-out prediction failed: ", as.character(pr), call. = FALSE)
    return(rep(NA_real_, nrow(xte)))
  }
  as.numeric(pr)
}

weighted_kappa <- function(actual, predicted, levels = 0:4) {
  actual <- factor(actual, levels = levels)
  predicted <- factor(predicted, levels = levels)
  ok <- !is.na(actual) & !is.na(predicted)
  if (sum(ok) < 20) return(NA_real_)
  tab <- table(actual[ok], predicted[ok])
  n <- sum(tab); k <- length(levels)
  if (n <= 0 || k < 2) return(NA_real_)
  i <- matrix(rep(seq_len(k), k), nrow = k)
  j <- t(i)
  w <- 1 - ((i - j) / (k - 1))^2
  obs <- sum(w * tab) / n
  exp_tab <- outer(rowSums(tab), colSums(tab)) / n
  exp_agree <- sum(w * exp_tab) / n
  if (!is.finite(exp_agree) || exp_agree >= 1) return(NA_real_)
  (obs - exp_agree) / (1 - exp_agree)
}

classification_metrics <- function(actual, predicted, levels = 0:4) {
  actual <- factor(actual, levels = levels)
  predicted <- factor(predicted, levels = levels)
  ok <- !is.na(actual) & !is.na(predicted)
  if (sum(ok) < 20) {
    return(list(balanced_accuracy = NA_real_, macro_F1 = NA_real_))
  }
  tab <- table(actual[ok], predicted[ok])
  recall <- diag(tab) / pmax(rowSums(tab), 1)
  precision <- diag(tab) / pmax(colSums(tab), 1)
  f1 <- 2 * precision * recall / pmax(precision + recall, 1e-12)
  present <- rowSums(tab) > 0
  list(
    balanced_accuracy = mean(recall[present], na.rm = TRUE),
    macro_F1 = mean(f1[present], na.rm = TRUE)
  )
}

make_biom_matrix <- function(dat, features, train_rows, test_rows, reference_train_rows = NULL) {
  idx <- c(train_rows, test_rows)
  ntr <- length(train_rows)
  X <- matrix(NA_real_, nrow = length(idx), ncol = length(features), dimnames = list(NULL, features))

  # Winsorization and imputation are learned from the full training set. For
  # prospective YY prediction, centering/scaling can be anchored to training
  # stage-0 non-event controls so that the matrix uses the same reference as
  # the fold-specific biomarker screening and YY evidence models.
  ref_rows <- intersect(as.integer(reference_train_rows %||% integer()), as.integer(train_rows))
  ref_pos <- match(ref_rows, idx)
  ref_pos <- ref_pos[is.finite(ref_pos) & ref_pos >= 1L & ref_pos <= ntr]

  for (j in seq_along(features)) {
    x <- safe_num(dat[[features[j]]][idx])
    tr_ok <- seq_len(ntr)[is.finite(x[seq_len(ntr)])]
    if (length(tr_ok) >= 20) {
      q <- quantile(x[tr_ok], c(.001, .999), na.rm = TRUE, names = FALSE, type = 8)
      x <- pmin(pmax(x, q[1]), q[2])
    }

    scale_pos <- ref_pos[is.finite(x[ref_pos])]
    if (length(scale_pos) < 100) scale_pos <- seq_len(ntr)[is.finite(x[seq_len(ntr)])]
    mu <- mean(x[scale_pos], na.rm = TRUE)
    ss <- sd(x[scale_pos], na.rm = TRUE)
    if (!is.finite(mu)) mu <- 0
    if (!is.finite(ss) || ss == 0) ss <- 1

    med <- median(x[seq_len(ntr)], na.rm = TRUE)
    if (!is.finite(med)) med <- mu
    x[!is.finite(x)] <- med
    X[, j] <- (x - mu) / ss
  }
  list(
    train = X[seq_len(ntr), , drop = FALSE],
    test = X[ntr + seq_along(test_rows), , drop = FALSE],
    all = X
  )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Biochemistry-based refinement analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_biom_refine <- function(dat, features, meta, marker_scan, clock_obj,
                                  fig_name = biom_fig_filename("FigS2.biom_refine.png"),
                                  book_name = biom_fig_filename("FigS2.biom_refine.out.xlsx"),
                                  out_prefix = "biom_refine",
                                  display_omic_label = omic_label) {
  pos_map <- get_clock_positions(clock_obj)
  eligible <- which(dat$ckm.stage %in% 0:3 & is.finite(unname(pos_map[as.character(dat$ckm.stage)])))
  if (length(eligible) < 1000) {
    warning("Too few participants for biochemistry refinement", call. = FALSE)
    return(list(skipped = TRUE))
  }
  candidates <- if (scan_all_biom) features else unique(c(
    marker_scan[seq_len(min(100, .N)), biomarker],
    resolve_display_biomarkers(features, meta, marker_scan, n = min(30, length(features)))
  ))
  candidates <- intersect(candidates, features)
  if (length(candidates) < 5) {
    warning("Too few candidate biomarkers for biochemistry refinement", call. = FALSE)
    return(list(skipped = TRUE))
  }
  covars <- get_demog_covars(dat)
  fold_id <- make_folds(dat$ckm.stage[eligible], K_outer, seed + 31)
  rows <- list(); weights_all <- list()
  progress_log(out_prefix, ": eligible rows=", fmt_int(length(eligible)),
               "; candidate biomarkers=", fmt_int(length(candidates)), "; folds=", K_outer)

  for (fold in seq_len(K_outer)) {
    t_fold <- Sys.time()
    progress_log(out_prefix, " fold ", fold, "/", K_outer, " started")
    test_rows <- eligible[fold_id == fold]
    train_rows <- eligible[fold_id != fold]
    ytr <- unname(pos_map[as.character(dat$ckm.stage[train_rows])])
    yte <- unname(pos_map[as.character(dat$ckm.stage[test_rows])])
    covd <- make_cov_design(dat, train_rows, test_rows, covars)

    screen <- rbindlist(map_features(candidates, function(feature_id) {
      feature_id <- as.character(feature_id)[1]
      x <- safe_num(dat[[feature_id]][train_rows])
      ok0 <- is.finite(x)
      if (sum(ok0) < 200) return(data.table(biomarker = feature_id, beta = NA_real_, p = NA_real_))
      q <- quantile(x[ok0], c(.001, .999), na.rm = TRUE, names = FALSE, type = 8)
      x <- pmin(pmax(x, q[1]), q[2])
      med <- median(x[ok0], na.rm = TRUE); x[!is.finite(x)] <- med
      mu <- mean(x); ss <- sd(x); if (!is.finite(ss) || ss == 0) return(data.table(biomarker = feature_id, beta = NA_real_, p = NA_real_))
      x <- (x - mu) / ss
      dd <- data.frame(y = ytr, x = x, covd$train, check.names = FALSE)
      ok <- complete.cases(dd) & is.finite(dd$y) & is.finite(dd$x)
      if (sum(ok) < 200) return(data.table(biomarker = feature_id, beta = NA_real_, p = NA_real_))
      fit <- try(lm(y ~ ., data = dd[ok, , drop = FALSE]), silent = TRUE)
      if (inherits(fit, "try-error")) return(data.table(biomarker = feature_id, beta = NA_real_, p = NA_real_))
      co <- summary(fit)$coef
      if (!"x" %in% rownames(co)) return(data.table(biomarker = feature_id, beta = NA_real_, p = NA_real_))
      data.table(biomarker = feature_id, beta = co["x", "Estimate"], p = co["x", "Pr(>|t|)"])
    }, label = paste0(out_prefix, " fold ", fold, " screen")), fill = TRUE)
    screen <- screen[is.finite(beta) & is.finite(p)][order(p)]
    keep <- screen[seq_len(min(pred_top_n, .N))]
    progress_log(out_prefix, " fold ", fold, "/", K_outer, ": selected biomarkers=", fmt_int(nrow(keep)))
    if (!nrow(keep)) {
      stop(out_prefix, " fold ", fold, " selected no usable biomarkers", call. = FALSE)
    }
    bm <- make_biom_matrix(dat, keep$biomarker, train_rows, test_rows)
    w <- keep$beta; names(w) <- keep$biomarker
    den <- sqrt(sum(w^2)); if (!is.finite(den) || den == 0) den <- 1
    score_tr <- as.numeric(bm$train[, keep$biomarker, drop = FALSE] %*% (w / den))
    score_te <- as.numeric(bm$test[, keep$biomarker, drop = FALSE] %*% (w / den))
    pred0 <- fit_lm_matrix_predict(ytr, covd$train, covd$test)
    pred1 <- fit_lm_matrix_predict(ytr, cbind(covd$train, biomarker_stage_score = score_tr), cbind(covd$test, biomarker_stage_score = score_te))
    if (any(!is.finite(pred0)) || any(!is.finite(pred1))) {
      stop(out_prefix, " fold ", fold, " failed to produce complete held-out risk-age predictions", call. = FALSE)
    }
    nearest <- function(x) {
      stops <- unname(pos_map[as.character(0:3)])
      vapply(x, function(v) if (!is.finite(v)) NA_integer_ else (which.min(abs(stops - v)) - 1L), integer(1))
    }
    null_pos <- rep(unname(pos_map["2"]), length(test_rows))
    rows[[fold]] <- rbindlist(list(
      data.table(eid = dat$eid[test_rows], fold, model = "Dnull_always_stage2", actual_stage = dat$ckm.stage[test_rows],
                 actual_position = yte, predicted_position = null_pos, predicted_stage = 2L),
      data.table(eid = dat$eid[test_rows], fold, model = "D0_demographic", actual_stage = dat$ckm.stage[test_rows],
                 actual_position = yte, predicted_position = pred0, predicted_stage = nearest(pred0)),
      data.table(eid = dat$eid[test_rows], fold, model = "D1_explainable_biomarker_stage", actual_stage = dat$ckm.stage[test_rows],
                 actual_position = yte, predicted_position = pred1, predicted_stage = nearest(pred1))
    ))
    weights_all[[fold]] <- keep[, `:=`(fold = fold, normalized_weight = beta / den)]
    progress_log(out_prefix, " fold ", fold, "/", K_outer, " done in ", elapsed_since(t_fold))
  }

  oof <- rbindlist(rows, fill = TRUE)
  weights <- rbindlist(weights_all, fill = TRUE)
  if (!nrow(oof)) return(list(skipped = TRUE))
  expected_oof_rows <- length(eligible) * 3L
  if (nrow(oof) != expected_oof_rows ||
      uniqueN(oof[model == "D1_explainable_biomarker_stage", eid]) != length(eligible)) {
    stop(
      out_prefix, " produced incomplete out-of-fold predictions: expected ",
      expected_oof_rows, " rows across three models but obtained ", nrow(oof),
      call. = FALSE
    )
  }
  metrics <- oof[, {
    cm <- classification_metrics(actual_stage, predicted_stage, levels = 0:3)
    list(
      N = .N,
      RMSE = sqrt(mean((predicted_position - actual_position)^2, na.rm = TRUE)),
      MAE = mean(abs(predicted_position - actual_position), na.rm = TRUE),
      Spearman = suppressWarnings(cor(predicted_position, actual_position, method = "spearman", use = "complete.obs")),
      stage_accuracy = mean(predicted_stage == actual_stage, na.rm = TRUE),
      balanced_accuracy = cm$balanced_accuracy,
      macro_F1 = cm$macro_F1,
      weighted_kappa = weighted_kappa(actual_stage, predicted_stage, levels = 0:3)
    )
  }, by = .(fold, model)]
  metrics_summary <- metrics[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = model,
                             .SDcols = c("RMSE", "MAE", "Spearman", "stage_accuracy", "balanced_accuracy", "macro_F1", "weighted_kappa")]
  metric_long <- data.table::melt(as.data.table(metrics_summary), id.vars = "model", variable.name = "metric", value.name = "value")
  metric_long[, model_plot := factor(model,
    levels = c("Dnull_always_stage2", "D0_demographic", "D1_explainable_biomarker_stage"),
    labels = c("Always stage 2", "Demographic", paste0("Demographic + ", display_omic_label)))]
  metric_plot <- metric_long[metric %in% c("RMSE", "Spearman", "balanced_accuracy", "macro_F1", "weighted_kappa")]
  metric_plot[, higher_better := !metric %in% c("RMSE")]
  metric_plot[, scaled := if (.N <= 1) 1 else {
    rng <- range(value, na.rm = TRUE)
    if (!all(is.finite(rng)) || diff(rng) == 0) rep(1, .N) else if (higher_better[1]) (value - rng[1]) / diff(rng) else (rng[2] - value) / diff(rng)
  }, by = metric]
  metric_plot[, metric_label := factor(metric, levels = c("RMSE", "Spearman", "balanced_accuracy", "macro_F1", "weighted_kappa"))]

  bio <- oof[model == "D1_explainable_biomarker_stage", .(eid, actual_stage, actual_position, predicted_position)]
  bio[, stage_advancement_raw := predicted_position - actual_position]
  bio[, stage_advancement := predicted_position - mean(predicted_position, na.rm = TRUE), by = actual_stage]
  bio[, stage_advancement_z := {
    ss <- sd(stage_advancement, na.rm = TRUE)
    if (!is.finite(ss) || ss == 0) rep(NA_real_, .N) else stage_advancement / ss
  }, by = actual_stage]
  dat2 <- merge(dat, bio[, .(eid, stage_advancement_raw, stage_advancement, stage_advancement_z)], by = "eid", all.x = TRUE, sort = FALSE)
  dat2$stage_f <- factor(dat2$ckm.stage, levels = 0:3, labels = unname(stage_short_labels[as.character(0:3)]))
  advancement_cox <- cox_one(dat2[dat2$yin & dat2$ckm.stage %in% 0:3, , drop = FALSE],
                             "stage_advancement_z", unique(c("stage_f", get_clock_covars(dat2, "le8_behavior"))), min_n = 500, min_event = 20)
  advancement_cox[, stage_group := "Overall"]
  stage_specific_hr <- rbindlist(lapply(0:3, function(s0) {
    rr <- cox_one(dat2[dat2$yin & dat2$ckm.stage == s0, , drop = FALSE], "stage_advancement_z", get_clock_covars(dat2, "le8_behavior"), min_n = 250, min_event = 12)
    rr[, stage_group := unname(stage_short_labels[as.character(s0)])]
    rr
  }), fill = TRUE)
  advancement_all <- rbindlist(list(advancement_cox, stage_specific_hr), fill = TRUE)
  advancement_all[, interpretation := paste0("HR per within-stage SD of cross-fitted ", display_omic_label, " residual risk-age")]
  bio_plot <- merge(bio, dat[, c("eid", "event", "time")], by = "eid", all.x = TRUE, sort = FALSE)
  bio_plot[, event_group := factor(ifelse(event == 1, "Future CKM stage 4", "No future event"), levels = c("No future event", "Future CKM stage 4"))]

  p1 <- ggplot(oof[model == "D1_explainable_biomarker_stage"], aes(actual_position, predicted_position, color = factor(actual_stage, levels = 0:3))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = .12, size = .7) + geom_smooth(method = "lm", se = FALSE, linewidth = .7) +
    labs(title = paste0("a. Out-of-fold ", display_omic_label, " reconstruction of CKM risk-age"), x = "Assigned CKM risk-age position", y = paste0(display_omic_label, "-predicted risk-age"), color = "Stage") +
    scale_color_manual(values = stage_palette[as.character(0:3)], labels = unname(stage_short_labels[as.character(0:3)]), drop = FALSE) + theme_ckm(10)
  p2 <- ggplot(bio_plot, aes(factor(actual_stage, levels = 0:3, labels = unname(stage_short_labels[as.character(0:3)])), stage_advancement_raw, fill = event_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + geom_boxplot(outlier.alpha = .08, position = position_dodge(width = .75)) +
    labs(title = paste0("b. Raw residual ", display_omic_label, " risk-age by stage"), subtitle = "Residual = predicted risk-age minus assigned stage position.\nFuture-event cases should sit higher if residual risk is informative.", x = "Conventional CKM stage", y = "Raw residual predicted risk-age", fill = NULL) + theme_ckm(10)
  p3 <- ggplot(metric_plot, aes(scaled, metric_label, color = model_plot, group = model_plot)) +
    geom_line(linewidth = .7) + geom_point(size = 2.5) +
    geom_text(aes(label = sprintf("%.3f", value)), nudge_x = 0.035, hjust = 0, size = 2.8, show.legend = FALSE) +
    scale_x_continuous(limits = c(-0.02, 1.16), breaks = c(0, .5, 1), labels = c("worst", "mid", "best")) +
    labs(title = "c. Held-out reconstruction summary", subtitle = "Metrics are scaled within measure for display.\nText labels show the original values.", x = "Within-metric relative performance", y = NULL, color = NULL) + theme_ckm(10) + theme(legend.position = "bottom")
  p4 <- ggplot(advancement_all[is.finite(HR)], aes(HR, factor(stage_group, levels = rev(c("Overall", unname(stage_short_labels[as.character(0:3)])))), xmin = lo, xmax = hi)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") + geom_errorbar(orientation = "y", width = .14) + geom_point(size = 2.8) +
    geom_text(aes(label = sprintf("HR %.2f", HR)), nudge_x = 0.01, hjust = 0, size = 3) +
    labs(title = paste0("d. Does within-stage ", display_omic_label, " residual risk-age predict future CKM stage 4?"), x = "Hazard ratio per within-stage SD", y = NULL) + theme_ckm(10)

  save_plot(((p1 | p2) / (p3 | p4)) + plot_layout(guides = "collect") & theme(legend.position = "bottom"), fig_name, 16, 12.2)
  write_tab(oof, paste0(out_prefix, ".oof.tsv.gz"))
  write_tab(metrics, paste0(out_prefix, ".metrics_by_fold.tsv"))
  write_tab(metrics_summary, paste0(out_prefix, ".metrics_summary.tsv"))
  write_tab(weights, paste0(out_prefix, ".weights_by_fold.tsv.gz"))
  write_tab(advancement_all, paste0(out_prefix, ".advancement_cox.tsv"))
  write_book(list(metrics_by_fold = metrics, metrics_summary = metrics_summary, advancement_cox = advancement_all,
                  weights_head = head(weights[order(fold, p)], 5000)), book_name)
  saveRDS(list(oof = oof, metrics = metrics, metrics_summary = metrics_summary, weights = weights, advancement_cox = advancement_all),
          file.path(rawdir, paste0(out_prefix, ".rds")), compress = "gzip")
  list(
    metrics = metrics,
    metrics_summary = metrics_summary,
    weights = weights,
    advancement_cox = advancement_all,
    oof_file = file.path(rawdir, paste0(out_prefix, ".oof.tsv.gz"))
  )
}

fold_cox_screen <- function(dat, features, train_yin, covars, label = "prediction Cox screen") {
  rbindlist(map_features(features, function(p) {
    cols <- unique(c("time", "event", "ckm.stage", p, covars))
    d <- dat[train_yin, intersect(cols, names(dat)), drop = FALSE]
    d$.z <- standardize_ref(d[[p]], d$ckm.stage == 0 & d$event == 0)
    rr <- cox_one(d, ".z", covars, min_n = 250, min_event = 10)
    data.table(biomarker = p, cox_beta = rr$beta, cox_p = rr$p, n = rr$n, events = rr$events)
  }, label = label), fill = TRUE)
}

learn_yy_weights <- function(dat, features, train_yin, yang_rows, covars, label = "prediction") {
  progress_log(label, ": learning YY weights from ", fmt_int(length(features)), " candidate biomarkers")
  screen <- fold_cox_screen(dat, features, train_yin, covars, label = paste0(label, " Cox screen"))
  screen <- screen[is.finite(cox_p)][order(cox_p)]
  screen <- screen[seq_len(min(pred_screen_n, .N))]
  progress_log(label, ": Cox screen retained ", fmt_int(nrow(screen)), " biomarker(s) for detailed weighting")
  if (!nrow(screen)) return(data.table())
  detailed <- rbindlist(map_features(screen$biomarker, function(p) {
    cols <- unique(c("ckm.stage", "event", "yang", "yin", "time", "b2e", p, covars))
    d <- dat[, intersect(cols, names(dat)), drop = FALSE]
    ref <- rep(FALSE, nrow(d)); ref[train_yin] <- d$ckm.stage[train_yin] == 0 & d$event[train_yin] == 0
    # Learn winsorization and scaling only from the training Yin risk set;
    # Yang participants are transformed by, but do not influence, that scale.
    train_scope <- rep(FALSE, nrow(d)); train_scope[train_yin] <- TRUE
    d$.z <- standardize_ref_train(d[[p]], ref, train_scope)
    dtr <- d[train_yin, , drop = FALSE]
    c0 <- screen[biomarker == p]
    lag5 <- cox_one(dtr, ".z", covars, lag = 5, min_n = 200, min_event = 8)
    ip <- interaction_p(dtr, ".z", covars)
    ctrl_rows <- train_yin[d$event[train_yin] == 0]
    cmp_rows <- unique(c(ctrl_rows, yang_rows))
    cd <- d[cmp_rows, , drop = FALSE]; cd$prevalent_yang <- as.integer(cd$yang)
    prev <- lm_one(cd, ".z", "prevalent_yang", covars, min_n = 200)
    yd <- d[yang_rows, , drop = FALSE]; yd <- yd[is.finite(yd$b2e) & yd$b2e >= -20 & yd$b2e <= 0, , drop = FALSE]
    yd$proximity10 <- yd$b2e / 10
    prox <- lm_one(yd, ".z", "proximity10", covars, min_n = 80)
    data.table(
      biomarker = p, cox_beta = c0$cox_beta[1], cox_p = c0$cox_p[1],
      lag5_beta = lag5$beta, lag5_p = lag5$p,
      prev_beta = prev$beta, prev_p = prev$p,
      proximity_beta = prox$beta, proximity_p = prox$p,
      interaction_p = ip
    )
  }, label = paste0(label, " YY detail")), fill = TRUE)
  if (!nrow(detailed)) return(data.table())
  detailed[, `:=`(
    lag5_q = p.adjust(lag5_p, method = "BH"),
    prev_q = p.adjust(prev_p, method = "BH"),
    proximity_q = p.adjust(proximity_p, method = "BH"),
    interaction_q = p.adjust(interaction_p, method = "BH"),
    lag_retention = abs(lag5_beta) / pmax(abs(cox_beta), 1e-8)
  )]
  detailed[, `:=`(
    lag_support = is.finite(lag5_q) & lag5_q < .05 &
      is.finite(cox_beta) & abs(cox_beta) >= marker_min_cox_loghr &
      is.finite(lag5_beta) & abs(lag5_beta) >= marker_min_cox_loghr &
      sign(lag5_beta) == sign(cox_beta) & lag_retention >= lag_min_retention,
    yang_support = is.finite(prev_q) & prev_q < .05 &
      is.finite(prev_beta) & abs(prev_beta) >= marker_min_prevalence_beta &
      sign(prev_beta) == sign(cox_beta),
    stage_support = is.finite(interaction_q) & interaction_q < .05
  )]
  # data.table evaluates all RHS expressions in one := before adding new columns.
  # Compute this separately because it depends on lag_support above.
  detailed[, reactive_penalty := is.finite(proximity_q) & proximity_q < .05 &
    is.finite(proximity_beta) & abs(proximity_beta) >= marker_min_proximity_beta &
    !(lag_support %in% TRUE)]
  detailed[, multiplier := pmax(.25, pmin(2,
    1 + .35 * lag_support + .25 * yang_support + .25 * stage_support - .45 * reactive_penalty
  ))]
  detailed[, filter_score := pmin(-log10(pmax(cox_p, 1e-300)), 50) +
             3 * lag_support + 2 * yang_support + 2 * stage_support - 2 * reactive_penalty]
  detailed[, weight := cox_beta * multiplier]
  detailed[, reason := fcase(
    lag_support & stage_support, "persistent + stage-specific",
    lag_support & yang_support, "persistent + Yang-supported",
    lag_support, "persistent prospective",
    stage_support, "stage-context-specific",
    yang_support, "Yang-supported",
    reactive_penalty, "reactive down-weighted",
    default = "Yin Cox"
  )]
  detailed[, `:=`(
    yy_rank = frank(-filter_score, ties.method = "first"),
    yin_rank = frank(cox_p, ties.method = "first")
  )]
  detailed[, `:=`(
    selected_yy = yy_rank <= min(pred_top_n, .N),
    selected_yin = yin_rank <= min(pred_top_n, .N)
  )]
  detailed[selected_yy | selected_yin][order(yy_rank, yin_rank)]
}


ensure_yy_weight_schema <- function(x) {
  if (is.null(x)) return(data.table())
  x <- as.data.table(copy(x))
  if (!nrow(x)) return(x)
  logical_cols <- c("lag_support", "yang_support", "stage_support", "reactive_penalty", "selected_yy", "selected_yin")
  numeric_cols <- c("weight", "cox_beta", "normalized_weight", "yin_normalized_weight", "filter_score")
  for (v in logical_cols) {
    if (!v %in% names(x)) x[, (v) := FALSE]
    x[is.na(get(v)), (v) := FALSE]
  }
  for (v in numeric_cols) if (!v %in% names(x)) x[, (v) := NA_real_]
  if (!"reason" %in% names(x)) x[, reason := "Yin Cox"]
  if (!"biomarker" %in% names(x)) stop("YY weight table lacks biomarker column", call. = FALSE)
  x[, biomarker := as.character(biomarker)]
  x
}

score_from_weights <- function(X, weights) {
  if (is.null(weights) || !nrow(weights)) {
    stop("YY score cannot be computed because no biomarker weights were learned", call. = FALSE)
  }
  w <- weights$weight; names(w) <- weights$biomarker
  w <- w[colnames(X)]; w[!is.finite(w)] <- 0
  den <- sqrt(sum(w^2))
  if (!is.finite(den) || den == 0) {
    stop("YY score cannot be computed because all learned weights are zero/non-finite", call. = FALSE)
  }
  as.numeric(X %*% (w / den))
}


score_from_weights_or_na <- function(X, weights) {
  if (is.null(weights) || !nrow(weights)) return(rep(NA_real_, nrow(X)))
  out <- try(score_from_weights(X, weights), silent = TRUE)
  if (inherits(out, "try-error")) rep(NA_real_, nrow(X)) else out
}

safe_fit_score_model <- function(ytr, etr, cov_tr, cov_te, score_tr, score_te, score_name = "score") {
  if (!length(score_tr) || !length(score_te) || all(!is.finite(score_tr)) || all(!is.finite(score_te))) {
    return(rep(NA_real_, nrow(cov_te)))
  }
  sd_tr <- sd(score_tr, na.rm = TRUE)
  sd_te <- sd(score_te, na.rm = TRUE)
  if (!is.finite(sd_tr) || sd_tr <= 1e-8 || !is.finite(sd_te) || sd_te <= 1e-8) {
    return(rep(NA_real_, nrow(cov_te)))
  }
  fit_cox_matrix_predict(ytr, etr,
                         cbind(cov_tr, setNames(data.frame(score_tr), score_name)),
                         cbind(cov_te, setNames(data.frame(score_te), score_name)))
}

pretty_prediction_model <- function(x) {
  recode(x,
         M0_base = "Base",
         M1_categorical_stage = "Base + categorical stage",
         M2_risk_age_stage = "Base + CKM risk-age",
         M2b_explainable_Yin_score = "Base + risk-age + Yin-only score",
         M3_risk_age_stage_YY_score = "Base + risk-age + full YY",
         M3a_persistent_YY = "Base + risk-age + persistent YY",
         M3b_stage_context_YY = "Base + risk-age + stage-context YY",
         M3c_yang_supported_YY = "Base + risk-age + Yang-supported YY",
         M4_Yin_glmnet_benchmark = "Yin glmnet benchmark",
         .default = x)
}

fit_glmnet_benchmark <- function(time_tr, event_tr, xtr, xte) {
  if (!benchmark_glmnet) return(rep(NA_real_, nrow(xte)))
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    install_try <- try(install.packages("glmnet", repos = "https://cloud.r-project.org"), silent = TRUE)
    if (inherits(install_try, "try-error") || !requireNamespace("glmnet", quietly = TRUE)) {
      warning("glmnet benchmark requested but the glmnet package is unavailable; benchmark values will be missing", call. = FALSE)
      return(rep(NA_real_, nrow(xte)))
    }
  }
  y <- Surv(time_tr, event_tr)
  foldid <- make_folds(event_tr, min(5, K_outer), seed + 100)
  fit <- try(
    glmnet::cv.glmnet(
      as.matrix(xtr), y, family = "cox", alpha = .5,
      foldid = foldid, standardize = FALSE
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    warning("glmnet benchmark failed; benchmark values will be missing: ", as.character(fit), call. = FALSE)
    return(rep(NA_real_, nrow(xte)))
  }
  pred <- try(predict(fit, newx = as.matrix(xte), s = fit$lambda.1se, type = "link"), silent = TRUE)
  if (inherits(pred, "try-error")) {
    warning("glmnet benchmark prediction failed; benchmark values will be missing: ", as.character(pred), call. = FALSE)
    return(rep(NA_real_, nrow(xte)))
  }
  as.numeric(pred)
}

run_prediction <- function(dat, dat_clock, features, meta, clock_obj) {
  t_pred <- Sys.time()
  eligible <- which(dat$yin & dat$ckm.stage %in% 0:3 & is.finite(dat$time) & dat$time > 0)
  events <- sum(dat$event[eligible])
  cohort <- data.table(trait, biom, N = length(eligible), events, Yang_N = sum(dat$yang), features = length(features),
					   base_mode = pred_base_mode,
                       clock_adjust_mode = clock_adjust_mode,
                       stage_clock_model = paste0("prospective_risk_age_", clock_adjust_mode))
  write_tab(cohort, "prediction.cohort.tsv")
  if (length(eligible) < pred_min_yin_n || events < pred_min_events) {
    msg <- paste0("Prediction skipped: Yin N/events = ", length(eligible), "/", events,
                  "; thresholds = ", pred_min_yin_n, "/", pred_min_events)
    warning(msg, call. = FALSE); writeLines(msg, file.path(rawdir, "prediction.SKIPPED.txt"))
    return(list(skipped = TRUE, cohort = cohort))
  }
  base_covars <- get_base_covars(dat)
  clock_covars <- get_clock_covars(dat_clock, clock_adjust_mode)
  fold_strata <- interaction(
    dat$ckm.stage[eligible],
    dat$event[eligible],
    drop = TRUE,
    lex.order = TRUE
  )
  folds <- make_folds(fold_strata, K_outer, seed)
  yang_rows <- which(dat$yang)
  pred_rows <- list(); weight_rows <- list(); fold_audit_rows <- list(); explain_rows <- list()
  progress_log("prediction: eligible Yin rows=", fmt_int(length(eligible)), "; events=", fmt_int(events),
               "; features=", fmt_int(length(features)), "; folds=", K_outer)

  for (fold in seq_len(K_outer)) {
    t_fold <- Sys.time()
    progress_log("Prediction fold ", fold, "/", K_outer, " started")
    test_rows <- eligible[folds == fold]
    train_rows <- eligible[folds != fold]
    test_eids <- dat$eid[test_rows]
    clock_train <- dat_clock[!dat_clock$eid %in% test_eids, , drop = FALSE]
    fold_road <- try(
      make_risk_age_stage_clock(
        clock_train,
        safe_covars(clock_train, clock_covars),
        reference_age = risk_age_reference_age
      ),
      silent = TRUE
    )
    if (inherits(fold_road, "try-error")) {
      stop(
        "Fold-specific CKM risk-age road failed in prediction fold ", fold,
        ". The pipeline will not mix risk-age and ordinal stage coordinates. Error: ",
        as.character(fold_road),
        call. = FALSE
      )
    }
    pos_map <- setNames(fold_road$road$risk_age_years, fold_road$road$stage)
    fold_clock_source <- "outer-training CKM risk-age road"
    stage_position <- unname(pos_map[as.character(dat$ckm.stage)])
    covd <- make_cov_design(dat, train_rows, test_rows, base_covars)
    stage_cat <- make_stage_design(c(dat$ckm.stage[train_rows], dat$ckm.stage[test_rows]), c(stage_position[train_rows], stage_position[test_rows]), length(train_rows), TRUE)
    stage_quant <- make_stage_design(c(dat$ckm.stage[train_rows], dat$ckm.stage[test_rows]), c(stage_position[train_rows], stage_position[test_rows]), length(train_rows), FALSE)
    ytr <- dat$time[train_rows]; etr <- dat$event[train_rows]
    yte <- dat$time[test_rows]; ete <- dat$event[test_rows]

    lp0 <- fit_cox_matrix_predict(ytr, etr, covd$train, covd$test)
    lp1 <- fit_cox_matrix_predict(ytr, etr, cbind(covd$train, stage_cat$train), cbind(covd$test, stage_cat$test))
    lp2 <- fit_cox_matrix_predict(ytr, etr, cbind(covd$train, stage_quant$train), cbind(covd$test, stage_quant$test))
    base_models <- list(M0_base = lp0, M1_categorical_stage = lp1, M2_risk_age_stage = lp2)
    for (nm in names(base_models)) {
      zz <- base_models[[nm]]
      if (length(zz) != length(test_rows) || any(!is.finite(zz))) {
        stop("Prediction fold ", fold, " failed to produce complete held-out predictions for ", nm, call. = FALSE)
      }
    }

    weights <- ensure_yy_weight_schema(learn_yy_weights(dat, features, train_rows, yang_rows, base_covars, label = paste0("prediction fold ", fold)))
    if (is.null(weights) || !nrow(weights)) {
      stop("Prediction fold ", fold, " learned no YY biomarker weights. This is a pipeline failure, not evidence of a null YY effect.", call. = FALSE)
    }
    yy_weights <- copy(weights[selected_yy %in% TRUE])
    yin_weights <- copy(weights[selected_yin %in% TRUE])
    if (!nrow(yy_weights) || !nrow(yin_weights)) {
      stop("Prediction fold ", fold, " did not retain both YY and Yin-only biomarker sets", call. = FALSE)
    }
    yy_norm <- sqrt(sum(yy_weights$weight^2, na.rm = TRUE))
    yin_norm <- sqrt(sum(yin_weights$cox_beta^2, na.rm = TRUE))
    if (!is.finite(yy_norm) || yy_norm <= 0 || !is.finite(yin_norm) || yin_norm <= 0) {
      stop("Prediction fold ", fold, " produced zero/non-finite explainable score weights", call. = FALSE)
    }
    weights[, `:=`(
      normalized_weight = NA_real_,
      yin_weight = fifelse(selected_yin %in% TRUE, cox_beta, NA_real_),
      yin_normalized_weight = NA_real_,
      fold = fold,
      stage_clock_source = fold_clock_source
    )]
    weights[selected_yy %in% TRUE, normalized_weight := weight / yy_norm]
    weights[selected_yin %in% TRUE, yin_normalized_weight := cox_beta / yin_norm]
    yy_weights <- copy(weights[selected_yy %in% TRUE])
    yin_weights <- copy(weights[selected_yin %in% TRUE])
    yin_weights[, weight := cox_beta]
    weight_rows[[fold]] <- weights

    selected <- unique(c(as.character(yy_weights$biomarker), as.character(yin_weights$biomarker)))
    selected <- intersect(selected, features)
    if (!length(selected)) stop("Prediction fold ", fold, " has no usable selected biomarkers", call. = FALSE)
    prediction_ref_rows <- train_rows[
      dat$ckm.stage[train_rows] == 0 & dat$event[train_rows] == 0
    ]
    bm <- make_biom_matrix(
      dat, selected, train_rows, test_rows,
      reference_train_rows = prediction_ref_rows
    )

    score_tr <- score_from_weights(bm$train, yy_weights)
    score_te <- score_from_weights(bm$test, yy_weights)
    lp3 <- safe_fit_score_model(
      ytr, etr,
      cbind(covd$train, stage_quant$train),
      cbind(covd$test, stage_quant$test),
      score_tr, score_te, "YY_score"
    )

    yin_score_tr <- score_from_weights(bm$train, yin_weights)
    yin_score_te <- score_from_weights(bm$test, yin_weights)
    lp_yin <- safe_fit_score_model(
      ytr, etr,
      cbind(covd$train, stage_quant$train),
      cbind(covd$test, stage_quant$test),
      yin_score_tr, yin_score_te, "Yin_score"
    )

    explain_rows[[length(explain_rows) + 1L]] <- data.table(
      fold = fold,
      model = "M2b_explainable_Yin_score",
      biomarkers = nrow(yin_weights),
      abs_weight = sum(abs(yin_weights$yin_normalized_weight), na.rm = TRUE),
      reasons = "Yin Cox ranking and weighting"
    )
    explain_rows[[length(explain_rows) + 1L]] <- data.table(
      fold = fold,
      model = "M3_risk_age_stage_YY_score",
      biomarkers = nrow(yy_weights),
      abs_weight = sum(abs(yy_weights$normalized_weight), na.rm = TRUE),
      reasons = paste(sort(unique(yy_weights$reason)), collapse = "; ")
    )

    component_defs <- list(
      M3a_persistent_YY = yy_weights[lag_support %in% TRUE],
      M3b_stage_context_YY = yy_weights[stage_support %in% TRUE],
      M3c_yang_supported_YY = yy_weights[yang_support %in% TRUE]
    )
    component_models <- lapply(names(component_defs), function(m) {
      subw <- component_defs[[m]]
      if (!nrow(subw)) {
        explain_rows[[length(explain_rows) + 1L]] <<- data.table(
          fold = fold, model = m, biomarkers = 0L, abs_weight = 0,
          reasons = "No biomarkers met this evidence channel in the training fold"
        )
        # Optional explanation channels must not terminate the primary model.
        return(lp2)
      }
      subtr <- score_from_weights_or_na(bm$train, subw)
      subte <- score_from_weights_or_na(bm$test, subw)
      if (all(!is.finite(subtr)) || all(!is.finite(subte)) ||
          sd(subtr, na.rm = TRUE) <= 1e-8 || sd(subte, na.rm = TRUE) <= 1e-8) {
        lp <- lp2
        reason0 <- "Evidence channel produced a constant/invalid score; risk-age model used as neutral fallback"
      } else {
        lp0 <- safe_fit_score_model(
          ytr, etr,
          cbind(covd$train, stage_quant$train),
          cbind(covd$test, stage_quant$test),
          subtr, subte, paste0(m, "_score")
        )
        lp <- if (all(!is.finite(lp0)) || sd(lp0, na.rm = TRUE) <= 1e-10) lp2 else lp0
        reason0 <- paste(sort(unique(subw$reason)), collapse = "; ")
      }
      explain_rows[[length(explain_rows) + 1L]] <<- data.table(
        fold = fold, model = m, biomarkers = nrow(subw),
        abs_weight = sum(abs(subw$normalized_weight), na.rm = TRUE),
        reasons = reason0
      )
      lp
    })
    names(component_models) <- names(component_defs)

    score_sd_tr <- sd(score_tr, na.rm = TRUE)
    score_sd_te <- sd(score_te, na.rm = TRUE)
    yin_score_sd_tr <- sd(yin_score_tr, na.rm = TRUE)
    yin_score_sd_te <- sd(yin_score_te, na.rm = TRUE)
    lp3_sd <- sd(lp3, na.rm = TRUE)
    lp_yin_sd <- sd(lp_yin, na.rm = TRUE)
    lp_delta_sd <- sd(lp3 - lp2, na.rm = TRUE)
    if (!is.finite(score_sd_tr) || score_sd_tr <= 1e-8 ||
        !is.finite(score_sd_te) || score_sd_te <= 1e-8) {
      stop("Prediction fold ", fold, " produced a constant/invalid YY score", call. = FALSE)
    }
    if (!is.finite(yin_score_sd_tr) || yin_score_sd_tr <= 1e-8 ||
        !is.finite(yin_score_sd_te) || yin_score_sd_te <= 1e-8) {
      stop("Prediction fold ", fold, " produced a constant/invalid Yin-only score", call. = FALSE)
    }
    if (!is.finite(lp3_sd) || lp3_sd <= 1e-10) {
      stop("Prediction fold ", fold, " failed to fit a non-constant full YY model", call. = FALSE)
    }
    if (!is.finite(lp_yin_sd) || lp_yin_sd <= 1e-10) {
      stop("Prediction fold ", fold, " failed to fit a non-constant explainable Yin-only model", call. = FALSE)
    }
    if (!is.finite(lp_delta_sd) || lp_delta_sd <= 1e-10) {
      stop("Prediction fold ", fold, " produced a full YY predictor identical to CKM risk-age", call. = FALSE)
    }
    fold_audit_rows[[fold]] <- data.table(
      fold = fold, train_N = length(train_rows), test_N = length(test_rows),
      train_events = sum(etr), test_events = sum(ete),
      selected_union = length(selected),
      selected_yy = nrow(yy_weights),
      selected_yin = nrow(yin_weights),
      selected_overlap = sum(weights$selected_yy %in% TRUE & weights$selected_yin %in% TRUE),
      yy_score_sd_train = score_sd_tr, yy_score_sd_test = score_sd_te,
      yin_score_sd_train = yin_score_sd_tr, yin_score_sd_test = yin_score_sd_te,
      lp3_sd = lp3_sd, lp_yin_sd = lp_yin_sd, lp3_minus_lp2_sd = lp_delta_sd,
      stage_clock_source = fold_clock_source
    )

    models <- c(
      list(
        M0_base = lp0,
        M1_categorical_stage = lp1,
        M2_risk_age_stage = lp2,
        M2b_explainable_Yin_score = lp_yin,
        M3_risk_age_stage_YY_score = lp3
      ),
      component_models
    )
    if (benchmark_glmnet && nrow(yin_weights)) {
      # Keep the optional elastic-net benchmark strictly Yin-only: its omic
      # inputs are the fold-specific biomarkers selected by the Yin Cox ranking,
      # never features added only through Yang/stage evidence.
      yin_glmnet_features <- intersect(as.character(yin_weights$biomarker), colnames(bm$train))
      xtr_g <- cbind(covd$train, stage_quant$train, bm$train[, yin_glmnet_features, drop = FALSE])
      xte_g <- cbind(covd$test, stage_quant$test, bm$test[, yin_glmnet_features, drop = FALSE])
      models$M4_Yin_glmnet_benchmark <- fit_glmnet_benchmark(ytr, etr, xtr_g, xte_g)
    }
    pred_rows[[fold]] <- rbindlist(lapply(names(models), function(m) data.table(
      eid = dat$eid[test_rows], fold, model = m, time = yte, event = ete, lp = models[[m]],
      stage_clock_source = fold_clock_source
    )), fill = TRUE)
    progress_log("Prediction fold ", fold, "/", K_outer, " done in ", elapsed_since(t_fold))
  }

  pred <- rbindlist(pred_rows, fill = TRUE)
  weights <- rbindlist(weight_rows, fill = TRUE)
  fold_audit <- rbindlist(fold_audit_rows, fill = TRUE)
  explain_folds <- rbindlist(explain_rows, fill = TRUE)
  core_prediction_models <- c(
    "M0_base", "M1_categorical_stage", "M2_risk_age_stage",
    "M2b_explainable_Yin_score", "M3_risk_age_stage_YY_score"
  )
  for (nm in core_prediction_models) {
    zz <- pred[model == nm]
    if (nrow(zz) != length(eligible) || uniqueN(zz$eid) != length(eligible) || any(!is.finite(zz$lp))) {
      stop("Incomplete out-of-fold predictions for core model ", nm, call. = FALSE)
    }
  }
  perf <- pred[is.finite(lp), .(
    n = .N, events = sum(event),
    cindex = calc_cindex(time, event, lp),
    calibration_slope = calc_calibration_slope(time, event, lp),
    AUC10 = { known <- time >= 10 | event == 1; calc_auc(as.integer(event[known] == 1 & time[known] <= 10), lp[known]) }
  ), by = .(fold, model)]
  base <- perf[model == "M0_base", .(fold, base_c = cindex, base_auc = AUC10)]
  quant <- perf[model == "M2_risk_age_stage", .(fold, quant_c = cindex, quant_auc = AUC10)]
  yinref <- perf[model == "M2b_explainable_Yin_score", .(fold, yin_c = cindex, yin_auc = AUC10)]
  perf <- merge(perf, base, by = "fold", all.x = TRUE)
  perf <- merge(perf, quant, by = "fold", all.x = TRUE)
  perf <- merge(perf, yinref, by = "fold", all.x = TRUE)
  perf[, `:=`(
    delta_vs_base = cindex - base_c,
    delta_vs_quant = cindex - quant_c,
    delta_vs_yin = cindex - yin_c,
    delta_auc_vs_base = AUC10 - base_auc,
    delta_auc_vs_quant = AUC10 - quant_auc,
    delta_auc_vs_yin = AUC10 - yin_auc
  )]
  perf_sum <- perf[, .(
    mean_cindex = mean(cindex, na.rm = TRUE),
    sd_cindex = sd(cindex, na.rm = TRUE),
    mean_delta_vs_base = mean(delta_vs_base, na.rm = TRUE),
    mean_delta_vs_quant = mean(delta_vs_quant, na.rm = TRUE),
    mean_delta_vs_yin = mean(delta_vs_yin, na.rm = TRUE),
    mean_AUC10 = mean(AUC10, na.rm = TRUE),
    mean_delta_auc_vs_yin = mean(delta_auc_vs_yin, na.rm = TRUE),
    mean_calibration_slope = mean(calibration_slope, na.rm = TRUE)
  ), by = model][order(-mean_cindex)]

  weights[, selection_class := fcase(
    selected_yy %in% TRUE & selected_yin %in% TRUE, "shared",
    selected_yy %in% TRUE, "YY-only",
    selected_yin %in% TRUE, "Yin-only",
    default = "not selected"
  )]
  weight_summary <- weights[, .(
    yy_select_freq = sum(selected_yy %in% TRUE),
    yin_select_freq = sum(selected_yin %in% TRUE),
    shared_freq = sum(selected_yy %in% TRUE & selected_yin %in% TRUE),
    mean_weight = mean(weight[selected_yy %in% TRUE], na.rm = TRUE),
    mean_abs_weight = mean(abs(weight[selected_yy %in% TRUE]), na.rm = TRUE),
    mean_norm_weight = mean(normalized_weight, na.rm = TRUE),
    abs_norm_weight = mean(abs(normalized_weight), na.rm = TRUE),
    mean_yin_norm_weight = mean(yin_normalized_weight, na.rm = TRUE),
    lag_support_rate = mean(lag_support[selected_yy %in% TRUE] %in% TRUE),
    stage_support_rate = mean(stage_support[selected_yy %in% TRUE] %in% TRUE),
    yang_support_rate = mean(yang_support[selected_yy %in% TRUE] %in% TRUE),
    reactive_penalty_rate = mean(reactive_penalty[selected_yy %in% TRUE] %in% TRUE),
    dominant_reason = {
      z <- reason[selected_yy %in% TRUE]
      z <- z[!is.na(z) & nzchar(z)]
      if (length(z)) names(sort(table(z), decreasing = TRUE))[1] else NA_character_
    }
  ), by = biomarker]
  for (v in c("mean_weight", "mean_abs_weight", "mean_norm_weight", "abs_norm_weight", "mean_yin_norm_weight",
              "lag_support_rate", "stage_support_rate", "yang_support_rate", "reactive_penalty_rate")) {
    set(weight_summary, which(!is.finite(weight_summary[[v]])), v, NA_real_)
  }
  meta0 <- as.data.table(copy(meta))
  meta_keep <- intersect(c("biomarker", "label", "category", "subgroup"), names(meta0))
  if ("biomarker" %in% meta_keep) {
    weight_summary <- merge(weight_summary, unique(meta0[, ..meta_keep], by = "biomarker"),
                            by = "biomarker", all.x = TRUE, sort = FALSE)
  }
  if (!"label" %in% names(weight_summary)) weight_summary[, label := biomarker]
  if (!"category" %in% names(weight_summary)) weight_summary[, category := NA_character_]
  if (!"subgroup" %in% names(weight_summary)) weight_summary[, subgroup := NA_character_]
  weight_summary[is.na(label), label := biomarker]
  top_weights <- head(weight_summary[yy_select_freq > 0][order(-yy_select_freq, -abs_norm_weight)], yy_top_show)
  top_weights[, label_plot := factor(label, levels = rev(unique(label)))]

  reason_summary <- weights[selected_yy %in% TRUE, .(
    abs_weight = sum(abs(normalized_weight), na.rm = TRUE),
    selections = .N
  ), by = reason][order(-abs_weight)]
  selection_counts <- weights[, .N, by = .(fold, selection_class)]
  selection_counts[, proportion := N / sum(N), by = fold]

  channel_perf <- perf[model %in% c(
    "M2b_explainable_Yin_score", "M3_risk_age_stage_YY_score",
    "M3a_persistent_YY", "M3b_stage_context_YY", "M3c_yang_supported_YY"
  )]
  channel_perf[, model_plot := factor(
    pretty_prediction_model(model),
    levels = pretty_prediction_model(c(
      "M2b_explainable_Yin_score", "M3_risk_age_stage_YY_score",
      "M3a_persistent_YY", "M3b_stage_context_YY", "M3c_yang_supported_YY"
    ))
  )]
  model_order <- c(
    "M0_base", "M1_categorical_stage", "M2_risk_age_stage",
    "M2b_explainable_Yin_score", "M3_risk_age_stage_YY_score", "M3a_persistent_YY",
    "M3b_stage_context_YY", "M3c_yang_supported_YY",
    "M4_Yin_glmnet_benchmark"
  )
  perf[, model_plot := factor(
    pretty_prediction_model(model),
    levels = pretty_prediction_model(model_order)
  )]
  perf_sum[, model_plot := pretty_prediction_model(model)]

  p1 <- ggplot(perf[model %in% c(
    "M0_base", "M1_categorical_stage", "M2_risk_age_stage",
    "M2b_explainable_Yin_score", "M3_risk_age_stage_YY_score",
    "M4_Yin_glmnet_benchmark"
  )],
               aes(cindex, fct_rev(model_plot), fill = model)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(height = .08), size = 1.5, alpha = .8) +
    labs(
      title = "a. Held-out Yin prediction",
      subtitle = "Full YY is compared with CKM stage models and an explainable Yin-only score.\nThe glmnet benchmark is shown only when requested.",
      x = "C-index", y = NULL, fill = NULL
    ) +
    scale_fill_manual(values = model_palette, guide = "none") + theme_ckm(10)
  p2 <- ggplot(channel_perf, aes(delta_vs_quant, fct_rev(model_plot), fill = model)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(height = .06), size = 1.6) +
    labs(
      title = "b. What information does YY add beyond CKM risk-age?",
      subtitle = "Each component is added separately to CKM risk-age.\nComponents may overlap and are not additive.",
      x = "Delta C-index vs CKM risk-age", y = NULL
    ) +
    scale_fill_manual(values = model_palette, guide = "none") + theme_ckm(10)
  p3 <- ggplot(top_weights, aes(mean_norm_weight, label_plot, color = dominant_reason, size = yy_select_freq)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    geom_segment(aes(x = 0, xend = mean_norm_weight, y = label_plot, yend = label_plot), color = "grey75") +
    geom_point() +
    labs(title = "c. Variables entering the explainable YY score", subtitle = "Point size: fold-selection frequency.\nPosition: average signed normalized weight.", x = "Average normalized YY weight", y = NULL, color = "Dominant reason", size = "Fold count") +
    theme_ckm(10)
  p4 <- ggplot(reason_summary, aes(abs_weight, fct_reorder(reason, abs_weight), fill = reason)) +
    geom_col(color = "white") + geom_text(aes(label = sprintf("%.1f", abs_weight)), hjust = -0.1, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, .2))) +
    labs(title = "d. Information decomposition of the YY score", subtitle = "Total absolute normalized weight by explanation class across folds.", x = "Total absolute normalized weight", y = NULL, fill = NULL) + theme_ckm(10) + theme(legend.position = "none")
  save_plot((p1 | p2) / (p3 | p4), biom_fig_filename("Fig6.prediction.png"), 16, 12.5)

  full_vs_yin <- perf[model == "M3_risk_age_stage_YY_score", .(
    fold,
    delta_cindex_full_vs_yin = delta_vs_yin,
    delta_auc10_full_vs_yin = delta_auc_vs_yin
  )]
  selection_total <- weights[, .N, by = selection_class][order(-N)]
  selection_total[, proportion := N / sum(N)]
  reweight_rows <- weights[
    selected_yy %in% TRUE & selected_yin %in% TRUE,
    .(
      fold, biomarker, normalized_weight, yin_normalized_weight,
      delta_normalized_weight = normalized_weight - yin_normalized_weight,
      reason
    )
  ]
  reweight_summary <- reweight_rows[, .(
    folds = .N,
    mean_yy_weight = mean(normalized_weight, na.rm = TRUE),
    mean_yin_weight = mean(yin_normalized_weight, na.rm = TRUE),
    mean_delta_weight = mean(delta_normalized_weight, na.rm = TRUE)
  ), by = biomarker]
  reweight_summary <- merge(
    reweight_summary,
    weight_summary[, .(biomarker, label, category, dominant_reason)],
    by = "biomarker", all.x = TRUE, sort = FALSE
  )
  reweight_summary[is.na(label), label := biomarker]
  reweight_top <- head(reweight_summary[order(-abs(mean_delta_weight))], 15)
  reweight_top[, label_plot := factor(label, levels = rev(unique(label)))]

  yy_only <- weights[selection_class == "YY-only", .(
    folds = .N,
    mean_weight = mean(normalized_weight, na.rm = TRUE),
    dominant_reason = {
      z <- reason[!is.na(reason)]
      if (length(z)) names(sort(table(z), decreasing = TRUE))[1] else NA_character_
    }
  ), by = biomarker]
  yy_only <- merge(
    yy_only,
    weight_summary[, .(biomarker, label, category)],
    by = "biomarker", all.x = TRUE, sort = FALSE
  )
  yy_only[is.na(label), label := biomarker]
  yy_only_top <- head(yy_only[order(-folds, -abs(mean_weight))], 15)
  if (nrow(yy_only_top)) yy_only_top[, label_plot := factor(label, levels = rev(unique(label)))]

  pS1 <- ggplot(full_vs_yin, aes(factor(fold), delta_cindex_full_vs_yin)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_col(fill = unname(model_palette["M3_risk_age_stage_YY_score"])) +
    geom_text(aes(
      label = sprintf("%+.4f", delta_cindex_full_vs_yin),
      vjust = ifelse(delta_cindex_full_vs_yin >= 0, -0.25, 1.2)
    ), size = 3) +
    scale_y_continuous(expand = expansion(mult = c(.12, .18))) +
    labs(
      title = "a. Increment of full YY over the explainable Yin-only score",
      x = "Outer fold", y = "Delta C-index: full YY minus Yin-only"
    ) + theme_ckm(10)
  pS2 <- ggplot(selection_total, aes(proportion, fct_reorder(selection_class, proportion), fill = selection_class)) +
    geom_col(color = "white") +
    geom_text(aes(label = paste0(N, " (", scales::percent(proportion, accuracy = .1), ")")), hjust = -0.08, size = 3) +
    scale_x_continuous(labels = percent, expand = expansion(mult = c(0, .22))) +
    labs(
      title = "b. Feature-selection overlap between Yin-only and full YY",
      x = "Proportion of fold-specific selected features", y = NULL, fill = NULL
    ) + theme_ckm(10) + theme(legend.position = "none")
  pS3 <- if (nrow(reweight_top)) {
    ggplot(reweight_top, aes(mean_delta_weight, label_plot, color = dominant_reason)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
      geom_segment(aes(x = 0, xend = mean_delta_weight, y = label_plot, yend = label_plot), color = "grey75") +
      geom_point(size = 2.8) +
      labs(
        title = "c. Shared variables most strongly reweighted by YY evidence",
        x = "Mean normalized weight change: full YY minus Yin-only", y = NULL, color = "YY evidence"
      ) + theme_ckm(10)
  } else ggplot() + annotate("text", x = 0, y = 0, label = "No shared selected features") + theme_void()
  pS4 <- if (nrow(yy_only_top)) {
    ggplot(yy_only_top, aes(mean_weight, label_plot, color = dominant_reason, size = folds)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
      geom_segment(aes(x = 0, xend = mean_weight, y = label_plot, yend = label_plot), color = "grey75") +
      geom_point() +
      labs(
        title = "d. Variables selected by full YY but not by Yin-only ranking",
        x = "Mean normalized YY weight", y = NULL, color = "YY evidence", size = "Fold count"
      ) + theme_ckm(10)
  } else ggplot() + annotate("text", x = 0, y = 0, label = "The selected feature sets were identical") + theme_void()
  save_plot(
    (pS1 | pS2) / (pS3 | pS4),
    biom_fig_filename("FigS4.YY_increment.png"),
    16, 12
  )

  yy_thresholds <- data.table(
    component = c("persistent YY", "Yang-supported YY", "stage-context YY", "reactive penalty"),
    definition = c(
      paste0("within-screen lag-5 FDR<0.05; main and lag-5 |log-HR|>=", sprintf("%.4f", marker_min_cox_loghr),
             "; same direction; lag/main retention>=", sprintf("%.2f", lag_min_retention)),
      paste0("within-screen Yang-vs-Yin FDR<0.05; |beta|>=", marker_min_prevalence_beta, "; same direction as Yin Cox beta"),
      "within-screen biomarker-by-stage interaction FDR<0.05",
      paste0("within-screen event-proximity FDR<0.05; |beta|>=", marker_min_proximity_beta, "; no persistent support")
    )
  )
  write_tab(yy_thresholds, "prediction.yy_explanation_thresholds.tsv")
  write_tab(full_vs_yin, "prediction.full_yy_vs_yin.tsv")
  write_tab(selection_total, "prediction.selection_overlap.tsv")
  write_tab(reweight_summary, "prediction.yy_reweighting.tsv")
  write_tab(yy_only, "prediction.yy_only_features.tsv")
  write_tab(pred, "prediction.rows.tsv.gz")
  write_tab(perf, "prediction.performance_by_fold.tsv")
  write_tab(perf_sum, "prediction.performance_summary.tsv")
  write_tab(weights, "prediction.yy_weights_by_fold.tsv.gz")
  write_tab(weight_summary, "prediction.yy_weight_summary.tsv")
  write_tab(reason_summary, "prediction.yy_reason_summary.tsv")
  write_tab(explain_folds, "prediction.yy_component_folds.tsv")
  write_tab(fold_audit, "prediction.fold_audit.tsv")
  write_book(list(
    cohort = cohort,
    fold_audit = fold_audit,
    performance_by_fold = perf,
    performance_summary = perf_sum,
    yy_explanation_thresholds = yy_thresholds,
    yy_component_folds = explain_folds,
    full_yy_vs_yin = full_vs_yin,
    selection_overlap = selection_total,
    yy_reweighting = reweight_summary,
    yy_only_features = yy_only,
    yy_weights_head = head(weights[order(fold, -filter_score)], 5000),
    yy_weight_summary = weight_summary,
    yy_reason_summary = reason_summary
  ), biom_fig_filename("Fig6.prediction.out.xlsx"))
  saveRDS(list(
    cohort = cohort,
    fold_audit = fold_audit,
    pred = pred,
    perf = perf,
    perf_sum = perf_sum,
    weights = weights,
    weight_summary = weight_summary,
    reason_summary = reason_summary,
    full_vs_yin = full_vs_yin,
    selection_overlap = selection_total,
    reweighting = reweight_summary,
    yy_only_features = yy_only
  ), file.path(rawdir, "prediction.rds"), compress = "gzip")
  progress_log("prediction: completed in ", elapsed_since(t_pred))
  list(
    cohort = cohort,
    fold_audit = fold_audit,
    perf = perf,
    perf_sum = perf_sum,
    weight_summary = weight_summary,
    reason_summary = reason_summary,
    prediction_rds = file.path(rawdir, "prediction.rds")
  )
}

run_met_strict_sensitivity <- function(dat, features, meta, marker_scan, clock_obj, stage_ref_full = NULL) {
  if (biom != "met" || !met_strict_exclude_lipid_block) return(invisible(NULL))
  meta2 <- as.data.table(copy(meta))
  if (!"strict_exclude_lipoprotein_block" %in% names(meta2)) return(invisible(NULL))
  strict_features <- unique(meta2[
    biomarker %in% features & !(strict_exclude_lipoprotein_block %in% TRUE),
    biomarker
  ])
  excluded <- unique(meta2[
    biomarker %in% features & strict_exclude_lipoprotein_block %in% TRUE
  ], by = "biomarker")
  if (nrow(excluded) < 20) {
    stop(
      "The strict metabolite sensitivity identified only ", nrow(excluded),
      " lipoprotein-lipid features. Verify CKM_MET_META and the met_name/group/subgroup columns in met.lst.",
      call. = FALSE
    )
  }
  if (length(strict_features) < 10) {
    warning("Met strict sensitivity skipped because too few non-lipoprotein features remain", call. = FALSE)
    return(invisible(NULL))
  }
  progress_log("met strict sensitivity: retained ", fmt_int(length(strict_features)), " non-lipoprotein features; excluded ", fmt_int(nrow(excluded)), " lipoprotein-lipid features")
  write_tab(excluded, "met_strict.excluded_lipoprotein_lipid_features.tsv")
  write_tab(meta2[biomarker %in% strict_features], "met_strict.retained_nonlipoprotein_features.tsv")
  strict_res <- make_biom_refine(dat, strict_features, meta2[biomarker %in% strict_features], marker_scan, clock_obj,
                                      fig_name = biom_fig_filename("FigS2.strict_biom_refine.png"),
                                      book_name = biom_fig_filename("FigS2.strict_biom_refine.out.xlsx"),
                                      out_prefix = "met_strict_biom_refine",
                                      display_omic_label = "non-lipoprotein metabolomic")
  if (is.null(strict_res) || isTRUE(strict_res$skipped)) return(invisible(NULL))

  all_hr <- if (!is.null(stage_ref_full) && !isTRUE(stage_ref_full$skipped)) as.data.table(copy(stage_ref_full$advancement_cox)) else data.table()
  strict_hr <- as.data.table(strict_res$advancement_cox)
  if (nrow(all_hr)) all_hr[, analysis := "All metabolites"]
  if (nrow(strict_hr)) strict_hr[, analysis := "Non-lipoprotein metabolites"]
  compare_hr <- rbindlist(list(all_hr, strict_hr), fill = TRUE)
  compare_metrics <- rbindlist(list(
    if (!is.null(stage_ref_full) && !isTRUE(stage_ref_full$skipped)) as.data.table(copy(stage_ref_full$metrics_summary))[, analysis := "All metabolites"] else data.table(),
    as.data.table(strict_res$metrics_summary)[, analysis := "Non-lipoprotein metabolites"]
  ), fill = TRUE)
  exclusion_counts <- excluded[, .N, by = .(category)][order(-N)]
  if (!nrow(exclusion_counts)) exclusion_counts <- data.table(category = "Excluded lipid block", N = nrow(excluded))
  metric_keep <- melt(
    compare_metrics[, .(analysis, model, RMSE, Spearman, balanced_accuracy, macro_F1, weighted_kappa)],
    id.vars = c("analysis", "model"), variable.name = "metric", value.name = "value"
  )
  metric_keep <- metric_keep[model == "D1_explainable_biomarker_stage"]
  metric_keep[, higher_better := metric != "RMSE"]
  metric_keep[, scaled := {
    rr <- range(value, na.rm = TRUE)
    if (!all(is.finite(rr)) || diff(rr) == 0) rep(.5, .N) else if (higher_better[1]) {
      (value - rr[1]) / diff(rr)
    } else {
      (rr[2] - value) / diff(rr)
    }
  }, by = metric]

  p1 <- ggplot(compare_hr[stage_group == "Overall" & is.finite(HR)], aes(HR, analysis, xmin = lo, xmax = hi, color = analysis)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
    geom_errorbar(orientation = "y", width = .15) + geom_point(size = 3) +
    labs(title = "a. Residual risk-age after removing the lipoprotein-lipid block", x = "HR per within-stage SD", y = NULL, color = NULL) + theme_ckm(10)
  metric_span <- metric_keep[, .(x_from = min(scaled, na.rm = TRUE), x_to = max(scaled, na.rm = TRUE)), by = metric]
  p2 <- ggplot(metric_keep, aes(scaled, metric, color = analysis)) +
    geom_segment(
      data = metric_span,
      aes(x = x_from, xend = x_to, y = metric, yend = metric),
      inherit.aes = FALSE, color = "grey70", linewidth = .6
    ) +
    geom_point(size = 2.8) +
    geom_text(aes(label = sprintf("%.3f", value)), nudge_x = .04, hjust = 0, size = 2.8, show.legend = FALSE) +
    scale_x_continuous(limits = c(-.03, 1.16), breaks = c(0, .5, 1), labels = c("worse", "mid", "better")) +
    labs(title = "b. Reconstruction performance with and without the lipid block", x = "Within-metric relative performance", y = NULL, color = NULL) + theme_ckm(10)
  p3 <- ggplot(exclusion_counts, aes(N, fct_reorder(category, N), fill = category)) +
    geom_col(color = "white") + geom_text(aes(label = N), hjust = -0.1, size = 3) + scale_x_continuous(expand = expansion(mult = c(0, .15))) +
    labs(title = "c. What was removed in the strict metabolite sensitivity?", x = "Excluded features", y = NULL, fill = NULL) + theme_ckm(10) + theme(legend.position = "none")
  top_strict <- head(as.data.table(strict_res$weights)[order(fold, p)][, .(mean_beta = mean(beta, na.rm = TRUE), select_freq = .N), by = biomarker][order(-select_freq, -abs(mean_beta))], 15)
  top_strict <- merge(top_strict, meta2[, .(biomarker, label, category)], by = "biomarker", all.x = TRUE, sort = FALSE)
  top_strict[is.na(label), label := biomarker]
  top_strict[, label_plot := factor(label, levels = rev(unique(label)))]
  p4 <- ggplot(top_strict, aes(mean_beta, label_plot, fill = category)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") + geom_col() +
    labs(title = "d. Top non-lipoprotein metabolites driving the residual risk-age", x = "Mean stage-refinement beta", y = NULL, fill = "Category") + theme_ckm(10)
  save_plot((p1 | p2) / (p3 | p4), biom_fig_filename("FigS3.strict_sensitivity.png"), 16, 12)
  write_book(list(hr_comparison = compare_hr, metrics_comparison = compare_metrics, excluded_features = excluded,
                  exclusion_counts = exclusion_counts, strict_weights_head = head(as.data.table(strict_res$weights), 5000)),
             biom_fig_filename("FigS3.strict_sensitivity.out.xlsx"))
  invisible(list(strict = strict_res, compare_hr = compare_hr, compare_metrics = compare_metrics, excluded = excluded))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Consolidation and pipeline execution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_data_prep_outputs <- function(main, dat_target, dat_clock) {
  stage_tab <- as.data.table(dat_target)[
    ,
    .(
      target_N = .N,
      target_Yang = sum(yang, na.rm = TRUE),
      target_incident = sum(event, na.rm = TRUE),
      target_Yin = sum(yin, na.rm = TRUE)
    ),
    by = ckm.stage
  ]
  clock_tab <- as.data.table(dat_clock)[
    ,
    .(
      clock_N = .N,
      ckm4_Yang = sum(yang, na.rm = TRUE),
      ckm4_incident = sum(event, na.rm = TRUE),
      ckm4_Yin = sum(yin, na.rm = TRUE)
    ),
    by = ckm.stage
  ]
  stage_tab <- merge(stage_tab, clock_tab, by = "ckm.stage", all = TRUE, sort = TRUE)[order(ckm.stage)]

  biom_qc <- data.table(
    biomarker = main$features,
    missing_rate = vapply(main$features, function(v) mean(is.na(dat_target[[v]])), numeric(1))
  )
  meta <- data.table(
    item = c(
      "CKM definition", "Stage-3 method", "fixed CKM analysis target", "stage-4 date-source definition",
      "RMST target", "quantitative-axis cohort",
      "blood-biochemistry source", "biomarker mode", "biomarker analysis N",
      "full clock cohort N", "features after QC", "follow-up end", "ethnic filter",
      "scan all biomarkers", "base adjustment",
      "RMST estimand", "primary quantitative-stage estimand", "risk-calibrated sensitivity estimand",
      "clock adjustment mode", "clock requested tau", "strict metabolite sensitivity"
    ),
    value = c(
      "2026_AHA", stage3_method, "incident clinical CKM stage 4 (first CHD/HF/stroke/PAD/AF)", stage4_source_primary,
      "ckm4 clinical-CVD proxy", "full UKB before omic restriction",
      "all.rds", biom, nrow(dat_target),
      nrow(dat_clock), length(main$features), as.character(date_follow_end_global), ethnic_keep,
      scan_all_biom, pred_base_mode,
      "standardized restricted mean time free of first clinical CVD",
      "prospective CKM risk-age among s0-s3; not time to next stage",
      "adjusted future clinical-CVD risk coordinate for s0-s3; s4 shown separately",
      clock_adjust_mode, clock_tau,
      ifelse(biom == "met", met_strict_exclude_lipid_block, NA)
    )
  )

  write_tab(stage_tab, "data_prep.stage_counts.tsv")
  write_tab(biom_qc, "data_prep.biomarker_qc.tsv.gz")
  write_tab(main$metadata, "data_prep.biomarker_metadata.tsv")
  write_tab(meta, "data_prep.run_metadata.tsv")
  write_book(
    list(
      run_metadata = meta,
      stage_counts = stage_tab,
      biomarker_qc = biom_qc,
      biomarker_metadata = main$metadata
    ),
    "data_prep.xlsx"
  )
  invisible(list(stage = stage_tab, biom_qc = biom_qc, meta = meta))
}

make_cross_omic_summary <- function() {
  layers <- c("prot", "met")
  layer_dirs <- file.path(global_outroot, layers)
  meta_paths <- file.path(layer_dirs, "raw", "data_prep.run_metadata.tsv")
  if (!all(file.exists(meta_paths))) return(invisible(NULL))
  metas <- lapply(meta_paths, function(f) try(fread(f, showProgress = FALSE), silent = TRUE))
  if (any(vapply(metas, inherits, logical(1), "try-error"))) return(invisible(NULL))
  definition_ok <- vapply(metas, function(x) any(x$item == "CKM definition" & x$value == "2026_AHA"), logical(1))
  if (!all(definition_ok)) {
    progress_log("Cross-omic summary skipped: prot/met outputs do not use the current CKM definition")
    return(invisible(NULL))
  }

  read_layer <- function(layer, file) {
    f <- file.path(global_outroot, layer, "raw", file)
    if (!file.exists(f)) return(data.table())
    x <- try(fread(f, showProgress = FALSE), silent = TRUE)
    if (inherits(x, "try-error")) data.table() else as.data.table(x)[, layer := layer]
  }
  marker <- rbindlist(lapply(layers, read_layer, file = "marker_map.tsv.gz"), fill = TRUE)
  triage <- rbindlist(lapply(layers, read_layer, file = "causal_triage.summary.tsv.gz"), fill = TRUE)
  advance <- rbindlist(lapply(layers, read_layer, file = "biom_refine.advancement_cox.tsv"), fill = TRUE)
  if (nrow(advance) && "stage_group" %in% names(advance)) advance <- advance[stage_group == "Overall"]
  pred <- rbindlist(lapply(layers, read_layer, file = "prediction.performance_summary.tsv"), fill = TRUE)
  if (!nrow(marker) || !nrow(advance)) return(invisible(NULL))

  marker_counts <- marker[, .N, by = .(layer, marker_class)]
  marker_counts[, proportion := N / sum(N), by = layer]
  triage_counts <- if (nrow(triage)) {
    triage[, .N, by = .(layer, evidence_class)][, proportion := N / sum(N), by = layer]
  } else data.table()
  pred_keep <- if (nrow(pred)) pred[model %in% c("M0_base", "M2_risk_age_stage", "M2b_explainable_Yin_score", "M3_risk_age_stage_YY_score")] else data.table()
  if (nrow(pred_keep)) pred_keep[, model_plot := factor(
    pretty_prediction_model(model),
    levels = rev(pretty_prediction_model(c("M0_base", "M2_risk_age_stage", "M2b_explainable_Yin_score", "M3_risk_age_stage_YY_score")))
  )]

  p1 <- ggplot(marker_counts, aes(layer, proportion, fill = marker_class)) +
    geom_col(position = "fill") + scale_y_continuous(labels = percent) +
    labs(title = "a. Marker-evidence composition", x = NULL, y = "Proportion", fill = NULL) +
    theme_ckm(10) + guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")

  advance[, label_x := pmax(HR, hi, na.rm = TRUE) + 0.01]
  p2 <- ggplot(advance, aes(HR, layer, xmin = lo, xmax = hi)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
    geom_errorbar(orientation = "y", width = .15) + geom_point(size = 3) +
    geom_text(aes(x = label_x, label = sprintf("HR %.2f", HR)), hjust = 0, size = 3) +
    labs(title = "b. Within-stage omic residual risk-age", x = "HR per SD for incident stage 4", y = NULL) +
    theme_ckm(10)

  p3 <- if (nrow(pred_keep)) {
    ggplot(pred_keep, aes(mean_cindex, model_plot, color = layer)) +
      geom_errorbar(aes(xmin = mean_cindex - sd_cindex, xmax = mean_cindex + sd_cindex),
                    orientation = "y", width = .14, position = position_dodge(width = .5)) +
      geom_point(size = 2.8, position = position_dodge(width = .5)) +
      geom_text(aes(label = sprintf("%.3f", mean_cindex)), hjust = -0.25,
                position = position_dodge(width = .5), size = 2.8, show.legend = FALSE) +
      scale_x_continuous(expand = expansion(mult = c(.03, .14))) +
      labs(title = "c. Held-out prediction", x = "Mean C-index across outer folds", y = NULL, color = NULL) +
      theme_ckm(9) + theme(legend.position = "bottom")
  } else ggplot() + annotate("text", x = 0, y = 0, label = "Prediction outputs unavailable") + theme_void()

  p4 <- if (nrow(triage_counts)) {
    ggplot(triage_counts, aes(layer, proportion, fill = evidence_class)) +
      geom_col(position = "fill") + scale_y_continuous(labels = percent) +
      labs(title = "d. Etiologic-prioritization composition", x = NULL, y = "Proportion", fill = NULL) +
      theme_ckm(9) + guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + theme(legend.position = "bottom")
  } else ggplot() + annotate("text", x = 0, y = 0, label = "Triage outputs unavailable") + theme_void()

  save_plot((p1 | p2) / (p3 | p4), "Fig7.cross_omic_comparison.png", 15, 10, dir = global_outroot)
  write_book(list(
    marker_class_counts = marker_counts,
    residual_risk_age = advance,
    prediction = pred_keep,
    etiologic_class_counts = triage_counts
  ), "Fig7.cross_omic_comparison.out.xlsx", dir = global_outroot)
  invisible(list(marker = marker_counts, advance = advance, prediction = pred_keep, triage = triage_counts))
}

make_consolidate <- function(results) {
  files <- data.table(file = list.files(outroot, recursive = TRUE, full.names = FALSE))
  files[, size_MB := round(file.info(file.path(outroot, file))$size / 1024^2, 3)]
  interpretation <- data.table(
    item = c(
      "Primary target", "RMST estimand", "Primary quantitative-stage estimand", "Yang use", "Yin use",
      "Causality interpretation", "Prediction design", "Major data limitation"
    ),
    value = c(
      trait,
      "Population-equivalent restricted mean time free of the clinical-CVD stage-4 proxy within the supported follow-up horizon",
      "Prospective CKM risk-age for s0-s3; s4 is retained as a separate Yang disease state",
      "Pre-baseline clinical-CVD Yang cases contribute disease-state, event-proximity, and fold-specific weighting evidence",
      "Baseline stage 0-3 participants provide prospective risk for the CKM/ckm4 target; RMST targets first CHD/HF/stroke/PAD/AF",
      "Lag-persistent prospective association is observational triage only; external cis-pQTL/mQTL MR and colocalization are required",
      "Base -> categorical stage -> CKM risk-age -> explainable Yin-only score -> full explainable YY score; optional Yin-only elastic-net benchmark",
      "No repeat stage 0-3 measurements; risk coordinates are not time, the s3-to-s4 anchor gap is not estimated, and no adjacent-stage years are claimed"
    )
  )
  write_book(list(interpretation = interpretation, output_index = files), "ckm.summary.xlsx")
  saveRDS(list(results = results, interpretation = interpretation, files = files), file.path(rawdir, "ckm.all_results.rds"), compress = "gzip")
  if (biom == "met") make_cross_omic_summary()
}

# CKM construction is an explicit first phase. Downstream defaults never create
# or replace ckm.rds implicitly.
if (rebuild_ckm) {
  build_ckm_data(force = rebuild_ckm)
} else if (run_step("build_ckm")) {
  build_ckm_data(force = FALSE)
} else if (!file.exists(ckm_file)) {
  stop(
    "ckm.rds is required for downstream analysis. First run ./ckm.sh --rebuild-ckm.",
    call. = FALSE
  )
} else {
  progress_log("Using existing compact ckm.rds: ", ckm_file)
}

main <- load_analysis_data()
dat_clock <- prepare_target_data(main$clock_dat, "ckm4")
dat_target <- if (analysis_requires_biom) prepare_target_data(main$dat, trait) else dat_clock
if (analysis_requires_biom && !length(main$features)) {
  stop("No biomarkers passed the >=80% nonmissing QC", call. = FALSE)
}
if (analysis_requires_biom) {
  progress_log(
    "Biomarker analysis N = ", fmt_int(nrow(dat_target)),
    "; biomarkers = ", fmt_int(length(main$features)),
    "; target = ", trait,
    "; target Yang = ", fmt_int(sum(dat_target$yang)),
    "; target eligible Yin/events = ", fmt_int(sum(dat_target$yin)), "/", fmt_int(sum(dat_target$event))
  )
}
progress_log(
  "Full-cohort CKM4 clock N = ", fmt_int(nrow(dat_clock)),
  "; Yin/events = ", fmt_int(sum(dat_clock$yin)), "/", fmt_int(sum(dat_clock$event)),
  "; Yang = ", fmt_int(sum(dat_clock$yang))
)

results <- list()
run_if("summarize_ckm", {
  results$summarize_ckm <- make_ckm_summary(main$clock_dat)
})

run_if("data_prep", {
  results$data_prep <- make_data_prep_outputs(main, dat_target, dat_clock)
})

clock_path <- file.path(global_rawdir, "stage_clock.rds")
if (run_step("stage_clock")) {
  if (run_shared_outputs || !file.exists(clock_path)) {
    t0 <- Sys.time()
    progress_log("========== RUN STEP: stage_clock ==========")
    results$stage_clock <- make_stage_clock(dat_clock, write_outputs = FALSE)
    save_rds_safe(results$stage_clock, clock_path, compress = FALSE)
    invisible(gc())
    progress_log("========== DONE STEP: stage_clock (", elapsed_since(t0), ") ==========")
  } else {
    progress_log("Skip shared stage_clock; using existing global output: ", clock_path)
  }
} else {
  progress_log("Skip step: stage_clock")
}
clock_needed <- any(vapply(
  c("stage_clock", "marker_map", "biom_refine", "yy_timeline",
    "stage_specific", "causal_triage", "prediction"),
  run_step,
  logical(1)
))
clock_obj <- NULL
if (clock_needed) {
  clock_obj <- if (!is.null(results$stage_clock)) {
    results$stage_clock
  } else if (file.exists(clock_path)) {
    readRDS(clock_path)
  } else {
    stop("stage_clock.rds is required; start from stage_clock", call. = FALSE)
  }
  if (is.null(clock_obj$risk_age) || !identical(clock_obj$stage3_method %||% "", stage3_method)) {
    stop(
      "Existing stage_clock.rds is missing or was generated for a different Stage-3 method. ",
      "Re-run with --start-step stage_clock and the intended --stage3-method.",
      call. = FALSE
    )
  }
}

marker_path <- file.path(rawdir, "marker_map.rds")
run_if("marker_map", {
  results$marker_map <- make_marker_map(dat_target, main$features, main$metadata, clock_obj)
})
marker_available <- if (!is.null(results$marker_map)) {
  results$marker_map
} else if (file.exists(marker_path)) {
  readRDS(marker_path)
} else {
  NULL
}
if (!is.null(marker_available)) {
  marker_available <- as.data.table(marker_available)
  required_marker_fields <- c("biomarker", "stage_p", "yin_cox_p", "yang_proximity_p")
  if (!all(required_marker_fields %in% names(marker_available)))
    stop("Existing marker_map.rds lacks required fields; re-run with --start-step marker_map.", call. = FALSE)
}

run_if("biom_refine", {
  if (is.null(marker_available)) stop("marker_map.rds is required; start from marker_map", call. = FALSE)
  results$biom_refine <- make_biom_refine(dat_target, main$features, main$metadata, marker_available, clock_obj)
  if (biom == "met" && met_strict_exclude_lipid_block) {
    results$met_strict_sensitivity <- run_met_strict_sensitivity(dat_target, main$features, main$metadata, marker_available, clock_obj, results$biom_refine)
  }
})

run_if("yy_timeline", {
  results$yy_timeline <- make_yy_timeline(dat_target, main$features, main$metadata, clock_obj, marker_available)
})

marker_scan <- NULL
if (run_step("stage_specific") || run_step("causal_triage")) {
  marker_scan <- marker_available
  if (is.null(marker_scan)) stop("marker_map.rds is required; start from marker_map", call. = FALSE)
}

stage_path <- file.path(rawdir, "stage_specific.rds")
run_if("stage_specific", {
  results$stage_specific <- make_stage_specific(dat_target, main$features, main$metadata, marker_scan)
})
stage_obj <- NULL
if (run_step("causal_triage")) {
  stage_obj <- if (!is.null(results$stage_specific)) {
    results$stage_specific
  } else if (file.exists(stage_path)) {
    readRDS(stage_path)
  } else {
    stop("stage_specific.rds is required; start from stage_specific", call. = FALSE)
  }
  if (
    is.null(stage_obj$summary) || !nrow(stage_obj$summary)
  ) {
    stop(
      "Existing stage_specific.rds is missing a valid summary table. ",
      "Re-run with --start-step stage_specific before causal_triage.",
      call. = FALSE
    )
  }
}

run_if("causal_triage", {
  results$causal_triage <- make_causal_triage(dat_target, main$features, main$metadata, marker_scan, stage_obj)
})

run_if("prediction", {
  if (skip_prediction) {
    progress_log("Prediction skipped because CKM_SKIP_PREDICTION=TRUE")
  } else {
    results$prediction <- run_prediction(dat_target, dat_clock, main$features, main$metadata, clock_obj)
  }
})

run_if("consolidate", {
  make_consolidate(results)
})

if (biom_outputs_active) {
  progress_log("Done. Shared CKM outputs written under: ", global_outroot,
               "; biomarker-specific outputs written under: ", outroot)
} else {
  progress_log("Done. Shared CKM outputs written under: ", global_outroot)
}
