#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Yin-Yang multi-omic biomarker pipeline v18: explainable YY-guided filtering/weighting for prediction
# v18 changes:
#   1) The primary YY prediction model is now an interpretable weighted score:
#      Yin training-fold Cox beta supplies direction, while Yang prevalence,
#      recurrence, and pseudo-time evidence filters variables and modifies
#      weights. This is not a second Yang glmnet or a larger black box.
#   2) Older YY score-block/adaptive-glmnet models are disabled by default and
#      only run with YYP_LEGACY_COMPLEX=TRUE / --pred-legacy-complex.
#   3) Fold-level YY-guided weights are written for manuscript tables.
# v17 changes:
#   1) Fig2 is an AOD/occurrence-state diagnostic instead of duplicating the
#      pseudo-time trajectory panel now carried by Fig5.
#   2) Fig4 emphasizes target-vs-off-target specificity and keeps representative
#      reactive, mirror/protective, and prospective-only candidates in the scan.
#   3) Fig6/Fig7 write explicit diagnostics for eta=0 model degeneracy and
#      mirror/protective outliers such as RAB6A.
# v16 changes:
#   1) prot,met is handled by yy.sh as two independent runs; prot+met is the only merged mode.
#   2) M2 uses the exact M1 biomarker set plus an unpenalized fold-learned YY score block.
#   3) Fig6 uses a fixed bottom-to-top M0-M6 order and directly shows paired YY-vs-M1 increments.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_arg <- function(key, default = NA_character_) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- match(key, args)
  if (!is.na(i) && i < length(args)) args[i + 1] else default
}

parse_vec <- function(x, default = character()) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  z <- trimws(unlist(strsplit(x, ",|;|\\s+")))
  z[nzchar(z)]
}

step_order <- c("data_prep", "root_qc", "aod_qc", "yy_timeline", "aod_recurrence", "adjustment", "temporal_decomp", "scan", "prediction", "consolidate")
start_step <- Sys.getenv("START_STEP", unset = get_arg("--start_step", get_arg("--start-step", "all")))
if (is.na(start_step) || start_step %in% c("", "all", "ALL")) start_step <- step_order[1]
if (start_step %in% c("help", "--help", "-h")) {
  cat("Step order:\n")
  cat(paste0(seq_along(step_order), ". ", step_order, collapse = "\n"), "\n")
  quit(save = "no")
}
if (!start_step %in% step_order) stop("Unknown START_STEP: ", start_step, "; available: ", paste(step_order, collapse = ", "), call. = FALSE)
run_step <- function(step) match(step, step_order) >= match(start_step, step_order)
run_if <- function(step, expr) {
  if (run_step(step)) {
    cat("\n========== RUN STEP:", step, "==========\n")
    eval(substitute(expr), envir = parent.frame())
    invisible(gc())
    cat("========== DONE STEP:", step, "==========\n")
  } else {
    cat("Skip step:", step, "\n")
    invisible(gc())
  }
}

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(pacman::p_load(
  data.table, tidyverse, survival, lubridate, splines, patchwork, scales, writexl, ggrepel, forcats, glmnet
))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Paths and global settings
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir0 <- Sys.getenv("DIR0", unset = "")
if (!nzchar(dir0)) dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
indir <- Sys.getenv("PHEDIR", unset = file.path(dir0, "data", "ukb", "phe"))
script_dir <- Sys.getenv("SCRIPT_DIR", unset = file.path(dir0, "scripts", "0f"))
biom_request <- tolower(gsub("[[:space:]]+", "", Sys.getenv("YY_BIOM", unset = get_arg("--biom", "prot"))))
if (grepl(",", biom_request, fixed = TRUE)) {
  stop("yy.R accepts one biomarker mode per process. Use yy.sh --biom prot,met for two independent runs.", call. = FALSE)
}
if (biom_request %in% c("met+prot", "prot.met", "met.prot", "protein+metabolite", "metabolite+protein")) biom_request <- "prot+met"
biom <- biom_request
if (!biom %in% c("prot", "met", "prot+met")) {
  stop("Unknown YY_BIOM/--biom: ", biom, "; use prot, met, or prot+met", call. = FALSE)
}
biom_layers <- if (biom == "prot+met") c("prot", "met") else biom
biom_label <- switch(biom, prot = "protein", met = "metabolite", `prot+met` = "multi-omic biomarker")
biom_label_plural <- switch(biom, prot = "proteins", met = "metabolites", `prot+met` = "multi-omic biomarkers")
outbase <- Sys.getenv("YY_OUTDIR", unset = file.path(dir0, "analysis", "yy"))
outroot <- file.path(outbase, biom)
rawroot <- file.path(outroot, "raw")
dir.create(outroot, recursive = TRUE, showWarnings = FALSE)
dir.create(rawroot, recursive = TRUE, showWarnings = FALSE)
setwd(outroot)

helper_file <- file.path(script_dir, "0phe.f.R")
if (!file.exists(helper_file)) stop("Cannot find 0phe.f.R: ", helper_file, call. = FALSE)
source(helper_file)

# Core adjustment variables are defined centrally in 0phe.f.R.  Do not infer
# alternative LE8 names here: the YY pipeline uses these vectors verbatim.
if (!exists("vars.basic", inherits = TRUE)) {
  stop("0phe.f.R must define vars.basic.", call. = FALSE)
}
if (!exists("vars.le8", inherits = TRUE)) {
  stop("0phe.f.R must define vars.le8.", call. = FALSE)
}
yy_vars_basic <- unique(as.character(vars.basic))
yy_vars_le8 <- unique(as.character(vars.le8))

Ys <- parse_vec(Sys.getenv("Ys", unset = get_arg("--traits", get_arg("--trait", get_arg("-Y", "cvd_cad")))), default = "cvd_cad")
Ys <- unique(Ys)
Y_first <- Ys[1]
date_source <- Sys.getenv("YY_DATE_SOURCE", unset = "icd10")
# AOD source handling:
#   This pipeline does NOT build raw AOD from icd10.code0/date0 on the fly.
#   It only uses aod_* / aodtxt_* columns already present in all.rds/fod.*.
#   If absent, recurrence panels are explicitly labeled as FOD fallback.
#   Use RUN_STEPS=icd,merge after updating get_icd_fast() to create aodtxt_* first.
require_real_aod <- toupper(Sys.getenv("YY_REQUIRE_AOD", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")
fig_standardize <- Sys.getenv("YY_FIG_STANDARDIZE", unset = Sys.getenv("YY_FIG1_STANDARDIZE", unset = "control"))
ethnic_keep <- Sys.getenv("YY_ETHNIC", unset = "White")
trend_p_show <- as.numeric(Sys.getenv("YY_TREND_P_SHOW", unset = "1e-4"))
min_bin_n <- as.integer(Sys.getenv("YY_MIN_BIN_N", unset = "20"))
scan_all_proteins <- toupper(Sys.getenv("YY_SCAN_ALL_PROTEINS", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")
max_scan_proteins <- as.integer(Sys.getenv("YY_MAX_PROTEINS", unset = "0"))

# Integrated YY-informed prediction settings. The prediction step uses Yang samples
# only inside each training fold to learn YY-derived protein weights, then evaluates
# risk prediction strictly in held-out baseline disease-free Yin samples.
skip_prediction <- toupper(Sys.getenv("YY_SKIP_PREDICTION", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")
K_outer <- as.integer(Sys.getenv("YYP_FOLDS", unset = "10"))
screen_n <- as.integer(Sys.getenv("YYP_SCREEN_N", unset = "700"))
top_n <- as.integer(Sys.getenv("YYP_TOP_N", unset = "300"))
lambda_choice <- Sys.getenv("YYP_LAMBDA", unset = "1se")
if (!lambda_choice %in% c("1se", "min")) lambda_choice <- "1se"
alpha_glmnet <- as.numeric(Sys.getenv("YYP_ALPHA", unset = "0.5"))
pred_max_proteins <- as.integer(Sys.getenv("YYP_MAX_PROTEINS", unset = "0"))
pred_min_bin_n <- as.integer(Sys.getenv("YYP_MIN_BIN_N", unset = as.character(min_bin_n)))
pred_min_yin_n <- as.integer(Sys.getenv("YYP_MIN_YIN_N", unset = "1000"))
pred_min_events <- as.integer(Sys.getenv("YYP_MIN_EVENTS", unset = "50"))
pred_residualize <- toupper(Sys.getenv("YYP_RESIDUALIZE", unset = "TRUE")) %in% c("TRUE", "T", "1", "YES", "Y")
pred_block_size <- as.integer(Sys.getenv("YYP_BLOCK_SIZE", unset = "250"))
pred_base_mode <- tolower(Sys.getenv("YYP_BASE_MODE", unset = "le8"))
# Backward compatibility with yy.v13.sh: the former name "clinical" now means LE8.
if (identical(pred_base_mode, "clinical")) pred_base_mode <- "le8"
if (!pred_base_mode %in% c("demog", "le8", "full")) pred_base_mode <- "le8"
pred_adaptive_eta <- suppressWarnings(as.numeric(parse_vec(Sys.getenv("YYP_ADAPTIVE_ETA", unset = "0,0.25,0.5,1"), default = c("0", "0.25", "0.5", "1"))))
pred_adaptive_eta <- unique(pred_adaptive_eta[is.finite(pred_adaptive_eta) & pred_adaptive_eta >= 0])
if (!length(pred_adaptive_eta)) pred_adaptive_eta <- c(0, 0.25, 0.5, 1)
pred_layer_balance <- tolower(Sys.getenv("YYP_LAYER_BALANCE", unset = "auto"))
if (pred_layer_balance == "auto") pred_layer_balance <- ifelse(biom == "prot+met", "sqrt", "none")
if (!pred_layer_balance %in% c("none", "proportional", "sqrt", "equal")) pred_layer_balance <- "none"
pred_layer_ablation <- toupper(Sys.getenv("YYP_LAYER_ABLATION", unset = "TRUE")) %in% c("TRUE", "T", "1", "YES", "Y")
pred_legacy_complex <- toupper(Sys.getenv("YYP_LEGACY_COMPLEX", unset = "FALSE")) %in% c("TRUE", "T", "1", "YES", "Y")
base_covars_env <- parse_vec(Sys.getenv("YYP_BASE_COVARS", unset = ""), default = character())
prot_meta_file <- Sys.getenv("YY_PROT_META", unset = file.path(dir0, "data.BIG", "gwas", "ppp", "ppp_3k_b38.bed"))
met_meta_file <- Sys.getenv("YY_MET_META", unset = file.path(dir0, "data", "ukb", "phe", "common", "met.lst"))

# When --scan-all-proteins is used, YY quick timeline plots still show only a
# small pWAS-prioritized subset rather than thousands of proteins.
yy_plot_top_n <- as.integer(Sys.getenv("YY_PLOT_TOP_N", unset = "6"))
yy_plot_always <- parse_vec(Sys.getenv("YY_PLOT_ALWAYS_PROTEINS", unset = "GDF15,PCSK9"), default = character())

timeline_adjust <- Sys.getenv("YY_TIMELINE_ADJUST", unset = "age_sex_PC_tdi")
context_protein_scope <- tolower(Sys.getenv("YY_CONTEXT_PROTEINS", unset = "all"))
root_balance_top_n <- as.integer(Sys.getenv("YY_ROOT_BALANCE_TOP_N", unset = "16"))
seed <- as.integer(Sys.getenv("SEED", unset = "2026"))
set.seed(seed)

vip_requested <- parse_vec(Sys.getenv("YY_BIOMARKERS", unset = Sys.getenv("YY_PROTEINS", unset = "")), default = character())
extra_requested <- parse_vec(Sys.getenv("YY_EXTRA_BIOMARKERS", unset = Sys.getenv("YY_EXTRA_PROTEINS", unset = "")), default = character())
vip_proteins <- character()
extra_proteins <- character()
all_show_proteins <- character()

protein_cols <- c(
  ANGPTL3 = "#D84A3A", FGA = "#008C95", APOE = "#9A7A00", GDF15 = "#2B70D6",
  PCSK9 = "#00A000", IL6 = "#D64FD3", APOB = "#E76F51", LDLR = "#8E6BBE",
  TNNI3 = "#6096BA", NTPROBNP = "#5E60CE", REN = "#2A9D8F", MMP12 = "#3A86FF"
)
BLUE <- "#2C7FB8"; GOLD <- "#D8A03D"; TEAL <- "#00A6B2"; ORANGE <- "#E07A3F"; GREEN <- "#4D9221"

# Broad, reproducible protein annotation by gene-symbol rules. This avoids a
# small hand-written "prior" list in the root context figure. If a protein is not
# matched by these conservative rules, it remains Other/unclassified. The output
# table can be manually refined later or joined to external annotation.
annotate_proteins_fallback <- function(proteins) {
  proteins <- unique(as.character(proteins))
  proteins <- proteins[!is.na(proteins) & nzchar(proteins)]
  p <- toupper(proteins)
  has <- function(pattern) grepl(pattern, p, perl = TRUE, ignore.case = TRUE)
  category <- rep("Other/unclassified", length(p))
  category[has("^(APO|PCSK9$|ANGPTL[0-9]|LDLR$|LPL$|LIPA$|LIPC$|CETP$|SORT1$|PLTP$|LCAT$|ABCA1$|ABCG[0-9]|NPC1L1$|OLR1$|SCARB1$|HMGCR$|MTTP$|LPA$|PON[0-9])")] <- "Lipid/lipoprotein"
  category[has("^(IL|TNF|TNFR|IFN|IFNGR|CXCL|CCL|CCR|CXCR|CSF|TGFB|TGFBR|CRP$|SAA|OSM$|LIF$|JAK|STAT)")] <- "Inflammation/cytokine"
  category[has("^(CD[0-9]|HLA|LILR|KIR|FCGR|FCER|IG[HKL]|MS4A|NKG|NCR|TLR|CLEC|SIGLEC|ITG|SELL$|SELE$|SELP$|ICAM|VCAM)")] <- "Immune/cell-surface"
  category[has("^(FGA$|FGB$|FGG$|F2$|F3$|F5$|F7$|F8$|F9$|F10$|F11$|F12$|F13|VWF$|SERPIN|PROC$|PROS1$|PLG$|C[1-9][A-Z]|CF[ABDHIP]|C1Q|MASP)")] <- "Coagulation/complement"
  category[has("^(COL|MMP|TIMP|LAM|LAMA|LAMB|LAMC|FN1$|ELN$|FBN|VCAN$|DCN$|LUM$|THBS|POSTN$|SPARC$|TNC$|ADAM|ADAMTS)")] <- "ECM/remodeling"
  category[has("^(TNN|TNNT|TNNI|MYH|MYL|ACTN|NPPB$|NPPA$|NT[._-]?PRO[._-]?BNP$|MYBPC|DES$|CKM$|MB$|FABP3$|MYO)")] <- "Cardiac/muscle injury"
  category[has("^(REN$|ACE$|ACE2$|AGT$|AGTR|UMOD$|HAVCR1$|CST3$|NPHS|AQP|SLC|ALDOB$)")] <- "Renal/RAAS/transport"
  category[has("^(INS$|IGF|IGFBP|FGF|LEP$|LEPR$|ADIPOQ$|RETN$|GCG$|GHRL$|GHR$|SHBG$|GDF15$|GDF[0-9]+$|MSTN$|INH|AMH$)")] <- "Metabolic/growth-stress"
  category[has("^(VEGF|FLT|KDR$|TEK$|ANGPT|PECAM|ENG$|EDN|NOS|SEMA|EPHB|EPHA|PDGF|PDGFR|FGFR)")] <- "Vascular/angiogenesis"
  category[has("^(ALB$|HP$|HPX$|TF$|TTR$|ORM|AHSG$|FETUB$|A2M$|CP$|APCS$|SERPINA|SERPINF|GC$|RBP4$)")] <- "Liver/acute-phase/carrier"
  category[has("^(CASP|BCL|BAX$|TP53|MDM|FAS$|FASLG$|TRAIL|TNFSF|TNFRSF|MICA$|MICB$|BIRC)")] <- "Apoptosis/tumor-immunity"
  category[p %in% c("PCSK9", "ANGPTL3", "APOB", "APOE", "LDLR")] <- "Lipid/lipoprotein"
  category[p %in% c("IL6")] <- "Inflammation/cytokine"
  category[p %in% c("GDF15")] <- "Metabolic/growth-stress"
  category[p %in% c("FGA")] <- "Coagulation/complement"
  category[p %in% c("TNNI3", "NTPROBNP", "NPPB", "NPPA")] <- "Cardiac/muscle injury"
  category[p %in% c("REN")] <- "Renal/RAAS/transport"
  category[p %in% c("MMP12")] <- "ECM/remodeling"
  prior <- ifelse(p %in% c("PCSK9", "ANGPTL3", "APOE", "APOB", "LDLR", "LPA", "HMGCR", "CETP", "LPL", "APOC3"),
                  "lipid/CAD-related prior",
                  ifelse(p %in% c("GDF15", "IL6", "NTPROBNP", "TNNI3", "CRP", "TNFRSF1A", "TNFRSF1B"),
                         "stress/injury marker prior", "other"))
  data.table(protein = proteins, protein_upper = p, category = category, prior = prior)
}


# Metadata from the official PPP BED and UKB metabolite list is preferred. The
# rule-based annotation above is only a fallback when metadata is unavailable.
biomarker_meta_global <- data.table(
  protein = character(), label = character(), layer = character(),
  category = character(), subgroup = character(), prior = character()
)

annotate_proteins <- function(proteins) {
  proteins <- unique(as.character(proteins))
  proteins <- proteins[!is.na(proteins) & nzchar(proteins)]
  fb <- annotate_proteins_fallback(proteins)
  if (!nrow(biomarker_meta_global)) {
    fb[, `:=`(label = protein, layer = ifelse(biom == "met", "met", "prot"), subgroup = NA_character_)]
    return(fb[])
  }
  mm <- unique(as.data.table(biomarker_meta_global)[protein %in% proteins])
  out <- merge(data.table(protein = proteins), mm, by = "protein", all.x = TRUE, sort = FALSE)
  out <- merge(out, fb[, .(protein, protein_upper, fallback_category = category, fallback_prior = prior)], by = "protein", all.x = TRUE, sort = FALSE)
  out[is.na(label) | !nzchar(label), label := protein]
  out[is.na(layer) | !nzchar(layer), layer := ifelse(grepl("__met$", protein), "met", "prot")]
  out[is.na(category) | !nzchar(category), category := fallback_category]
  out[is.na(category) | !nzchar(category), category := "Other/unclassified"]
  out[is.na(prior) | !nzchar(prior), prior := fallback_prior]
  out[is.na(prior) | !nzchar(prior), prior := "other"]
  out[, c("fallback_category", "fallback_prior") := NULL]
  out[]
}

category_palette <- function(categories) {
  categories <- unique(as.character(categories))
  categories <- categories[!is.na(categories) & nzchar(categories)]
  fixed <- category_cols[intersect(categories, names(category_cols))]
  miss <- setdiff(categories, names(fixed))
  if (length(miss)) {
    add <- scales::hue_pal()(length(miss)); names(add) <- miss
    fixed <- c(fixed, add)
  }
  fixed[categories]
}


protein_colors <- function(proteins) {
  proteins <- as.character(unique(proteins))
  known <- protein_cols[intersect(proteins, names(protein_cols))]
  miss <- setdiff(proteins, names(known))
  if (length(miss)) {
    add <- scales::hue_pal()(length(miss))
    names(add) <- miss
    known <- c(known, add)
  }
  known[proteins]
}

protein_prior <- annotate_proteins(names(protein_cols))[, .(protein, prior, category)]
prior_cols <- c(`lipid/CAD-related prior` = "#E76F51", `stress/injury marker prior` = "#6096FF", `metabolic prior` = "#9A7A00", other = "#22AA55")
category_cols <- c(
  `Lipid/lipoprotein` = "#E76F51", `Inflammation/cytokine` = "#D64FD3", `Immune/cell-surface` = "#6A4C93",
  `Coagulation/complement` = "#008C95", `ECM/remodeling` = "#3A86FF", `Cardiac/muscle injury` = "#6096BA",
  `Renal/RAAS/transport` = "#2A9D8F", `Metabolic/growth-stress` = "#2B70D6", `Vascular/angiogenesis` = "#00A6B2",
  `Liver/acute-phase/carrier` = "#9A7A00", `Apoptosis/tumor-immunity` = "#A23B72", `Other/unclassified` = "grey70"
)

message("YY Yin-Yang multi-omic pipeline v18")
message("Ys = ", paste(Ys, collapse = ", "))
message("biom = ", biom, " (", biom_label_plural, ")")
message("outbase = ", outbase)
message("outroot = ", outroot)
message("date_source = ", date_source)
message("AOD handling = precomputed aod_* / aodtxt_* only; require_real_aod = ", require_real_aod)
message("pipeline features = death single-event handling; adjusted timelines; AOD/FOD decomposition; lagged Cox; recurrence and treatment triage")
message("timeline_adjust = ", timeline_adjust)
message("context_protein_scope = ", context_protein_scope)
message("fig_standardize = ", fig_standardize)
message("prediction step = ", ifelse(skip_prediction, "skipped", "enabled"), "; folds = ", K_outer, "; screen_n = ", screen_n, "; top_n = ", top_n, "; min Yin/events = ", pred_min_yin_n, "/", pred_min_events, "; residualize = ", pred_residualize)
message("prediction design = base_mode: ", pred_base_mode, "; adaptive eta: ", paste(pred_adaptive_eta, collapse = ","), "; layer balance: ", pred_layer_balance, "; layer ablation: ", pred_layer_ablation, "; legacy complex YY models: ", pred_legacy_complex)
message("metadata: prot = ", prot_meta_file, "; met = ", met_meta_file)
message("vars.basic from 0phe.f.R = ", paste(yy_vars_basic, collapse = ", "))
message("vars.le8 from 0phe.f.R = ", paste(yy_vars_le8, collapse = ", "))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# General helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`%||%` <- function(a, b) if (!is.null(a)) a else b
bt <- function(x) paste0("`", gsub("`", "", x), "`")
Y_safe <- function(Y) gsub("[^A-Za-z0-9_.-]+", "_", Y)

set_scope <- function(scope = c("root", "trait"), Y = NULL) {
  scope <- match.arg(scope)
  if (scope == "root") {
    outdir <<- outroot
    rawdir <<- rawroot
  } else {
    ydir <- file.path(outroot, Y_safe(Y))
    outdir <<- ydir
    rawdir <<- file.path(ydir, "raw")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
  setwd(outdir)
  invisible(NULL)
}

cleanup_legacy_trait_outputs <- function(Y) {
  set_scope("trait", Y)
  legacy_stems <- c(
    "Fig1.yy_plot", "Fig2.yy_timeline", "Fig3.aod_recurrence",
    "Fig4.adjustment", "Fig5.prospective_specificity",
    "Fig6.temporal_decomp", "Fig7.performance", "Fig7.prediction",
    "Fig8.yy_parameter_map"
  )
  ff <- list.files(outdir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  if (!length(ff)) return(invisible(character()))
  pat <- paste0("^", Y_safe(Y), "\\.(", paste(gsub("\\.", "\\\\.", legacy_stems), collapse = "|"), ")")
  hit <- ff[grepl(pat, basename(ff))]
  if (length(hit)) {
    unlink(hit, recursive = FALSE, force = TRUE)
    message("Removed ", length(hit), " legacy disease-figure file(s) for ", Y)
  }
  invisible(hit)
}

fig_file <- function(n, stem, ext = ".png", Y = NULL) {
  stem <- gsub("^\\.+|\\.+$", "", stem)
  parts <- if (is.null(Y)) c(paste0("Fig", n), stem) else c(Y_safe(Y), paste0("Fig", n), stem)
  paste0(paste(parts[nzchar(parts)], collapse = "."), ext)
}

fig_data_file <- function(n, stem, Y = NULL) {
  fig_file(n, stem, ext = "", Y = Y)
}

summary_file <- function(Y, ext = ".xlsx") {
  paste0("summary.", Y_safe(Y), ext)
}

as_date2 <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "IDate")) return(as.Date(x))
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))
  if (is.character(x)) return(suppressWarnings(as.Date(x)))
  y <- suppressWarnings(as.numeric(unclass(x)))
  y[!is.finite(y)] <- NA_real_
  as.Date(y, origin = "1970-01-01")
}

pmin_date2 <- function(...) {
  lst <- list(...)
  if (length(lst) == 1 && is.data.frame(lst[[1]])) lst <- as.list(lst[[1]])
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(as.Date(rep(NA, 0)))
  lst <- lapply(lst, as_date2)
  z <- do.call(pmin, c(lst, list(na.rm = TRUE)))
  z <- as_date2(z)
  z[!is.finite(as.numeric(z))] <- as.Date(NA)
  z
}

first_present <- function(x, candidates) {
  hit <- intersect(candidates, names(x))
  if (length(hit)) hit[1] else NA_character_
}
present_cols <- function(x, candidates) intersect(candidates, names(x))

n_aod_dates <- function(x) {
  out <- integer(length(x))
  ok <- !is.na(x) & nzchar(as.character(x))
  if (any(ok)) out[ok] <- lengths(strsplit(as.character(x[ok]), "|", fixed = TRUE))
  out
}

cap_count_group <- function(n) {
  factor(ifelse(is.na(n), NA_character_, ifelse(n <= 0, "0", ifelse(n == 1, "1", ifelse(n == 2, "2", ifelse(n == 3, "3", "4+"))))),
         levels = c("0", "1", "2", "3", "4+"))
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_covars <- function(dat, covars) {
  covars <- unique(covars)
  covars <- covars[covars %in% names(dat)]
  covars[vapply(covars, function(v) {
    z <- dat[[v]]
    if (all(is.na(z))) return(FALSE)
    length(unique(z[!is.na(z)])) > 1
  }, logical(1))]
}

# Use the centrally defined vectors from 0phe.f.R exactly.  Missing variables
# are treated as a data/configuration error rather than being guessed from
# similarly named columns.
get_le8_covars <- function(dat) {
  missing <- setdiff(yy_vars_le8, names(dat))
  if (length(missing)) {
    stop("Variables listed in vars.le8 are missing from the analysis data: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  safe_covars(dat, yy_vars_le8)
}

default_demog_covars <- function(dat, include_tdi = TRUE) {
  basic <- yy_vars_basic
  if (!include_tdi) basic <- setdiff(basic, "tdi")
  missing <- setdiff(basic, names(dat))
  if (length(missing)) {
    stop("Variables listed in vars.basic are missing from the analysis data: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  safe_covars(dat, basic)
}

get_timeline_covars <- function(dat, mode = timeline_adjust) {
  mode <- tolower(mode %||% "raw")
  if (identical(mode, "clinical")) mode <- "le8"  # backward compatibility
  basic_no_tdi <- default_demog_covars(dat, include_tdi = FALSE)
  demog <- default_demog_covars(dat, include_tdi = TRUE)
  le8 <- get_le8_covars(dat)
  drugs <- c("drug.lipid", "drug.dm", "drug.htn", "drug.aspirin")
  genetic <- grep("cad.*prs|cad[._]6m|apoe|PCSK9|ANGPTL3|APOB|LDLR|rare|lof|plof", names(dat), value = TRUE, ignore.case = TRUE)
  genetic <- setdiff(genetic[!grepl("__", genetic)], intersect(c(vip_proteins, extra_proteins), names(dat)))
  covs <- switch(mode,
    raw = character(), none = character(), unadjusted = character(),
    age_sex_pc = basic_no_tdi,
    age_sex_pc_tdi = demog,
    demog = demog,
    le8 = c(demog, le8),
    drugs = c(demog, le8, drugs),
    genetic = c(demog, le8, genetic),
    full = c(demog, le8, drugs, genetic),
    full_available = c(demog, le8, drugs, genetic),
    demog
  )
  safe_covars(dat, covs)
}

zscore <- function(x) {
  x <- safe_num(x)
  s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

winsorize <- function(x, p = c(0.001, 0.999)) {
  x <- safe_num(x)
  if (sum(is.finite(x)) < 10) return(x)
  q <- stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8)
  pmin(pmax(x, q[1]), q[2])
}

standardize_by_ref <- function(x, ref_idx) {
  x <- safe_num(x)
  ref_idx <- ref_idx %in% TRUE & is.finite(x)
  mu <- mean(x[ref_idx], na.rm = TRUE)
  s <- sd(x[ref_idx], na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mu) / s)
}

safe_mean <- function(x) { x <- safe_num(x); x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
safe_sd <- function(x) { x <- safe_num(x); x <- x[is.finite(x)]; if (length(x) > 1) sd(x) else NA_real_ }
safe_median <- function(x) { x <- safe_num(x); x <- x[is.finite(x)]; if (length(x)) median(x) else NA_real_ }

row_any_nonmissing <- function(dat, cols, block_size = 100) {
  cols <- intersect(cols, names(dat))
  out <- rep(FALSE, nrow(dat))
  if (!length(cols)) return(out)
  blocks <- split(cols, ceiling(seq_along(cols) / block_size))
  for (bb in blocks) out <- out | (rowSums(!is.na(dat[, bb, drop = FALSE])) > 0)
  out
}


write_tab <- function(x, fn) {
  f <- file.path(rawdir, fn)
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(as.data.table(x), f)
  cat("Wrote:", f, "\n")
  invisible(f)
}

write_book <- function(x, fn) {
  x <- x[vapply(x, function(z) is.data.frame(z) || data.table::is.data.table(z), logical(1))]
  x <- x[vapply(x, function(z) nrow(as.data.frame(z)) > 0, logical(1))]
  if (!length(x)) x <- list(empty = data.frame(note = "No rows"))
  names(x) <- substr(make.unique(gsub("[^A-Za-z0-9_]+", "_", names(x))), 1, 31)
  f <- file.path(outdir, fn)
  writexl::write_xlsx(lapply(x, as.data.frame), f)
  cat("Wrote:", f, "\n")
  invisible(f)
}

trim_result_payload <- function(x) {
  if (inherits(x, c("gg", "ggplot", "patchwork"))) return(NULL)
  if (!is.list(x) || is.data.frame(x) || data.table::is.data.table(x)) return(x)
  x$dat <- NULL
  x$plot <- NULL
  out <- lapply(x, trim_result_payload)
  out[!vapply(out, is.null, logical(1))]
}

save_plot <- function(p, fn, width = 12, height = 8, dpi = 300) {
  f <- file.path(outdir, fn)
  ggplot2::ggsave(f, p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  cat("Wrote:", f, "\n")
  invisible(f)
}

theme_yy <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", hjust = 0),
          plot.subtitle = element_text(color = "grey30"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "grey15"),
          strip.background = element_rect(fill = "grey92", color = NA),
          strip.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Protein transforms and model helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_protein_z <- function(dat, proteins, covars = character(), ref_idx = rep(TRUE, nrow(dat)), suffix = "z", winsor_p = c(0.001, 0.999)) {
  proteins <- intersect(proteins, names(dat))
  covars <- safe_covars(dat, covars)
  out <- dat
  ref_idx <- ref_idx %in% TRUE
  if (sum(ref_idx, na.rm = TRUE) < 200) ref_idx <- rep(TRUE, nrow(dat))
  for (p in proteins) {
    y <- winsorize(out[[p]], winsor_p)
    z <- rep(NA_real_, nrow(out))
    if (!length(covars)) {
      z <- standardize_by_ref(y, ref_idx)
    } else {
      dfit <- out[, covars, drop = FALSE]
      dfit$y <- y
      ok_ref <- ref_idx & complete.cases(dfit)
      ok_new <- complete.cases(dfit[, covars, drop = FALSE]) & is.finite(y)
      if (sum(ok_ref) >= 200) {
        fit <- try(lm(reformulate(covars, response = "y"), data = dfit[ok_ref, , drop = FALSE]), silent = TRUE)
        if (!inherits(fit, "try-error")) {
          pred <- rep(NA_real_, nrow(out))
          pred[ok_new] <- as.numeric(predict(fit, newdata = dfit[ok_new, , drop = FALSE]))
          r <- y - pred
          z <- standardize_by_ref(r, ref_idx)
        }
      }
      if (all(is.na(z))) z <- standardize_by_ref(y, ref_idx)
    }
    out[[paste0(p, "__", suffix)]] <- z
  }
  out
}

trend_bin_lm <- function(dd, x_col = "mid", y_col = "mean_z", w_col = "n") {
  dd <- dd %>% filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]), is.finite(.data[[w_col]]), .data[[w_col]] > 0)
  if (nrow(dd) < 4) return(tibble(beta10 = NA_real_, se10 = NA_real_, p = NA_real_))
  fit <- lm(as.formula(paste0(y_col, " ~ I(", x_col, "/10)")), data = dd, weights = dd[[w_col]])
  co <- summary(fit)$coef
  tibble(beta10 = co[2, 1], se10 = co[2, 2], p = co[2, 4])
}

bin_protein_means <- function(dat, xvar, z_cols, proteins, breaks = seq(-16, 16, 2), skip = c(1, 0), min_n = min_bin_n) {
  x <- safe_num(dat[[xvar]])
  mids <- head(breaks, -1) + diff(breaks) / 2
  bin <- cut(x, breaks = breaks, include.lowest = TRUE, right = TRUE)
  lv <- levels(cut(mids, breaks = breaks, include.lowest = TRUE, right = TRUE))
  out <- purrr::map2_dfr(proteins, z_cols, function(protein_name, zc) {
    y <- safe_num(dat[[zc]])
    tibble(bin = bin, y = y) %>%
      filter(!is.na(bin), is.finite(y)) %>%
      group_by(bin) %>%
      summarise(mean_z = mean(y, na.rm = TRUE), sd_z = sd(y, na.rm = TRUE), n = n(), .groups = "drop") %>%
      right_join(tibble(bin = factor(lv, levels = lv), mid = mids), by = "bin") %>%
      mutate(protein = protein_name, mean_z = if_else(n < min_n, NA_real_, mean_z))
  })
  if (sum(skip) > 0) {
    drop_i <- unique(c(seq_len(skip[1]), seq(length(mids) - skip[2] + 1, length(mids))))
    drop_i <- drop_i[drop_i >= 1 & drop_i <= length(mids)]
    out <- out %>% mutate(mean_z = if_else(match(mid, mids) %in% drop_i, NA_real_, mean_z))
  }
  out
}

predict_trend_lines <- function(bm) {
  bm %>% group_by(protein) %>% group_modify(function(.x, .g) {
    dd <- .x %>% filter(is.finite(mid), is.finite(mean_z), is.finite(n), n > 0)
    if (nrow(dd) < 4) return(tibble(mid = .x$mid, yhat = NA_real_))
    fit <- lm(mean_z ~ I(mid / 10), data = dd, weights = dd$n)
    tibble(mid = .x$mid, yhat = as.numeric(predict(fit, newdata = .x)))
  }) %>% ungroup()
}

cox_one <- function(dat, time_col, event_col, x_col, covars = character(), min_n = 500, min_event = 20) {
  covars <- safe_covars(dat, covars)
  need <- unique(c(time_col, event_col, x_col, covars))
  if (!all(need %in% names(dat))) return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, HR = NA_real_, lo = NA_real_, hi = NA_real_, n = 0, event = 0))
  dd <- dat[, need, drop = FALSE]
  dd <- dd[complete.cases(dd), , drop = FALSE]
  dd <- dd[is.finite(dd[[time_col]]) & dd[[time_col]] > 0, , drop = FALSE]
  if (nrow(dd) < min_n || sum(dd[[event_col]] == 1, na.rm = TRUE) < min_event) {
    return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, HR = NA_real_, lo = NA_real_, hi = NA_real_, n = nrow(dd), event = sum(dd[[event_col]] == 1, na.rm = TRUE)))
  }
  fm <- as.formula(paste0("survival::Surv(", bt(time_col), ", ", bt(event_col), ") ~ ", paste(bt(c(x_col, covars)), collapse = " + ")))
  fit <- try(survival::coxph(fm, data = dd), silent = TRUE)
  if (inherits(fit, "try-error")) return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, HR = NA_real_, lo = NA_real_, hi = NA_real_, n = nrow(dd), event = sum(dd[[event_col]] == 1, na.rm = TRUE)))
  sm <- summary(fit)
  co <- sm$coef[1, ]
  ci <- sm$conf.int[1, ]
  tibble(beta = co[["coef"]], se = co[["se(coef)"]], p = co[["Pr(>|z|)"]], HR = ci[["exp(coef)"]], lo = ci[["lower .95"]], hi = ci[["upper .95"]], n = nrow(dd), event = sum(dd[[event_col]] == 1, na.rm = TRUE))
}

cox_lag_one <- function(dat, time_col, event_col, x_col, covars = character(), lag_year = 0, min_n = 500, min_event = 20) {
  if (!all(c(time_col, event_col, x_col) %in% names(dat))) {
    return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, HR = NA_real_, lo = NA_real_, hi = NA_real_, n = 0, event = 0, lag_year = lag_year))
  }
  covars <- safe_covars(dat, covars)
  need <- unique(c(time_col, event_col, x_col, covars))
  if (!all(need %in% names(dat))) {
    return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, HR = NA_real_, lo = NA_real_, hi = NA_real_, n = 0, event = 0, lag_year = lag_year))
  }
  tt <- safe_num(dat[[time_col]])
  ee <- as.integer(safe_num(dat[[event_col]]) == 1)
  keep <- is.finite(tt) & tt > lag_year & !(ee == 1 & tt <= lag_year)
  dd <- dat[keep, need, drop = FALSE]
  if (!nrow(dd)) return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, HR = NA_real_, lo = NA_real_, hi = NA_real_, n = 0, event = 0, lag_year = lag_year))
  dd$.yy_lag_time <- safe_num(dd[[time_col]]) - lag_year
  dd$.yy_lag_event <- as.integer(safe_num(dd[[event_col]]) == 1)
  out <- cox_one(dd, ".yy_lag_time", ".yy_lag_event", x_col, covars = covars, min_n = min_n, min_event = min_event)
  out$lag_year <- lag_year
  out
}


lm_contrast_one <- function(dat, y_col, x_col, covars = character(), min_n = 500) {
  covars <- safe_covars(dat, covars)
  need <- unique(c(y_col, x_col, covars))
  if (!all(need %in% names(dat))) return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, n = 0))
  dd <- dat[, need, drop = FALSE]
  dd <- dd[complete.cases(dd), , drop = FALSE]
  if (nrow(dd) < min_n || length(unique(dd[[x_col]])) < 2) return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, n = nrow(dd)))
  fit <- try(lm(as.formula(paste0(bt(y_col), " ~ ", paste(bt(c(x_col, covars)), collapse = " + "))), data = dd), silent = TRUE)
  if (inherits(fit, "try-error")) return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, n = nrow(dd)))
  co <- summary(fit)$coef
  if (!x_col %in% rownames(co)) return(tibble(beta = NA_real_, se = NA_real_, p = NA_real_, n = nrow(dd)))
  tibble(beta = co[x_col, "Estimate"], se = co[x_col, "Std. Error"], p = co[x_col, "Pr(>|t|)"], n = nrow(dd))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Outcome date and AOD helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

candidate_date_cols <- function(dat, Y, source = date_source) {
  if (Y == "death") return(present_cols(dat, c("date_death", "fod_icd10_death")))
  cols_icd10 <- present_cols(dat, c(paste0("fod_icd10_", Y), paste0("date_", Y)))
  cols_all <- present_cols(dat, c(paste0("fod_icd10_", Y), paste0("fod_icd9_", Y), paste0("fod_srd_", Y), paste0("fod_gp_", Y), paste0("fod_ref_", Y), paste0("date_", Y), paste0(Y, "_date")))
  if (Y %in% c("cvd_cad", "mi", "ihd")) cols_all <- unique(c(cols_all, present_cols(dat, c("date_mi", "date_mace3"))))
  if (tolower(source) == "icd10" && length(cols_icd10)) cols_icd10 else cols_all
}

prepare_trait_data <- function(dat, Y) {
  datY <- dat
  y_cols_old <- grep(paste0("^", Y, "\\.(Yt2e|Yr2e|t2e|r2e|b2e)$"), names(datY), value = TRUE)
  if (length(y_cols_old)) datY <- datY[, setdiff(names(datY), y_cols_old), drop = FALSE]
  cols <- candidate_date_cols(datY, Y, date_source)
  if (!length(cols)) stop("No date columns found for outcome ", Y, ".", call. = FALSE)
  y_date_col <- paste0(Y, ".yy_date")
  datY[[y_date_col]] <- pmin_date2(datY[, cols, drop = FALSE])
  datY <- t2e(datY, NA, y_date_col, "birth_date", "date_attend", "date_lost", "date_death", date_follow_end, Y, "year")
  datY$b2e <- safe_num(datY[[paste0(Y, ".b2e")]])
  datY$t2e_y <- safe_num(datY[[paste0(Y, ".t2e")]])

  # Robust 0/1 flags.  Some t2e() helper versions leave Yr2e/Yt2e as NA
  # for non-events; treating those NA values literally can make the Yin
  # disease-free-at-baseline prediction subset nearly empty.  Here NA means
  # "no such event recorded" and is set to 0.
  yt_raw <- safe_num(datY[[paste0(Y, ".Yt2e")]])
  yr_raw <- safe_num(datY[[paste0(Y, ".Yr2e")]])
  datY$Yt2e_y <- fifelse(is.finite(yt_raw) & yt_raw == 1, 1L, 0L)
  datY$Yr2e_y <- fifelse(is.finite(yr_raw) & yr_raw == 1, 1L, 0L)

  # Yang: event recorded at or before baseline.  Yin: no pre-baseline event
  # and positive follow-up time; includes both incident cases and non-cases.
  # Same-day baseline events are conservatively treated as Yang.
  datY$yang_y <- (is.finite(datY$b2e) & datY$b2e <= 0) | datY$Yr2e_y == 1L
  datY$yin_y <- !datY$yang_y & is.finite(datY$t2e_y) & datY$t2e_y > 0
  datY$never_y <- datY$yin_y & datY$Yt2e_y == 0L
  attr(datY, "y_date_cols") <- cols
  attr(datY, "y_date_col") <- y_date_col
  make_occurrence_features(datY, Y)
}

parse_aod_one <- function(x) {
  if (inherits(x, "Date")) return(sort(unique(x[!is.na(x)])))
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  if (!length(x) || all(is.na(x))) return(as.Date(character()))
  x <- as.character(x)
  x <- unlist(strsplit(x, "\\|", fixed = FALSE), use.names = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x) & !is.na(x)]
  d <- suppressWarnings(as.Date(x))
  sort(unique(d[!is.na(d)]))
}

# Raw-code/date AOD construction is not performed here. Use precomputed aod_* or aodtxt_* columns.

get_aod_list <- function(dat, Y) {
  # Death is a single-event endpoint, not a recurrent disease phenotype. It has
  # no meaningful AOD list; --require-aod should therefore not stop death runs.
  if (identical(Y, "death")) {
    death_col <- first_present(dat, c("date_death", "fod_icd10_death", "death_date"))
    y_date_col <- attr(dat, "y_date_col")
    if (is.na(death_col) && !is.null(y_date_col) && y_date_col %in% names(dat)) death_col <- y_date_col
    if (!is.na(death_col)) {
      out <- as_date2(dat[[death_col]])
      attr(out, "source") <- paste0(death_col, "_single_event")
      attr(out, "is_real_aod") <- FALSE
      attr(out, "is_single_event") <- TRUE
      attr(out, "n_with_occurrence") <- sum(!is.na(out))
      return(out)
    }
  }
  list_col <- first_present(dat, c(paste0("aod_icd10_", Y), paste0("aod_", Y)))
  text_col <- first_present(dat, c(paste0("aodtxt_icd10_", Y), paste0("aodtxt_", Y)))
  if (!is.na(list_col)) {
    z <- dat[[list_col]]
    if (is.list(z)) {
      out <- lapply(z, parse_aod_one)
      attr(out, "source") <- list_col
      attr(out, "is_real_aod") <- TRUE
      attr(out, "n_with_occurrence") <- sum(lengths(out) > 0)
      return(out)
    }
  }
  if (!is.na(text_col)) {
    out <- lapply(dat[[text_col]], parse_aod_one)
    attr(out, "source") <- text_col
    attr(out, "is_real_aod") <- TRUE
    attr(out, "n_with_occurrence") <- sum(lengths(out) > 0)
    return(out)
  }

  if (isTRUE(require_real_aod)) {
    stop("No real AOD columns found for ", Y,
         ". Rebuild fod.* and all.rds with aodtxt_* columns first, or run without YY_REQUIRE_AOD.", call. = FALSE)
  }

  y_date_col <- attr(dat, "y_date_col")
  if (!is.null(y_date_col) && y_date_col %in% names(dat)) {
    out <- as_date2(dat[[y_date_col]])
    attr(out, "source") <- "FOD_fallback_only"
    attr(out, "is_real_aod") <- FALSE
    attr(out, "n_with_occurrence") <- sum(!is.na(out))
    warning("No real AOD source available for ", Y, "; occurrence diagnostics are FOD fallback only, not true AOD.", call. = FALSE)
    return(out)
  }
  out <- vector("list", nrow(dat))
  attr(out, "source") <- "none"
  attr(out, "is_real_aod") <- FALSE
  attr(out, "n_with_occurrence") <- 0L
  out
}

make_occurrence_features <- function(dat, Y, baseline_col = "date_attend") {
  aods <- get_aod_list(dat, Y)
  aod_source <- attr(aods, "source") %||% "unknown"
  aod_real <- isTRUE(attr(aods, "is_real_aod"))
  aod_single_event <- isTRUE(attr(aods, "is_single_event"))
  aod_n_with_occurrence <- attr(aods, "n_with_occurrence") %||% if (inherits(aods, "Date")) sum(!is.na(aods)) else sum(lengths(aods) > 0)
  message("Occurrence source for ", Y, ": ", aod_source, " | real AOD = ", aod_real, " | single-event = ", aod_single_event, " | N with occurrence = ", aod_n_with_occurrence)
  b0 <- as_date2(dat[[baseline_col]])
  n <- nrow(dat)
  out <- tibble(
    occ_n_all = integer(n), occ_n_pre = integer(n), occ_n_post = integer(n),
    occ_n_pre_2y = integer(n), occ_n_post_2y = integer(n),
    occ_first_date = as.Date(rep(NA, n)), occ_last_date = as.Date(rep(NA, n)),
    occ_last_pre_date = as.Date(rep(NA, n)), occ_first_post_date = as.Date(rep(NA, n)), occ_nearest_date = as.Date(rep(NA, n)),
    occ_first_b2e = NA_real_, occ_last_b2e = NA_real_,
    occ_last_pre_b2e = NA_real_, occ_first_post_b2e = NA_real_, occ_nearest_b2e = NA_real_,
    occ_span_year = NA_real_
  )
  if (inherits(aods, "Date")) {
    d <- as_date2(aods)
    ok <- !is.na(d) & !is.na(b0)
    dt <- rep(NA_real_, n)
    dt[ok] <- as.numeric(d[ok] - b0[ok]) / 365.25
    pre <- ok & dt < 0
    post <- ok & dt > 0
    out$occ_n_all[ok] <- 1L
    out$occ_n_pre[pre] <- 1L
    out$occ_n_post[post] <- 1L
    out$occ_n_pre_2y[pre & dt >= -2] <- 1L
    out$occ_n_post_2y[post & dt <= 2] <- 1L
    out$occ_first_date[ok] <- d[ok]
    out$occ_last_date[ok] <- d[ok]
    out$occ_nearest_date[ok] <- d[ok]
    out$occ_first_b2e[ok] <- dt[ok]
    out$occ_last_b2e[ok] <- dt[ok]
    out$occ_nearest_b2e[ok] <- dt[ok]
    out$occ_span_year[ok] <- 0
    out$occ_last_pre_date[pre] <- d[pre]
    out$occ_last_pre_b2e[pre] <- dt[pre]
    out$occ_first_post_date[post] <- d[post]
    out$occ_first_post_b2e[post] <- dt[post]
  } else {
    for (i in seq_len(n)) {
      d <- as.Date(aods[[i]])
      d <- sort(unique(d[!is.na(d)]))
      if (!length(d) || is.na(b0[i])) next
      dt <- as.numeric(d - b0[i]) / 365.25
      pre <- dt < 0; post <- dt > 0
      out$occ_n_all[i] <- length(d)
      out$occ_n_pre[i] <- sum(pre)
      out$occ_n_post[i] <- sum(post)
      out$occ_n_pre_2y[i] <- sum(dt >= -2 & dt < 0)
      out$occ_n_post_2y[i] <- sum(dt > 0 & dt <= 2)
      out$occ_first_date[i] <- d[1]
      out$occ_last_date[i]  <- d[length(d)]
      out$occ_first_b2e[i] <- dt[1]
      out$occ_last_b2e[i]  <- dt[length(dt)]
      out$occ_span_year[i] <- if (length(d) > 1) as.numeric(max(d) - min(d)) / 365.25 else 0
      if (any(pre)) {
        j <- which.max(dt[pre]); out$occ_last_pre_date[i] <- d[pre][j]; out$occ_last_pre_b2e[i] <- dt[pre][j]
      }
      if (any(post)) {
        j <- which.min(dt[post]); out$occ_first_post_date[i] <- d[post][j]; out$occ_first_post_b2e[i] <- dt[post][j]
      }
      j <- which.min(abs(dt)); out$occ_nearest_date[i] <- d[j]; out$occ_nearest_b2e[i] <- dt[j]
    }
  }
  res <- bind_cols(dat, out) %>% mutate(
    occ_source = aod_source,
    occ_is_real_aod = aod_real,
    occ_is_single_event = aod_single_event,
    occ_has_post = occ_n_post > 0,
    occ_has_pre = occ_n_pre > 0,
    occ_recent_pre = occ_n_pre_2y > 0,
    occ_recent_post = occ_n_post_2y > 0,
    occ_recurrent = occ_n_all > 1,
    occ_count_group = cap_count_group(occ_n_all),
    occ_prepost_pattern = case_when(
      occ_n_all == 0 ~ "none",
      occ_n_pre > 0 & occ_n_post > 0 ~ "pre + post",
      occ_n_pre > 0 & occ_n_post == 0 ~ "pre only",
      occ_n_post > 0 & occ_n_pre == 0 ~ "post only",
      TRUE ~ "same-day/other"
    ),
    occ_prepost_pattern = factor(occ_prepost_pattern, levels = c("none", "pre only", "post only", "pre + post", "same-day/other")),
    yy_group = case_when(
      is.na(b2e) ~ "Never disease",
      b2e > 5 ~ "Incident >5y",
      b2e > 0 & b2e <= 5 ~ "Incident 0-5y",
      b2e <= 0 & (occ_recent_pre | occ_recent_post | b2e >= -2) ~ "Prevalent recent/active",
      b2e < -2 & occ_has_post ~ "Prevalent remote + post-AOD",
      b2e < -2 & !occ_has_post ~ "Prevalent remote only",
      TRUE ~ "Other"
    ),
    yy_group = factor(yy_group, levels = c("Never disease", "Incident >5y", "Incident 0-5y", "Prevalent remote only", "Prevalent remote + post-AOD", "Prevalent recent/active", "Other")),
    occ_state = case_when(
      is.na(b2e) ~ "Never disease",
      b2e > 0 & occ_n_all <= 1 ~ "Incident single-date",
      b2e > 0 & occ_n_all > 1 ~ "Incident multi-date",
      b2e < 0 & occ_n_post > 0 ~ "Prevalent + post-baseline recurrence",
      b2e < 0 & (occ_recent_pre | occ_recent_post | b2e >= -2) ~ "Prevalent recent/active",
      b2e < 0 & occ_n_all > 1 ~ "Prevalent multi-date remote",
      b2e < 0 & occ_n_all <= 1 ~ "Prevalent single-date remote",
      TRUE ~ "Other"
    ),
    occ_state = factor(occ_state, levels = c("Never disease", "Incident single-date", "Incident multi-date", "Prevalent single-date remote", "Prevalent multi-date remote", "Prevalent recent/active", "Prevalent + post-baseline recurrence", "Other"))
  )
  attr(res, "aod_source") <- aod_source
  attr(res, "aod_is_real") <- aod_real
  res
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data loading
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

read_layer_metadata <- function(layer, raw_features, rename_map = NULL) {
  raw_features <- as.character(raw_features)
  if (is.null(rename_map)) rename_map <- setNames(raw_features, raw_features)
  if (layer == "prot") {
    meta <- data.table(label = raw_features, category = NA_character_, subgroup = NA_character_)
    if (file.exists(prot_meta_file)) {
      mm <- try(data.table::fread(prot_meta_file, header = FALSE, select = c(4, 5),
                                  col.names = c("label", "category"), data.table = TRUE,
                                  showProgress = FALSE), silent = TRUE)
      if (!inherits(mm, "try-error") && nrow(mm)) {
        mm[, label := toupper(as.character(label))]
        mm[, category := as.character(category)]
        mm <- mm[!is.na(label) & nzchar(label), .(
          category = {z <- category[!is.na(category) & nzchar(category)]; if (length(z)) z[1] else NA_character_}
        ), by = label]
        meta <- merge(data.table(label = raw_features), mm, by = "label", all.x = TRUE, sort = FALSE)
        meta[, subgroup := NA_character_]
      }
    }
    fb <- annotate_proteins_fallback(raw_features)[, .(label = protein, fallback_category = category, fallback_prior = prior)]
    meta <- merge(meta, fb, by = "label", all.x = TRUE, sort = FALSE)
    meta[is.na(category) | !nzchar(category), category := fallback_category]
    meta[, prior := fallback_prior]
    meta[, fallback_category := NULL]
  } else {
    meta <- data.table(label = raw_features, category = NA_character_, subgroup = NA_character_, prior = "other")
    if (file.exists(met_meta_file)) {
      mm <- try(data.table::fread(met_meta_file, header = TRUE, fill = TRUE, data.table = TRUE,
                                  showProgress = FALSE), silent = TRUE)
      if (!inherits(mm, "try-error") && nrow(mm)) {
        nn <- tolower(names(mm))
        ncol_name <- match("met_name", nn)
        gcol_name <- match("group", nn)
        sgcol_name <- match("subgroup", nn)
        if (is.na(ncol_name) && ncol(mm) >= 2) ncol_name <- 2L
        if (is.na(gcol_name) && ncol(mm) >= 4) gcol_name <- 4L
        if (!is.na(ncol_name)) {
          mm2 <- data.table(label = as.character(mm[[ncol_name]]))
          mm2[, category := if (!is.na(gcol_name)) as.character(mm[[gcol_name]]) else NA_character_]
          mm2[, subgroup := if (!is.na(sgcol_name)) as.character(mm[[sgcol_name]]) else NA_character_]
          mm2[, category := gsub("_", " ", category, fixed = TRUE)]
          mm2 <- mm2[!is.na(label) & nzchar(label)]
          meta <- merge(data.table(label = raw_features), unique(mm2, by = "label"), by = "label", all.x = TRUE, sort = FALSE)
          meta[, prior := fifelse(grepl("lipid|lipoprotein|cholesterol|fatty|apolipoprotein", category, ignore.case = TRUE),
                                  "metabolic prior", "other")]
        }
      }
    }
    meta[is.na(category) | !nzchar(category), category := "Metabolite/unclassified"]
  }
  meta[, layer := layer]
  meta[, protein := unname(rename_map[label])]
  meta[is.na(protein) | !nzchar(protein), protein := label]
  meta[, category := trimws(as.character(category))]
  meta[is.na(category) | !nzchar(category), category := "Other/unclassified"]
  meta[is.na(prior) | !nzchar(prior), prior := "other"]
  unique(meta[, .(protein, label, layer, category, subgroup, prior)], by = "protein")
}

read_main_dataset <- function() {
  dat0 <- readRDS(file.path(indir, "Rdata", "all.rds")) %>% as.data.frame()
  dat0$eid <- as.character(dat0$eid)
  all_n <- nrow(dat0)
  if (nzchar(ethnic_keep) && !tolower(ethnic_keep) %in% c("all", "none") && "ethnic.c" %in% names(dat0)) {
    dat0 <- dat0 %>% filter(as.character(ethnic.c) == ethnic_keep)
  }
  all_filtered_n <- nrow(dat0)

  layer_dat <- list(); raw_features <- list(); layer_counts <- data.table()
  for (layer in biom_layers) {
    f <- file.path(indir, "Rdata", paste0(layer, ".rds"))
    if (!file.exists(f)) stop("Cannot find biomarker data: ", f, call. = FALSE)
    dd <- readRDS(f) %>% as.data.frame()
    dd$eid <- as.character(dd$eid)
    dd <- dd[!duplicated(dd$eid), , drop = FALSE]
    if (layer == "prot") names(dd)[-1] <- toupper(names(dd)[-1])
    layer_dat[[layer]] <- dd
    raw_features[[layer]] <- setdiff(names(dd), "eid")
    layer_counts <- rbind(layer_counts, data.table(layer = layer, rds_N = nrow(dd), raw_features = length(raw_features[[layer]])), fill = TRUE)
  }

  # Exact name collisions are explicitly suffixed only in the combined mode.
  rename_maps <- lapply(raw_features, function(v) setNames(v, v))
  if (length(biom_layers) > 1) {
    collisions <- Reduce(intersect, raw_features[biom_layers])
    if (length(collisions)) {
      for (layer in biom_layers) {
        rename_maps[[layer]][collisions] <- paste0(collisions, "__", layer)
        old <- collisions
        new <- unname(rename_maps[[layer]][collisions])
        idx <- match(old, names(layer_dat[[layer]]))
        names(layer_dat[[layer]])[idx] <- new
      }
      message("Renamed ", length(collisions), " exact prot/met name collisions with __prot/__met suffixes.")
    }
  }

  # Also protect against a biomarker name colliding with an all.rds covariate.
  covariate_names <- setdiff(names(dat0), "eid")
  for (layer in biom_layers) {
    current <- unname(rename_maps[[layer]][raw_features[[layer]]])
    hit_current <- intersect(current, covariate_names)
    if (length(hit_current)) {
      for (nm in hit_current) {
        original <- names(rename_maps[[layer]])[match(nm, rename_maps[[layer]])]
        new_nm <- paste0(nm, "__", layer)
        idx <- match(nm, names(layer_dat[[layer]]))
        if (!is.na(idx)) names(layer_dat[[layer]])[idx] <- new_nm
        rename_maps[[layer]][original] <- new_nm
      }
      message("Renamed ", length(hit_current), " biomarker/covariate name collisions in ", layer, ".")
    }
  }

  layer_eids <- lapply(layer_dat, function(dd) intersect(dat0$eid, dd$eid))
  cohort_eids <- if (length(layer_eids) == 1) layer_eids[[1]] else Reduce(intersect, layer_eids)
  dat <- dat0[dat0$eid %in% cohort_eids, , drop = FALSE]
  for (layer in biom_layers) {
    dd <- layer_dat[[layer]][layer_dat[[layer]]$eid %in% cohort_eids, , drop = FALSE]
    dat <- merge(dat, dd, by = "eid", all.x = TRUE, sort = FALSE)
  }
  for (cc in intersect(c("birth_date", "date_attend", "date_lost", "date_death"), names(dat))) dat[[cc]] <- as_date2(dat[[cc]])

  vars_by_layer <- lapply(biom_layers, function(layer) unname(rename_maps[[layer]][raw_features[[layer]]]))
  names(vars_by_layer) <- biom_layers
  vars_biom <- unique(unlist(vars_by_layer, use.names = FALSE))
  meta <- rbindlist(lapply(biom_layers, function(layer) read_layer_metadata(layer, raw_features[[layer]], rename_maps[[layer]])), fill = TRUE)
  layer_counts[, analysis_N := vapply(layer, function(z) sum(dat$eid %in% layer_dat[[z]]$eid), numeric(1))]
  list(dat = dat, vars_biom = vars_biom, vars_by_layer = vars_by_layer, metadata = meta,
       biom_n = nrow(dat), all_n = all_n, all_filtered_n = all_filtered_n, layer_counts = layer_counts)
}

resolve_feature_names <- function(requested, meta, available) {
  requested <- unique(as.character(requested)); requested <- requested[nzchar(requested)]
  if (!length(requested)) return(character())
  mm <- as.data.table(meta)[protein %in% available]
  ans <- character()
  for (x in requested) {
    hit <- mm[tolower(protein) == tolower(x) | tolower(label) == tolower(x), protein]
    if (length(hit)) ans <- c(ans, hit[1])
  }
  unique(ans)
}

set_default_display_features <- function(available, meta) {
  prot_main <- c("ANGPTL3", "FGA", "APOE", "GDF15", "PCSK9", "IL6")
  prot_extra <- c("APOB", "LDLR", "TNNI3", "NTPROBNP", "REN", "MMP12")
  met_main <- c("ApoB", "ApoA1", "GlycA", "LDL_C", "LDL.C", "HDL_C", "HDL.C", "TG", "Glucose", "Gln", "Gly", "Phe", "Lactate")
  vip <- resolve_feature_names(vip_requested, meta, available)
  extra <- resolve_feature_names(extra_requested, meta, available)
  if (!length(vip)) {
    defaults <- if (biom == "prot") prot_main else if (biom == "met") met_main else c(prot_main, met_main)
    vip <- resolve_feature_names(defaults, meta, available)
  }
  if (!length(extra)) {
    defaults <- if (biom == "prot") prot_extra else if (biom == "met") character() else prot_extra
    extra <- resolve_feature_names(defaults, meta, available)
  }
  if (length(vip) < min(6L, length(available))) vip <- unique(c(vip, head(setdiff(available, vip), min(6L, length(available)) - length(vip))))
  vip_proteins <<- head(vip, if (biom == "prot+met") 12L else 6L)
  extra_proteins <<- head(setdiff(extra, vip_proteins), 12L)
  all_show_proteins <<- unique(c(vip_proteins, extra_proteins))
  invisible(NULL)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Root-level figures: not disease-specific
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

assoc_covariate_heatmap <- function(dat, proteins) {
  proteins <- intersect(proteins, names(dat))
  if (!length(proteins)) return(data.table())
  ref_idx <- rep(TRUE, nrow(dat))
  dat2 <- make_protein_z(dat, proteins, covars = character(), ref_idx = ref_idx, suffix = "root")
  covars <- present_cols(dat2, c("age", "tdi", "bmi", "nonhdl", "bb_LDL", "bb_HDL", "bb_TG", "bb_TC", "bb_CRE", "egfr", "bb_HBA1C", "bb_GLU"))
  adj0 <- default_demog_covars(dat2, include_tdi = FALSE)
  rbindlist(lapply(proteins, function(p) {
    zc <- paste0(p, "__root")
    rbindlist(lapply(covars, function(cv) {
      adj <- setdiff(adj0, cv)
      dd <- dat2[, unique(c(zc, cv, adj)), drop = FALSE]
      dd <- dd[complete.cases(dd), , drop = FALSE]
      if (nrow(dd) < 500 || safe_sd(dd[[cv]]) == 0) return(data.table())
      fit <- try(lm(as.formula(paste0(bt(zc), " ~ ", paste(bt(c(cv, adj)), collapse = " + "))), data = dd), silent = TRUE)
      if (inherits(fit, "try-error")) return(data.table())
      co <- summary(fit)$coef
      if (!cv %in% rownames(co)) return(data.table())
      data.table(protein = p, covariate = cv, beta = co[cv, "Estimate"], se = co[cv, "Std. Error"], p = co[cv, "Pr(>|t|)"], n = nrow(dd))
    }), fill = TRUE)
  }), fill = TRUE)[, signed_logp := sign(beta) * pmin(-log10(pmax(p, 1e-300)), 20)][]
}

make_root_figures <- function(dat, vars_biom, biom_n, all_n, all_filtered_n = all_n, layer_counts = data.table()) {
  set_scope("root")
  proteins <- intersect(vip_proteins, names(dat))
  extra_present <- intersect(extra_proteins, names(dat))
  show <- unique(c(proteins, extra_present))
  all_proteins <- intersect(vars_biom, names(dat))
  context_proteins <- if (context_protein_scope %in% c("all", "all_ppp", "ppp")) all_proteins else show
  if (!length(context_proteins)) context_proteins <- show
  any_protein_idx <- if (length(all_proteins)) row_any_nonmissing(dat, all_proteins) else rep(FALSE, nrow(dat))
  if (sum(any_protein_idx, na.rm = TRUE) < 100 && length(show)) any_protein_idx <- rowSums(!is.na(dat[, show, drop = FALSE])) > 0
  nonmiss_vip <- if (length(proteins)) sum(rowSums(!is.na(dat[, proteins, drop = FALSE])) == length(proteins)) else 0
  nonmiss_vip_ppp <- if (length(proteins)) sum(any_protein_idx & rowSums(!is.na(dat[, proteins, drop = FALSE])) == length(proteins)) else 0
  sample_tbl <- rbindlist(list(
    data.table(item = c("all.rds participants", "all.rds after demographic filter"), N = c(all_n, all_filtered_n)),
    if (nrow(layer_counts)) layer_counts[, .(item = paste0(layer, ".rds participants"), N = rds_N)] else data.table(),
    data.table(item = paste0("analysis cohort: ", biom), N = nrow(dat)),
    data.table(item = paste0("complete selected ", biom_label_plural), N = nonmiss_vip)
  ), fill = TRUE)
  miss_tbl <- data.table(
    protein = show,
    denominator = paste0("analysis cohort: ", biom),
    missing_rate = sapply(show, function(protein_name) mean(is.na(dat[[protein_name]])))
  )

  value_dist <- rbindlist(lapply(show, function(protein_name) {
    z <- standardize_by_ref(winsorize(dat[[protein_name]]), any_protein_idx)
    zz <- z[any_protein_idx & is.finite(z)]
    if (!length(zz)) return(data.table())
    qs <- as.numeric(stats::quantile(zz, probs = c(.01, .25, .50, .75, .99), na.rm = TRUE, names = FALSE, type = 8))
    data.table(protein = protein_name, n = length(zz), p01 = qs[1], q25 = qs[2], median = qs[3], q75 = qs[4], p99 = qs[5], mean = mean(zz), sd = sd(zz))
  }), fill = TRUE)

  corr_tbl <- data.table()
  if (length(show) >= 2) {
    cm <- cor(dat[, show, drop = FALSE], use = "pairwise.complete.obs")
    corr_tbl <- as.data.table(as.table(cm)); names(corr_tbl) <- c("protein1", "protein2", "cor")
  }
  cov_tbl <- assoc_covariate_heatmap(dat, show)

  pA <- ggplot(sample_tbl, aes(reorder(item, N), N, fill = item)) +
    geom_col(width = .68, color = "white") + geom_text(aes(label = comma(N)), hjust = -0.08, fontface = "bold", size = 3.4) +
    coord_flip(clip = "off") + scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .25))) +
    scale_fill_manual(values = setNames(scales::hue_pal()(nrow(sample_tbl)), sample_tbl$item), guide = "none") +
    labs(title = "a. Sample flow", x = NULL, y = "Participants") + theme_yy(10)

  pB <- if (nrow(value_dist)) ggplot(value_dist, aes(y = reorder(protein, median), color = protein)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(aes(x = p01, xend = p99, yend = reorder(protein, median)), linewidth = .50, alpha = .65) +
    geom_segment(aes(x = q25, xend = q75, yend = reorder(protein, median)), linewidth = 1.6) +
    geom_point(aes(x = median), size = 2.1) +
    scale_color_manual(values = protein_colors(unique(value_dist$protein)), guide = "none") +
    labs(title = paste0("b. Selected-", biom_label, " distributions within ", toupper(biom)), subtitle = "Points are medians; thick bars are IQR; thin bars are 1st-99th percentiles after winsorized z-scoring", x = paste0(tools::toTitleCase(biom_label), " z-score"), y = NULL) + theme_yy(10)
  else ggplot() + labs(title = paste0("b. ", tools::toTitleCase(biom_label), " distributions unavailable")) + theme_yy(10)

  pC <- if (nrow(corr_tbl)) ggplot(corr_tbl, aes(protein1, protein2, fill = cor)) +
    geom_tile(color = "white") + geom_text(aes(label = sprintf("%.2f", cor)), size = 2.6) +
    scale_fill_gradient2(low = BLUE, mid = "white", high = ORANGE, midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(title = paste0("c. ", tools::toTitleCase(biom_label), " correlation"), x = NULL, y = NULL) + theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) else ggplot() + labs(title = "c. Correlation unavailable") + theme_yy(10)

  pD <- if (nrow(cov_tbl)) ggplot(cov_tbl, aes(covariate, protein, fill = signed_logp)) +
    geom_tile(color = "white") + geom_text(aes(label = sprintf("%.2f", beta)), size = 2.4) +
    scale_fill_gradient2(low = BLUE, mid = "white", high = ORANGE, midpoint = 0, name = "signed\n-log10(P)") +
    labs(title = "d. Baseline covariate association", subtitle = "Numbers are adjusted beta estimates", x = NULL, y = NULL) + theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1)) else ggplot() + labs(title = "d. Covariate association unavailable") + theme_yy(10)

  Fig1 <- (pA | pB) / (pC | pD)
  save_plot(Fig1, fig_file(1, paste0(biom, ".QC")), 15.5, 10)
  write_book(list(sample_flow = sample_tbl, biom_missingness = miss_tbl, biom_value_distribution = value_dist, biom_correlation = corr_tbl, covariate_associations = cov_tbl), fig_file(1, paste0(biom, ".QC"), ".out.xlsx"))
  write_tab(sample_tbl, fig_data_file(1, paste0(biom, ".sample_flow.tsv")))
  write_tab(miss_tbl, fig_data_file(1, paste0(biom, ".missingness.tsv")))
  if (nrow(value_dist)) write_tab(value_dist, fig_data_file(1, paste0(biom, ".value_distribution.tsv")))
  if (nrow(cov_tbl)) write_tab(cov_tbl, fig_data_file(1, paste0(biom, ".covariate_associations.tsv.gz")))

  # Root biomarker context: category composition over all available biomarkers, plus selected-biomarker covariate/correlation context.
  anno_all <- annotate_proteins(context_proteins)
  anno_show <- annotate_proteins(show)
  cat_count <- anno_all[, .N, by = category][order(-N)]
  p2A <- ggplot(cat_count, aes(N, reorder(category, N), fill = category)) + geom_col(width = .72, color = "white") +
    geom_text(aes(label = N), hjust = -0.08, fontface = "bold", size = 3.0) +
    scale_x_continuous(expand = expansion(mult = c(0, .18))) + scale_fill_manual(values = category_palette(cat_count$category), guide = "none") +
    labs(title = paste0("a. ", tools::toTitleCase(biom_label), " category composition"), subtitle = paste0("Metadata-based annotation of ", length(context_proteins), " ", toupper(biom), " biomarkers"), x = paste0("N ", biom_label_plural), y = NULL) + theme_yy(10)
  p2B <- if (nrow(cov_tbl)) {
    dd <- cov_tbl[covariate %in% c("age", "bb_CRE", "egfr", "nonhdl", "bb_LDL", "bmi", "tdi")]
    dd <- merge(dd, anno_show[, .(protein, category)], by = "protein", all.x = TRUE, sort = FALSE)
    ggplot(dd, aes(covariate, protein, size = -log10(pmax(p, 1e-300)), color = category)) +
      geom_point(alpha = .85) + scale_size_continuous(range = c(.8, 5), name = "-log10(P)") +
      scale_color_manual(values = category_palette(dd$category), name = "Category") +
      labs(title = "b. Covariate fingerprint", x = NULL, y = NULL) + theme_yy(10) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  } else ggplot() + labs(title = "b. Covariate fingerprint unavailable") + theme_yy(10)
  p2C <- if (nrow(corr_tbl)) {
    dd <- corr_tbl[as.character(protein1) < as.character(protein2)]
    dd <- dd[order(-abs(cor))]
    dd <- dd[seq_len(min(20, nrow(dd)))]
    dd[, pair := paste(protein1, protein2, sep = " / ")]
    ggplot(dd, aes(abs(cor), reorder(pair, abs(cor)), fill = cor > 0)) + geom_col(width = .68, color = "white") +
      scale_fill_manual(values = c(`TRUE` = ORANGE, `FALSE` = BLUE), guide = "none") +
      labs(title = paste0("c. Strongest selected-", biom_label, " correlations"), x = "Absolute correlation", y = NULL) + theme_yy(10)
  } else ggplot() + labs(title = "c. Correlation pairs unavailable") + theme_yy(10)
  p2D <- if (nrow(cov_tbl)) {
    dd <- cov_tbl[, .(mean_abs_logp = mean(abs(signed_logp), na.rm = TRUE)), by = protein]
    dd <- merge(dd, anno_show[, .(protein, category, prior)], by = "protein", all.x = TRUE, sort = FALSE)
    ggplot(dd, aes(mean_abs_logp, reorder(protein, mean_abs_logp), fill = category)) + geom_col(width = .68, color = "white") +
      scale_fill_manual(values = category_palette(dd$category)) + labs(title = "d. Overall baseline-covariate sensitivity", x = "Mean absolute signed -log10(P)", y = NULL, fill = "Category") + theme_yy(10)
  } else ggplot() + labs(title = "d. Sensitivity unavailable") + theme_yy(10)
  Fig2 <- (p2A | p2B) / (p2C | p2D)
  save_plot(Fig2, fig_file(2, paste0(biom, ".context")), 15.5, 10)
  write_book(list(biom_categories = anno_all, selected_biom_categories = anno_show, category_counts = cat_count, covariate_associations = cov_tbl, correlation = corr_tbl), fig_file(2, paste0(biom, ".context"), ".out.xlsx"))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Root-level AOD text QC: not disease-specific
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aodtxt_col_for_trait <- function(dat, trait, src = "icd10") {
  first_present(dat, c(paste0("aodtxt_", src, "_", trait), paste0("aodtxt_", trait)))
}

single_event_date_col_for_trait <- function(dat, trait) {
  if (identical(trait, "death")) return(first_present(dat, c("date_death", "fod_icd10_death", "death_date")))
  NA_character_
}

aod_dates_from_text_vector <- function(x) {
  lapply(as.character(x), parse_aod_one)
}

audit_aodtxt_trait <- function(dat, trait, src = "icd10") {
  single_col <- single_event_date_col_for_trait(dat, trait)
  if (!is.na(single_col)) {
    d <- as_date2(dat[[single_col]])
    return(data.table(trait = trait, src = src, aod_col = single_col, aod_type = "single_event", n_total = nrow(dat),
                      n_any = sum(!is.na(d)), n_single = sum(!is.na(d)), n_multi = 0L, pct_multi_among_any = 0,
                      max_dates = ifelse(sum(!is.na(d)) > 0, 1L, 0L), bad_fod = 0L, bad_lod = 0L,
                      bad_n_dates_gt_cnt = 0L, same_day_duplicate_code_records = 0L,
                      fod_col = single_col, lod_col = single_col, cnt_col = NA_character_))
  }
  aod_col <- aodtxt_col_for_trait(dat, trait, src)
  if (is.na(aod_col)) return(data.table(trait = trait, src = src, aod_col = NA_character_, aod_type = "missing", n_total = nrow(dat), n_any = 0L, n_single = 0L, n_multi = 0L, pct_multi_among_any = NA_real_, max_dates = NA_integer_, bad_fod = NA_integer_, bad_lod = NA_integer_, bad_n_dates_gt_cnt = NA_integer_, same_day_duplicate_code_records = NA_integer_, fod_col = NA_character_, lod_col = NA_character_, cnt_col = NA_character_))
  x <- as.character(dat[[aod_col]])
  nd <- n_aod_dates(x)
  ok <- which(nd > 0)
  fod_col <- first_present(dat, c(paste0("fod_", src, "_", trait), paste0("fod_", trait)))
  lod_col <- first_present(dat, c(paste0("lod_", src, "_", trait), paste0("lod_", trait)))
  cnt_col <- first_present(dat, c(paste0("cnt_", src, "_", trait), paste0("cnt_", trait)))
  bad_fod <- bad_lod <- integer(0)
  if (length(ok) && !is.na(fod_col)) {
    parts <- strsplit(x[ok], "|", fixed = TRUE)
    first_date <- vapply(parts, function(z) z[1], character(1))
    fod_chr <- format(as.Date(dat[[fod_col]][ok]), "%Y-%m-%d")
    bad_fod <- ok[is.na(fod_chr) | first_date != fod_chr]
  }
  if (length(ok) && !is.na(lod_col)) {
    parts <- strsplit(x[ok], "|", fixed = TRUE)
    last_date <- vapply(parts, function(z) z[length(z)], character(1))
    lod_chr <- format(as.Date(dat[[lod_col]][ok]), "%Y-%m-%d")
    bad_lod <- ok[is.na(lod_chr) | last_date != lod_chr]
  }
  bad_cnt <- same_day_dup <- integer(0)
  if (!is.na(cnt_col)) {
    cnt0 <- safe_num(dat[[cnt_col]]); cnt0[is.na(cnt0)] <- 0
    bad_cnt <- which(nd > cnt0)
    same_day_dup <- which(nd > 0 & nd < cnt0)
  }
  data.table(trait = trait, src = src, aod_col = aod_col, aod_type = "aodtxt", n_total = nrow(dat), n_any = sum(nd > 0), n_single = sum(nd == 1), n_multi = sum(nd > 1), pct_multi_among_any = ifelse(sum(nd > 0) > 0, sum(nd > 1) / sum(nd > 0), NA_real_), max_dates = max(nd), bad_fod = length(bad_fod), bad_lod = length(bad_lod), bad_n_dates_gt_cnt = length(bad_cnt), same_day_duplicate_code_records = length(same_day_dup), fod_col = fod_col, lod_col = lod_col, cnt_col = cnt_col)
}

aodtxt_distribution_trait <- function(dat, trait, src = "icd10") {
  single_col <- single_event_date_col_for_trait(dat, trait)
  if (!is.na(single_col)) {
    nd <- as.integer(!is.na(as_date2(dat[[single_col]])))
    return(data.table(trait = trait, n_dates = nd)[, .N, by = .(trait, n_dates)][order(trait, n_dates)])
  }
  aod_col <- aodtxt_col_for_trait(dat, trait, src)
  if (is.na(aod_col)) return(data.table())
  nd <- n_aod_dates(dat[[aod_col]])
  data.table(trait = trait, n_dates = nd)[, .N, by = .(trait, n_dates)][order(trait, n_dates)]
}

aodtxt_prepost_trait <- function(dat, trait, src = "icd10") {
  if (!"date_attend" %in% names(dat)) return(data.table())
  b0 <- as_date2(dat$date_attend)
  single_col <- single_event_date_col_for_trait(dat, trait)
  if (!is.na(single_col)) {
    d <- as_date2(dat[[single_col]])
    dt <- as.numeric(d - b0) / 365.25
    pattern <- ifelse(is.na(d) | is.na(b0), NA_character_, ifelse(dt < 0, "pre only", ifelse(dt > 0, "post only", "baseline day")))
    out <- data.table(trait = trait, pattern = pattern)
    return(out[!is.na(pattern), .N, by = .(trait, pattern)])
  }
  aod_col <- aodtxt_col_for_trait(dat, trait, src)
  if (is.na(aod_col)) return(data.table())
  aods <- aod_dates_from_text_vector(dat[[aod_col]])
  pattern <- rep(NA_character_, nrow(dat))
  for (i in seq_len(nrow(dat))) {
    d <- aods[[i]]
    if (!length(d) || is.na(b0[i])) next
    dt <- as.numeric(d - b0[i]) / 365.25
    has_pre <- any(dt < 0, na.rm = TRUE)
    has_post <- any(dt > 0, na.rm = TRUE)
    has_base <- any(dt == 0, na.rm = TRUE)
    pattern[i] <- if (has_pre && has_post) "pre + post" else if (has_pre) "pre only" else if (has_post) "post only" else if (has_base) "baseline day" else "unknown"
  }
  out <- data.table(trait = trait, pattern = pattern)
  out[!is.na(pattern), .N, by = .(trait, pattern)]
}

aodtxt_span_trait <- function(dat, trait, src = "icd10") {
  aod_col <- aodtxt_col_for_trait(dat, trait, src)
  if (is.na(aod_col)) return(data.table())
  aods <- aod_dates_from_text_vector(dat[[aod_col]])
  span <- vapply(aods, function(d) if (length(d) > 1) as.numeric(max(d) - min(d)) / 365.25 else NA_real_, numeric(1))
  nd <- lengths(aods)
  data.table(trait = trait, n_dates = nd, span_year = span)[is.finite(span_year)]
}

aodtxt_examples_trait <- function(dat, trait, src = "icd10", n = 30) {
  aod_col <- aodtxt_col_for_trait(dat, trait, src)
  if (is.na(aod_col)) return(data.table())
  nd <- n_aod_dates(dat[[aod_col]])
  idx <- which(nd > 1)
  if (!length(idx)) return(data.table())
  fod_col <- first_present(dat, c(paste0("fod_", src, "_", trait), paste0("fod_", trait)))
  lod_col <- first_present(dat, c(paste0("lod_", src, "_", trait), paste0("lod_", trait)))
  cnt_col <- first_present(dat, c(paste0("cnt_", src, "_", trait), paste0("cnt_", trait)))
  cols <- intersect(c("eid", fod_col, lod_col, cnt_col, aod_col), names(dat))
  out <- as.data.table(dat[idx, cols, drop = FALSE])
  out[, trait := trait]
  out[, n_aod_dates := nd[idx]]
  out[order(-n_aod_dates)][seq_len(min(n, .N))]
}

make_root_aod_qc_figure <- function() {
  set_scope("root")
  dat0 <- readRDS(file.path(indir, "Rdata", "all.rds")) %>% as.data.frame()
  if ("date_attend" %in% names(dat0)) dat0$date_attend <- as_date2(dat0$date_attend)
  if ("date_death" %in% names(dat0)) dat0$date_death <- as_date2(dat0$date_death)
  traits <- parse_vec(Sys.getenv("YY_AOD_QC_TRAITS", unset = paste(Ys, "mi,cvd_cad,cvd_hfail,cvd_afib,t2dm,ckd,death", collapse = ",")), default = c("mi", "cvd_cad", "cvd_hfail", "cvd_afib", "t2dm", "ckd", "death"))
  traits <- unique(traits)
  src <- Sys.getenv("YY_AOD_QC_SRC", unset = "icd10")
  audit <- rbindlist(lapply(traits, function(tr) audit_aodtxt_trait(dat0, tr, src)), fill = TRUE)
  dist <- rbindlist(lapply(traits, function(tr) aodtxt_distribution_trait(dat0, tr, src)), fill = TRUE)
  prepost <- rbindlist(lapply(traits, function(tr) aodtxt_prepost_trait(dat0, tr, src)), fill = TRUE)
  spans <- rbindlist(lapply(traits, function(tr) aodtxt_span_trait(dat0, tr, src)), fill = TRUE)
  examples <- rbindlist(lapply(traits, function(tr) aodtxt_examples_trait(dat0, tr, src, n = 20)), fill = TRUE)
  if (!nrow(audit) || all(audit$n_any == 0, na.rm = TRUE)) {
    warning("No aodtxt or single-event date columns found for root AOD QC.", call. = FALSE)
    write_book(list(aod_audit = audit), fig_file(3, paste0(biom, ".aodtxt_QC"), ".out.xlsx"))
    return(invisible(list(audit = audit, distribution = dist, prepost = prepost, spans = spans, examples = examples)))
  }
  audit[, trait := factor(trait, levels = traits)]
  if (nrow(dist)) {
    dist[, trait := factor(trait, levels = traits)]
    dist[, n_dates_group := factor(ifelse(n_dates >= 6, "6+", as.character(n_dates)), levels = c("0", "1", "2", "3", "4", "5", "6+"))]
    dist2 <- dist[, .(N = sum(N)), by = .(trait, n_dates_group)]
  } else dist2 <- data.table(trait = factor(character(), levels = traits), n_dates_group = factor(character(), levels = c("0", "1", "2", "3", "4", "5", "6+")), N = integer())
  if (nrow(prepost)) {
    prepost[, trait := factor(trait, levels = traits)]
    prepost[, pattern := factor(pattern, levels = c("pre only", "post only", "pre + post", "baseline day", "unknown"))]
  }
  pA <- ggplot(audit, aes(trait, n_any, fill = trait)) + geom_col(width = .68, color = "white") + geom_text(aes(label = comma(n_any)), vjust = -.30, fontface = "bold", size = 3.0) + scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .18))) + labs(title = "a. Participants with any occurrence date", subtitle = paste0("Source: aodtxt_", src, "_*; death is treated as a single-event endpoint"), x = NULL, y = "N with any date") + theme_yy(10) + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
  pB <- ggplot(audit, aes(trait, n_multi, fill = trait)) + geom_col(width = .68, color = "white") + geom_text(aes(label = comma(n_multi)), vjust = -.30, fontface = "bold", size = 3.0) + scale_y_continuous(labels = comma, expand = expansion(mult = c(0, .18))) + labs(title = "b. Multi-date AOD records", subtitle = "Single-event endpoints should have zero multi-date AOD", x = NULL, y = "N with >1 unique date") + theme_yy(10) + theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
  pC <- ggplot(dist2[n_dates_group != "0"], aes(n_dates_group, N, fill = n_dates_group)) + geom_col(width = .72, color = "white") + facet_wrap(~ trait, scales = "free_y") + scale_y_continuous(labels = comma) + labs(title = "c. Distribution of unique occurrence dates among cases", x = "Unique AOD dates", y = "Participants") + theme_yy(9) + theme(legend.position = "none")
  pD <- if (nrow(prepost)) ggplot(prepost, aes(trait, N, fill = pattern)) + geom_col(width = .72, color = "white", position = "fill") +
    scale_y_continuous(labels = percent) + labs(title = "d. Timing of occurrence dates relative to baseline", subtitle = "Composition among AOD-positive records; this replaces the less informative same-day duplicate-code bar", x = NULL, y = "Proportion among cases", fill = NULL) +
    theme_yy(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom") else ggplot() + labs(title = "d. Pre/post timing unavailable") + theme_yy(9)
  Fig <- (pA | pB) / (pC | pD)
  save_plot(Fig, fig_file(3, paste0(biom, ".aodtxt_QC")), 16, 11)
  audit_long <- melt(audit[, .(trait, bad_fod, bad_lod, bad_n_dates_gt_cnt, same_day_duplicate_code_records)], id.vars = "trait", variable.name = "check", value.name = "N")
  write_tab(audit, fig_data_file(3, paste0(biom, ".aodtxt_QC.audit.tsv")))
  if (nrow(dist)) write_tab(dist, fig_data_file(3, paste0(biom, ".aodtxt_QC.distribution.tsv")))
  if (nrow(prepost)) write_tab(prepost, fig_data_file(3, paste0(biom, ".aodtxt_QC.prepost_timing.tsv")))
  if (nrow(spans)) write_tab(spans, fig_data_file(3, paste0(biom, ".aodtxt_QC.multi_date_spans.tsv.gz")))
  if (nrow(examples)) write_tab(examples, fig_data_file(3, paste0(biom, ".aodtxt_QC.multi_date_examples.tsv.gz")))
  write_book(list(aod_audit = audit, aod_distribution = dist, prepost_timing = prepost, multi_date_spans = spans, internal_consistency = audit_long, multi_date_examples = examples), fig_file(3, paste0(biom, ".aodtxt_QC"), ".out.xlsx"))
  invisible(list(audit = audit, distribution = dist, prepost = prepost, spans = spans, examples = examples))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Disease-specific figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_yy_timeline_panel <- function(datY, Y, proteins = vip_proteins, z_suffix = "fig") {
  proteins <- intersect(proteins, names(datY))
  z_cols <- paste0(proteins, "__", z_suffix)
  keep <- z_cols %in% names(datY)
  proteins <- proteins[keep]; z_cols <- z_cols[keep]
  breaks <- seq(-16, 16, 2)
  cases <- datY %>% filter(is.finite(b2e), b2e >= min(breaks), b2e <= max(breaks))
  if (!nrow(cases) || !length(proteins)) return(list(plot = ggplot() + labs(title = "YY timeline unavailable") + theme_yy(), bin = data.table(), trend = data.table()))
  mids <- head(breaks, -1) + diff(breaks) / 2
  hist_df <- as.data.table(hist(cases$b2e, breaks = breaks, plot = FALSE)[c("mids", "counts")])
  setnames(hist_df, c("mid", "count"))
  hist_df[, period := fifelse(mid < 0, "Pre-baseline", "Post-baseline")]
  bm <- bin_protein_means(cases, "b2e", z_cols, proteins, breaks = breaks, skip = c(1, 0), min_n = min_bin_n)
  tr <- bm %>% group_by(protein) %>% group_modify(~ trend_bin_lm(.x)) %>% ungroup()
  yh <- predict_trend_lines(bm) %>% left_join(tr %>% dplyr::select(protein, p), by = "protein") %>% filter(is.finite(p), p < trend_p_show)
  yr <- range(bm$mean_z, yh$yhat, 0, na.rm = TRUE)
  if (!all(is.finite(yr))) yr <- c(-1, 1)
  yr <- range(pretty(yr + c(-.08, .08), n = 5))
  plot_max <- max(hist_df$count, na.rm = TRUE) * 1.12
  scale_factor <- plot_max / diff(yr)
  z_to_count <- function(z) (z - yr[1]) * scale_factor
  bm <- bm %>% mutate(y_plot = z_to_count(mean_z))
  yh <- yh %>% mutate(y_plot = z_to_count(yhat))
  tr_lab <- tr %>% filter(is.finite(p), p < trend_p_show) %>% arrange(p) %>% slice_head(n = 5) %>%
    mutate(label = sprintf("%s\nbeta = %.3f\np = %.2E", protein, beta10, p),
           x = min(breaks) + 0.9,
           y_z = seq(yr[2] - .08 * diff(yr), yr[2] - .70 * diff(yr), length.out = n()),
           y_plot = z_to_count(y_z))

  p <- ggplot() +
    geom_col(data = hist_df, aes(mid, count, fill = period), width = diff(breaks)[1] * .98, color = "white", alpha = .55) +
    geom_hline(yintercept = z_to_count(0), linetype = "dashed", linewidth = .35, color = "grey35") +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = .45, color = "grey25") +
    geom_line(data = bm, aes(mid, y_plot, color = protein, group = protein), linewidth = .75, na.rm = TRUE) +
    geom_point(data = bm, aes(mid, y_plot, color = protein), size = 2.1, na.rm = TRUE) +
    geom_line(data = yh, aes(mid, y_plot, color = protein, group = protein), linewidth = .55, linetype = "dashed", na.rm = TRUE) +
    geom_text(data = tr_lab, aes(x, y_plot, label = label, color = protein), hjust = 0, vjust = 1, size = 3.1, lineheight = .95, show.legend = FALSE) +
    scale_fill_manual(values = c(`Pre-baseline` = "grey70", `Post-baseline` = "#8EC5FF"), guide = "none") +
    scale_color_manual(values = protein_colors(proteins), breaks = proteins) +
    scale_x_continuous(limits = range(breaks), breaks = seq(-15, 15, 5)) +
    scale_y_continuous(name = "Patient count", limits = c(0, plot_max), sec.axis = sec_axis(~ . / scale_factor + yr[1], name = "Biomarker level (standardized)")) +
    labs(title = paste0("a. Yin-Yang ", biom_label, " timeline for ", Y), subtitle = paste0("Biomarker z-score mode: ", timeline_adjust, "; reference: ", fig_standardize), x = "Time from baseline to diagnosis (year)", color = NULL) +
    theme_yy(11) + theme(legend.position = c(.84, .86), legend.background = element_rect(fill = alpha("white", .70), color = "grey80"), legend.text = element_text(face = "bold"))
  list(plot = p, bin = as.data.table(bm), trend = as.data.table(tr))
}

plot_trend_forest <- function(tr, title = "b. YY trend test") {
  if (!nrow(tr)) return(ggplot() + labs(title = title) + theme_yy())
  tr <- as.data.table(tr)
  tr[, `:=`(lo = beta10 - 1.96 * se10, hi = beta10 + 1.96 * se10)]
  tr[, protein := factor(protein, levels = rev(protein[order(beta10)]))]
  ggplot(tr, aes(beta10, protein, color = protein)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .18, linewidth = .7, na.rm = TRUE) +
    geom_point(size = 2.5, na.rm = TRUE) +
    scale_color_manual(values = protein_colors(as.character(tr$protein)), guide = "none") +
    labs(title = title, subtitle = "Beta per 10-year later diagnosis date", x = "YY trend beta", y = NULL) + theme_yy(10)
}

plot_counts_panel <- function(datY, Y) {
  breaks <- seq(-16, 16, 2)
  d <- datY %>% filter(is.finite(b2e), b2e >= min(breaks), b2e <= max(breaks)) %>% mutate(bin = cut(b2e, breaks, include.lowest = TRUE), mid = head(breaks, -1)[as.integer(bin)] + 1)
  if (!nrow(d)) return(ggplot() + labs(title = "c. Event counts unavailable") + theme_yy())
  tab <- d %>% count(bin, name = "N") %>% mutate(mid = head(breaks, -1)[as.integer(bin)] + 1, period = if_else(mid < 0, "Pre-baseline", "Post-baseline"))
  ggplot(tab, aes(mid, N, fill = period)) + geom_col(width = 1.85, color = "white") +
    geom_vline(xintercept = 0, linetype = "dotted") + scale_fill_manual(values = c(`Pre-baseline` = "grey70", `Post-baseline` = "#8EC5FF"), guide = "none") +
    scale_y_continuous(labels = comma) + labs(title = paste0("c. ", Y, " cases by timeline bin"), x = "Time from baseline to diagnosis (year)", y = "Cases") + theme_yy(10)
}

plot_key_pair_panel <- function(bin_df) {
  d <- as.data.table(bin_df)
  focus <- intersect(c("PCSK9", "GDF15", "ANGPTL3", "APOE"), unique(d$protein))
  if (length(focus) < 2) focus <- head(unique(d$protein), min(4L, uniqueN(d$protein)))
  d <- d[protein %in% focus]
  if (!nrow(d)) return(ggplot() + labs(title = "d. Key biomarker trajectories unavailable") + theme_yy())
  ggplot(d, aes(mid, mean_z, color = protein, group = protein)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + geom_vline(xintercept = 0, linetype = "dotted", color = "grey35") +
    geom_line(linewidth = .8, na.rm = TRUE) + geom_point(size = 2, na.rm = TRUE) +
    scale_color_manual(values = protein_colors(unique(d$protein))) +
    labs(title = "d. Key trajectories on biomarker scale", x = "Time from baseline to diagnosis (year)", y = "Biomarker z-score", color = NULL) + theme_yy(10) + theme(legend.position = "top")
}


choose_yy_plot_proteins <- function(datY, Y) {
  set_scope("trait", Y)
  # Default behavior: for small candidate mode, use the curated proteins. In full
  # scan mode, rank all available PPP proteins by incident pWAS P value and show
  # only a small number in Fig1/Fig2 timeline panels.
  if (!isTRUE(scan_all_proteins)) return(intersect(vip_proteins, names(datY)))
  candidates <- intersect(vars_biom_global, names(datY))
  if (max_scan_proteins > 0 && length(candidates) > max_scan_proteins) candidates <- head(candidates, max_scan_proteins)
  if (!length(candidates)) return(intersect(vip_proteins, names(datY)))
  base_covs <- default_demog_covars(datY, include_tdi = TRUE)
  pwas <- rbindlist(lapply(candidates, function(p) {
    out <- cox_one(datY, "t2e_y", "Yt2e_y", p, base_covs, min_n = 500, min_event = 20)
    data.table(protein = p, beta = out$beta, se = out$se, p = out$p, HR = out$HR, lo = out$lo, hi = out$hi, n = out$n, event = out$event)
  }), fill = TRUE)
  if (nrow(pwas)) {
    pwas <- merge(pwas, annotate_proteins(pwas$protein), by = "protein", all.x = TRUE, sort = FALSE)
    pwas <- pwas[order(p, na.last = TRUE)]
    write_tab(pwas, fig_data_file(0, "incident_pwas_for_yy_plot.tsv.gz", Y = Y))
  }
  top <- pwas[is.finite(p)][order(p)][seq_len(min(yy_plot_top_n, .N))]$protein
  top <- unique(c(intersect(yy_plot_always, candidates), top))
  top <- head(top, yy_plot_top_n)
  if (length(top) < 2) top <- intersect(vip_proteins, names(datY))
  top
}

make_yy_timeline_figure <- function(datY, Y, plot_proteins = NULL) {
  set_scope("trait", Y)
  if (is.null(plot_proteins) || !length(plot_proteins)) plot_proteins <- vip_proteins
  proteins <- intersect(plot_proteins, names(datY))
  cases <- datY %>% filter(is.finite(b2e), b2e >= -16, b2e <= 16)
  ref_idx <- if (tolower(fig_standardize) == "control") datY$never_y else rep(FALSE, nrow(datY))
  if (tolower(fig_standardize) == "case") ref_idx[match(cases$eid, datY$eid)] <- TRUE
  timeline_covars <- get_timeline_covars(datY, timeline_adjust)
  datY <- make_protein_z(datY, proteins, covars = timeline_covars, ref_idx = ref_idx, suffix = "fig")
  yy <- make_yy_timeline_panel(datY, Y, proteins = proteins, z_suffix = "fig")
  pB <- plot_trend_forest(yy$trend, "b. Weighted bin-level trend test")
  pC <- plot_counts_panel(datY, Y)
  pD <- plot_key_pair_panel(yy$bin)
  Fig <- yy$plot / (pB | pC | pD) + plot_layout(heights = c(1.4, 1))
  save_plot(Fig, fig_file(1, "yy_timeline", Y = Y), 16, 12)
  write_tab(yy$bin, fig_data_file(1, "yy_timeline.bins.tsv.gz", Y = Y))
  write_tab(yy$trend, fig_data_file(1, "yy_timeline.trend.tsv", Y = Y))
  write_book(list(binned_means = yy$bin, trend = yy$trend, timeline_covariates = data.frame(covariate = timeline_covars)), fig_file(1, "yy_timeline", ".out.xlsx", Y = Y))
  invisible(list(dat = datY, yy = yy))
}

make_aod_recurrence_figure <- function(datY, Y) {
  set_scope("trait", Y)
  proteins <- intersect(vip_proteins, names(datY))
  datY <- make_protein_z(datY, proteins, covars = character(), ref_idx = datY$never_y, suffix = "ctrl")
  z_cols <- paste0(proteins, "__ctrl")
  aod_source <- unique(as.character(datY$occ_source))[1] %||% "unknown"
  aod_real <- any(datY$occ_is_real_aod %in% TRUE, na.rm = TRUE)
  aod_single_event <- any(datY$occ_is_single_event %in% TRUE, na.rm = TRUE)
  xvars <- if (isTRUE(aod_single_event)) {
    c(b2e = "Single-event date alignment")
  } else if (isTRUE(aod_real)) {
    c(b2e = "FOD alignment", occ_last_pre_b2e = "Last pre-baseline occurrence", occ_nearest_b2e = "Nearest occurrence", occ_first_post_b2e = "First post-baseline occurrence")
  } else {
    c(b2e = "FOD alignment", occ_nearest_b2e = "FOD fallback: nearest=FOD", occ_last_pre_b2e = "FOD fallback: pre-baseline only")
  }
  bmall <- rbindlist(lapply(names(xvars), function(xv) {
    if (!xv %in% names(datY)) return(data.table())
    dd <- datY %>% filter(is.finite(.data[[xv]]), .data[[xv]] >= -16, .data[[xv]] <= 16)
    as.data.table(bin_protein_means(dd, xv, z_cols, proteins, breaks = seq(-16, 16, 2), skip = c(1, 0), min_n = min_bin_n))[, alignment := xvars[[xv]]][]
  }), fill = TRUE)
  tr <- bmall %>% group_by(alignment, protein) %>% group_modify(~ trend_bin_lm(.x)) %>% ungroup() %>% as.data.table()
  occ_count_tbl <- as.data.table(datY)[, .N, by = .(occ_count_group, occ_prepost_pattern)][order(occ_count_group, occ_prepost_pattern)]
  occ_state_tbl <- as.data.table(datY)[, .N, by = occ_state][order(occ_state)]
  group_sm <- rbindlist(lapply(proteins, function(protein_name) {
    zc <- paste0(protein_name, "__ctrl")
    datY %>% transmute(protein = .env$protein_name, occ_state, z = .data[[zc]]) %>% filter(!is.na(occ_state), is.finite(z)) %>%
      group_by(protein, occ_state) %>% summarise(n = n(), mean_z = mean(z), se = sd(z) / sqrt(n), lo = mean_z - 1.96 * se, hi = mean_z + 1.96 * se, .groups = "drop") %>% as.data.table()
  }), fill = TRUE)
  active_contrast <- rbindlist(lapply(proteins, function(protein_name) {
    zc <- paste0(protein_name, "__ctrl")
    covs0 <- default_demog_covars(datY, include_tdi = TRUE)
    rbindlist(lapply(c("Prevalent single-date remote", "Prevalent multi-date remote", "Prevalent recent/active", "Prevalent + post-baseline recurrence", "Incident multi-date"), function(grp) {
      dd <- datY %>% filter(occ_state %in% c("Never disease", grp), is.finite(.data[[zc]])) %>% mutate(x = as.integer(occ_state == grp))
      lm_contrast_one(dd, zc, "x", covs0) %>% mutate(protein = .env$protein_name, contrast = paste0(grp, " vs never"))
    }), fill = TRUE)
  }), fill = TRUE)
  burden_contrast <- rbindlist(lapply(proteins, function(protein_name) {
    zc <- paste0(protein_name, "__ctrl")
    dd <- datY %>% filter(!is.na(b2e), is.finite(.data[[zc]]), occ_n_all > 0) %>% mutate(occ_n_log1p = log1p(occ_n_all), occ_span_log1p = log1p(pmax(occ_span_year, 0)))
    burden_covs <- unique(c(default_demog_covars(dd, include_tdi = TRUE), "b2e"))
    r1 <- lm_contrast_one(dd, zc, "occ_n_log1p", burden_covs) %>% mutate(protein = .env$protein_name, contrast = "per log1p(unique AOD dates)")
    r2 <- lm_contrast_one(dd, zc, "occ_span_log1p", burden_covs) %>% mutate(protein = .env$protein_name, contrast = "per log1p(AOD span years)")
    bind_rows(r1, r2)
  }), fill = TRUE)

  aod_source_table <- data.frame(Y = Y, occ_source = aod_source, occ_is_real_aod = aod_real, occ_is_single_event = aod_single_event, n_any = sum(datY$occ_n_all > 0), n_multi = sum(datY$occ_n_all > 1), max_dates = max(datY$occ_n_all), stringsAsFactors = FALSE)
  source_summary <- data.table(
    metric = c("source", "real AOD", "any occurrence", "multi-date", "max dates/person"),
    value = c(aod_source, ifelse(aod_real, "yes", "no"), comma(aod_source_table$n_any), comma(aod_source_table$n_multi), comma(aod_source_table$max_dates))
  )
  source_summary[, metric := factor(metric, levels = rev(metric))]
  pA <- ggplot(occ_count_tbl, aes(occ_count_group, N, fill = occ_prepost_pattern)) +
    geom_col(color = "white", width = .72) +
    scale_y_continuous(labels = comma) +
    labs(title = "a. Occurrence-date burden",
         subtitle = "Number of unique dates and whether dates are before/after baseline",
         x = "Unique occurrence dates", y = "Participants", fill = NULL) +
    theme_yy(9)
  pB <- ggplot(occ_state_tbl, aes(N, fct_reorder(as.character(occ_state), N))) +
    geom_col(width = .72, fill = "#6AA5FF", color = "white") +
    geom_text(aes(label = comma(N)), hjust = -0.08, size = 2.8) +
    scale_x_continuous(labels = comma, expand = expansion(mult = c(0, .16))) +
    labs(title = "b. Yin/Yang occurrence states", x = "Participants", y = NULL) +
    theme_yy(9)
  pC <- ggplot(source_summary, aes(1, metric)) +
    geom_tile(fill = "#F2F4F7", color = "white", width = .98, height = .82) +
    geom_text(aes(label = paste0(as.character(metric), ": ", value)), hjust = 0, nudge_x = -.44, size = 3.1) +
    coord_cartesian(xlim = c(.52, 1.48), clip = "off") +
    labs(title = ifelse(aod_single_event, "c. Single-event endpoint", "c. AOD source summary"), x = NULL, y = NULL) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0))
  tr[, `:=`(lo = beta10 - 1.96 * se10, hi = beta10 + 1.96 * se10)]
  pD <- ggplot(group_sm, aes(occ_state, mean_z, ymin = lo, ymax = hi, color = protein)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3) + geom_pointrange(position = position_dodge(.65), linewidth = .45, size = 1.25, na.rm = TRUE) +
    facet_wrap(~ protein, ncol = 3, scales = "free_y") + scale_color_manual(values = protein_colors(proteins), guide = "none") +
    labs(title = "d. Baseline protein levels by occurrence-state", x = NULL, y = "Protein z-score") + theme_yy(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7.5))
  active_contrast[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]
  active_contrast[, contrast_short := fcase(
    grepl("^Incident multi-date", contrast), "Incident multi-date",
    grepl("^Prevalent \\+ post-baseline recurrence", contrast), "Prevalent + recurrence",
    grepl("^Prevalent multi-date remote", contrast), "Prevalent multi-date remote",
    grepl("^Prevalent recent/active", contrast), "Prevalent recent/active",
    grepl("^Prevalent single-date remote", contrast), "Prevalent single-date remote",
    default = contrast
  )]
  pE <- ggplot(active_contrast, aes(beta, protein, xmin = lo, xmax = hi, color = contrast_short)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") + geom_errorbar(orientation = "y", width = .18, linewidth = .55, na.rm = TRUE) + geom_point(size = 1.8, na.rm = TRUE) +
    labs(title = "e. Occurrence-state contrasts", subtitle = "Adjusted linear contrast against never disease", x = "Adjusted difference in protein z-score", y = NULL, color = NULL) + theme_yy(9) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7.2))
  burden_contrast[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]
  pF <- ggplot(burden_contrast, aes(beta, protein, xmin = lo, xmax = hi, color = contrast)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") + geom_errorbar(orientation = "y", width = .18, linewidth = .55, na.rm = TRUE) + geom_point(size = 1.8, na.rm = TRUE) +
    labs(title = "f. Within-case AOD burden", subtitle = "Adjusted for b2e and baseline covariates", x = "Adjusted protein difference", y = NULL, color = NULL) + theme_yy(9) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7.5))
  Fig <- (pA | pB | pC) / (pD | pE) / pF + plot_layout(heights = c(.72, 1.15, .78))
  save_plot(Fig, fig_file(2, "aod_recurrence", Y = Y), 17, 13.5)
  write_tab(bmall, fig_data_file(2, "aod_alignment.bins.tsv.gz", Y = Y))
  write_tab(tr, fig_data_file(2, "aod_alignment.trend.tsv", Y = Y))
  write_tab(occ_count_tbl, fig_data_file(2, "occurrence_count.tsv", Y = Y))
  write_tab(occ_state_tbl, fig_data_file(2, "occurrence_state_counts.tsv", Y = Y))
  write_tab(group_sm, fig_data_file(2, "occurrence_state_protein.tsv", Y = Y))
  write_tab(active_contrast, fig_data_file(2, "occurrence_state_contrast.tsv", Y = Y))
  write_tab(burden_contrast, fig_data_file(2, "aod_burden_contrast.tsv", Y = Y))
  write_tab(aod_source_table, fig_data_file(2, "aod_source.tsv", Y = Y))
  write_book(list(aod_source = aod_source_table, aod_bins = bmall, trend = tr, occurrence_counts = occ_count_tbl, occurrence_state_counts = occ_state_tbl, occurrence_state_protein = group_sm, occurrence_state_contrast = active_contrast, aod_burden_contrast = burden_contrast), fig_file(2, "aod_recurrence", ".out.xlsx", Y = Y))
  invisible(list(dat = datY, aod_bins = bmall, trend = tr, group = group_sm, active_contrast = active_contrast, burden_contrast = burden_contrast))
}

make_adjustment_figure <- function(datY, Y) {
  set_scope("trait", Y)
  proteins <- intersect(vip_proteins, names(datY))
  basic_no_tdi <- default_demog_covars(datY, include_tdi = FALSE)
  demog <- default_demog_covars(datY, include_tdi = TRUE)
  le8 <- get_le8_covars(datY)
  drugs <- c("drug.lipid", "drug.dm", "drug.htn", "drug.aspirin")
  genetic <- grep("cad.*prs|cad\\.6m|apoe|PCSK9|ANGPTL3|APOB|LDLR|rare|lof|pLoF", names(datY), value = TRUE, ignore.case = TRUE)
  genetic <- setdiff(genetic[!grepl("__", genetic)], proteins)
  adj_sets <- list(
    raw = character(),
    age_sex_PC = basic_no_tdi,
    plus_le8 = c(demog, le8),
    plus_drugs = c(demog, le8, drugs),
    plus_genetic = c(demog, le8, genetic),
    full_available = c(demog, le8, drugs, genetic)
  )
  ref_idx <- datY$never_y
  trend_res <- data.table(); cox_res <- data.table(); dat_adj <- datY
  for (nm in names(adj_sets)) {
    dat_adj <- make_protein_z(dat_adj, proteins, covars = adj_sets[[nm]], ref_idx = ref_idx, suffix = paste0("adj_", nm))
    z_cols <- paste0(proteins, "__adj_", nm)
    cases <- dat_adj %>% filter(is.finite(b2e), b2e >= -16, b2e <= 16)
    bm <- bin_protein_means(cases, "b2e", z_cols, proteins, breaks = seq(-16, 16, 2), skip = c(1, 0), min_n = min_bin_n)
    tt <- bm %>% group_by(protein) %>% group_modify(~ trend_bin_lm(.x)) %>% ungroup() %>% mutate(adjustment = nm) %>% as.data.table()
    trend_res <- rbind(trend_res, tt, fill = TRUE)
    cc <- rbindlist(lapply(seq_along(proteins), function(i) {
      cox_one(dat_adj, "t2e_y", "Yt2e_y", z_cols[i], safe_covars(dat_adj, adj_sets[[nm]])) %>% mutate(protein = proteins[i], adjustment = nm)
    }), fill = TRUE)
    cox_res <- rbind(cox_res, cc, fill = TRUE)
  }
  trend_res[, `:=`(lo = beta10 - 1.96 * se10, hi = beta10 + 1.96 * se10, adjustment = factor(adjustment, levels = names(adj_sets)))]
  cox_res[, `:=`(lo_beta = beta - 1.96 * se, hi_beta = beta + 1.96 * se, adjustment = factor(adjustment, levels = names(adj_sets)))]
  pA <- ggplot(trend_res, aes(adjustment, beta10, ymin = lo, ymax = hi, color = protein, group = protein)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + geom_pointrange(position = position_dodge(.65), linewidth = .45, size = 1.6, na.rm = TRUE) +
    facet_wrap(~ protein, ncol = 3, scales = "free_y") + scale_color_manual(values = protein_colors(proteins), guide = "none") + labs(title = "a. YY trend after sequential adjustment", x = NULL, y = "Trend beta") + theme_yy(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  pB <- ggplot(cox_res, aes(adjustment, beta, ymin = lo_beta, ymax = hi_beta, color = protein, group = protein)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + geom_pointrange(position = position_dodge(.65), linewidth = .45, size = 1.6, na.rm = TRUE) +
    facet_wrap(~ protein, ncol = 3, scales = "free_y") + scale_color_manual(values = protein_colors(proteins), guide = "none") + labs(title = "b. Incident-disease Cox beta after adjustment", x = NULL, y = "Cox beta") + theme_yy(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  drug_tbl <- data.table()
  if ("drug.lipid" %in% names(dat_adj)) {
    drug_tbl <- rbindlist(lapply(intersect(c("PCSK9", "GDF15", "ANGPTL3", "APOE"), proteins), function(p) {
      zc <- paste0(p, "__adj_raw")
      dat_adj %>% filter(yy_group %in% c("Never disease", "Incident 0-5y", "Prevalent remote only", "Prevalent recent/active"), is.finite(.data[[zc]])) %>%
        mutate(drug_lipid = factor(if_else(safe_num(drug.lipid) == 1, "lipid drug", "no lipid drug"), levels = c("no lipid drug", "lipid drug"))) %>%
        group_by(protein = .env$p, yy_group, drug_lipid) %>% summarise(n = n(), mean_z = mean(.data[[zc]], na.rm = TRUE), .groups = "drop") %>% as.data.table()
    }), fill = TRUE)
  }
  pC <- if (nrow(drug_tbl)) ggplot(drug_tbl, aes(yy_group, mean_z, fill = drug_lipid)) + geom_col(position = position_dodge(.72), width = .68, color = "white") +
    facet_wrap(~ protein, ncol = 2, scales = "free_y") + scale_fill_manual(values = c(`no lipid drug` = "grey70", `lipid drug` = ORANGE)) + labs(title = "c. Lipid-lowering drug stratification", x = NULL, y = "Mean protein z-score", fill = NULL) + theme_yy(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1)) else ggplot() + labs(title = "c. drug.lipid not available") + theme_yy(9)
  avail <- data.table(set = names(adj_sets), n_cov = sapply(adj_sets, function(v) length(safe_covars(datY, v))))
  pD <- ggplot(avail, aes(factor(set, levels = names(adj_sets)), n_cov, fill = set)) + geom_col(width = .68, color = "white") + geom_text(aes(label = n_cov), vjust = -.35, fontface = "bold") +
    scale_fill_manual(values = setNames(scales::hue_pal()(length(adj_sets)), names(adj_sets)), guide = "none") + labs(title = "d. Available covariates in adjustment sets", x = NULL, y = "N covariates") + theme_yy(9) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  Fig <- (pA / pB) / (pC | pD) + plot_layout(heights = c(1.2, 1.2, .9))
  save_plot(Fig, fig_file(3, "adjustment", Y = Y), 17, 15)
  write_tab(trend_res, fig_data_file(3, "adjustment_trend.tsv", Y = Y))
  write_tab(cox_res, fig_data_file(3, "adjustment_cox.tsv", Y = Y))
  if (nrow(drug_tbl)) write_tab(drug_tbl, fig_data_file(3, "drug_stratified.tsv", Y = Y))
  write_book(list(trend_adjustment = trend_res, cox_adjustment = cox_res, drug_stratified = drug_tbl, covariate_sets = avail), fig_file(3, "adjustment", ".out.xlsx", Y = Y))
  invisible(list(dat = dat_adj, trend = trend_res, cox = cox_res))
}

make_common_outcomes <- function(dat) {
  cand <- parse_vec(Sys.getenv("YY_SPECIFICITY_OUTCOMES", unset = "cvd_cad,t2dm,ckd,cvd_hfail,death"), default = character())
  cand[vapply(cand, function(Y) length(candidate_date_cols(dat, Y, "all")) > 0, logical(1))]
}

prepare_t2e_simple <- function(dat, Y) {
  if (all(c(paste0(Y, ".t2e"), paste0(Y, ".Yt2e")) %in% names(dat))) return(dat)
  cols <- candidate_date_cols(dat, Y, "all")
  if (!length(cols)) return(dat)
  ydate <- paste0(Y, ".yy_specificity_date")
  dat[[ydate]] <- pmin_date2(dat[, cols, drop = FALSE])
  old <- grep(paste0("^", Y, "\\.(Yt2e|Yr2e|t2e|r2e|b2e)$"), names(dat), value = TRUE)
  if (length(old)) dat <- dat[, setdiff(names(dat), old), drop = FALSE]
  t2e(dat, NA, ydate, "birth_date", "date_attend", "date_lost", "date_death", date_follow_end, Y, "year")
}


force_protein_character <- function(x) {
  x <- as.data.table(x)
  if (!"protein" %in% names(x)) x[, protein := character()]
  x[, protein := as.character(protein)]
  x
}

safe_merge_protein <- function(x, y, all = TRUE) {
  x <- force_protein_character(x)
  y <- force_protein_character(y)
  merge(x, y, by = "protein", all = all, sort = FALSE)
}

make_scan_figure <- function(datY, dat0_full, Y) {
  set_scope("trait", Y)
  scan_proteins <- if (scan_all_proteins) intersect(names(datY), vars_biom_global) else unique(c(vip_proteins, extra_proteins))
  scan_proteins <- intersect(scan_proteins, names(datY))
  if (max_scan_proteins > 0 && length(scan_proteins) > max_scan_proteins) scan_proteins <- head(scan_proteins, max_scan_proteins)
  if (!length(scan_proteins)) {
    warning("No proteins available for scan step in ", Y, call. = FALSE)
    return(list(scan = data.table(), specificity = data.table()))
  }
  datY <- make_protein_z(datY, scan_proteins, covars = character(), ref_idx = datY$never_y, suffix = "scan")
  z_cols <- paste0(scan_proteins, "__scan")
  cases <- datY %>% filter(is.finite(b2e), b2e >= -16, b2e <= 16)
  yy_res <- bin_protein_means(cases, "b2e", z_cols, scan_proteins, breaks = seq(-16, 16, 2), skip = c(1, 0), min_n = min_bin_n) %>%
    group_by(protein) %>% group_modify(~ trend_bin_lm(.x)) %>% ungroup() %>% transmute(protein, yy_beta10 = beta10, yy_se10 = se10, yy_p = p) %>% as.data.table()
  yy_res <- force_protein_character(yy_res)
  base_covs <- default_demog_covars(datY, include_tdi = TRUE)
  cox_res <- rbindlist(lapply(seq_along(scan_proteins), function(i) {
    cox_one(datY, "t2e_y", "Yt2e_y", z_cols[i], base_covs) %>% mutate(protein = as.character(scan_proteins[i]))
  }), fill = TRUE) %>% as.data.table()
  cox_res <- force_protein_character(cox_res)
  setnames(cox_res, c("beta", "se", "p", "HR", "lo", "hi", "n", "event"), c("cox_beta", "cox_se", "cox_p", "cox_HR", "cox_lo", "cox_hi", "cox_n", "cox_event"), skip_absent = TRUE)
  scan_res <- safe_merge_protein(yy_res, cox_res, all = TRUE)
  scan_anno <- annotate_proteins(scan_res$protein)[, .(protein, prior, category)]
  scan_res <- safe_merge_protein(scan_res, scan_anno, all = TRUE)
  scan_res[is.na(prior), prior := "other"]
  scan_res[is.na(category), category := "Other/unclassified"]

  scan_res[, specificity_rank_p := pmin(fifelse(is.finite(cox_p), cox_p, 1), fifelse(is.finite(yy_p), yy_p, 1))]
  scan_res[, quadrant := fcase(
    is.finite(yy_beta10) & yy_beta10 < 0 & is.finite(cox_beta) & cox_beta > 0 & yy_p < .05 & cox_p < .05, "reactive/risk pattern",
    is.finite(yy_beta10) & yy_beta10 > 0 & is.finite(cox_beta) & cox_beta < 0 & yy_p < .05 & cox_p < .05, "protective/mirror pattern",
    is.finite(cox_beta) & cox_p < .05 & (!is.finite(yy_p) | yy_p >= .05), "prospective-specific pattern",
    is.finite(yy_beta10) & yy_p < .05 & (!is.finite(cox_p) | cox_p >= .05), "timeline-only",
    is.finite(cox_beta) & cox_p < .05, "prospective-only",
    default = "mixed/uncertain"
  )]

  outcomes <- make_common_outcomes(dat0_full)
  top_proteins <- function(dt, n, by = "specificity_rank_p") {
    dt <- as.data.table(dt)
    if (!nrow(dt) || n < 1) return(character())
    by <- intersect(by, names(dt))
    if (length(by)) data.table::setorderv(dt, by, na.last = TRUE)
    head(dt$protein, n)
  }
  always_specificity <- intersect(unique(c(vip_proteins, yy_plot_always, "PCSK9", "GDF15", "RAB6A")), scan_res$protein)
  spec_features <- unique(c(
    always_specificity,
    top_proteins(scan_res[quadrant == "reactive/risk pattern"], 18),
    top_proteins(scan_res[quadrant == "protective/mirror pattern"], 10),
    top_proteins(scan_res[quadrant %in% c("prospective-specific pattern", "prospective-only")], 12),
    top_proteins(scan_res[quadrant == "timeline-only"], 8),
    top_proteins(scan_res, 24)
  ))
  spec_features <- head(spec_features, 80)
  dat_spec <- make_protein_z(dat0_full, intersect(spec_features, names(dat0_full)), covars = character(), ref_idx = rep(TRUE, nrow(dat0_full)), suffix = "spec")
  spec_res <- rbindlist(lapply(outcomes, function(O) {
    dO <- prepare_t2e_simple(dat_spec, O)
    tcol <- paste0(O, ".t2e"); ecol <- paste0(O, ".Yt2e")
    if (!all(c(tcol, ecol) %in% names(dO))) return(data.table())
    rbindlist(lapply(intersect(spec_features, names(dO)), function(protein_name) {
      cox_one(dO, tcol, ecol, paste0(protein_name, "__spec"), default_demog_covars(dO, include_tdi = TRUE)) %>% mutate(outcome = .env$O, protein = .env$protein_name)
    }), fill = TRUE)
  }), fill = TRUE)
  if (nrow(spec_res)) {
    spec_res <- force_protein_character(spec_res)
    spec_res[, signed_logp := sign(beta) * pmin(-log10(pmax(p, 1e-300)), 20)]
  }

  breadth <- if (nrow(spec_res)) spec_res[, .(outcome_breadth_p05 = sum(is.finite(p) & p < .05), outcome_breadth_fdrlike = sum(is.finite(p) & p < .05 / max(1, length(unique(outcome))))), by = protein] else data.table(protein = character())
  breadth <- force_protein_character(breadth)
  scan_res <- safe_merge_protein(scan_res, breadth, all = TRUE)
  spec_summary <- data.table(protein = character())
  if (nrow(spec_res)) {
    target_spec <- spec_res[outcome == Y, .(
      protein,
      target_spec_beta = beta,
      target_spec_p = p,
      target_spec_logp = pmin(-log10(pmax(p, 1e-300)), 50)
    )]
    off_spec <- spec_res[outcome != Y, .(
      offtarget_max_abs_logp = if (any(is.finite(signed_logp))) max(abs(signed_logp), na.rm = TRUE) else NA_real_,
      offtarget_mean_abs_logp = if (any(is.finite(signed_logp))) mean(abs(signed_logp), na.rm = TRUE) else NA_real_
    ), by = protein]
    spec_summary <- safe_merge_protein(target_spec, off_spec, all = TRUE)
  }
  scan_res <- safe_merge_protein(scan_res, spec_summary, all = TRUE)
  scan_res[!is.finite(target_spec_logp), target_spec_logp := pmin(-log10(pmax(cox_p, 1e-300)), 50)]
  scan_res[, specificity_delta := target_spec_logp - offtarget_max_abs_logp]
  scan_res[, specificity_class := fcase(
    quadrant == "protective/mirror pattern", "Mirror/protective",
    prior == "lipid/CAD-related prior", "Lipid/CAD prior",
    is.finite(outcome_breadth_p05) & outcome_breadth_p05 >= 4, "Broad/systemic",
    is.finite(specificity_delta) & specificity_delta > 1, paste0(Y, " enriched"),
    quadrant == "reactive/risk pattern", "Reactive/risk",
    quadrant %in% c("prospective-specific pattern", "prospective-only"), "Prospective-only",
    quadrant == "timeline-only", "Timeline-only",
    default = "Mixed/uncertain"
  )]

  quadrant_cols <- c(
    `reactive/risk pattern` = "#D95F02",
    `protective/mirror pattern` = "#0072B2",
    `prospective-specific pattern` = "#009E73",
    `prospective-only` = "#E69F00",
    `timeline-only` = "#CC79A7",
    `mixed/uncertain` = "grey60"
  )
  label_keep <- unique(c(
    always_specificity,
    top_proteins(scan_res[quadrant == "reactive/risk pattern"], 18),
    top_proteins(scan_res[quadrant == "protective/mirror pattern"], 8),
    top_proteins(scan_res[quadrant %in% c("prospective-specific pattern", "prospective-only")], 8)
  ))
  scan_res[, scan_label := fifelse(protein %in% label_keep, protein, "")]

  pA <- ggplot(scan_res, aes(yy_beta10, cox_beta, color = quadrant, shape = prior, label = scan_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.15, alpha = .78, na.rm = TRUE) +
    ggrepel::geom_text_repel(size = 2.65, max.overlaps = 34, seed = 9, na.rm = TRUE) +
    scale_color_manual(values = quadrant_cols, na.value = "grey70") +
    labs(title = "a. YY trend vs prospective disease association",
         x = "YY trend beta per 10-year later diagnosis", y = "Prospective Cox beta",
         color = "Pattern", shape = "Prior") +
    theme_yy(10) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7.5))

  spec_plot <- copy(scan_res[protein %in% spec_features & is.finite(target_spec_logp)])
  if (nrow(spec_plot)) {
    spec_plot[, offtarget_plot := fifelse(is.finite(offtarget_max_abs_logp), offtarget_max_abs_logp, 0)]
    spec_label_keep <- unique(c(always_specificity, top_proteins(spec_plot, 14)))
    spec_plot[, spec_label := fifelse(protein %in% spec_label_keep, protein, "")]
  }
  pB <- if (nrow(spec_plot)) {
    ggplot(spec_plot, aes(offtarget_plot, target_spec_logp, color = specificity_class, label = spec_label)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
      geom_point(size = 2.4, alpha = .82, na.rm = TRUE) +
      ggrepel::geom_text_repel(size = 2.7, max.overlaps = 26, seed = 12, na.rm = TRUE) +
      labs(title = "b. Target evidence vs off-target breadth",
           subtitle = "Points above the diagonal are more target-enriched; broad/systemic markers move rightward",
           x = "Max absolute -log10(P) in non-target outcomes", y = paste0(Y, " -log10(P)"),
           color = "Class") +
      theme_yy(10) + theme(legend.position = "bottom")
  } else ggplot() + labs(title = "b. Specificity unavailable") + theme_yy(10)

  quad <- scan_res[, .N, by = .(quadrant, specificity_class)]
  quad[, quadrant := factor(quadrant, levels = names(quadrant_cols))]
  pC <- ggplot(quad, aes(N, fct_rev(quadrant), fill = specificity_class)) +
    geom_col(color = "white") +
    scale_x_continuous(labels = comma) +
    labs(title = "c. Pattern classification", x = "N proteins", y = NULL, fill = "Class") +
    theme_yy(10) + theme(legend.position = "bottom")

  keep_prot <- unique(c(
    always_specificity,
    top_proteins(scan_res[quadrant == "reactive/risk pattern"], 9),
    top_proteins(scan_res[quadrant == "protective/mirror pattern"], 6),
    top_proteins(scan_res[quadrant %in% c("prospective-specific pattern", "prospective-only")], 8),
    top_proteins(scan_res[quadrant == "timeline-only"], 5)
  ))
  triage <- copy(scan_res[protein %in% keep_prot])
  triage[, triage_order := pmin(-log10(pmax(cox_p, 1e-300)), 50) + 0.5 * pmin(-log10(pmax(yy_p, 1e-300)), 50)]
  keep_prot <- triage[order(-triage_order)]$protein
  score_long <- triage[, .(
    protein,
    specificity_class,
    `YY reactive trend` = -sign(yy_beta10) * pmin(-log10(pmax(yy_p, 1e-300)), 10),
    `Incident Cox` = sign(cox_beta) * pmin(-log10(pmax(cox_p, 1e-300)), 10),
    `Target enrichment` = pmax(pmin(specificity_delta, 10), -10),
    `Outcome breadth` = -pmin(outcome_breadth_p05, 10)
  )]
  score_long <- melt(score_long, id.vars = c("protein", "specificity_class"), variable.name = "axis", value.name = "score")
  score_long[, protein := factor(protein, levels = rev(keep_prot))]
  pD <- ggplot(score_long, aes(axis, protein, fill = score)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.finite(score), sprintf("%.1f", score), "")), size = 2.4) +
    scale_fill_gradient2(low = BLUE, mid = "white", high = ORANGE, midpoint = 0, limits = c(-10, 10), oob = scales::squish, na.value = "grey92", name = "score") +
    labs(title = "d. Representative candidate triage",
         subtitle = "Negative breadth means less outcome specificity; negative YY score marks mirror/protective direction",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
  Fig <- (pA | pB) / (pC | pD)
  save_plot(Fig, fig_file(4, "prospective_specificity", Y = Y), 16, 11)
  write_tab(scan_res, fig_data_file(4, "yy_vs_cox_scan.tsv.gz", Y = Y))
  if (nrow(spec_res)) write_tab(spec_res, fig_data_file(4, "outcome_specificity.tsv.gz", Y = Y))
  write_book(list(yy_vs_cox = scan_res, outcome_specificity = spec_res, specificity_summary = spec_summary, pattern_counts = quad), fig_file(4, "prospective_specificity", ".out.xlsx", Y = Y))
  invisible(list(scan = scan_res, specificity = spec_res))
}



fit_aligned_trajectory_spline <- function(dd, grid_step = 0.25, spline_df = 4L) {
  dd <- as.data.table(dd)[is.finite(mid) & is.finite(mean_z) & is.finite(n) & n > 0]
  if (nrow(dd) < 5) return(list(curve = data.table(), features = data.table()))
  df_use <- min(as.integer(spline_df), nrow(dd) - 1L)
  fit_lin <- try(lm(mean_z ~ I(mid / 10), data = dd, weights = n), silent = TRUE)
  fit_spl <- try(lm(mean_z ~ splines::ns(mid, df = df_use), data = dd, weights = n), silent = TRUE)
  if (inherits(fit_spl, "try-error")) return(list(curve = data.table(), features = data.table()))
  grid <- seq(min(dd$mid), max(dd$mid), by = grid_step)
  yhat <- as.numeric(predict(fit_spl, newdata = data.frame(mid = grid)))
  dy <- c(NA_real_, diff(yhat) / diff(grid))
  nonlinear_p <- NA_real_
  if (!inherits(fit_lin, "try-error")) {
    aa <- try(anova(fit_lin, fit_spl), silent = TRUE)
    if (!inherits(aa, "try-error") && nrow(aa) >= 2) nonlinear_p <- aa$`Pr(>F)`[2]
  }
  near0 <- which.min(abs(grid))
  pre_idx <- which(grid < 0); post_idx <- which(grid > 0)
  feat <- data.table(
    amplitude = max(yhat, na.rm = TRUE) - min(yhat, na.rm = TRUE),
    peak_time = grid[which.max(yhat)], trough_time = grid[which.min(yhat)],
    value_near0 = yhat[near0], nonlinear_p = nonlinear_p,
    mean_pre_slope = if (length(pre_idx)) mean(dy[pre_idx], na.rm = TRUE) else NA_real_,
    mean_post_slope = if (length(post_idx)) mean(dy[post_idx], na.rm = TRUE) else NA_real_,
    max_abs_slope = max(abs(dy), na.rm = TRUE)
  )
  list(curve = data.table(mid = grid, fitted_z = yhat, derivative = dy), features = feat)
}

make_aligned_trajectory_splines <- function(case_split_bins) {
  d <- as.data.table(case_split_bins)
  if (!nrow(d)) return(list(curves = data.table(), features = data.table()))
  keys <- unique(d[, .(alignment, protein)])
  curves <- data.table(); features <- data.table()
  for (i in seq_len(nrow(keys))) {
    al <- keys$alignment[i]; pr <- keys$protein[i]
    ff <- fit_aligned_trajectory_spline(d[alignment == al & protein == pr])
    if (nrow(ff$curve)) curves <- rbind(curves, ff$curve[, `:=`(alignment = al, protein = pr)], fill = TRUE)
    if (nrow(ff$features)) features <- rbind(features, ff$features[, `:=`(alignment = al, protein = pr)], fill = TRUE)
  }
  list(curves = curves, features = features)
}

make_temporal_decomp_figure <- function(datY, Y) {
  set_scope("trait", Y)
  proteins <- intersect(vip_proteins, names(datY))
  if (!length(proteins)) {
    warning("No selected proteins available for temporal_decomp in ", Y, call. = FALSE)
    return(list(case_split = data.table(), lag_cox = data.table(), occurrence_contrast = data.table(), evidence = data.table()))
  }
  base_covs <- default_demog_covars(datY, include_tdi = TRUE)
  td_cols <- unique(c("eid", "never_y", "b2e", "t2e_y", "Yt2e_y", "occ_last_pre_b2e", "occ_state", "drug.lipid", proteins, base_covs))
  dat_td <- datY[, td_cols[td_cols %in% names(datY)], drop = FALSE]
  datZ <- make_protein_z(dat_td, proteins, covars = character(), ref_idx = dat_td$never_y, suffix = "td_raw")
  datZ <- make_protein_z(datZ, proteins, covars = base_covs, ref_idx = datZ$never_y, suffix = "td_adj")
  z_raw <- paste0(proteins, "__td_raw")
  z_adj <- paste0(proteins, "__td_adj")

  make_alignment_bins <- function(label, xcol, keep_fun) {
    if (!xcol %in% names(datZ)) return(data.table())
    dd <- keep_fun(datZ[, unique(c(xcol, z_raw)), drop = FALSE])
    if (!nrow(dd)) return(data.table())
    dd$time_align <- safe_num(dd[[xcol]])
    dd <- dd %>% filter(is.finite(time_align), time_align >= -16, time_align <= 16)
    if (!nrow(dd)) return(data.table())
    as.data.table(bin_protein_means(dd, "time_align", z_raw, proteins, breaks = seq(-16, 16, 2), skip = c(1, 0), min_n = min_bin_n))[, alignment := label][]
  }
  case_split_bins <- rbindlist(list(
    make_alignment_bins("Incident: first post-baseline date", "b2e", function(d) d %>% filter(is.finite(b2e), b2e > 0)),
    make_alignment_bins("Prevalent: first/FOD before baseline", "b2e", function(d) d %>% filter(is.finite(b2e), b2e < 0)),
    make_alignment_bins("Prevalent: last pre-baseline AOD", "occ_last_pre_b2e", function(d) d %>% filter(is.finite(occ_last_pre_b2e), occ_last_pre_b2e < 0))
  ), fill = TRUE)
  case_split_trend <- if (nrow(case_split_bins)) {
    case_split_bins %>% group_by(alignment, protein) %>% group_modify(~ trend_bin_lm(.x)) %>% ungroup() %>% as.data.table()
  } else data.table()
  trajectory_splines <- make_aligned_trajectory_splines(case_split_bins)
  trajectory_curves <- trajectory_splines$curves
  trajectory_features <- trajectory_splines$features

  lags <- safe_num(parse_vec(Sys.getenv("YY_LAG_YEARS", unset = "0,2,5,10"), default = c("0", "2", "5", "10")))
  lags <- unique(lags[is.finite(lags) & lags >= 0])
  lag_res <- rbindlist(lapply(seq_along(proteins), function(i) {
    protein_name <- proteins[i]
    zc <- z_raw[i]
    rbindlist(lapply(lags, function(lag0) {
      cox_lag_one(datZ, "t2e_y", "Yt2e_y", zc, covars = base_covs, lag_year = lag0) %>% mutate(protein = .env$protein_name)
    }), fill = TRUE)
  }), fill = TRUE)
  if (nrow(lag_res)) lag_res[, `:=`(lo_beta = beta - 1.96 * se, hi_beta = beta + 1.96 * se)]

  occurrence_groups <- c("Prevalent single-date remote", "Prevalent multi-date remote", "Prevalent recent/active", "Prevalent + post-baseline recurrence", "Incident multi-date")
  occurrence_contrast <- rbindlist(lapply(seq_along(proteins), function(i) {
    protein_name <- proteins[i]
    zc <- z_adj[i]
    rbindlist(lapply(occurrence_groups, function(grp) {
      dd <- datZ[, c("occ_state", zc), drop = FALSE] %>%
        filter(occ_state %in% c("Never disease", grp), is.finite(.data[[zc]])) %>%
        mutate(x = as.integer(occ_state == grp))
      lm_contrast_one(dd, zc, "x", covars = character()) %>% mutate(protein = .env$protein_name, contrast = paste0(grp, " vs never"))
    }), fill = TRUE)
  }), fill = TRUE)
  if (nrow(occurrence_contrast)) occurrence_contrast[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]

  drug_delta <- data.table()
  if ("drug.lipid" %in% names(datZ)) {
    drug_delta <- rbindlist(lapply(seq_along(proteins), function(i) {
      protein_name <- proteins[i]
      zc <- z_raw[i]
      keep_cols <- unique(c("never_y", "drug.lipid", zc, base_covs))
      dd <- datZ[, keep_cols[keep_cols %in% names(datZ)], drop = FALSE] %>%
        filter(never_y, is.finite(.data[[zc]]), !is.na(drug.lipid)) %>%
        mutate(lipid_drug = as.integer(safe_num(drug.lipid) == 1))
      lm_contrast_one(dd, zc, "lipid_drug", covars = base_covs) %>% mutate(protein = .env$protein_name, contrast = "lipid drug vs no lipid drug among never disease")
    }), fill = TRUE)
    if (nrow(drug_delta)) drug_delta[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]
  }

  trend_all <- if (nrow(case_split_bins)) {
    trend_dat <- datZ[, unique(c("b2e", z_raw)), drop = FALSE]
    all_bins <- bin_protein_means(trend_dat %>% filter(is.finite(b2e), b2e >= -16, b2e <= 16), "b2e", z_raw, proteins, breaks = seq(-16, 16, 2), skip = c(1, 0), min_n = min_bin_n)
    all_bins %>% group_by(protein) %>% group_modify(~ trend_bin_lm(.x)) %>% ungroup() %>% mutate(axis = "YY all-case trend", beta = beta10, se = se10) %>% as.data.table()
  } else data.table()
  evidence <- rbindlist(list(
    if (nrow(trend_all)) trend_all[, .(protein, axis, beta, se, p)] else data.table(),
    if (nrow(lag_res)) lag_res[lag_year %in% c(0, 5), .(protein, axis = paste0("Incident Cox lag", lag_year, "y"), beta, se, p)] else data.table(),
    if (nrow(occurrence_contrast)) occurrence_contrast[contrast %in% c("Prevalent recent/active vs never", "Prevalent + post-baseline recurrence vs never"), .(protein, axis = contrast, beta, se, p)] else data.table(),
    if (nrow(drug_delta)) drug_delta[, .(protein, axis = contrast, beta, se, p)] else data.table()
  ), fill = TRUE)
  if (nrow(evidence)) {
    evidence <- merge(evidence, protein_prior, by = "protein", all.x = TRUE, sort = FALSE)
    evidence[is.na(prior), prior := "other"]
    evidence[, signed_logp := sign(beta) * pmin(-log10(pmax(p, 1e-300)), 20)]
  }

  keep_plot <- intersect(c("PCSK9", "GDF15", "IL6", "ANGPTL3", "APOE", "FGA"), proteins)
  if (length(keep_plot) < 2) keep_plot <- head(proteins, min(6L, length(proteins)))
  pA <- if (nrow(case_split_bins)) {
    ggplot(case_split_bins[protein %in% keep_plot], aes(mid, mean_z, color = protein, group = protein)) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = .3, color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dotted", linewidth = .3, color = "grey35") +
      geom_point(size = 1.35, alpha = .75, na.rm = TRUE) +
      geom_line(data = trajectory_curves[protein %in% keep_plot], aes(mid, fitted_z, color = protein, group = protein), linewidth = .78, na.rm = TRUE) +
      facet_grid(protein ~ alignment) + scale_color_manual(values = protein_colors(proteins), guide = "none") +
      labs(title = "a. Integrated Yin-Yang/AOD pseudo-time trajectories",
           subtitle = "Single baseline biomarker measurement aligned to incident, FOD, and last pre-baseline occurrence dates",
           x = "Time from baseline to occurrence (year)", y = "Biomarker z-score vs never disease") + theme_yy(8.5)
  } else ggplot() + labs(title = "a. Split timelines unavailable") + theme_yy(9)
  pB <- if (nrow(lag_res)) {
    ggplot(lag_res, aes(beta, protein, xmin = lo_beta, xmax = hi_beta, color = protein)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") + geom_errorbar(orientation = "y", width = .18, linewidth = .55, na.rm = TRUE) + geom_point(size = 1.8, na.rm = TRUE) +
      facet_wrap(~ paste0("Exclude first ", lag_year, " y"), nrow = 1) + scale_color_manual(values = protein_colors(proteins), guide = "none") +
      labs(title = "b. Prospective Cox after lag exclusion", subtitle = "Participants with incident events within the lag are excluded; model uses raw protein z-score + baseline covariates", x = "Cox beta", y = NULL) + theme_yy(9)
  } else ggplot() + labs(title = "b. Lagged Cox unavailable") + theme_yy(9)
  pC <- if (nrow(occurrence_contrast)) {
    ggplot(occurrence_contrast, aes(beta, protein, xmin = lo, xmax = hi, color = contrast)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") + geom_errorbar(orientation = "y", width = .18, linewidth = .55, na.rm = TRUE) + geom_point(size = 1.7, na.rm = TRUE) +
      labs(title = "c. Occurrence-state contrasts after covariate residualization", x = "Adjusted protein z-score difference", y = NULL, color = NULL) + theme_yy(9) + theme(legend.position = "bottom")
  } else ggplot() + labs(title = "c. Occurrence-state contrasts unavailable") + theme_yy(9)
  pD <- if (nrow(evidence)) {
    ggplot(evidence, aes(axis, protein, fill = signed_logp)) + geom_tile(color = "white") + geom_text(aes(label = ifelse(is.finite(beta), sprintf("%.2f", beta), "")), size = 2.5) +
      scale_fill_gradient2(low = BLUE, mid = "white", high = ORANGE, midpoint = 0, name = "signed\n-log10(P)") +
      labs(title = "d. Triage grid: timeline, incident risk, recurrence, and treatment", subtitle = "Numbers are beta estimates; color is signed evidence strength", x = NULL, y = NULL) + theme_minimal(base_size = 9) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1), strip.text = element_text(face = "bold"))
  } else ggplot() + labs(title = "d. Evidence grid unavailable") + theme_yy(9)
  Fig <- (pA / pB) / (pC | pD) + plot_layout(heights = c(1.25, .75, 1.05))
  save_plot(Fig, fig_file(5, "temporal_decomp", Y = Y), 17, 17)
  write_tab(case_split_bins, fig_data_file(5, "case_split_bins.tsv.gz", Y = Y))
  if (nrow(case_split_trend)) write_tab(case_split_trend, fig_data_file(5, "case_split_trend.tsv", Y = Y))
  if (nrow(trajectory_curves)) write_tab(trajectory_curves, fig_data_file(5, "trajectory_spline_curves.tsv.gz", Y = Y))
  if (nrow(trajectory_features)) write_tab(trajectory_features, fig_data_file(5, "trajectory_spline_features.tsv", Y = Y))
  if (nrow(lag_res)) write_tab(lag_res, fig_data_file(5, "lagged_cox.tsv", Y = Y))
  if (nrow(occurrence_contrast)) write_tab(occurrence_contrast, fig_data_file(5, "occurrence_contrast.tsv", Y = Y))
  if (nrow(drug_delta)) write_tab(drug_delta, fig_data_file(5, "drug_delta.tsv", Y = Y))
  if (nrow(evidence)) write_tab(evidence, fig_data_file(5, "triage_grid.tsv", Y = Y))
  write_book(list(case_split_bins = case_split_bins, case_split_trend = case_split_trend, trajectory_spline_features = trajectory_features, trajectory_spline_curves = trajectory_curves, lagged_cox = lag_res, occurrence_contrast = occurrence_contrast, drug_delta = drug_delta, triage_grid = evidence), fig_file(5, "temporal_decomp", ".out.xlsx", Y = Y))
  invisible(list(case_split = case_split_bins, trajectory_features = trajectory_features, lag_cox = lag_res, occurrence_contrast = occurrence_contrast, drug_delta = drug_delta, evidence = evidence))
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Integrated YY-informed prediction helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pred_make_dir <- function(...) { d <- file.path(...); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
pred_write_tab <- function(x, path) { data.table::fwrite(as.data.table(x), path); message("Wrote: ", path); invisible(path) }
pred_write_book <- function(x, path) {
  x <- x[vapply(x, function(z) is.data.frame(z) || data.table::is.data.table(z), logical(1))]
  if (!length(x)) x <- list(empty = data.frame(note = "No rows"))
  names(x) <- substr(make.unique(gsub("[^A-Za-z0-9_]+", "_", names(x))), 1, 31)
  writexl::write_xlsx(lapply(x, as.data.frame), path)
  message("Wrote: ", path)
  invisible(path)
}
pred_save_plot <- function(p, path, width = 9, height = 6, dpi = 300) {
  ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  message("Wrote: ", path)
  invisible(path)
}

pred_winsorize_train <- function(x, train_idx, p = c(0.001, 0.999)) {
  x <- safe_num(x)
  tr <- train_idx & is.finite(x)
  if (sum(tr) < 20) return(x)
  q <- stats::quantile(x[tr], p, na.rm = TRUE, names = FALSE, type = 8)
  pmin(pmax(x, q[1]), q[2])
}

pred_scale_by_train <- function(x, train_idx, winsor = TRUE) {
  x <- safe_num(x)
  if (winsor) x <- pred_winsorize_train(x, train_idx)
  tr <- train_idx & is.finite(x)
  mu <- mean(x[tr], na.rm = TRUE)
  s <- sd(x[tr], na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- sd(x[is.finite(x)], na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mu) / s
}

pred_impute_train_median <- function(x, train_idx) {
  x <- safe_num(x)
  m <- median(x[train_idx & is.finite(x)], na.rm = TRUE)
  if (!is.finite(m)) m <- median(x[is.finite(x)], na.rm = TRUE)
  if (!is.finite(m)) m <- 0
  x[!is.finite(x)] <- m
  x
}

# Model helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pred_make_folds <- function(event, k = 10, seed = 2026) {
  set.seed(seed)
  fold <- integer(length(event))
  for (e in sort(unique(event))) {
    ii <- which(event == e & !is.na(event))
    fold[ii] <- sample(rep(seq_len(k), length.out = length(ii)))
  }
  fold[fold == 0] <- sample(rep(seq_len(k), length.out = sum(fold == 0)))
  fold
}

pred_pick_base_covars <- function(dat) {
  if (length(base_covars_env)) return(safe_covars(dat, base_covars_env))
  demog <- default_demog_covars(dat, include_tdi = TRUE)
  le8 <- get_le8_covars(dat)
  genetic <- grep("cad.*prs|cad[._]6m|apoe|rare|lof|plof", names(dat), value = TRUE, ignore.case = TRUE)
  out <- switch(pred_base_mode,
                demog = demog,
                le8 = c(demog, le8),
                full = c(demog, le8, genetic),
                c(demog, le8))
  safe_covars(dat, out)
}

pred_select_top_balanced <- function(params, score_col = "cox_p", n_total = 700L, decreasing = FALSE) {
  x <- as.data.table(params)
  if (!nrow(x) || !score_col %in% names(x)) return(character())
  x <- x[is.finite(get(score_col))]
  if (!nrow(x)) return(character())
  ord_fun <- function(dd) {
    if (decreasing) dd[order(-get(score_col))] else dd[order(get(score_col))]
  }
  if (pred_layer_balance == "none" || uniqueN(x$layer[!is.na(x$layer)]) < 2) {
    return(ord_fun(x)[seq_len(min(as.integer(n_total), .N))]$protein)
  }
  cnt <- x[, .N, by = layer]
  if (pred_layer_balance == "equal") cnt[, weight := 1]
  if (pred_layer_balance == "proportional") cnt[, weight := N]
  if (pred_layer_balance == "sqrt") cnt[, weight := sqrt(N)]
  cnt[, quota := pmax(1L, floor(as.integer(n_total) * weight / sum(weight)))]
  while (sum(cnt$quota) > n_total) {
    j <- which.max(cnt$quota); cnt$quota[j] <- cnt$quota[j] - 1L
  }
  keep <- rbindlist(lapply(seq_len(nrow(cnt)), function(i) {
    dd <- ord_fun(x[layer == cnt$layer[i]])
    dd[seq_len(min(cnt$quota[i], .N))]
  }), fill = TRUE)
  if (nrow(keep) < n_total) {
    extra <- ord_fun(x[!protein %in% keep$protein])
    keep <- rbind(keep, extra[seq_len(min(n_total - nrow(keep), .N))], fill = TRUE)
  }
  unique(keep$protein)[seq_len(min(n_total, uniqueN(keep$protein)))]
}

pred_prep_model_frame <- function(dat, rows, covars) {
  dd <- dat[rows, covars, drop = FALSE]
  for (v in names(dd)) {
    if (is.numeric(dd[[v]]) || is.integer(dd[[v]])) {
      dd[[v]] <- pred_impute_train_median(dd[[v]], rep(TRUE, nrow(dd)))
    } else {
      dd[[v]] <- as.character(dd[[v]])
      mode_v <- names(sort(table(dd[[v]][!is.na(dd[[v]]) & nzchar(dd[[v]])]), decreasing = TRUE))[1]
      if (is.na(mode_v) || !length(mode_v)) mode_v <- "Missing"
      dd[[v]][is.na(dd[[v]]) | !nzchar(dd[[v]])] <- mode_v
      dd[[v]] <- factor(dd[[v]])
    }
  }
  dd
}

pred_make_cov_design <- function(dat, train_rows, test_rows, covars) {
  covars <- safe_covars(dat, covars)
  if (!length(covars)) {
    return(list(x_train = matrix(nrow = length(train_rows), ncol = 0), x_test = matrix(nrow = length(test_rows), ncol = 0), terms = character()))
  }
  train_df <- dat[train_rows, covars, drop = FALSE]
  all_df <- dat[c(train_rows, test_rows), covars, drop = FALSE]
  for (v in covars) {
    if (is.numeric(all_df[[v]]) || is.integer(all_df[[v]])) {
      x <- safe_num(all_df[[v]])
      train_idx <- seq_along(train_rows)
      x <- pred_impute_train_median(x, seq_along(x) %in% train_idx)
      mu <- mean(x[train_idx], na.rm = TRUE); s <- sd(x[train_idx], na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1
      all_df[[v]] <- (x - mu) / s
    } else {
      tr <- as.character(train_df[[v]])
      tr[is.na(tr) | !nzchar(tr)] <- "Missing"
      lev <- names(sort(table(tr), decreasing = TRUE))
      if (!length(lev)) lev <- "Missing"
      z <- as.character(all_df[[v]])
      z[is.na(z) | !nzchar(z)] <- "Missing"
      z[!z %in% lev] <- "Other"
      lev <- unique(c(lev, "Other"))
      all_df[[v]] <- factor(z, levels = lev)
    }
  }
  mm <- model.matrix(~ . , data = all_df)
  if (ncol(mm) > 1) mm <- mm[, -1, drop = FALSE] else mm <- matrix(nrow = nrow(all_df), ncol = 0)
  list(x_train = mm[seq_along(train_rows), , drop = FALSE],
       x_test = mm[length(train_rows) + seq_along(test_rows), , drop = FALSE],
       terms = colnames(mm))
}

pred_make_cov_design_all <- function(dat, train_rows, covars) {
  covars <- safe_covars(dat, covars)
  if (!length(covars)) return(matrix(1, nrow = nrow(dat), ncol = 1, dimnames = list(NULL, "(Intercept)")))
  all_df <- dat[, covars, drop = FALSE]
  for (v in covars) {
    if (is.numeric(all_df[[v]]) || is.integer(all_df[[v]])) {
      x <- safe_num(all_df[[v]])
      x <- pred_impute_train_median(x, seq_len(nrow(dat)) %in% train_rows)
      mu <- mean(x[train_rows], na.rm = TRUE); s <- sd(x[train_rows], na.rm = TRUE)
      if (!is.finite(s) || s == 0) s <- 1
      all_df[[v]] <- (x - mu) / s
    } else {
      tr <- as.character(all_df[[v]][train_rows])
      tr[is.na(tr) | !nzchar(tr)] <- "Missing"
      lev <- names(sort(table(tr), decreasing = TRUE)); if (!length(lev)) lev <- "Missing"
      z <- as.character(all_df[[v]])
      z[is.na(z) | !nzchar(z)] <- "Missing"
      z[!z %in% lev] <- "Other"
      all_df[[v]] <- factor(z, levels = unique(c(lev, "Other")))
    }
  }
  model.matrix(~ ., data = all_df)
}

pred_residualize_matrix <- function(z_all, x_all, train_rows, block_size = 250L) {
  if (!ncol(z_all) || !length(train_rows)) return(z_all)
  x_all <- as.matrix(x_all)
  x_tr <- x_all[train_rows, , drop = FALSE]
  qx <- qr(x_tr, LAPACK = FALSE)
  blocks <- split(seq_len(ncol(z_all)), ceiling(seq_len(ncol(z_all)) / max(1L, block_size)))
  for (idx in blocks) {
    yy <- z_all[train_rows, idx, drop = FALSE]
    bb <- try(qr.coef(qx, yy), silent = TRUE)
    if (inherits(bb, "try-error")) next
    bb[!is.finite(bb)] <- 0
    rr <- z_all[, idx, drop = FALSE] - x_all %*% bb
    mu <- colMeans(rr[train_rows, , drop = FALSE], na.rm = TRUE)
    ss <- apply(rr[train_rows, , drop = FALSE], 2, sd, na.rm = TRUE)
    ss[!is.finite(ss) | ss == 0] <- 1
    rr <- sweep(rr, 2, mu, "-")
    rr <- sweep(rr, 2, ss, "/")
    z_all[, idx] <- rr
  }
  z_all
}

pred_make_protein_matrix <- function(dat, proteins, train_rows, test_rows) {
  proteins <- intersect(proteins, names(dat))
  if (!length(proteins)) return(list(x_train = matrix(nrow = length(train_rows), ncol = 0), x_test = matrix(nrow = length(test_rows), ncol = 0), proteins = character(), z_all = matrix(nrow = nrow(dat), ncol = 0)))
  train_idx_all <- rep(FALSE, nrow(dat)); train_idx_all[train_rows] <- TRUE
  z_all <- matrix(NA_real_, nrow = nrow(dat), ncol = length(proteins), dimnames = list(NULL, proteins))
  for (j in seq_along(proteins)) {
    x <- pred_scale_by_train(dat[[proteins[j]]], train_idx_all, winsor = TRUE)
    m <- median(x[train_idx_all & is.finite(x)], na.rm = TRUE); if (!is.finite(m)) m <- 0
    x[!is.finite(x)] <- m
    z_all[, j] <- x
  }
  list(x_train = z_all[train_rows, , drop = FALSE], x_test = z_all[test_rows, , drop = FALSE], proteins = proteins, z_all = z_all)
}

pred_cox_screen_one <- function(dat, y_time, y_event, x, covars_df, min_event = 20) {
  dd <- covars_df
  dd$time <- y_time; dd$event <- y_event; dd$x <- x
  ok <- complete.cases(dd) & is.finite(dd$time) & dd$time > 0 & is.finite(dd$x)
  dd <- dd[ok, , drop = FALSE]
  if (nrow(dd) < 500 || sum(dd$event == 1, na.rm = TRUE) < min_event || safe_sd(dd$x) == 0) {
    return(c(beta = NA_real_, se = NA_real_, p = NA_real_, n = nrow(dd), event = sum(dd$event == 1, na.rm = TRUE)))
  }
  vars <- setdiff(names(dd), c("time", "event"))
  fit <- try(survival::coxph(as.formula(paste0("Surv(time, event) ~ ", paste(bt(vars), collapse = " + "))), data = dd), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(beta = NA_real_, se = NA_real_, p = NA_real_, n = nrow(dd), event = sum(dd$event == 1, na.rm = TRUE)))
  co <- summary(fit)$coef
  if (!"x" %in% rownames(co)) return(c(beta = NA_real_, se = NA_real_, p = NA_real_, n = nrow(dd), event = sum(dd$event == 1, na.rm = TRUE)))
  c(beta = co["x", "coef"], se = co["x", "se(coef)"], p = co["x", "Pr(>|z|)"], n = nrow(dd), event = sum(dd$event == 1, na.rm = TRUE))
}

pred_calc_auc <- function(y, score) {
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]; score <- score[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(score)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

pred_calc_cindex <- function(time, event, lp) {
  ok <- is.finite(time) & is.finite(event) & is.finite(lp) & time > 0
  if (sum(ok) < 10 || sum(event[ok] == 1) < 2) return(NA_real_)
  # Cox linear predictors increase with risk, whereas concordance() defaults to
  # larger predictor = longer survival. reverse=TRUE is therefore essential.
  fit <- try(survival::concordance(survival::Surv(time[ok], event[ok]) ~ lp[ok], reverse = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  as.numeric(fit$concordance)
}

pred_fit_coxph_predict <- function(time_tr, event_tr, x_tr, x_te) {
  x_tr <- as.data.frame(x_tr); x_te <- as.data.frame(x_te)
  names(x_tr) <- make.names(names(x_tr), unique = TRUE); names(x_te) <- names(x_tr)
  dd <- x_tr; dd$time <- time_tr; dd$event <- event_tr
  ok <- complete.cases(dd) & is.finite(dd$time) & dd$time > 0
  if (sum(ok) < 500 || sum(dd$event[ok] == 1, na.rm = TRUE) < 20 || ncol(x_tr) == 0) return(rep(0, nrow(x_te)))
  fit <- try(survival::coxph(as.formula(paste0("Surv(time, event) ~ ", paste(bt(names(x_tr)), collapse = " + "))), data = dd[ok, , drop = FALSE], x = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(0, nrow(x_te)))
  as.numeric(predict(fit, newdata = x_te, type = "lp"))
}

pred_fit_glmnet_cox_predict <- function(time_tr, event_tr, x_tr, x_te, lambda_choice = "1se", alpha = 0.5,
                                           nfolds = 10, penalty_factor = NULL, seed_inner = seed) {
  x_tr <- as.matrix(x_tr); x_te <- as.matrix(x_te)
  ok <- is.finite(time_tr) & time_tr > 0 & is.finite(event_tr) & complete.cases(x_tr)
  x_tr <- x_tr[ok, , drop = FALSE]
  time_tr <- time_tr[ok]; event_tr <- event_tr[ok]
  if (nrow(x_tr) < 500 || sum(event_tr == 1, na.rm = TRUE) < 20 || ncol(x_tr) < 1) {
    return(list(lp = rep(0, nrow(x_te)), selected = character(), fit = NULL))
  }
  if (is.null(penalty_factor)) penalty_factor <- rep(1, ncol(x_tr))
  penalty_factor <- as.numeric(penalty_factor)
  if (length(penalty_factor) != ncol(x_tr)) stop("penalty_factor length does not match x_tr columns", call. = FALSE)
  y <- survival::Surv(time_tr, event_tr)
  nfolds <- min(nfolds, max(3, sum(event_tr == 1, na.rm = TRUE)))
  foldid <- pred_make_folds(event_tr, nfolds, seed_inner)
  fit <- try(glmnet::cv.glmnet(x_tr, y, family = "cox", alpha = alpha, foldid = foldid,
                               standardize = FALSE, penalty.factor = penalty_factor), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(lp = rep(0, nrow(x_te)), selected = character(), fit = NULL))
  lam <- if (lambda_choice == "min") fit$lambda.min else fit$lambda.1se
  lp <- as.numeric(predict(fit, newx = x_te, s = lam, type = "link"))
  cf <- as.matrix(coef(fit, s = lam))
  selected <- rownames(cf)[as.numeric(cf[, 1]) != 0]
  list(lp = lp, selected = selected, fit = fit)
}


pred_fit_adaptive_glmnet_cox_predict <- function(time_tr, event_tr, x_tr, x_te, n_base,
                                                   prior, eta_grid = pred_adaptive_eta,
                                                   lambda_choice = "1se", alpha = 0.5,
                                                   nfolds = 10, seed_inner = seed) {
  x_tr <- as.matrix(x_tr); x_te <- as.matrix(x_te)
  ok <- is.finite(time_tr) & time_tr > 0 & is.finite(event_tr) & complete.cases(x_tr)
  x_fit <- x_tr[ok, , drop = FALSE]
  time_fit <- time_tr[ok]; event_fit <- event_tr[ok]
  if (nrow(x_fit) < 500 || sum(event_fit == 1, na.rm = TRUE) < 20 || ncol(x_fit) < 1) {
    return(list(lp = rep(0, nrow(x_te)), selected = character(), fit = NULL,
                eta = NA_real_, cv_error = NA_real_, tuning = data.table()))
  }
  n_biom <- ncol(x_fit) - n_base
  if (n_biom < 1) {
    return(list(lp = rep(0, nrow(x_te)), selected = character(), fit = NULL,
                eta = NA_real_, cv_error = NA_real_, tuning = data.table()))
  }
  prior <- safe_num(prior)
  if (length(prior) != n_biom) stop("Adaptive prior length does not match biomarker columns", call. = FALSE)
  prior[!is.finite(prior)] <- median(prior[is.finite(prior)], na.rm = TRUE)
  if (!any(is.finite(prior))) prior[] <- 0
  med <- median(prior, na.rm = TRUE)
  sc <- IQR(prior, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- sd(prior, na.rm = TRUE)
  if (!is.finite(sc) || sc == 0) sc <- 1
  zprior <- pmax(-3, pmin(3, (prior - med) / sc))
  eta_grid <- unique(c(0, eta_grid[is.finite(eta_grid) & eta_grid >= 0]))
  y <- survival::Surv(time_fit, event_fit)
  nfolds <- min(nfolds, max(3, sum(event_fit == 1, na.rm = TRUE)))
  foldid <- pred_make_folds(event_fit, nfolds, seed_inner)
  fits <- vector("list", length(eta_grid)); tune <- data.table()
  for (i in seq_along(eta_grid)) {
    eta <- eta_grid[i]
    pf_biom <- pmax(0.20, pmin(5, exp(-eta * zprior)))
    penalty <- c(rep(0, n_base), pf_biom)
    fit <- try(glmnet::cv.glmnet(
      x_fit, y, family = "cox", alpha = alpha, foldid = foldid,
      standardize = FALSE, penalty.factor = penalty
    ), silent = TRUE)
    if (inherits(fit, "try-error")) next
    lam <- if (lambda_choice == "min") fit$lambda.min else fit$lambda.1se
    j <- which.min(abs(fit$lambda - lam))
    err <- fit$cvm[j]
    fits[[i]] <- fit
    tune <- rbind(tune, data.table(eta = eta, lambda = lam, cv_error = err), fill = TRUE)
  }
  if (!nrow(tune)) {
    return(list(lp = rep(0, nrow(x_te)), selected = character(), fit = NULL,
                eta = NA_real_, cv_error = NA_real_, tuning = tune))
  }
  best_row <- which.min(tune$cv_error)
  best_eta <- tune$eta[best_row]
  fit_i <- which.min(abs(eta_grid - best_eta))
  fit <- fits[[fit_i]]
  lam <- tune$lambda[best_row]
  lp <- as.numeric(predict(fit, newx = x_te, s = lam, type = "link"))
  cf <- as.matrix(coef(fit, s = lam))
  selected <- rownames(cf)[as.numeric(cf[, 1]) != 0]
  list(lp = lp, selected = selected, fit = fit, eta = best_eta,
       cv_error = tune$cv_error[best_row], tuning = tune)
}

pred_calc_segment_trend_one <- function(x_time, z, lo, hi, min_n = 60L) {
  ok <- is.finite(x_time) & is.finite(z) & x_time >= lo & x_time < hi
  if (sum(ok) < min_n) return(c(beta10 = NA_real_, se10 = NA_real_, p = NA_real_, n = sum(ok)))
  tt <- x_time[ok] / 10
  fit <- try(lm(z[ok] ~ tt), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(beta10 = NA_real_, se10 = NA_real_, p = NA_real_, n = sum(ok)))
  co <- summary(fit)$coef
  c(beta10 = co[2, 1], se10 = co[2, 2], p = co[2, 4], n = sum(ok))
}

pred_calc_boundary_jump <- function(x_time, z, width = 2) {
  pre <- is.finite(x_time) & is.finite(z) & x_time >= -width & x_time < 0
  post <- is.finite(x_time) & is.finite(z) & x_time > 0 & x_time <= width
  pred_calc_mean_diff(z, post, pre)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# YY parameter learning inside training fold
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pred_trend_bin_lm <- function(mid, mean_z, n) {
  ok <- is.finite(mid) & is.finite(mean_z) & is.finite(n) & n > 0
  if (sum(ok) < 4) return(c(beta10 = NA_real_, se10 = NA_real_, p = NA_real_))
  fit <- try(lm(mean_z[ok] ~ I(mid[ok] / 10), weights = n[ok]), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(beta10 = NA_real_, se10 = NA_real_, p = NA_real_))
  co <- summary(fit)$coef
  c(beta10 = co[2, 1], se10 = co[2, 2], p = co[2, 4])
}

pred_calc_yy_trend_one <- function(x_time, z, breaks = seq(-16, 16, 2), min_n = 20) {
  ok <- is.finite(x_time) & is.finite(z) & x_time >= min(breaks) & x_time <= max(breaks)
  if (sum(ok) < min_n * 4) return(c(beta10 = NA_real_, se10 = NA_real_, p = NA_real_, bins = 0))
  bin <- cut(x_time[ok], breaks = breaks, include.lowest = TRUE, right = TRUE)
  mids <- head(breaks, -1) + diff(breaks) / 2
  lv <- levels(cut(mids, breaks = breaks, include.lowest = TRUE, right = TRUE))
  n <- as.numeric(table(factor(bin, levels = lv)))
  mean_z <- tapply(z[ok], factor(bin, levels = lv), mean, na.rm = TRUE)
  mean_z <- as.numeric(mean_z)
  mean_z[n < min_n] <- NA_real_
  rr <- pred_trend_bin_lm(mids, mean_z, n)
  c(rr, bins = sum(is.finite(mean_z)))
}

pred_calc_mean_diff <- function(z, idx1, idx0) {
  idx1 <- idx1 & is.finite(z); idx0 <- idx0 & is.finite(z)
  if (sum(idx1) < 20 || sum(idx0) < 50) return(c(beta = NA_real_, se = NA_real_, p = NA_real_, n1 = sum(idx1), n0 = sum(idx0)))
  tt <- try(t.test(z[idx1], z[idx0]), silent = TRUE)
  if (inherits(tt, "try-error")) return(c(beta = NA_real_, se = NA_real_, p = NA_real_, n1 = sum(idx1), n0 = sum(idx0)))
  c(beta = mean(z[idx1], na.rm = TRUE) - mean(z[idx0], na.rm = TRUE), se = unname(diff(tt$conf.int)) / (2 * 1.96), p = tt$p.value, n1 = sum(idx1), n0 = sum(idx0))
}

pred_learn_yy_parameters <- function(datY, proteins, z_all, train_rows, base_covars) {
  train_flag <- rep(FALSE, nrow(datY)); train_flag[train_rows] <- TRUE
  yin_train <- train_flag & datY$yin_y
  yang_train <- train_flag & datY$yang_y
  control_train <- yin_train & datY$Yt2e_y == 0
  incident_train <- yin_train & datY$Yt2e_y == 1
  incident_early <- incident_train & datY$t2e_y <= 5
  cases_train <- train_flag & is.finite(datY$b2e) & datY$b2e >= -16 & datY$b2e <= 16
  prevalent_recent <- yang_train & datY$occ_state %in% c("Prevalent recent/active", "Prevalent + post-baseline recurrence")
  prevalent_remote <- yang_train & datY$occ_state %in% c("Prevalent single-date remote", "Prevalent multi-date remote")
  recurrent <- yang_train & datY$occ_recurrent %in% TRUE
  nonrecurrent <- yang_train & !datY$occ_recurrent %in% TRUE

  covars <- safe_covars(datY, base_covars)
  cov_df <- pred_prep_model_frame(datY, which(yin_train), covars)
  time_y <- datY$t2e_y[yin_train]
  event_y <- datY$Yt2e_y[yin_train]

  out <- vector("list", length(proteins))
  for (j in seq_along(proteins)) {
    z <- z_all[, j]
    inc <- pred_cox_screen_one(datY[which(yin_train), , drop = FALSE], time_y, event_y, z[yin_train], cov_df)
    tr_all <- pred_calc_yy_trend_one(datY$b2e[cases_train], z[cases_train], min_n = pred_min_bin_n)
    tr_inc <- pred_calc_yy_trend_one(datY$b2e[incident_train], z[incident_train], breaks = seq(0, 16, 2), min_n = pred_min_bin_n)
    tr_prev_fod <- pred_calc_yy_trend_one(datY$b2e[yang_train], z[yang_train], breaks = seq(-16, 0, 2), min_n = pred_min_bin_n)
    tr_prev_last <- pred_calc_yy_trend_one(datY$occ_last_pre_b2e[yang_train], z[yang_train], breaks = seq(-16, 0, 2), min_n = pred_min_bin_n)
    prev_remote <- pred_calc_segment_trend_one(datY$b2e[yang_train], z[yang_train], -16, -8)
    prev_near <- pred_calc_segment_trend_one(datY$b2e[yang_train], z[yang_train], -8, 0)
    incident_early_slope <- pred_calc_segment_trend_one(datY$b2e[incident_train], z[incident_train], 0, 5)
    incident_late_slope <- pred_calc_segment_trend_one(datY$b2e[incident_train], z[incident_train], 5, 16)
    boundary_jump <- pred_calc_boundary_jump(datY$b2e[cases_train], z[cases_train], width = 2)
    prev_vs_never <- pred_calc_mean_diff(z, yang_train, control_train)
    early_vs_never <- pred_calc_mean_diff(z, incident_early, control_train)
    recent_vs_remote <- pred_calc_mean_diff(z, prevalent_recent, prevalent_remote)
    recurrent_vs_non <- pred_calc_mean_diff(z, recurrent, nonrecurrent)
    out[[j]] <- data.table(
      protein = proteins[j],
      cox_beta = inc["beta"], cox_se = inc["se"], cox_p = inc["p"], cox_n = inc["n"], cox_event = inc["event"],
      yy_beta10 = tr_all["beta10"], yy_se10 = tr_all["se10"], yy_p = tr_all["p"], yy_bins = tr_all["bins"],
      yy_incident_beta10 = tr_inc["beta10"], yy_incident_p = tr_inc["p"], yy_incident_bins = tr_inc["bins"],
      yy_prevalent_fod_beta10 = tr_prev_fod["beta10"], yy_prevalent_fod_p = tr_prev_fod["p"], yy_prevalent_fod_bins = tr_prev_fod["bins"],
      yy_prevalent_lastAOD_beta10 = tr_prev_last["beta10"], yy_prevalent_lastAOD_p = tr_prev_last["p"], yy_prevalent_lastAOD_bins = tr_prev_last["bins"],
      prev_remote_beta10 = prev_remote["beta10"], prev_remote_p = prev_remote["p"],
      prev_near_beta10 = prev_near["beta10"], prev_near_p = prev_near["p"],
      incident_early_beta10 = incident_early_slope["beta10"], incident_early_slope_p = incident_early_slope["p"],
      incident_late_beta10 = incident_late_slope["beta10"], incident_late_slope_p = incident_late_slope["p"],
      boundary_jump_beta = boundary_jump["beta"], boundary_jump_p = boundary_jump["p"],
      prev_vs_never_beta = prev_vs_never["beta"], prev_vs_never_p = prev_vs_never["p"],
      early_incident_vs_never_beta = early_vs_never["beta"], early_incident_vs_never_p = early_vs_never["p"],
      recent_vs_remote_beta = recent_vs_remote["beta"], recent_vs_remote_p = recent_vs_remote["p"],
      recurrent_vs_nonrecurrent_beta = recurrent_vs_non["beta"], recurrent_vs_nonrecurrent_p = recurrent_vs_non["p"]
    )
  }
  res <- rbindlist(out, fill = TRUE)
  res <- merge(res, annotate_proteins(res$protein), by = "protein", all.x = TRUE, sort = FALSE)
  ev20 <- function(p) pmin(-log10(pmax(p, 1e-300)), 20)
  res[, incident_evidence := pmin(-log10(pmax(cox_p, 1e-300)), 50)]
  res[, yy_evidence := pmin(-log10(pmax(yy_p, 1e-300)), 50)]
  res[, yy_strength := abs(yy_beta10) * pmin(yy_evidence, 20)]
  res[, early_strength := abs(early_incident_vs_never_beta) * ev20(early_incident_vs_never_p)]
  res[, prevalent_strength := abs(prev_vs_never_beta) * ev20(prev_vs_never_p)]
  res[, yang_state_strength :=
        abs(prev_vs_never_beta) * ev20(prev_vs_never_p) +
        abs(recent_vs_remote_beta) * ev20(recent_vs_remote_p) +
        abs(recurrent_vs_nonrecurrent_beta) * ev20(recurrent_vs_nonrecurrent_p)]
  res[, trajectory_strength :=
        abs(yy_prevalent_fod_beta10) * ev20(yy_prevalent_fod_p) +
        abs(yy_prevalent_lastAOD_beta10) * ev20(yy_prevalent_lastAOD_p) +
        abs(prev_remote_beta10) * ev20(prev_remote_p) +
        abs(prev_near_beta10) * ev20(prev_near_p) +
        abs(boundary_jump_beta) * ev20(boundary_jump_p)]
  res[, yang_prior_score := log1p(pmax(0, yang_state_strength) + pmax(0, trajectory_strength))]
  res[, trajectory_prior_score := log1p(pmax(0, trajectory_strength))]
  res[, reactive_flag := is.finite(cox_beta) & is.finite(yy_beta10) & cox_beta > 0 & yy_beta10 < 0]
  yy_cut <- if (any(is.finite(res$yy_beta10))) quantile(abs(res$yy_beta10), 0.35, na.rm = TRUE) else NA_real_
  res[, causal_like_flag := is.finite(cox_beta) & is.finite(yy_beta10) & is.finite(yy_cut) & abs(yy_beta10) <= yy_cut & cox_p < 0.05]
  res[, yy_rank_score := incident_evidence + 0.8 * yy_evidence + 0.6 * pmin(prevalent_strength, 20) + 0.6 * pmin(early_strength, 20)]
  res[, yy_rank_score := ifelse(is.finite(yy_rank_score), yy_rank_score, incident_evidence)]
  res[order(-yy_rank_score)]
}

pred_make_weighted_scores <- function(z_mat, params, proteins) {
  params <- as.data.table(params)
  params <- params[protein %in% proteins]
  if (!nrow(params)) return(matrix(0, nrow = nrow(z_mat), ncol = 0))
  ev <- function(p) pmin(-log10(pmax(p, 1e-300)), 20) / 20
  params[, e_prev := fifelse(is.finite(prev_vs_never_p), ev(prev_vs_never_p), 0)]
  params[, e_recent := fifelse(is.finite(recent_vs_remote_p), ev(recent_vs_remote_p), 0)]
  params[, e_recur := fifelse(is.finite(recurrent_vs_nonrecurrent_p), ev(recurrent_vs_nonrecurrent_p), 0)]
  params[, e_prev_fod := fifelse(is.finite(yy_prevalent_fod_p), ev(yy_prevalent_fod_p), 0)]
  params[, e_prev_last := fifelse(is.finite(yy_prevalent_lastAOD_p), ev(yy_prevalent_lastAOD_p), 0)]
  params[, e_boundary := fifelse(is.finite(boundary_jump_p), ev(boundary_jump_p), 0)]
  params[, w_yin_risk := fifelse(is.finite(cox_beta), cox_beta, 0)]
  params[, w_yang_prevalence := fifelse(is.finite(prev_vs_never_beta), prev_vs_never_beta * e_prev, 0)]
  params[, w_yang_activity := fifelse(is.finite(recent_vs_remote_beta), recent_vs_remote_beta * e_recent, 0)]
  params[, w_yang_recurrence := fifelse(is.finite(recurrent_vs_nonrecurrent_beta), recurrent_vs_nonrecurrent_beta * e_recur, 0)]
  params[, w_yang_proximity :=
           fifelse(is.finite(yy_prevalent_fod_beta10), -yy_prevalent_fod_beta10 * e_prev_fod, 0) +
           fifelse(is.finite(yy_prevalent_lastAOD_beta10), -yy_prevalent_lastAOD_beta10 * e_prev_last, 0)]
  params[, w_yy_boundary := fifelse(is.finite(boundary_jump_beta), boundary_jump_beta * e_boundary, 0)]
  prior_den <- median(params$yang_prior_score[is.finite(params$yang_prior_score)], na.rm = TRUE)
  if (!is.finite(prior_den) || prior_den <= 0) prior_den <- 1
  params[, w_yy_joint := w_yin_risk * (1 + pmin(2, yang_prior_score / prior_den))]
  wcols <- c("w_yin_risk", "w_yang_prevalence", "w_yang_activity", "w_yang_recurrence", "w_yang_proximity", "w_yy_boundary", "w_yy_joint")
  if (length(unique(params$layer[!is.na(params$layer)])) > 1) {
    params[, w_yang_prot := fifelse(layer == "prot", w_yang_prevalence + w_yang_activity + w_yang_recurrence + w_yang_proximity, 0)]
    params[, w_yang_met := fifelse(layer == "met", w_yang_prevalence + w_yang_activity + w_yang_recurrence + w_yang_proximity, 0)]
    wcols <- c(wcols, "w_yang_prot", "w_yang_met")
  }
  scores <- matrix(0, nrow = nrow(z_mat), ncol = length(wcols),
                   dimnames = list(NULL, sub("^w_", "score_", wcols)))
  for (k in seq_along(wcols)) {
    w <- params[[wcols[k]]]; names(w) <- params$protein
    w <- w[colnames(z_mat)]; w[!is.finite(w)] <- 0
    denom <- sqrt(sum(w^2, na.rm = TRUE)); if (!is.finite(denom) || denom == 0) denom <- 1
    scores[, k] <- as.numeric(z_mat %*% (w / denom))
  }
  scores
}

pred_make_yy_guided_weights <- function(params, n_total = top_n) {
  p <- copy(as.data.table(params))
  p <- p[is.finite(cox_beta) & is.finite(cox_p)]
  if (!nrow(p)) return(data.table())
  ev <- function(x, cap = 20) pmin(-log10(pmax(x, 1e-300)), cap)
  sgn <- function(x) fifelse(is.finite(x) & x > 0, 1, fifelse(is.finite(x) & x < 0, -1, 0))
  p[, `:=`(
    cox_ev = ev(cox_p, 50),
    yy_ev = ev(yy_p, 20),
    prev_ev = ev(prev_vs_never_p, 20),
    early_ev = ev(early_incident_vs_never_p, 20),
    recur_ev = ev(recurrent_vs_nonrecurrent_p, 20),
    cox_sign = sgn(cox_beta),
    yy_sign = sgn(yy_beta10),
    prev_sign = sgn(prev_vs_never_beta),
    early_sign = sgn(early_incident_vs_never_beta)
  )]
  p[, `:=`(
    yy_trajectory_support = is.finite(yy_beta10) & yy_sign == -cox_sign & yy_p < 0.05,
    yang_prevalence_support = is.finite(prev_vs_never_beta) & prev_sign == cox_sign & prev_vs_never_p < 0.05,
    early_incident_support = is.finite(early_incident_vs_never_beta) & early_sign == cox_sign & early_incident_vs_never_p < 0.05
  )]
  p[, mechanism_modifier := fcase(
    prior == "lipid/CAD-related prior", 1.25,
    prior == "stress/injury marker prior", 0.70,
    prior == "metabolic prior", 1.05,
    grepl("lipid|lipoprotein|cardiometabolic", category, ignore.case = TRUE), 1.10,
    default = 1.00
  )]
  p[, yy_support_score :=
      0.35 * fifelse(yy_trajectory_support, pmin(yy_ev, 10) / 10, 0) +
      0.35 * fifelse(yang_prevalence_support, pmin(prev_ev, 10) / 10, 0) +
      0.20 * fifelse(early_incident_support, pmin(early_ev, 10) / 10, 0) +
      0.10 * fifelse(is.finite(recurrent_vs_nonrecurrent_beta) & sgn(recurrent_vs_nonrecurrent_beta) == cox_sign & recurrent_vs_nonrecurrent_p < 0.05, pmin(recur_ev, 10) / 10, 0)]
  p[, prior_bonus := fcase(
    prior == "lipid/CAD-related prior", 1.5,
    prior == "stress/injury marker prior", -0.75,
    default = 0
  )]
  p[, yy_filter_score := cox_ev + 8 * yy_support_score + prior_bonus]
  p[, keep_by_yy := cox_p < 0.05 &
      (yy_trajectory_support | yang_prevalence_support | early_incident_support |
         prior == "lipid/CAD-related prior" | cox_p < 0.05 / max(1, nrow(p)))]
  keep <- p[keep_by_yy %in% TRUE]
  if (nrow(keep) < max(5L, min(20L, n_total))) {
    keep <- unique(rbindlist(list(keep, p[order(cox_p)][seq_len(min(max(20L, n_total), .N))]), fill = TRUE), by = "protein")
  }
  keep <- keep[order(-yy_filter_score, cox_p)]
  keep <- keep[seq_len(min(n_total, .N))]
  if (!nrow(keep)) return(data.table())
  keep[, weight_multiplier := pmax(0.25, pmin(2.00, mechanism_modifier * (0.75 + yy_support_score)))]
  keep[, yy_guided_weight := cox_beta * weight_multiplier]
  keep[!is.finite(yy_guided_weight), yy_guided_weight := 0]
  keep[, yy_keep_reason := fcase(
    yy_trajectory_support & yang_prevalence_support, "Yin Cox + Yang trajectory/prevalence",
    yy_trajectory_support, "Yin Cox + Yang trajectory",
    yang_prevalence_support, "Yin Cox + Yang prevalence",
    early_incident_support, "Yin Cox + early incident support",
    prior == "lipid/CAD-related prior", "Yin Cox + lipid/CAD prior",
    default = "Yin Cox fallback"
  )]
  keep[, yy_pattern := fcase(
    cox_beta > 0 & yy_beta10 < 0 & cox_p < 0.05 & yy_p < 0.05, "reactive/risk",
    cox_beta < 0 & yy_beta10 > 0 & cox_p < 0.05 & yy_p < 0.05, "mirror/protective",
    cox_p < 0.05 & (!is.finite(yy_p) | yy_p >= 0.05), "Yin-only",
    default = "mixed"
  )]
  keep[, .(protein, label, layer, category, prior,
           cox_beta, cox_p, yy_beta10, yy_p, prev_vs_never_beta, prev_vs_never_p,
           early_incident_vs_never_beta, early_incident_vs_never_p,
           yy_filter_score, yy_support_score, mechanism_modifier, weight_multiplier,
           yy_guided_weight, yy_keep_reason, yy_pattern)]
}

pred_linear_score_from_weights <- function(z_mat, weights, weight_col = "yy_guided_weight") {
  weights <- as.data.table(weights)
  if (!nrow(weights) || !weight_col %in% names(weights)) return(rep(0, nrow(z_mat)))
  w <- weights[[weight_col]]
  names(w) <- weights$protein
  w <- w[colnames(z_mat)]
  w[!is.finite(w)] <- 0
  denom <- sqrt(sum(w^2, na.rm = TRUE))
  if (!is.finite(denom) || denom <= 0) return(rep(0, nrow(z_mat)))
  as.numeric(z_mat %*% (w / denom))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# One-trait prediction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_yy_prediction <- function(datY, vars_biom, Y) {
  set_scope("trait", Y)
  message("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  message("# YY prediction for outcome: ", Y)
  message("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

  features <- intersect(vars_biom, names(datY))
  all_N0 <- nrow(datY)
  has_any <- row_any_nonmissing(datY, features)
  cohort_keep <- has_any
  cohort_N0 <- sum(cohort_keep, na.rm = TRUE)
  if (cohort_N0 < 1) stop("No participants have biomarker measurements for prediction.", call. = FALSE)
  features <- features[colMeans(!is.na(datY[cohort_keep, features, drop = FALSE])) >= 0.80]
  datY <- datY[cohort_keep, , drop = FALSE]
  if (pred_max_proteins > 0 && length(features) > pred_max_proteins) features <- features[seq_len(pred_max_proteins)]
  if (length(features) < 10) stop("Too few biomarkers available for prediction: ", length(features), call. = FALSE)

  base_covars <- pred_pick_base_covars(datY)
  yin <- which(datY$yin_y %in% TRUE & is.finite(datY$t2e_y) & datY$t2e_y > 0 & !is.na(datY$Yt2e_y))
  yang <- which(datY$yang_y %in% TRUE)
  yin_events <- sum(datY$Yt2e_y[yin] == 1, na.rm = TRUE)
  meta_pred <- annotate_proteins(features)
  layer_feature_counts <- meta_pred[, .N, by = layer]

  cohort <- data.table(
    outcome = Y, biom = biom, all_N = all_N0, biomarker_cohort_N = cohort_N0,
    features = length(features), yin_N = length(yin), yin_events = yin_events,
    yang_N = length(yang), residualized = pred_residualize,
    base_mode = pred_base_mode,
    le8_covars = paste(get_le8_covars(datY), collapse = ";"),
    base_covars = paste(base_covars, collapse = ";"),
    features_by_layer = paste(layer_feature_counts$layer, layer_feature_counts$N, sep = ":", collapse = ";")
  )
  pred_write_tab(cohort, file.path(rawdir, fig_data_file(6, "prediction.cohort.tsv", Y = Y)))
  message("Prediction cohort for ", Y, ": N=", cohort_N0,
          "; Yin_N=", length(yin), "; Yin_events=", yin_events,
          "; Yang_N=", length(yang), "; biomarkers=", length(features),
          "; residualized=", pred_residualize)

  if (length(yin) < pred_min_yin_n || yin_events < pred_min_events) {
    msg <- paste0("Prediction skipped for ", Y, ": too few Yin samples/events after biomarker filtering. ",
                  "Yin_N=", length(yin), ", Yin_events=", yin_events,
                  ", thresholds=", pred_min_yin_n, "/", pred_min_events, ".")
    warning(msg, call. = FALSE)
    writeLines(msg, file.path(rawdir, fig_data_file(6, "prediction.SKIPPED.txt", Y = Y)))
    pred_write_book(list(cohort = cohort, note = data.table(note = msg)), file.path(outdir, fig_file(6, "prediction", ".out.xlsx", Y = Y)))
    return(list(skipped = TRUE, cohort = cohort, message = msg))
  }

  fold_yin <- pred_make_folds(datY$Yt2e_y[yin], K_outer, seed)
  pred_rows <- list(); params_rows <- list(); selected_rows <- list(); weights_rows <- list(); adaptive_rows <- list(); yy_guided_rows <- list()

  for (fold in seq_len(K_outer)) {
    message("Fold ", fold, "/", K_outer)
    test_rows <- yin[fold_yin == fold]
    train_yin_rows <- yin[fold_yin != fold]
    # Yang is an auxiliary disease-state training set and is never part of the
    # held-out Yin evaluation set; all available Yang samples can therefore be
    # reused in each outer fold without leaking held-out Yin outcomes.
    train_learn_rows <- unique(c(train_yin_rows, yang))

    pm_all <- pred_make_protein_matrix(datY, features, train_yin_rows, integer())
    z_all <- pm_all$z_all
    if (pred_residualize && length(base_covars)) {
      x_cov_all <- pred_make_cov_design_all(datY, train_yin_rows, base_covars)
      z_all <- pred_residualize_matrix(z_all, x_cov_all, train_yin_rows, pred_block_size)
      rm(x_cov_all); invisible(gc())
    }

    yy_params <- pred_learn_yy_parameters(datY, features, z_all, train_learn_rows, base_covars)
    yy_params[, fold := fold]
    params_rows[[fold]] <- yy_params

    inc_set <- pred_select_top_balanced(yy_params, "cox_p", screen_n, decreasing = FALSE)
    yy_guided_weights <- pred_make_yy_guided_weights(yy_params, n_total = top_n)
    yy_guided_weights[, fold := fold]
    yy_guided_rows[[fold]] <- yy_guided_weights

    model_sets <- list(M1_Yin_glmnet = unique(inc_set))
    if (isTRUE(pred_legacy_complex)) {
      reactive_pool <- yy_params[reactive_flag %in% TRUE & is.finite(yy_rank_score)]
      reactive_set <- pred_select_top_balanced(reactive_pool, "yy_rank_score", top_n, decreasing = TRUE)
      if (length(reactive_set) < 10) reactive_set <- pred_select_top_balanced(yy_params, "yy_rank_score", top_n, decreasing = TRUE)
      model_sets <- c(model_sets, list(
        M2_legacy_Yin_plus_Yang_scores = unique(inc_set),
        M3_legacy_YY_adaptive_glmnet = unique(inc_set),
        M4_legacy_YY_trajectory_glmnet = unique(inc_set),
        M5_legacy_YY_reactive_glmnet = unique(reactive_set)
      ))
    }
    if (biom == "prot+met" && isTRUE(pred_layer_ablation)) {
      prot_set <- yy_params[layer == "prot" & is.finite(cox_p)][order(cox_p)][seq_len(min(screen_n, .N))]$protein
      met_set <- yy_params[layer == "met" & is.finite(cox_p)][order(cox_p)][seq_len(min(screen_n, .N))]$protein
      model_sets$M1a_Yin_prot <- unique(prot_set)
      model_sets$M1b_Yin_met <- unique(met_set)
    }
    selected_rows[[paste0("candidate_", fold)]] <- rbindlist(lapply(names(model_sets), function(nm) {
      data.table(fold = fold, model = nm, selection_type = "candidate", protein = model_sets[[nm]])
    }), fill = TRUE)
    if (nrow(yy_guided_weights)) {
      selected_rows[[paste0("yy_guided_", fold)]] <- yy_guided_weights[, .(
        fold, model = "M2_YY_guided_score", selection_type = "YY_guided_weighted_score",
        protein, yy_guided_weight, weight_multiplier, yy_keep_reason, yy_pattern
      )]
    }

    cov_des <- pred_make_cov_design(datY, train_yin_rows, test_rows, base_covars)
    y_time_tr <- datY$t2e_y[train_yin_rows]
    y_event_tr <- datY$Yt2e_y[train_yin_rows]
    y_time_te <- datY$t2e_y[test_rows]
    y_event_te <- datY$Yt2e_y[test_rows]
    n_base <- ncol(cov_des$x_train)

    get_z_subset <- function(rows, feature_set) {
      feature_set <- intersect(feature_set, colnames(z_all))
      if (!length(feature_set)) return(matrix(nrow = length(rows), ncol = 0))
      z_all[rows, feature_set, drop = FALSE]
    }
    lp0 <- pred_fit_coxph_predict(y_time_tr, y_event_tr, cov_des$x_train, cov_des$x_test)
    pred_rows[[paste(fold, "M0", sep = "_")]] <- data.table(
      row_id = test_rows, eid = datY$eid[test_rows], fold = fold, model = "M0_base",
      time = y_time_te, event = y_event_te, lp = lp0)

    yy_score_train <- pred_linear_score_from_weights(z_all[train_yin_rows, , drop = FALSE], yy_guided_weights)
    yy_score_test <- pred_linear_score_from_weights(z_all[test_rows, , drop = FALSE], yy_guided_weights)
    xtr_yy <- cbind(cov_des$x_train, YY_guided_protein_score = yy_score_train)
    xte_yy <- cbind(cov_des$x_test, YY_guided_protein_score = yy_score_test)
    lp2 <- pred_fit_coxph_predict(y_time_tr, y_event_tr, xtr_yy, xte_yy)
    pred_rows[[paste(fold, "M2_YY_guided", sep = "_")]] <- data.table(
      row_id = test_rows, eid = datY$eid[test_rows], fold = fold, model = "M2_YY_guided_score",
      time = y_time_te, event = y_event_te, lp = lp2)

    scores_train_all <- scores_test_all <- scores_train_m2 <- scores_test_m2 <- matrix(nrow = 0, ncol = 0)
    yang_score_cols <- yy_aug_score_cols <- character()
    if (isTRUE(pred_legacy_complex)) {
      scores_train_all <- pred_make_weighted_scores(z_all[train_yin_rows, , drop = FALSE], yy_params, features)
      scores_test_all <- pred_make_weighted_scores(z_all[test_rows, , drop = FALSE], yy_params, features)
      scores_train_m2 <- pred_make_weighted_scores(z_all[train_yin_rows, , drop = FALSE], yy_params, inc_set)
      scores_test_m2 <- pred_make_weighted_scores(z_all[test_rows, , drop = FALSE], yy_params, inc_set)
      yang_score_cols <- grep("^score_yang_", colnames(scores_train_all), value = TRUE)
      yy_aug_score_cols <- unique(c(
        grep("^score_yang_", colnames(scores_train_m2), value = TRUE),
        intersect(c("score_yy_boundary", "score_yy_joint"), colnames(scores_train_m2))
      ))
      keep_variable_scores <- function(score_mat, score_cols) {
        score_cols[vapply(score_cols, function(v) {
          z <- safe_num(score_mat[, v]); is.finite(sd(z, na.rm = TRUE)) && sd(z, na.rm = TRUE) > 1e-8
        }, logical(1))]
      }
      yang_score_cols <- keep_variable_scores(scores_train_all, yang_score_cols)
      yy_aug_score_cols <- keep_variable_scores(scores_train_m2, yy_aug_score_cols)
      xtr_scores <- cbind(cov_des$x_train, scores_train_all[, yang_score_cols, drop = FALSE])
      xte_scores <- cbind(cov_des$x_test, scores_test_all[, yang_score_cols, drop = FALSE])
      lp6 <- pred_fit_coxph_predict(y_time_tr, y_event_tr, xtr_scores, xte_scores)
      pred_rows[[paste(fold, "M6", sep = "_")]] <- data.table(
        row_id = test_rows, eid = datY$eid[test_rows], fold = fold, model = "M6_legacy_Yang_scores_only",
        time = y_time_te, event = y_event_te, lp = lp6)
    }

    for (m in names(model_sets)) {
      feature_set <- model_sets[[m]]
      xtr_p <- get_z_subset(train_yin_rows, feature_set)
      xte_p <- get_z_subset(test_rows, feature_set)
      xtr <- cbind(cov_des$x_train, xtr_p)
      xte <- cbind(cov_des$x_test, xte_p)
      if (m == "M2_legacy_Yin_plus_Yang_scores") {
        xtr <- cbind(xtr, scores_train_m2[, yy_aug_score_cols, drop = FALSE])
        xte <- cbind(xte, scores_test_m2[, yy_aug_score_cols, drop = FALSE])
      }
      if (m %in% c("M3_legacy_YY_adaptive_glmnet", "M4_legacy_YY_trajectory_glmnet")) {
        prior_col <- ifelse(m == "M3_legacy_YY_adaptive_glmnet", "yang_prior_score", "trajectory_prior_score")
        prior_vec <- yy_params[[prior_col]][match(feature_set, yy_params$protein)]
        fitp <- pred_fit_adaptive_glmnet_cox_predict(
          y_time_tr, y_event_tr, xtr, xte, n_base = n_base, prior = prior_vec,
          eta_grid = pred_adaptive_eta, lambda_choice = lambda_choice,
          alpha = alpha_glmnet, seed_inner = seed + fold)
        if (nrow(fitp$tuning)) {
          adaptive_rows[[paste(fold, m, sep = "_")]] <- fitp$tuning[, `:=`(fold = fold, model = m, selected_eta = eta == fitp$eta)]
        }
      } else {
        # Base covariates are always unpenalized. In M2, the low-dimensional YY
        # score block is also unpenalized so that the test is genuinely
        # "M1 plus Yang information" rather than allowing glmnet to discard the
        # Yang block at the 1-SE lambda. Raw biomarkers remain penalized.
        if (m == "M2_legacy_Yin_plus_Yang_scores") {
          penalty <- c(rep(0, n_base), rep(1, ncol(xtr_p)), rep(0, length(yy_aug_score_cols)))
        } else {
          penalty <- c(rep(0, n_base), rep(1, ncol(xtr) - n_base))
        }
        fitp <- pred_fit_glmnet_cox_predict(
          y_time_tr, y_event_tr, xtr, xte, lambda_choice = lambda_choice,
          alpha = alpha_glmnet, penalty_factor = penalty, seed_inner = seed + fold)
      }
      pred_rows[[paste(fold, m, sep = "_")]] <- data.table(
        row_id = test_rows, eid = datY$eid[test_rows], fold = fold, model = m,
        time = y_time_te, event = y_event_te, lp = fitp$lp)
      if (length(fitp$selected)) {
        selected_rows[[paste(fold, m, "glmnet", sep = "_")]] <- data.table(
          fold = fold, model = m, selection_type = "glmnet_nonzero", protein = fitp$selected)
      }
    }

    weights_rows[[fold]] <- yy_params[, .(fold, protein, label, layer, category, prior,
                                           cox_beta, cox_p, yy_beta10, yy_p, yy_rank_score,
                                           yang_prior_score, trajectory_prior_score,
                                           prev_remote_beta10, prev_near_beta10,
                                           incident_early_beta10, incident_late_beta10,
                                           boundary_jump_beta, reactive_flag, causal_like_flag)]
    rm(z_all, pm_all, scores_train_all, scores_test_all, scores_train_m2, scores_test_m2); invisible(gc())
  }

  pred <- rbindlist(pred_rows, fill = TRUE)
  yy_params_all <- rbindlist(params_rows, fill = TRUE)
  selected <- rbindlist(selected_rows, fill = TRUE)
  if (nrow(selected)) {
    selected[, feature_type := fcase(
      protein %in% features, "biomarker",
      grepl("^score_", protein), "YY_score",
      default = "base_covariate"
    )]
  }
  weights_all <- rbindlist(weights_rows, fill = TRUE)
  adaptive_tuning <- rbindlist(adaptive_rows, fill = TRUE)
  yy_guided_all <- rbindlist(yy_guided_rows, fill = TRUE)

  perf <- pred[, .(
    n = .N, events = sum(event == 1, na.rm = TRUE),
    cindex = pred_calc_cindex(time, event, lp)
  ), by = .(fold, model)]
  auc_rows <- pred[, {
    elig <- time >= 10 | event == 1
    .(AUC10 = pred_calc_auc(as.integer(event[elig] == 1 & time[elig] <= 10), lp[elig]),
      n10 = sum(elig), events10 = sum(event[elig] == 1 & time[elig] <= 10, na.rm = TRUE))
  }, by = .(fold, model)]
  perf <- merge(perf, auc_rows, by = c("fold", "model"), all.x = TRUE)
  base <- perf[model == "M0_base", .(fold, base_cindex = cindex, base_AUC10 = AUC10)]
  perf <- merge(perf, base, by = "fold", all.x = TRUE)
  perf[, `:=`(delta_cindex = cindex - base_cindex, delta_AUC10 = AUC10 - base_AUC10)]
  perf_sum <- perf[, .(
    mean_cindex = mean(cindex, na.rm = TRUE), sd_cindex = sd(cindex, na.rm = TRUE),
    mean_delta_cindex = mean(delta_cindex, na.rm = TRUE), sd_delta_cindex = sd(delta_cindex, na.rm = TRUE),
    mean_AUC10 = mean(AUC10, na.rm = TRUE), sd_AUC10 = sd(AUC10, na.rm = TRUE),
    mean_delta_AUC10 = mean(delta_AUC10, na.rm = TRUE), sd_delta_AUC10 = sd(delta_AUC10, na.rm = TRUE)
  ), by = model][order(-mean_cindex)]

  contrast_pairs <- data.table(
    contrast = c("YY-guided explainable score vs traditional Yin glmnet",
                 "YY-guided explainable score vs base"),
    model_a = c("M2_YY_guided_score", "M2_YY_guided_score"),
    model_b = c("M1_Yin_glmnet", "M0_base")
  )
  if (isTRUE(pred_legacy_complex)) {
    contrast_pairs <- rbind(contrast_pairs, data.table(
      contrast = c("Legacy Yang score block vs traditional Yin",
                   "Legacy Yang-informed adaptive penalty vs traditional Yin",
                   "Legacy trajectory-informed adaptive penalty vs traditional Yin",
                   "Legacy reactive-feature glmnet vs traditional Yin",
                   "Legacy Yang-only scores vs base"),
      model_a = c("M2_legacy_Yin_plus_Yang_scores", "M3_legacy_YY_adaptive_glmnet", "M4_legacy_YY_trajectory_glmnet",
                  "M5_legacy_YY_reactive_glmnet", "M6_legacy_Yang_scores_only"),
      model_b = c("M1_Yin_glmnet", "M1_Yin_glmnet", "M1_Yin_glmnet", "M1_Yin_glmnet", "M0_base")
    ), fill = TRUE)
  }
  if (all(c("M1a_Yin_prot", "M1_Yin_glmnet") %in% perf$model)) {
    contrast_pairs <- rbind(contrast_pairs, data.table(contrast = "Integrated Yin vs protein-only Yin on the same cohort", model_a = "M1_Yin_glmnet", model_b = "M1a_Yin_prot"))
  }
  if (all(c("M1b_Yin_met", "M1_Yin_glmnet") %in% perf$model)) {
    contrast_pairs <- rbind(contrast_pairs, data.table(contrast = "Integrated Yin vs metabolite-only Yin on the same cohort", model_a = "M1_Yin_glmnet", model_b = "M1b_Yin_met"))
  }
  perf_contrast_by_fold <- rbindlist(lapply(seq_len(nrow(contrast_pairs)), function(i) {
    a <- perf[model == contrast_pairs$model_a[i], .(fold, c_a = cindex, a_a = AUC10)]
    b <- perf[model == contrast_pairs$model_b[i], .(fold, c_b = cindex, a_b = AUC10)]
    d <- merge(a, b, by = "fold")
    if (!nrow(d)) return(data.table())
    d[, `:=`(contrast = contrast_pairs$contrast[i], model_a = contrast_pairs$model_a[i], model_b = contrast_pairs$model_b[i],
             delta_cindex = c_a - c_b, delta_AUC10 = a_a - a_b)]
    d[]
  }), fill = TRUE)
  perf_contrasts <- perf_contrast_by_fold[, {
    dc <- delta_cindex[is.finite(delta_cindex)]; da <- delta_AUC10[is.finite(delta_AUC10)]
    tc <- if (length(dc) >= 2) try(t.test(dc), silent = TRUE) else NULL
    ta <- if (length(da) >= 2) try(t.test(da), silent = TRUE) else NULL
    .(n_folds = .N,
      mean_delta_cindex = mean(dc, na.rm = TRUE), sd_delta_cindex = sd(dc, na.rm = TRUE),
      lo_delta_cindex = if (!is.null(tc) && !inherits(tc, "try-error")) tc$conf.int[1] else NA_real_,
      hi_delta_cindex = if (!is.null(tc) && !inherits(tc, "try-error")) tc$conf.int[2] else NA_real_,
      p_delta_cindex = if (!is.null(tc) && !inherits(tc, "try-error")) tc$p.value else NA_real_,
      mean_delta_AUC10 = mean(da, na.rm = TRUE), sd_delta_AUC10 = sd(da, na.rm = TRUE),
      lo_delta_AUC10 = if (!is.null(ta) && !inherits(ta, "try-error")) ta$conf.int[1] else NA_real_,
      hi_delta_AUC10 = if (!is.null(ta) && !inherits(ta, "try-error")) ta$conf.int[2] else NA_real_,
      p_delta_AUC10 = if (!is.null(ta) && !inherits(ta, "try-error")) ta$p.value else NA_real_)
  }, by = .(contrast, model_a, model_b)]

  pred_write_tab(pred, file.path(rawdir, fig_data_file(6, "prediction_rows.tsv.gz", Y = Y)))
  pred_write_tab(perf, file.path(rawdir, fig_data_file(6, "performance_by_fold.tsv", Y = Y)))
  pred_write_tab(perf_sum, file.path(rawdir, fig_data_file(6, "performance_summary.tsv", Y = Y)))
  pred_write_tab(perf_contrasts, file.path(rawdir, fig_data_file(6, "performance_contrasts.tsv", Y = Y)))
  pred_write_tab(perf_contrast_by_fold, file.path(rawdir, fig_data_file(6, "performance_contrasts_by_fold.tsv", Y = Y)))
  if (nrow(adaptive_tuning)) pred_write_tab(adaptive_tuning, file.path(rawdir, fig_data_file(6, "adaptive_penalty_tuning.tsv", Y = Y)))
  pred_write_tab(yy_params_all, file.path(rawdir, fig_data_file(7, "yy_parameters_by_fold.tsv.gz", Y = Y)))
  pred_write_tab(selected, file.path(rawdir, fig_data_file(6, "selected_features_by_fold.tsv.gz", Y = Y)))
  pred_write_tab(weights_all, file.path(rawdir, fig_data_file(7, "weights_by_fold.tsv.gz", Y = Y)))
  if (nrow(yy_guided_all)) pred_write_tab(yy_guided_all, file.path(rawdir, fig_data_file(6, "yy_guided_weights_by_fold.tsv.gz", Y = Y)))

  best_model <- perf_sum[model != "M0_base"][order(-mean_cindex)][1, model]
  if (length(best_model) && !is.na(best_model)) {
    risk <- copy(pred[model == best_model])
    risk[, lp_z := {s <- sd(lp, na.rm = TRUE); if (!is.finite(s) || s == 0) rep(0, .N) else (lp - mean(lp, na.rm = TRUE)) / s}, by = fold]
    risk[, decile := pmin(10L, pmax(1L, ceiling(10 * frank(lp_z, ties.method = "average") / .N))), by = fold]
    risk_hr <- tryCatch({
      d <- as.data.frame(risk[is.finite(decile)])
      d$decile <- relevel(factor(d$decile), ref = "1")
      fit <- coxph(Surv(time, event) ~ decile + strata(fold), data = d)
      sm <- summary(fit)
      data.table(term = rownames(sm$coef), HR = sm$conf.int[, "exp(coef)"],
                 lo = sm$conf.int[, "lower .95"], hi = sm$conf.int[, "upper .95"],
                 p = sm$coef[, "Pr(>|z|)"])
    }, error = function(e) data.table())
    pred_write_tab(risk, file.path(rawdir, fig_data_file(6, paste0("risk_deciles.", best_model, ".tsv.gz"), Y = Y)))
    if (nrow(risk_hr)) pred_write_tab(risk_hr, file.path(rawdir, fig_data_file(6, paste0("risk_decile_HR.", best_model, ".tsv"), Y = Y)))
  }

  model_desc <- data.table(
    model = c("M0_base", "M1_Yin_glmnet", "M2_YY_guided_score"),
    description = c(
      paste0("Base covariates only; mode=", pred_base_mode),
      paste0("Yin-only incident-risk elastic net after training-fold Cox screening; integrated candidate balance=", pred_layer_balance),
      "Explainable YY-guided protein score: Yin training-fold Cox beta provides the base direction/weight; Yang prevalence, recurrence, and pseudo-time evidence filters variables and modifies weights; no biomarker glmnet is fitted for this score"
    )
  )
  if (isTRUE(pred_legacy_complex)) {
    model_desc <- rbind(model_desc, data.table(
      model = c("M2_legacy_Yin_plus_Yang_scores", "M3_legacy_YY_adaptive_glmnet",
                "M4_legacy_YY_trajectory_glmnet", "M5_legacy_YY_reactive_glmnet",
                "M6_legacy_Yang_scores_only"),
      description = c(
        "Legacy diagnostic: identical Yin-screened raw feature set as M1 plus an unpenalized fold-learned YY score block",
        "Legacy diagnostic: identical Yin-screened raw feature set as M1 with Yang-state and Yang-trajectory evidence changing feature-specific elastic-net penalties",
        "Legacy diagnostic: identical Yin-screened raw feature set as M1 with only aligned pseudo-time trajectory evidence changing penalties",
        "Legacy diagnostic: raw biomarkers selected for concordant prospective-risk and disease-proximity/reactive patterns, fitted by glmnet",
        "Legacy diagnostic: base covariates plus Yang-only compact scores; no Yin incident-risk score"
      )
    ), fill = TRUE)
  }
  if (biom == "prot+met" && isTRUE(pred_layer_ablation)) {
    model_desc <- rbind(model_desc,
                        data.table(model = c("M1a_Yin_prot", "M1b_Yin_met"),
                                   description = c("Protein-only traditional Yin model evaluated in the same integrated cohort", "Metabolite-only traditional Yin model evaluated in the same integrated cohort")))
  }
  pred_write_tab(model_desc, file.path(rawdir, fig_data_file(6, "model_descriptions.tsv", Y = Y)))
  pred_write_book(list(cohort = cohort, performance_by_fold = perf, performance_summary = perf_sum,
                       Yang_increment_contrasts = perf_contrasts, contrasts_by_fold = perf_contrast_by_fold,
                       adaptive_penalty_tuning = adaptive_tuning, model_descriptions = model_desc,
                       yy_guided_weights_head = head(yy_guided_all[order(fold, -yy_filter_score)], 3000),
                       selected_head = head(selected, 3000), yy_params_head = head(yy_params_all[order(fold, cox_p)], 3000)),
                  file.path(outdir, fig_file(6, "prediction", ".out.xlsx", Y = Y)))
  saveRDS(list(cohort = cohort, pred = pred, perf = perf, perf_sum = perf_sum,
               perf_contrasts = perf_contrasts, perf_contrast_by_fold = perf_contrast_by_fold,
               adaptive_tuning = adaptive_tuning, yy_params = yy_params_all,
               yy_guided_weights = yy_guided_all, selected = selected, model_desc = model_desc),
          file.path(rawdir, fig_data_file(6, "prediction.rds", Y = Y)), compress = "xz")

  core_models <- c("M0_base", "M1_Yin_glmnet", "M2_YY_guided_score")
  if (isTRUE(pred_legacy_complex)) {
    core_models <- c(core_models,
      "M2_legacy_Yin_plus_Yang_scores",
      "M3_legacy_YY_adaptive_glmnet",
      "M4_legacy_YY_trajectory_glmnet",
      "M5_legacy_YY_reactive_glmnet",
      "M6_legacy_Yang_scores_only"
    )
  }
  model_diagnostics <- data.table()
  if (nrow(perf_contrast_by_fold)) {
    md_source <- perf_contrast_by_fold[model_b == "M1_Yin_glmnet" & model_a %in% c("M3_legacy_YY_adaptive_glmnet", "M4_legacy_YY_trajectory_glmnet")]
    if (nrow(md_source)) {
      model_diagnostics <- md_source[
        ,
        .(n_folds = .N,
          mean_delta_cindex = mean(delta_cindex, na.rm = TRUE),
          sd_delta_cindex = sd(delta_cindex, na.rm = TRUE),
          max_abs_delta_cindex = max(abs(delta_cindex), na.rm = TRUE),
          all_cindex_delta_zero = all(is.finite(delta_cindex) & abs(delta_cindex) < 1e-12),
          mean_delta_AUC10 = mean(delta_AUC10, na.rm = TRUE),
          max_abs_delta_AUC10 = max(abs(delta_AUC10), na.rm = TRUE),
          all_AUC10_delta_zero = all(is.finite(delta_AUC10) & abs(delta_AUC10) < 1e-12)),
        by = .(model = model_a)
      ]
    }
  }
  if (nrow(adaptive_tuning)) {
    eta_diag <- adaptive_tuning[selected_eta %in% TRUE, .(
      selected_eta_values = paste(sort(unique(eta)), collapse = ","),
      eta0_selected_folds = sum(eta == 0, na.rm = TRUE),
      nonzero_eta_selected_folds = sum(eta > 0, na.rm = TRUE),
      selected_eta_all_zero = all(eta == 0)
    ), by = model]
    model_diagnostics <- merge(model_diagnostics, eta_diag, by = "model", all = TRUE, sort = FALSE)
  }
  if (nrow(model_diagnostics)) {
    pred_write_tab(model_diagnostics, file.path(rawdir, fig_data_file(6, "model_diagnostics.tsv", Y = Y)))
  }

  core_labels <- c(
    M0_base = "M0  Base",
    M1_Yin_glmnet = "M1  Yin glmnet",
    M2_YY_guided_score = "M2  YY-guided score",
    M2_legacy_Yin_plus_Yang_scores = "L2  legacy Yin + Yang scores",
    M3_legacy_YY_adaptive_glmnet = "L3  legacy YY adaptive glmnet",
    M4_legacy_YY_trajectory_glmnet = "L4  legacy YY trajectory glmnet",
    M5_legacy_YY_reactive_glmnet = "L5  legacy YY reactive glmnet",
    M6_legacy_Yang_scores_only = "L6  legacy Yang scores only"
  )
  perf_plot <- copy(perf[model %in% core_models])
  perf_plot[, model_order := factor(model, levels = core_models, labels = unname(core_labels[core_models]))]

  # Fixed factor levels make the vertical order bottom-to-top M0, M1, ..., M6.
  p1 <- ggplot(perf_plot, aes(cindex, model_order, fill = model_order)) +
    geom_boxplot(alpha = .80, outlier.shape = NA, show.legend = FALSE) +
    geom_point(position = position_jitter(height = .08, width = 0), alpha = .80, size = 2, show.legend = FALSE) +
    labs(title = paste0("a. Cross-validated C-index for ", Y), x = "C-index", y = NULL) + theme_yy(11)

  p2_dat <- perf_plot[model != "M0_base"]
  p2_dat[, model_order := factor(model, levels = core_models[-1], labels = unname(core_labels[core_models[-1]]))]
  p2 <- ggplot(p2_dat, aes(model_order, delta_cindex, fill = model_order)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(alpha = .80, outlier.shape = NA, show.legend = FALSE) +
    geom_point(position = position_jitter(width = .06), alpha = .80, size = 2, show.legend = FALSE) +
    labs(title = "b. Incremental C-index vs M0 base", x = NULL, y = "Delta C-index") + theme_yy(11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  p3 <- ggplot(perf_plot, aes(AUC10, model_order, fill = model_order)) +
    geom_boxplot(alpha = .80, outlier.shape = NA, show.legend = FALSE) +
    geom_point(position = position_jitter(height = .08, width = 0), alpha = .80, size = 2, show.legend = FALSE) +
    labs(title = "c. 10-year known-status AUC",
         subtitle = "Participants censored before 10 years without an event are excluded; C-index remains the primary metric",
         x = "AUC", y = NULL) + theme_yy(11)

  yy_compare_models <- "M2_YY_guided_score"
  if (isTRUE(pred_legacy_complex)) {
    yy_compare_models <- c(
      yy_compare_models,
      "M2_legacy_Yin_plus_Yang_scores", "M3_legacy_YY_adaptive_glmnet",
      "M4_legacy_YY_trajectory_glmnet", "M5_legacy_YY_reactive_glmnet"
    )
  }
  yy_vs_yin <- copy(perf_contrast_by_fold[model_b == "M1_Yin_glmnet" & model_a %in% yy_compare_models])
  yy_vs_yin[, model_order := factor(model_a, levels = yy_compare_models, labels = unname(core_labels[yy_compare_models]))]
  if (nrow(yy_vs_yin)) {
    yy_mean <- yy_vs_yin[, .(mean_delta = mean(delta_cindex, na.rm = TRUE)), by = model_order]
    zero_models <- character()
    zero_eta_models <- character()
    if (nrow(model_diagnostics)) {
      zero_models <- model_diagnostics[all_cindex_delta_zero %in% TRUE & selected_eta_all_zero %in% TRUE, model]
      zero_eta_models <- unname(core_labels[intersect(zero_models, names(core_labels))])
      zero_models <- model_diagnostics[all_cindex_delta_zero %in% TRUE, model]
      zero_models <- unname(core_labels[intersect(zero_models, names(core_labels))])
    }
    p4_subtitle <- "Each point is one outer-CV fold; positive values favor the YY-informed model"
    if (length(zero_models)) {
      note <- paste0(paste(zero_models, collapse = " and "), " produced predictions identical to M1")
      if (length(zero_eta_models)) {
        note <- paste0(note, " (eta=0 in all folds for ", paste(zero_eta_models, collapse = " and "), ")")
      }
      p4_subtitle <- paste0(p4_subtitle, "; ", note)
    }
    p4 <- ggplot(yy_vs_yin, aes(model_order, delta_cindex, fill = model_order)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
      geom_boxplot(alpha = .80, outlier.shape = NA, show.legend = FALSE) +
      geom_point(aes(group = fold), position = position_jitter(width = .035), alpha = .85, size = 2, show.legend = FALSE) +
      geom_text(data = yy_mean, aes(model_order, mean_delta, label = sprintf("mean %.4f", mean_delta)),
                vjust = ifelse(yy_mean$mean_delta >= 0, -1.0, 1.4), size = 3.0, inherit.aes = FALSE) +
      labs(title = "d. Paired YY increment over M1 Yin model",
           subtitle = p4_subtitle,
           x = NULL, y = "Delta C-index vs M1") + theme_yy(10) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1))
  } else {
    p4 <- ggplot() + labs(title = "d. YY-vs-M1 contrasts unavailable") + theme_yy(10)
  }
  pred_save_plot((p1 | p2) / (p3 | p4), file.path(outdir, fig_file(6, "performance", Y = Y)), 15.5, 11)

  # Selection stability remains available as a separate diagnostic rather than
  # crowding the main performance figure.
  sel_top <- selected[
    model %in% core_models & selection_type == "glmnet_nonzero" & protein %in% features,
    .(folds_selected = uniqueN(fold)), by = .(model, protein)
  ]
  if (nrow(sel_top)) {
    sel_top <- sel_top[order(model, -folds_selected)]
    sel_top <- sel_top[, head(.SD, 10), by = model]
    sel_top[, model_order := factor(model, levels = core_models, labels = unname(core_labels[core_models]))]
    p_sel <- ggplot(sel_top, aes(folds_selected, reorder(protein, folds_selected), fill = model_order)) +
      geom_col(width = .70, show.legend = FALSE) +
      facet_wrap(~ model_order, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = seq_len(K_outer), limits = c(0, K_outer)) +
      labs(title = "Biomarker selection stability across outer folds",
           subtitle = "Top 10 nonzero biomarkers within each model",
           x = "Folds selected", y = NULL) + theme_yy(9)
    pred_save_plot(p_sel, file.path(outdir, fig_file(6, "selection_stability", Y = Y)), 11, 12)
    pred_write_tab(sel_top, file.path(rawdir, fig_data_file(6, "selection_stability.tsv", Y = Y)))
  }

  robust_z <- function(x) {
    med <- median(x, na.rm = TRUE)
    sc <- mad(x, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(sc) || sc <= 0) sc <- sd(x, na.rm = TRUE)
    if (!is.finite(sc) || sc <= 0) sc <- 1
    (x - med) / sc
  }
  yy_map_all <- yy_params_all[, .(
    n_folds = uniqueN(fold),
    cox_beta = median(cox_beta, na.rm = TRUE),
    cox_p = median(cox_p, na.rm = TRUE),
    yy_beta10 = median(yy_beta10, na.rm = TRUE),
    yy_p = median(yy_p, na.rm = TRUE),
    prev_vs_never_beta = median(prev_vs_never_beta, na.rm = TRUE),
    prev_vs_never_p = median(prev_vs_never_p, na.rm = TRUE),
    recurrent_vs_nonrecurrent_beta = median(recurrent_vs_nonrecurrent_beta, na.rm = TRUE),
    recurrent_vs_nonrecurrent_p = median(recurrent_vs_nonrecurrent_p, na.rm = TRUE),
    yy_rank_score = median(yy_rank_score, na.rm = TRUE),
    category = {zz <- category[!is.na(category)]; if (length(zz)) zz[1] else "Other/unclassified"},
    prior = {zz <- prior[!is.na(prior)]; if (length(zz)) zz[1] else "other"}
  ), by = protein][is.finite(cox_p) | is.finite(yy_p)]
  yy_map_all[, quadrant := fcase(
    is.finite(yy_beta10) & yy_beta10 < 0 & is.finite(cox_beta) & cox_beta > 0 & yy_p < .05 & cox_p < .05, "Reactive/risk",
    is.finite(yy_beta10) & yy_beta10 > 0 & is.finite(cox_beta) & cox_beta < 0 & yy_p < .05 & cox_p < .05, "Mirror/protective",
    is.finite(cox_beta) & cox_p < .05 & (!is.finite(yy_p) | yy_p >= .05), "Prospective-only",
    is.finite(yy_beta10) & yy_p < .05 & (!is.finite(cox_p) | cox_p >= .05), "Timeline-only",
    default = "Mixed/uncertain"
  )]
  yy_map_all[, `:=`(
    yy_robust_z = robust_z(yy_beta10),
    cox_robust_z = robust_z(cox_beta)
  )]
  yy_map_all[, outlier_score := sqrt(yy_robust_z^2 + cox_robust_z^2)]
  yy_map_all[, mirror_score := fifelse(quadrant == "Mirror/protective",
                                       pmin(-log10(pmax(cox_p, 1e-300)), 50) + pmin(-log10(pmax(yy_p, 1e-300)), 50),
                                       NA_real_)]
  value_features <- intersect(yy_map_all$protein, names(datY))
  value_diag <- rbindlist(lapply(value_features, function(p) {
    x <- safe_num(datY[[p]])
    ok <- is.finite(x)
    if (!any(ok)) {
      return(data.table(protein = p, value_nonmissing = 0L, value_missing_rate = mean(!ok),
                        value_mean = NA_real_, value_sd = NA_real_,
                        value_q001 = NA_real_, value_q01 = NA_real_, value_median = NA_real_,
                        value_q99 = NA_real_, value_q999 = NA_real_,
                        value_max_abs_robust_z = NA_real_, value_frac_abs_robust_z_gt5 = NA_real_))
    }
    qs <- quantile(x[ok], probs = c(.001, .01, .5, .99, .999), na.rm = TRUE, names = FALSE)
    rz <- robust_z(x)
    data.table(
      protein = p,
      value_nonmissing = sum(ok),
      value_missing_rate = mean(!ok),
      value_mean = mean(x[ok], na.rm = TRUE),
      value_sd = sd(x[ok], na.rm = TRUE),
      value_q001 = qs[1],
      value_q01 = qs[2],
      value_median = qs[3],
      value_q99 = qs[4],
      value_q999 = qs[5],
      value_max_abs_robust_z = if (any(is.finite(rz))) max(abs(rz), na.rm = TRUE) else NA_real_,
      value_frac_abs_robust_z_gt5 = mean(abs(rz[ok]) > 5, na.rm = TRUE)
    )
  }), fill = TRUE)
  if (nrow(value_diag)) yy_map_all <- merge(yy_map_all, value_diag, by = "protein", all.x = TRUE, sort = FALSE)
  yy_map_all <- yy_map_all[order(cox_p, yy_p)]
  outlier_diag <- yy_map_all[is.finite(outlier_score)][order(-outlier_score)][seq_len(min(50L, .N))]
  if (nrow(outlier_diag)) pred_write_tab(outlier_diag, file.path(rawdir, fig_data_file(7, "yy_parameter_outliers.tsv", Y = Y)))
  pred_write_tab(yy_map_all, file.path(rawdir, fig_data_file(7, "yy_parameter_map_summary.tsv", Y = Y)))

  keep_map <- unique(c(
    "PCSK9", "GDF15", "RAB6A",
    yy_map_all[order(cox_p, yy_p)][seq_len(min(55L, .N))]$protein,
    yy_map_all[quadrant == "Mirror/protective"][order(-mirror_score)][seq_len(min(12L, .N))]$protein,
    yy_map_all[quadrant == "Reactive/risk"][order(cox_p, yy_p)][seq_len(min(18L, .N))]$protein,
    outlier_diag[seq_len(min(10L, .N)), protein]
  ))
  yy_show <- yy_map_all[protein %in% keep_map][order(cox_p, yy_p)][seq_len(min(90L, .N))]
  label_map <- unique(c("PCSK9", "GDF15", "RAB6A",
                        yy_show[quadrant == "Mirror/protective"][order(-mirror_score)][seq_len(min(10L, .N))]$protein,
                        yy_show[order(cox_p, yy_p)][seq_len(min(35L, .N))]$protein,
                        outlier_diag[seq_len(min(8L, .N)), protein]))
  yy_show[, map_label := fifelse(protein %in% label_map, protein, "")]
  quadrant_cols2 <- c(
    `Reactive/risk` = "#D95F02",
    `Mirror/protective` = "#0072B2",
    `Prospective-only` = "#009E73",
    `Timeline-only` = "#CC79A7",
    `Mixed/uncertain` = "grey65"
  )
  p5 <- ggplot(yy_show, aes(yy_beta10, cox_beta, color = quadrant, label = map_label, size = yy_rank_score)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = .86) + ggrepel::geom_text_repel(size = 2.6, max.overlaps = 38, seed = 7) +
    scale_color_manual(values = quadrant_cols2, na.value = "grey70") +
    labs(title = paste0("YY modelling parameters learned across folds: ", Y),
         subtitle = "Mirror/protective proteins have positive YY trend but negative incident Cox beta; see raw Fig7 outlier diagnostics",
         x = "Median YY trend beta per 10-year later diagnosis", y = "Median incident Cox beta",
         color = "Pattern", size = "YY rank") + theme_yy(10) + theme(legend.position = "bottom")
  pred_save_plot(p5, file.path(outdir, fig_file(7, "yy_parameter_map", Y = Y)), 10.5, 8)

  invisible(list(cohort = cohort, performance = perf, performance_summary = perf_sum,
                 performance_contrasts = perf_contrasts, yy_params = yy_params_all, selected = selected))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

make_consolidate <- function(results, Y) {
  set_scope("trait", Y)
  files <- data.table(file = list.files(outdir, recursive = TRUE, full.names = FALSE))
  files[, size_MB := round(file.info(file.path(outdir, file))$size / 1024^2, 3)]
  summary <- data.table(
    item = c("Outcome", "Biomarker mode", "Date source", "Fig standardization", "Trend p shown", "Minimum bin N", "Prediction residualization", "Prediction base mode", "Available LE8 covariates", "Adaptive eta grid", "Layer balance", "Main interpretation"),
    value = c(Y, biom, date_source, fig_standardize, trend_p_show, min_bin_n, pred_residualize, pred_base_mode, paste(get_le8_covars(dat), collapse = ";"), paste(pred_adaptive_eta, collapse = ","), pred_layer_balance,
              "Root figures summarize the measured biomarker cohort and metadata groups; disease Fig1-Fig7 separate aligned pseudo-time trajectories, AOD recurrence, adjustment, specificity, temporal decomposition, prediction, and fold-learned YY parameters. Prediction is evaluated only in held-out baseline disease-free Yin participants. The primary YY model is an explainable score: Yin training-fold Cox beta provides direction, while Yang prevalence/recurrence/pseudo-time evidence filters variables and modifies weights. Legacy complex YY glmnet diagnostics are only run when YYP_LEGACY_COMPLEX=TRUE.")
  )
  write_book(list(run_summary = summary, output_index = files), summary_file(Y, ".xlsx"))
  saveRDS(list(summary = summary, files = files, results = results), file.path(rawdir, summary_file(Y, ".rds")), compress = "xz")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run pipeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

main <- read_main_dataset()
dat <- main$dat
vars_biom_global <- main$vars_biom
vars_biom_by_layer_global <- main$vars_by_layer
biomarker_meta_global <- as.data.table(main$metadata)
set_default_display_features(vars_biom_global, biomarker_meta_global)
protein_prior <- annotate_proteins(vars_biom_global)[, .(protein, prior, category)]
le8_covars_global <- get_le8_covars(dat)
message("Analysis cohort N = ", nrow(dat), "; features = ", length(vars_biom_global),
        "; selected display features = ", paste(vip_proteins, collapse = ", "))
message("LE8 covariates found (", length(le8_covars_global), ") = ",
        ifelse(length(le8_covars_global), paste(le8_covars_global, collapse = ", "), "none"))
if (pred_base_mode %in% c("le8", "full") && length(le8_covars_global) < 8) {
  warning("Only ", length(le8_covars_global), " LE8 component/total variables were found: ",
          paste(le8_covars_global, collapse = ", "), call. = FALSE)
}

if (run_step("data_prep")) {
  set_scope("root")
  data_meta <- rbindlist(list(
    data.table(item = c("all_rds_N", "all_rds_after_filter_N", "analysis_cohort_N", "ethnic_filter", "biomarker_columns_N", "selected_features_present"),
               value = as.character(c(main$all_n, main$all_filtered_n, nrow(dat), ethnic_keep, length(vars_biom_global), paste(intersect(vip_proteins, names(dat)), collapse = ",")))),
    main$layer_counts[, .(item = paste0(layer, c("_rds_N")), value = as.character(rds_N))]
  ), fill = TRUE)
  write_tab(data_meta, paste0("yy_data_meta.", biom, ".tsv"))
  write_tab(biomarker_meta_global, paste0("yy_biomarker_metadata.", biom, ".tsv"))
  saveRDS(list(meta = data_meta, vars_biom = vars_biom_global, vars_by_layer = vars_biom_by_layer_global,
               biomarker_metadata = biomarker_meta_global),
          file.path(rawroot, paste0("yy_data_meta.", biom, ".rds")), compress = "xz")
}

run_if("root_qc", {
  make_root_figures(dat, vars_biom_global, biom_n = main$biom_n, all_n = main$all_n,
                    all_filtered_n = main$all_filtered_n, layer_counts = main$layer_counts)
})

run_if("aod_qc", {
  make_root_aod_qc_figure()
})

all_results <- list()
for (Y in Ys) {
  message("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  message("# RUN OUTCOME: ", Y)
  message("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  cleanup_legacy_trait_outputs(Y)
  datY <- prepare_trait_data(dat, Y)
  message("Outcome ", Y, ": cases within [-16,16] = ", sum(is.finite(datY$b2e) & datY$b2e >= -16 & datY$b2e <= 16, na.rm = TRUE),
          "; pre-baseline = ", sum(is.finite(datY$b2e) & datY$b2e < 0 & datY$b2e >= -16, na.rm = TRUE),
          "; post-baseline = ", sum(is.finite(datY$b2e) & datY$b2e > 0 & datY$b2e <= 16, na.rm = TRUE))
  invisible(gc())
  resY <- list()
  run_if("yy_timeline", {
    plot_proteins <- choose_yy_plot_proteins(datY, Y)
    message("YY timeline display proteins for ", Y, ": ", paste(plot_proteins, collapse = ", "))
    yy_res <- make_yy_timeline_figure(datY, Y, plot_proteins = plot_proteins)
    if (!is.null(yy_res$dat)) datY <- yy_res$dat
    yy_res$dat <- NULL
    resY$yy_timeline <- yy_res
  })
  run_if("aod_recurrence", { resY$aod_recurrence <- trim_result_payload(make_aod_recurrence_figure(datY, Y)) })
  run_if("adjustment", { resY$adjustment <- trim_result_payload(make_adjustment_figure(datY, Y)) })
  run_if("temporal_decomp", { resY$temporal_decomp <- make_temporal_decomp_figure(datY, Y) })
  run_if("scan", { resY$scan <- make_scan_figure(datY, dat, Y) })
  run_if("prediction", {
    if (isTRUE(skip_prediction)) {
      message("Prediction step skipped because YY_SKIP_PREDICTION=TRUE")
    } else {
      resY$prediction <- trim_result_payload(run_yy_prediction(datY, vars_biom_global, Y))
    }
  })
  resY <- trim_result_payload(resY)
  run_if("consolidate", { make_consolidate(resY, Y) })
  all_results[[Y]] <- resY
  rm(datY, resY)
  invisible(gc())
}

set_scope("root")
root_files <- data.table(file = list.files(outroot, recursive = TRUE, full.names = FALSE))
root_files[, size_MB := round(file.info(file.path(outroot, file))$size / 1024^2, 3)]
write_book(list(outputs = root_files), paste0("yy_all_outputs.", biom, ".xlsx"))
saveRDS(all_results, file.path(rawroot, paste0("yy_all_results.", biom, ".rds")), compress = "xz")
message("Done. Outputs written under: ", outroot)
