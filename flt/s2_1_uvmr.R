# Bidirectional two-sample MR from local, harmonised-format GWAS files.
# Run with Rscript s2_1_uvmr.R, source(), or RStudio Source.
# Required: data.table, R.utils, TwoSampleMR, openxlsx; PLINK 1.9 + EUR LD panel.
# The two directions select instruments independently from their own exposure.
# Output includes five-method estimates and combined sensitivity tables for each direction.

# 1. Settings ---------------------------------------------------------------
win <- Sys.info()[["sysname"]] == "Windows"
dir0 <- if (win) "D:" else "/mnt/d"	# "D:/用户/sph-xwy" else "/mnt/d/用户/sph-xwy"
data_dirs <- c(
  flt = file.path(dir0, "data/gwas/flt/common"),
  covs = file.path(dir0, "data/gwas/covs/common"),
  brain = file.path(dir0, "data/gwas/brain/common")
)
outdir <- file.path(dir0, "projects/flt/analysis/02uvmr")
instrument_cache_dir <- file.path(outdir, "clump")

# A group, one ID, or a vector of IDs; group/ID resolves duplicated IDs.
trait1 <- "flt"                    # e.g. "flt", "flt/G755", c("flt/G755", "flt/G756")
trait2 <- "brain"                   # e.g. "brain/ad"
direction <- "both"                 # "forward", "reverse", or "both"
# forward = trait1 -> trait2; reverse = trait2 -> trait1.
fdr_scope <- "family"               # "family": each trait2 separately; "direction": all pairs

plink_path <- if (win) "software/windows/plink/plink.exe" else "software/linux/plink/plink"
plink <- Sys.getenv("PLINK_BIN", file.path(dir0, plink_path))
ld_prefix <- Sys.getenv("MR_LD_PREFIX", file.path(dir0, "data/ldref/opengwas/1kg.v3/EUR"))
p_cut <- 5e-8
clump_kb <- 10000
clump_r2 <- 0.001
min_f <- 10                         # Per-SNP F=(BETA/SE)^2, after clumping; 0 disables filtering.
reuse_clump <- TRUE
harmonise_action <- 2               # Allele-frequency-assisted palindromic SNP handling.
seed <- 12345L

run_bwmr <- TRUE                    # Retains the additional method in the original script.
run_presso <- FALSE                 # Optional and computationally expensive.
run_loo <- FALSE
save_harmonised <- FALSE            # TRUE writes one combined RDS, not one file per pair.
bwmr_file <- file.path(instrument_cache_dir, "BWMR_updated.R")
bwmr_url <- "https://raw.githubusercontent.com/jiazhao97/BWMR/master/updates/BWMR_updated.R"

# Used only to supply metadata; missing metadata does not block the main MR.
brain_info_file <- file.path(dir0, "data/gwas/brain/gwas_brain.xlsx")
brain_info_start_row <- 2L
# Optional CSV with columns key,type,N,ncase,ncontrol,prevalence.
# Example key: brain/ad or covs/smoking. Overrides the Excel/default metadata.
trait_info_file <- file.path(outdir, "trait_info.csv")
# Defaults: flt=continuous, brain=binary, covs=continuous. Explicitly override
brain_ids <- c("dementia_r13", "anystroke_eur", "sleep_disorder_r13", "mdd_2025", "anxiety_2026")

# 2. Trait selection and metadata -------------------------------------------
discover_traits <- function(dirs) {
  rbindlist(lapply(names(dirs), function(g) {
    ids <- sort(list.dirs(dirs[[g]], recursive = FALSE, full.names = FALSE))
    files <- file.path(dirs[[g]], ids, "gwas", paste0(ids, ".gz"))
    keep <- file.exists(files)
    d <- data.table(group = rep(g, sum(keep)), id = ids[keep], file = files[keep])
    d[, key := if (.N) paste(group, id, sep = "/") else character()]
    d
  }), use.names = TRUE)
}

select_traits <- function(selection, registry) {
  if (!length(selection) || anyNA(selection)) stop("trait1/trait2 must contain a group or trait ID.")
  keys <- unlist(lapply(selection, function(s) {
    if (s %in% names(data_dirs)) {
      hit <- registry[group == s]
      if (s == "brain") {
        missing <- setdiff(brain_ids, hit$id)
        if (length(missing)) message("brain group: unavailable IDs omitted: ", paste(missing, collapse = ", "))
        hit <- hit[id %in% brain_ids]
      }
    } else if (grepl("/", s, fixed = TRUE)) {
      hit <- registry[key == s]
    } else {
      hit <- registry[id == s]
      if (nrow(hit) > 1L) stop("Ambiguous ID; use group/ID: ", s)
    }
    if (!nrow(hit)) stop("No GWAS found for selection: ", s)
    hit$key
  }), use.names = FALSE)
  registry[match(unique(keys), key)]
}

make_pairs <- function(a, b, mode = direction) {
  if (!mode %in% c("forward", "reverse", "both")) stop("direction must be forward, reverse, or both.")
  if (length(intersect(a$key, b$key))) stop("trait1 and trait2 overlap; remove shared traits.")
  base <- CJ(trait1_key = a$key, trait2_key = b$key, unique = TRUE)
  forward <- copy(base)[, `:=`(direction = "forward", exposure_key = trait1_key, outcome_key = trait2_key)]
  reverse <- copy(base)[, `:=`(direction = "reverse", exposure_key = trait2_key, outcome_key = trait1_key)]
  z <- switch(mode, forward = forward, reverse = reverse, both = rbind(forward, reverse))
  reg <- rbind(a, b)
  z[, `:=`(pair_id = seq_len(.N), family = trait2_key,
    stage = ifelse(direction == "forward", "uvmr", "reverse"),
    exposure = reg$id[match(exposure_key, reg$key)], outcome = reg$id[match(outcome_key, reg$key)],
    exposure_group = reg$group[match(exposure_key, reg$key)], outcome_group = reg$group[match(outcome_key, reg$key)])]
  z[]
}

number <- function(x) suppressWarnings(as.numeric(gsub("[,[:space:]]", "", as.character(x))))
population_prevalence <- function(x) {
  is_percent <- grepl("%", as.character(x), fixed = TRUE)
  z <- number(gsub("%", "", as.character(x), fixed = TRUE))
  z[is_percent & !is.na(is_percent)] <- z[is_percent & !is.na(is_percent)] / 100
  z
}

load_trait_info <- function(traits) {
  info <- copy(traits[, .(key, group, id)])
  info[, `:=`(type = ifelse(group == "brain", "binary", "continuous"),
    N = NA_real_, ncase = NA_real_, ncontrol = NA_real_, prevalence = NA_real_, source = "group default")]
  if (any(info$group == "brain") && file.exists(brain_info_file)) {
    meta <- tryCatch(as.data.table(openxlsx::read.xlsx(brain_info_file, startRow = brain_info_start_row)),
      error = function(e) {message("Brain metadata skipped: ", conditionMessage(e)); NULL})
    if (!is.null(meta)) {
      canonical <- gsub("[^a-z0-9]", "", tolower(names(meta)))
      column <- function(label) {
        # Prefer the literal header: removing Chinese text can turn a
        # descriptive column such as prevalence.依据类型 into prevalence too.
        j <- which(tolower(trimws(names(meta))) == label)
        if (!length(j)) j <- which(canonical == label)
        if (length(j) == 1L) meta[[j]] else rep(NA_character_, nrow(meta))
      }
      ids <- trimws(as.character(column("outcomeshort")))
      aliases <- c("Alzheimer's disease" = "ad", "Alzheimer disease" = "ad", Dementia = "dementia",
        Dystonia = "dystonia", Encephalitis = "encephalitis", "Multiple sclerosis" = "ms",
        Epilepsy = "epilepsy", "Sleep disorder" = "sleep", "Sleep disorders" = "sleep",
        Schizophrenia = "scz", "Bipolar disorder" = "bd", "Anxiety disorders" = "anxiety",
        "Major depressive disorder" = "mdd", "Parkinson's disease" = "pd")
      blank <- is.na(ids) | ids == ""
      ids[blank] <- unname(aliases[trimws(as.character(column("outcome")))[blank]])
      # Do not guess that a non-EUR stroke row describes the *_eur GWAS.
      counts <- strsplit(as.character(column("ncasencontrol")), "/", fixed = TRUE)
      cases <- vapply(counts, function(x) if (length(x) == 2L) number(x[1L]) else NA_real_, numeric(1))
      controls <- vapply(counts, function(x) if (length(x) == 2L) number(x[2L]) else NA_real_, numeric(1))
      m <- data.table(trait_key = paste0("brain/", ids), N = number(column("samplesizen")),
        ncase = cases, ncontrol = controls, prevalence = population_prevalence(column("prevalence")))
      setnames(m, "trait_key", "key")
      m <- m[key %in% info$key]
      dup <- unique(m$key[duplicated(m$key)])
      if (length(dup)) {
        message("Ambiguous Excel metadata skipped (set exact rows in trait_info.csv): ", paste(dup, collapse = ", "))
        m <- m[!key %in% dup]
      }
      for (j in seq_len(nrow(m))) {
        i <- match(m$key[j], info$key)
        for (v in c("N", "ncase", "ncontrol", "prevalence")) set(info, i, v, m[[v]][j])
        set(info, i, "source", "brain Excel")
      }
    }
  } else if (any(info$group == "brain")) {
    message("Brain metadata Excel not found; main MR will continue: ", brain_info_file)
  }
  if (file.exists(trait_info_file)) {
    m <- fread(trait_info_file, na.strings = c("", "NA"), showProgress = FALSE)
    if (!"key" %in% names(m)) {
      if (!all(c("group", "id") %in% names(m))) stop("trait_info.csv needs key (group/ID), or group and id columns.")
      m[, key := paste(group, id, sep = "/")]
    }
    m[, key := trimws(as.character(key))]
    m <- m[key %in% info$key]
    if (anyDuplicated(m$key)) stop("Duplicate selected keys in ", trait_info_file)
    for (j in seq_len(nrow(m))) {
      i <- match(m$key[j], info$key)
      for (v in intersect(c("type", "N", "ncase", "ncontrol", "prevalence"), names(m))) {
        val <- if (v == "type") tolower(trimws(as.character(m[[v]][j]))) else if (v == "prevalence")
          population_prevalence(m[[v]][j]) else number(m[[v]][j])
        if (!is.na(val)) set(info, i, v, val)
      }
      set(info, i, "source", "trait_info.csv override")
    }
  }
  if (any(!info$type %in% c("continuous", "binary"))) stop("Metadata type must be continuous or binary.")
  for (v in c("N", "ncase", "ncontrol")) info[!is.finite(get(v)) | get(v) <= 0, (v) := NA_real_]
  bad_prev <- info$type == "binary" & !is.na(info$prevalence) &
    (!is.finite(info$prevalence) | info$prevalence <= 0 | info$prevalence >= 1)
  if (any(bad_prev)) {
    message("Invalid population prevalence set to NA: ", paste(info$key[bad_prev], collapse = ", "))
    info[bad_prev, prevalence := NA_real_]
  }
  info[is.na(N) & is.finite(ncase) & is.finite(ncontrol), N := ncase + ncontrol]
  if (any(info$group == "covs" & info$source == "group default"))
    message("covs defaults to continuous; set binary covariates explicitly in ", trait_info_file)
  info[]
}

# 3. Read one GWAS at a time; prepare and cache independent instruments ------
read_gwas <- function(file, snps = NULL, p_max = NULL) {
  required <- c("SNP", "CHR", "EA", "NEA", "BETA", "SE", "P")
  con <- gzfile(file, open = "rt")
  header <- tryCatch(readLines(con, n = 1L, warn = FALSE), finally = close(con))
  if (!length(header)) stop("Empty GWAS: ", file)
  header_cols <- names(fread(text = paste0(header, "\n"), nrows = 0L, showProgress = FALSE))
  missing <- setdiff(required, header_cols)
  if (length(missing)) stop("GWAS missing columns: ", paste(missing, collapse = ", "), "; ", file)
  d <- fread(file, select = c(required, intersect(c("EAF", "N"), header_cols)), showProgress = FALSE)
  d[, SNP := as.character(SNP)]
  if (!is.null(snps)) d <- d[SNP %chin% snps]
  for (v in c("EAF", "N")) if (!v %in% names(d)) d[, (v) := NA_real_]
  for (v in c("BETA", "SE", "P", "EAF", "N")) set(d, j = v, value = number(d[[v]]))
  d[, `:=`(CHR = suppressWarnings(as.integer(sub("^chr", "", as.character(CHR), ignore.case = TRUE))),
           EA = toupper(as.character(EA)), NEA = toupper(as.character(NEA)))]
  d <- d[!is.na(SNP) & grepl("^rs[0-9]+$", SNP) & CHR %in% 1:22 &
    EA %in% c("A", "C", "G", "T") & NEA %in% c("A", "C", "G", "T") & EA != NEA &
    is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P >= 0 & P <= 1]
  # Drop ambiguous duplicated rsIDs before selecting significant exposure rows.
  d <- d[!duplicated(SNP) & !duplicated(SNP, fromLast = TRUE)]
  d[!is.finite(EAF) | EAF <= 0 | EAF >= 1, EAF := NA_real_]
  d[!is.finite(N) | N <= 2, N := NA_real_]
  if (!is.null(p_max)) d <- d[P < p_max]
  d[, P := pmax(P, 1e-300)]
  setcolorder(d, c("SNP", "CHR", "EA", "NEA", "BETA", "SE", "P", "EAF", "N"))
  d[]
}

instrument_stamp <- function(file) {
  paths <- normalizePath(c(file, plink, paste0(ld_prefix, c(".bed", ".bim", ".fam"))), winslash = "/", mustWork = TRUE)
  stat <- file.info(paths)
  stamp <- list(version = 3L, paths = paths, size = stat$size, mtime = as.numeric(stat$mtime),
                p_cut = p_cut, clump_kb = clump_kb, clump_r2 = clump_r2)
  paste0("# inputs: ", paste(deparse(stamp, width.cutoff = 500L, control = "all"), collapse = " "))
}

clump <- function(d) {
  if (!nrow(d)) return(character())
  work_dir <- tempfile("uvmr-plink-")
  if (!dir.create(work_dir)) stop("Cannot create PLINK temporary directory.")
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  prefix <- file.path(work_dir, "plink")
  input <- paste0(prefix, ".tsv")
  console <- paste0(prefix, ".console.txt")
  fwrite(d[order(P), .(SNP, P)], input, sep = "\t")
  status <- system2(plink, c("--bfile", shQuote(ld_prefix), "--clump", shQuote(input),
    "--clump-p1", as.character(p_cut), "--clump-p2", "1", "--clump-r2", as.character(clump_r2),
    "--clump-kb", as.character(clump_kb), "--out", shQuote(prefix)), stdout = console, stderr = console)
  log <- if (file.exists(console)) readLines(console, warn = FALSE) else character()
  if (status != 0L) stop("PLINK failed (", status, "): ", paste(tail(log, 12L), collapse = "\n"))
  output <- paste0(prefix, ".clumped")
  if (!file.exists(output)) {
    if (any(grepl("No significant --clump results|No clumps", log, ignore.case = TRUE))) return(character())
    stop("PLINK returned no .clumped file: ", paste(tail(log, 12L), collapse = "\n"))
  }
  hits <- fread(output, showProgress = FALSE)
  if (!"SNP" %in% names(hits)) stop("PLINK output has no SNP column; use PLINK 1.9.")
  unique(as.character(hits$SNP))
}

get_instruments <- function(trait) {
  cache_file <- file.path(instrument_cache_dir, paste0(trait$group, "__", trait$id, "_clump.txt"))
  stamp <- instrument_stamp(trait$file)
  z <- NULL
  if (reuse_clump && file.exists(cache_file)) {
    z <- tryCatch({
      meta <- readLines(cache_file, n = 3L, warn = FALSE)
      if (!identical(meta[1L], stamp)) stop("stale cache")
      n_sig <- as.integer(sub("^# n_significant=", "", meta[2L]))
      n_clumped <- as.integer(sub("^# n_clumped=", "", meta[3L]))
      d <- fread(cache_file, skip = 3L, sep = "\t", na.strings = "NA", showProgress = FALSE,
        colClasses = list(character = c("SNP", "EA", "NEA"), integer = "CHR",
                          numeric = c("BETA", "SE", "P", "EAF", "N")))
      if (!setequal(names(d), c("SNP", "CHR", "EA", "NEA", "BETA", "SE", "P", "EAF", "N")) ||
          anyNA(c(n_sig, n_clumped)) || n_clumped < 0 || n_sig < n_clumped ||
          nrow(d) != n_clumped || anyNA(d$SNP) || anyDuplicated(d$SNP)) stop("invalid cache")
      list(n_significant = n_sig, n_clumped = n_clumped, data = d)
    }, error = function(e) NULL)
    if (!is.null(z)) message("  Instrument cache: ", trait$key)
  }
  if (is.null(z)) {
    message("  Select and clump: ", trait$key)
    significant <- read_gwas(trait$file, p_max = p_cut)
    snps <- clump(significant)
    d <- significant[SNP %chin% snps]
    z <- list(n_significant = nrow(significant), n_clumped = nrow(d), data = d)
    tmp <- tempfile(".iv-", tmpdir = instrument_cache_dir)
    on.exit(unlink(tmp), add = TRUE)
    writeLines(c(stamp, paste0("# n_significant=", z$n_significant), paste0("# n_clumped=", z$n_clumped)), tmp)
    fwrite(d, tmp, sep = "\t", na = "NA", append = TRUE, col.names = TRUE)
    if (!suppressWarnings(file.rename(tmp, cache_file)) && !file.copy(tmp, cache_file, overwrite = TRUE))
      warning("Could not save instrument cache: ", cache_file)
  }
  # Filtering after cache loading keeps min_f changes independent of clumping.
  z$data <- z$data[is.finite((BETA / SE)^2) & (BETA / SE)^2 >= min_f]
  z$n_fpass <- nrow(z$data)
  z
}

# 4. MR methods, harmonisation and diagnostics --------------------------------
format_mr_data <- function(d, key, side, info) {
  if (!nrow(d)) stop("No matching SNPs: ", key)
  z <- as.data.frame(copy(d))
  z$trait <- z$id <- key
  fallback <- info$N[match(key, info$key)]
  z$N[!is.finite(z$N) | z$N <= 2] <- fallback
  suppressMessages(TwoSampleMR::format_data(z, type = side, phenotype_col = "trait", id_col = "id",
    snp_col = "SNP", beta_col = "BETA", se_col = "SE", effect_allele_col = "EA",
    other_allele_col = "NEA", eaf_col = "EAF", pval_col = "P", samplesize_col = "N", min_pval = 1e-300))
}

load_bwmr <- function() {
  if (!run_bwmr) return(list(fun = NULL, status = "disabled"))
  tryCatch({
    if (!file.exists(bwmr_file)) {
      tmp <- tempfile(".bwmr-", tmpdir = instrument_cache_dir)
      on.exit(unlink(tmp), add = TRUE)
      utils::download.file(bwmr_url, tmp, mode = "wb", quiet = TRUE)
      if (!file.copy(tmp, bwmr_file, overwrite = TRUE)) stop("Cannot save BWMR source.")
    }
    env <- new.env(parent = globalenv())
    sys.source(bwmr_file, envir = env)
    if (!is.function(env$BWMR)) stop("BWMR source does not define BWMR().")
    list(fun = env$BWMR, status = "enabled")
  }, error = function(e) {
    message("BWMR unavailable; other MR methods continue: ", conditionMessage(e))
    list(fun = NULL, status = paste("unavailable:", conditionMessage(e)))
  })
}

add_context <- function(d, pair) {
  if (is.null(d) || !NROW(d)) return(data.table())
  z <- as.data.table(copy(d))
  for (v in names(pair)) set(z, j = v, value = rep(pair[[v]][1L], nrow(z)))
  setcolorder(z, c(names(pair), setdiff(names(z), names(pair))))
  z[]
}

run_steiger <- function(h, x, y, info) {
  z <- as.data.frame(h)
  for (side in c("exposure", "outcome")) {
    trait_key <- if (side == "exposure") x else y
    meta <- info[match(trait_key, key)]
    b <- z[[paste0("beta.", side)]]
    se <- z[[paste0("se.", side)]]
    n <- z[[paste0("samplesize.", side)]]
    if (is.null(n)) n <- rep(NA_real_, nrow(z))
    if (meta$type == "binary") {
      if (!all(is.finite(c(meta$ncase, meta$ncontrol, meta$prevalence))) ||
          meta$ncase <= 0 || meta$ncontrol <= 0 || meta$prevalence <= 0 || meta$prevalence >= 1)
        stop("Missing valid ncase/ncontrol/population prevalence for ", trait_key)
      af <- z[[paste0("eaf.", side)]]
      valid <- is.finite(af) & af > 0 & af < 1
      r <- rep(NA_real_, nrow(z))
      if (any(valid)) r[valid] <- TwoSampleMR::get_r_from_lor(b[valid], af[valid], meta$ncase, meta$ncontrol, meta$prevalence)
      n <- rep(meta$ncase + meta$ncontrol, nrow(z))
    } else {
      r <- rep(NA_real_, nrow(z))
      valid <- is.finite(n) & n > 2
      if (any(valid)) r[valid] <- TwoSampleMR::get_r_from_bsen(b[valid], se[valid], n[valid])
    }
    z[[paste0("r.", side)]] <- r
    z[[paste0("samplesize.", side)]] <- n
  }
  ok <- is.finite(z$r.exposure) & is.finite(z$r.outcome) &
    is.finite(z$samplesize.exposure) & z$samplesize.exposure > 2 &
    is.finite(z$samplesize.outcome) & z$samplesize.outcome > 2
  z <- z[ok, , drop = FALSE]
  if (!nrow(z)) stop("Steiger has no complete SNPs; supply correct N/EAF and binary-trait metadata.")
  if (sum(z$r.exposure^2) >= 1 || sum(z$r.outcome^2) >= 1) stop("Steiger total R2 >= 1; check N, effect units and LD.")
  ans <- as.data.table(TwoSampleMR::directionality_test(z))
  if (nrow(ans) != 1L || !is.finite(ans$steiger_pval)) stop("Steiger did not return a finite P value.")
  ans[, `:=`(nsnp = nrow(z), nsnp_mr = nrow(h), status = "ok", message = "",
             supports_direction = correct_causal_direction & steiger_pval < 0.05)]
  ans[]
}

run_pair <- function(pair, instruments, outcome_data, info, bwmr_fun = NULL, outcome_error = NULL) {
  issues <- list()
  counts <- data.table(n_significant = NA_integer_, n_clumped = NA_integer_, n_fpass = NA_integer_,
    n_matched = NA_integer_, n_mr = NA_integer_, status = "failed", message = "")
  optional <- function(label, expr) {
    tryCatch(expr, error = function(e) {
      issues[[length(issues) + 1L]] <<- data.table(component = label, message = conditionMessage(e))
      NULL
    })
  }
  fit <- tryCatch({
    if (!is.null(instruments$error)) stop(instruments$error)
    counts[, `:=`(n_significant = instruments$n_significant, n_clumped = instruments$n_clumped,
                  n_fpass = instruments$n_fpass)]
    if (!counts$n_significant) stop("No significant exposure SNPs at the configured P threshold.")
    if (!counts$n_clumped) stop("No clumped exposure instruments found in the LD panel.")
    if (!counts$n_fpass) stop("No exposure instruments pass min_f.")
    if (!is.null(outcome_error)) stop(outcome_error)
    exposure <- instruments$data
    outcome <- outcome_data[SNP %chin% exposure$SNP]
    counts[, n_matched := nrow(outcome)]
    if (!nrow(outcome)) stop("No exposure instruments matched in the outcome GWAS.")
    h_all <- suppressMessages(TwoSampleMR::harmonise_data(
      format_mr_data(exposure, pair$exposure_key, "exposure", info),
      format_mr_data(outcome, pair$outcome_key, "outcome", info), action = harmonise_action))
    h <- h_all[h_all$mr_keep %in% TRUE, , drop = FALSE]
    counts[, n_mr := nrow(h)]
    if (!nrow(h)) stop("No SNPs remain after allele harmonisation.")
    n <- nrow(h)
    primary_method <- if (n == 1L) "mr_wald_ratio" else "mr_ivw"
    compute <- function(method_code) {
      a <- as.data.table(suppressMessages(TwoSampleMR::mr(h, method_list = method_code)))
      if (nrow(a) != 1L || !all(is.finite(c(a$b, a$se, a$pval))) ||
          a$se <= 0 || a$pval < 0 || a$pval > 1) stop("Invalid estimate from ", method_code)
      a[, `:=`(method_id = method_code, primary = method_code == primary_method)]
      a
    }
    estimates <- list(compute(primary_method))
    if (n >= 3L) for (method in c("mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
      estimates[[length(estimates) + 1L]] <- optional(method, compute(method))
    if (n >= 2L && is.function(bwmr_fun)) {
      bwmr <- optional("BWMR", {
        f <- bwmr_fun(h$beta.exposure, h$beta.outcome, h$se.exposure, h$se.outcome)
        if (!all(is.finite(c(f$beta, f$se_beta, f$P_value))) || f$se_beta <= 0 || f$P_value < 0 || f$P_value > 1)
          stop("BWMR returned an invalid estimate.")
        data.table(id.exposure = pair$exposure_key, id.outcome = pair$outcome_key,
          exposure = pair$exposure_key, outcome = pair$outcome_key, method = "BWMR", method_id = "BWMR",
          primary = FALSE, nsnp = n, b = f$beta, se = f$se_beta, pval = f$P_value)
      })
      estimates[[length(estimates) + 1L]] <- bwmr
    }
    res <- rbindlist(estimates, fill = TRUE)
    x_type <- info$type[match(pair$exposure_key, info$key)]
    y_type <- info$type[match(pair$outcome_key, info$key)]
    f_values <- (h$beta.exposure / h$se.exposure)^2
    res[, `:=`(exposure_type = x_type, outcome_type = y_type, F_min = min(f_values), F_mean = mean(f_values),
      n_F_below_10 = sum(f_values < 10), beta_low = b - qnorm(0.975) * se, beta_high = b + qnorm(0.975) * se,
      OR = NA_real_, OR_low = NA_real_, OR_high = NA_real_,
      exposure_unit = if (x_type == "binary") "1 log-odds of exposure" else "1 exposure GWAS unit",
      outcome_unit = if (y_type == "binary") "log-odds" else "outcome GWAS unit")]
    if (y_type == "binary") res[, `:=`(OR = exp(b), OR_low = exp(beta_low), OR_high = exp(beta_high))]
    heterogeneity <- if (n >= 2L) optional("heterogeneity", {
      ivw <- TwoSampleMR::mr_ivw(h$beta.exposure, h$beta.outcome, h$se.exposure, h$se.outcome)
      q <- data.table(method = "Inverse variance weighted", Q = ivw$Q, Q_df = ivw$Q_df, Q_pval = ivw$Q_pval)
      if (n >= 3L) q <- rbindlist(list(q, as.data.table(TwoSampleMR::mr_heterogeneity(h,
        method_list = "mr_egger_regression"))), fill = TRUE)
      q
    }) else NULL
    egger <- if (n >= 3L) optional("egger_intercept", TwoSampleMR::mr_pleiotropy_test(h)) else NULL
    steiger <- tryCatch(run_steiger(h, pair$exposure_key, pair$outcome_key, info),
      error = function(e) data.table(status = "skipped", message = conditionMessage(e),
        nsnp = NA_integer_, nsnp_mr = n, supports_direction = NA))
    loo <- if (run_loo && n >= 3L) optional("leave_one_out", TwoSampleMR::mr_leaveoneout(h, method = TwoSampleMR::mr_ivw)) else NULL
    presso <- if (run_presso && n >= 4L) optional("MR_PRESSO", MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure",
      data = as.data.frame(h, row.names = h$SNP), OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
      SignifThreshold = 0.05, NbDistribution = max(1000L, 100L * n), seed = seed)) else NULL
    counts[, `:=`(status = "ok", message = "")]
    list(mr = res, heterogeneity = heterogeneity, egger = egger, steiger = steiger,
      leave_one_out = loo, presso = presso,
      harmonised = if (save_harmonised) h_all else NULL)
  }, error = function(e) {
    counts[, message := conditionMessage(e)]
    list()
  })
  fit$counts <- counts
  fit$issues <- rbindlist(issues, fill = TRUE)
  for (v in c("mr", "counts", "heterogeneity", "egger", "steiger", "leave_one_out", "issues", "harmonised"))
    fit[[v]] <- add_context(fit[[v]], pair)
  fit
}

# 5. Multiple testing and a direct forward/reverse comparison ---------------
adjust_results <- function(result, pairs, scope = fdr_scope) {
  if (!scope %in% c("family", "direction")) stop("fdr_scope must be family or direction.")
  if (!nrow(result)) return(data.table(direction = character(), exposure = character(), outcome = character(),
    primary = logical(), method_id = character(), method = character(), b = numeric(), se = numeric(), pval = numeric(),
    FDR = numeric(), FDR_method = numeric()))
  p <- copy(pairs)[, fdr_group := if (scope == "family") family else "all"]
  planned <- p[, .(n_tests = .N), by = .(direction, fdr_group)]
  z <- copy(result)[, fdr_group := if (scope == "family") family else "all"]
  z <- merge(z, planned, by = c("direction", "fdr_group"), all.x = TRUE, sort = FALSE)
  if (anyDuplicated(z[primary == TRUE, pair_id])) stop("Internal error: multiple primary tests for one pair.")
  z[, FDR := NA_real_]
  # Only IVW P-values enter BH. All planned pairs remain in n, including pairs
  # without an IVW estimate. Wald ratio and the other four methods are unadjusted.
  z[method_id == "mr_ivw", FDR := p.adjust(pval, method = "BH", n = n_tests[1L]), by = .(direction, fdr_group)]
  # The report column mirrors the IVW-only FDR; other methods remain NA.
  z[, FDR_method := FDR]
  flags <- z[, {
    m <- .SD[primary == TRUE]
    s <- .SD[primary == FALSE]
    support <- s[sign(b) == sign(m$b) & m$b != 0 & pval < 0.05, method]
    list(primary_significant = m$FDR < 0.05, sensitivity_n_available = nrow(s),
      sensitivity_n_supporting = if (nrow(s)) length(support) else NA_integer_,
      sensitivity_support = if (nrow(s)) length(support) > 0 else NA,
      supporting_methods = paste(support, collapse = "; "))
  }, by = pair_id]
  z <- merge(z, flags, by = "pair_id", all.x = TRUE, sort = FALSE)
  setorder(z, pair_id, -primary, method)
  z[]
}

compare_directions <- function(result, counts, pairs) {
  ans <- unique(pairs[, .(trait1_key, trait2_key)])
  for (d in c("forward", "reverse")) {
    selected <- pairs[direction == d, .(pair_id, trait1_key, trait2_key)]
    z <- merge(selected, counts[, .(pair_id, status, message)], by = "pair_id", all.x = TRUE, sort = FALSE)
    if (nrow(result)) {
      m <- result[primary == TRUE & direction == d,
        .(pair_id, method, nsnp, b, se, beta_low, beta_high, OR, OR_low, OR_high, pval, FDR)]
      z <- merge(z, m, by = "pair_id", all.x = TRUE, sort = FALSE)
    }
    z[, pair_id := NULL]
    cols <- setdiff(names(z), c("trait1_key", "trait2_key"))
    setnames(z, cols, paste0(d, "_", cols))
    ans <- merge(ans, z, by = c("trait1_key", "trait2_key"), all.x = TRUE, sort = FALSE)
    col <- paste0(d, "_status")
    ans[is.na(get(col)), (col) := "not_requested"]
  }
  ans[]
}

# 6. One instrument-preparation pass, then one read per outcome --------------
run_batch <- function(pairs, traits, info, bwmr_fun = NULL) {
  exposures <- unique(pairs$exposure_key)
  instruments <- setNames(vector("list", length(exposures)), exposures)
  message("Preparing instruments for ", length(exposures), " exposure traits.")
  for (exp_key in exposures) instruments[[exp_key]] <- tryCatch(get_instruments(traits[match(exp_key, traits$key)]),
    error = function(e) list(error = conditionMessage(e)))
  fits <- vector("list", nrow(pairs))
  outcomes <- unique(pairs$outcome_key)
  completed <- 0L
  for (y in outcomes) {
    indexes <- which(pairs$outcome_key == y)
    xs <- unique(pairs$exposure_key[indexes])
    snps <- unique(unlist(lapply(xs, function(x) instruments[[x]]$data$SNP), use.names = FALSE))
    outcome_error <- NULL
    outcome_data <- data.table(SNP = character())
    if (length(snps)) {
      message("Read outcome once: ", y, " (", length(snps), " requested SNPs; ", length(indexes), " pairs)")
      outcome_data <- tryCatch(read_gwas(traits$file[match(y, traits$key)], snps = snps),
        error = function(e) {outcome_error <<- conditionMessage(e); data.table(SNP = character())})
    }
    for (i in indexes) {
      pair <- pairs[i]
      # Stable per-pair seed: forward-only and both give the same stochastic fit.
      pair_text <- paste(pair$direction, pair$exposure_key, pair$outcome_key, sep = "|")
      code <- utf8ToInt(pair_text)
      pair_seed <- as.integer((as.double(seed) + sum(as.double(code) * seq_along(code))) %% 2147483646 + 1)
      set.seed(pair_seed)
      fits[[i]] <- run_pair(pair, instruments[[pair$exposure_key]], outcome_data, info, bwmr_fun, outcome_error)
      completed <- completed + 1L
      message("[", completed, "/", nrow(pairs), "] ", pair$direction, ": ", pair$exposure_key, " -> ",
        pair$outcome_key, " : ", fits[[i]]$counts$status, if (fits[[i]]$counts$status != "ok")
          paste0(" (", fits[[i]]$counts$message, ")") else "")
    }
    rm(outcome_data)
  }
  collect <- function(name) rbindlist(lapply(fits, `[[`, name), fill = TRUE)
  counts <- collect("counts")
  results <- adjust_results(collect("mr"), pairs)
  list(mr = results, counts = counts, hits = results[method_id == "mr_ivw" & is.finite(FDR) & FDR < 0.05],
    comparison = compare_directions(results, counts, pairs), failures = counts[status != "ok"],
    heterogeneity = collect("heterogeneity"), egger = collect("egger"), steiger = collect("steiger"),
    leave_one_out = collect("leave_one_out"), issues = collect("issues"),
    harmonised = if (save_harmonised) collect("harmonised") else NULL,
    presso = if (run_presso) setNames(lapply(fits, `[[`, "presso"),
      paste(pairs$direction, pairs$exposure_key, pairs$outcome_key, sep = "|")) else NULL)
}

# 7. Five-method main table and pair-level sensitivity table -----------------
make_report_tables <- function(result, pairs, info, bwmr_status = "enabled") {
  base <- merge(pairs[, .(pair_id, direction, exposure, outcome, exposure_key, outcome_key)],
    result$counts[, .(pair_id, n_mr, analysis_status = status, analysis_message = message)],
    by = "pair_id", all.x = TRUE, sort = FALSE)
  base[, `:=`(outcome_type = info$type[match(outcome_key, info$key)], F_mean = NA_real_, F_min = NA_real_)]
  # Use qualified labels only when a bare ID would be ambiguous in this run.
  ambiguous <- info[, .(n_keys = uniqueN(key)), by = id][n_keys > 1L, id]
  base[exposure %in% ambiguous, exposure := exposure_key]
  base[outcome %in% ambiguous, outcome := outcome_key]
  if (nrow(result$mr)) {
    f <- result$mr[primary == TRUE, .(pair_id, F_mean, F_min)]
    at <- match(f$pair_id, base$pair_id)
    set(base, at, "F_mean", f$F_mean)
    set(base, at, "F_min", f$F_min)
  }
  methods <- data.table(
    method_id = c("mr_ivw", "BWMR", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression", "mr_wald_ratio"),
    method = c("IVW", "Bayesian weighted (BWMR)", "Weighted median", "Weighted mode", "MR-Egger", "Wald ratio"),
    min_snps = c(2L, 2L, 3L, 3L, 3L, 1L), method_order = 1:6)
  main <- base[rep(seq_len(nrow(base)), each = nrow(methods))]
  for (v in names(methods)) set(main, j = v, value = rep(methods[[v]], times = nrow(base)))
  # A single harmonised SNP yields one Wald-ratio row, including its failure
  # reason if estimation failed. Other pairs retain the five-method structure.
  main <- main[(n_mr == 1L & method_id == "mr_wald_ratio") |
    ((is.na(n_mr) | n_mr != 1L) & method_id != "mr_wald_ratio")]
  main[, `:=`(nsnp = n_mr, b = NA_real_, se = NA_real_, beta_low = NA_real_, beta_high = NA_real_,
    OR = NA_real_, OR_low = NA_real_, OR_high = NA_real_, pval = NA_real_, FDR_method = NA_real_)]
  if (nrow(result$mr)) {
    estimates <- result$mr[, .(pair_id, method_id, nsnp, b, se, beta_low, beta_high, OR, OR_low, OR_high, pval, FDR_method)]
    # pair_id identifies both trait IDs AND direction; methods cannot cross-join.
    index <- match(paste(estimates$pair_id, estimates$method_id), paste(main$pair_id, main$method_id))
    keep <- !is.na(index)
    for (v in setdiff(names(estimates), c("pair_id", "method_id")))
      set(main, i = index[keep], j = v, value = estimates[[v]][keep])
  }
  main[, `:=`(status = "method_failed", reason = "Estimate unavailable; see optional_method_issues.")]
  main[analysis_status != "ok", `:=`(status = "pair_failed", reason = analysis_message)]
  main[analysis_status == "ok" & n_mr < min_snps,
    `:=`(status = "insufficient_snps", reason = paste0("Requires >=", min_snps, " harmonised SNPs; available: ", n_mr, "."))]
  if (!identical(bwmr_status, "enabled")) main[analysis_status == "ok" & method_id == "BWMR" & n_mr >= min_snps,
    `:=`(status = "unavailable", reason = paste("BWMR", bwmr_status))]
  if (nrow(result$issues)) {
    issues <- result$issues[component %in% methods$method_id]
    index <- match(paste(issues$pair_id, issues$component), paste(main$pair_id, main$method_id))
    for (j in which(!is.na(index))) set(main, index[j], "reason", issues$message[j])
  }
  main[is.finite(b) & is.finite(se) & se > 0 & is.finite(pval), `:=`(status = "ok", reason = "")]
  main[, `:=`(effect_type = ifelse(outcome_type == "binary", "OR", "Beta"), effect_ci = NA_character_)]
  fmt <- function(x) {
    out <- sprintf("%.3f", x)
    small <- is.finite(x) & x != 0 & abs(x) < 0.001
    out[small] <- sprintf("%.2e", x[small])
    out[is.na(x)] <- NA_character_
    out
  }
  ci <- function(x, lo, hi) paste0(fmt(x), " (", fmt(lo), ", ", fmt(hi), ")")
  main[status == "ok" & effect_type == "OR", effect_ci := ci(OR, OR_low, OR_high)]
  main[status == "ok" & effect_type == "Beta", effect_ci := ci(b, beta_low, beta_high)]
  setorder(main, pair_id, method_order)
  main <- main[, .(pair_id, direction, exposure, outcome, method_id, method, nsnp, F_mean, F_min,
    effect_type, effect_ci, pval, FDR_method, status, reason)]

  sensitivity <- copy(base)[, `:=`(egger_intercept = NA_real_, egger_se = NA_real_, egger_pval = NA_real_,
    ivw_Q = NA_real_, ivw_Q_df = NA_real_, ivw_Q_pval = NA_real_, steiger_pval = NA_real_,
    correct_causal_direction = NA, supports_direction = NA, steiger_nsnp = NA_integer_,
    steiger_status = "not_run", steiger_message = analysis_message)]
  fill <- function(source, map) {
    if (!nrow(source)) return(invisible(NULL))
    if (anyDuplicated(source$pair_id)) stop("Duplicate diagnostic rows for a pair in report assembly.")
    index <- match(source$pair_id, sensitivity$pair_id)
    keep <- !is.na(index)
    for (v in names(map)) if (map[[v]] %in% names(source))
      set(sensitivity, index[keep], v, source[[map[[v]]]][keep])
    invisible(NULL)
  }
  fill(result$egger, c(egger_intercept = "egger_intercept", egger_se = "se", egger_pval = "pval"))
  if (nrow(result$heterogeneity)) fill(result$heterogeneity[method == "Inverse variance weighted"],
    c(ivw_Q = "Q", ivw_Q_df = "Q_df", ivw_Q_pval = "Q_pval"))
  fill(result$steiger, c(steiger_pval = "steiger_pval", correct_causal_direction = "correct_causal_direction",
    supports_direction = "supports_direction", steiger_nsnp = "nsnp", steiger_status = "status", steiger_message = "message"))
  # Numerical failures (for example a degenerate Egger fit) are unavailable
  # diagnostics, not Excel error cells. Original diagnostic values remain in raw outputs.
  for (v in c("egger_intercept", "egger_se", "egger_pval", "ivw_Q", "ivw_Q_df", "ivw_Q_pval", "steiger_pval"))
    sensitivity[!is.finite(get(v)), (v) := NA_real_]
  setorder(sensitivity, pair_id)
  notes <- main[status != "ok", .(direction, exposure, outcome, component = method, reason)]
  for (component_name in c("Egger intercept", "IVW heterogeneity", "Steiger")) {
    d <- switch(component_name,
      "Egger intercept" = sensitivity[!is.finite(egger_pval)],
      "IVW heterogeneity" = sensitivity[!is.finite(ivw_Q_pval)],
      "Steiger" = sensitivity[steiger_status != "ok"])
    if (!nrow(d)) next
    reason <- if (component_name == "Steiger") d$steiger_message else ifelse(
      d$analysis_status != "ok", d$analysis_message, ifelse(
        d$n_mr < if (component_name == "Egger intercept") 3L else 2L,
        "Too few harmonised SNPs for this diagnostic.", "Diagnostic unavailable; see the detailed diagnostic and optional-method sheets."))
    notes <- rbindlist(list(notes, d[, .(direction, exposure, outcome, component = component_name)][, reason := reason]), fill = TRUE)
  }
  sensitivity <- sensitivity[, .(pair_id, direction, exposure, outcome, n_mr,
    egger_intercept, egger_se, egger_pval, ivw_Q, ivw_Q_df, ivw_Q_pval, steiger_pval,
    correct_causal_direction, supports_direction, steiger_nsnp, steiger_status, steiger_message, analysis_status)]
  list(main = main, sensitivity = sensitivity, notes = unique(notes))
}

# Publication tables plus unmodified analysis sheets in one openxlsx workbook.
# report contains main, sensitivity and notes data.frames/data.tables.
write_results_workbook <- function(sheets, report, file) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("openxlsx is required.")
  wb <- openxlsx::createWorkbook(creator = "Bidirectional MR analysis")
  openxlsx::modifyBaseFont(wb, fontName = "Calibri", fontSize = 11, fontColour = "#172B3A")
  st <- list(
    title = openxlsx::createStyle(fontSize = 15, textDecoration = "bold", fontColour = "#173B50"),
    header = openxlsx::createStyle(fgFill = "#173B50", fontColour = "#FFFFFF",
      textDecoration = "bold", halign = "center", valign = "center", wrapText = TRUE),
    group = openxlsx::createStyle(fgFill = "#DCEAF0", fontColour = "#173B50",
      textDecoration = "bold", halign = "center", valign = "center"),
    body = openxlsx::createStyle(valign = "center", wrapText = TRUE),
    centered = openxlsx::createStyle(halign = "center"),
    band = openxlsx::createStyle(fgFill = "#F3F7F9"),
    boundary = openxlsx::createStyle(border = "top", borderColour = "#CBD9E1"),
    integer = openxlsx::createStyle(numFmt = "0", halign = "right"),
    decimal = openxlsx::createStyle(numFmt = "0.00", halign = "right"),
    precise = openxlsx::createStyle(numFmt = "0.0000", halign = "right"),
    probability = openxlsx::createStyle(numFmt = "0.00E+00", halign = "right"),
    footnote = openxlsx::createStyle(fontSize = 10, fontColour = "#526875",
      wrapText = TRUE, valign = "top")
  )
  apply_style <- function(sheet, style, rows, cols) {
    if (length(rows) && length(cols)) openxlsx::addStyle(wb, sheet, style,
      rows = rows, cols = cols, gridExpand = TRUE, stack = TRUE)
  }
  new_report_sheet <- function(sheet, title, width, header_row) {
    # openxlsx 4.2.8 can emit invalid XML when tabColour and fitToWidth are combined.
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE, zoom = 90)
    openxlsx::writeData(wb, sheet, title, startRow = 1, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = seq_along(width), rows = 1)
    apply_style(sheet, st$title, 1, seq_along(width))
    openxlsx::setRowHeights(wb, sheet, rows = 1, heights = 27)
    openxlsx::setColWidths(wb, sheet, cols = seq_along(width), widths = width)
    openxlsx::freezePane(wb, sheet, firstActiveRow = header_row + 1L, firstActiveCol = 3)
    openxlsx::pageSetup(wb, sheet, orientation = "landscape", paperSize = 9,
      fitToWidth = TRUE, fitToHeight = FALSE, printTitleRows = seq_len(header_row))
  }
  write_table <- function(sheet, x, header_row, pair_id = NULL) {
    openxlsx::writeData(wb, sheet, x, startRow = header_row, headerStyle = st$header,
      withFilter = nrow(x) > 0L, keepNA = TRUE, na.string = "NA")
    openxlsx::setRowHeights(wb, sheet, rows = header_row, heights = 34)
    rows <- seq_len(nrow(x)) + header_row
    apply_style(sheet, st$body, rows, seq_len(ncol(x)))
    if (length(rows)) {
      openxlsx::setRowHeights(wb, sheet, rows = rows, heights = 29)
      block <- if (is.null(pair_id)) seq_len(nrow(x)) else match(pair_id, unique(pair_id))
      apply_style(sheet, st$band, rows[block %% 2L == 0L], seq_len(ncol(x)))
      starts <- c(TRUE, diff(block) != 0L)
      apply_style(sheet, st$boundary, rows[starts], seq_len(ncol(x)))
    }
    rows
  }
  footnotes <- function(sheet, text, start_row, ncol) {
    for (i in seq_along(text)) {
      row <- start_row + i - 1L
      openxlsx::writeData(wb, sheet, text[i], startRow = row, colNames = FALSE)
      openxlsx::mergeCells(wb, sheet, cols = seq_len(ncol), rows = row)
      apply_style(sheet, st$footnote, row, seq_len(ncol))
      openxlsx::setRowHeights(wb, sheet, rows = row, heights = 29)
    }
  }
  main <- as.data.frame(report$main, stringsAsFactors = FALSE)
  sensitivity <- as.data.frame(report$sensitivity, stringsAsFactors = FALSE)
  directions <- c("forward", "reverse")
  direction_label <- c(forward = "Forward MR: trait1 to trait2", reverse = "Reverse MR: trait2 to trait1")
  # Keep complete method rows, including unavailable estimates supplied as NA.
  for (d in directions) {
    x <- main[main$direction %in% d, , drop = FALSE]
    if (!nrow(x)) next
    sheet <- paste0("main_", d)
    type <- unique(x$effect_type[!is.na(x$effect_type) & nzchar(x$effect_type)])
    mixed <- length(type) != 1L
    label <- if (mixed) "Effect (95% CI)\u1d43" else paste0(type, " (95% CI)\u1d43")
    ci <- as.character(x$effect_ci)
    if (mixed) {
      valid <- !is.na(ci) & ci != "NA" & !is.na(x$effect_type)
      ci[valid] <- paste(x$effect_type[valid], ci[valid])
    }
    table <- data.frame(x$exposure, x$outcome, x$method, x$nsnp, x$F_mean,
      ci, x$pval, x$FDR_method, stringsAsFactors = FALSE, check.names = FALSE)
    names(table) <- c("Exposure", "Outcome", "Method", "No. of SNP", "F statistics",
      label, "P-value", "P-FDR")
    new_report_sheet(sheet, paste(if (d == "forward") "Forward" else "Reverse", "MR estimates"),
      c(24, 24, 28, 12, 14, 32, 14, 14), 3L)
    rows <- write_table(sheet, table, 3L, x$pair_id)
    apply_style(sheet, st$integer, rows, 4)
    apply_style(sheet, st$decimal, rows, 5)
    apply_style(sheet, st$centered, rows, 6)
    apply_style(sheet, st$probability, rows, 7:8)
    footnotes(sheet, c(
      "No. of SNP is the number of harmonised instruments used for this exposure-outcome pair; F statistics is mean((beta_exposure / SE_exposure)^2).",
      "a. Effect sizes were estimated using IVW, Bayesian weighted, weighted median, weighted mode and MR-Egger, as applicable. Wald ratio was used when only one SNP was available.",
      "For disease outcomes, effect sizes are ORs with 95% CIs; for continuous outcomes, they are Beta estimates with 95% CIs. CIs use beta +/- 1.96 x SE, exponentiated for ORs.",
      "For a binary exposure, estimates are per 1-unit increase in its genetically predicted log odds, not a contrast between diagnosed cases and controls.",
      "Only IVW P-values were adjusted using the Benjamini-Hochberg procedure, within the configured test family and separately by direction. Planned pairs remain in the correction denominator.",
      "IVW P-FDR < 0.05 was considered statistically significant. The other methods, including Wald ratio, report unadjusted P-values; their P-FDR is not applicable (NA).",
      "NA indicates an unavailable or inapplicable estimate; see report_notes for the reason. Detailed values and all methods are retained in the raw analysis sheets."
    ), nrow(table) + 6L, 8L)
  }
  for (d in directions) {
    x <- sensitivity[sensitivity$direction %in% d, , drop = FALSE]
    if (!nrow(x)) next
    sheet <- paste0("sens_", d)
    correct <- ifelse(is.na(x$correct_causal_direction), "NA",
      ifelse(x$correct_causal_direction, "Yes", "No"))
    table <- data.frame(x$exposure, x$outcome, x$egger_intercept, x$egger_se,
      x$egger_pval, x$ivw_Q, x$ivw_Q_df, x$ivw_Q_pval, x$steiger_pval, correct,
      stringsAsFactors = FALSE, check.names = FALSE)
    names(table) <- c("Exposure", "Outcome", "Egger intercept", "Standard error", "P value",
      "Q statistic", "Q df", "P value", "P value", "Correct direction")
    new_report_sheet(sheet, paste(if (d == "forward") "Forward" else "Reverse", "MR sensitivity analyses"),
      c(24, 24, 18, 17, 14, 17, 11, 14, 14, 18), 4L)
    groups <- list("Directional pleiotropy test" = 3:5, "Heterogeneity test" = 6:8, "Steiger test" = 9:10)
    for (name in names(groups)) {
      cols <- groups[[name]]
      openxlsx::writeData(wb, sheet, name, startCol = min(cols), startRow = 3, colNames = FALSE)
      openxlsx::mergeCells(wb, sheet, cols = cols, rows = 3)
      apply_style(sheet, st$group, 3, cols)
    }
    openxlsx::setRowHeights(wb, sheet, rows = 3, heights = 25)
    rows <- write_table(sheet, table, 4L, x$pair_id)
    apply_style(sheet, st$precise, rows, 3:4)
    apply_style(sheet, st$decimal, rows, 6)
    apply_style(sheet, st$integer, rows, 7)
    apply_style(sheet, st$probability, rows, c(5, 8, 9))
    apply_style(sheet, st$centered, rows, 10)
    footnotes(sheet, c(
      "Directional pleiotropy is assessed by the MR-Egger intercept test. This requires at least three harmonised instruments.",
      "Heterogeneity columns contain IVW Cochran's Q, its degrees of freedom and P value. MR-Egger Q remains in the raw heterogeneity sheet.",
      "Correct direction is the Steiger correct_causal_direction flag (Yes/No), not an inference from the P value alone. Yes with P < 0.05 supports the direction tested.",
      "Steiger may use fewer instruments than MR because N, EAF or binary-trait metadata are missing. Its SNP count and status are retained in the raw Steiger sheet.",
      "NA indicates an unavailable or inapplicable test; see report_notes. A nonsignificant sensitivity test does not establish that the corresponding bias is absent."
    ), nrow(table) + 7L, 10L)
  }
  notes <- as.data.frame(report$notes, stringsAsFactors = FALSE)
  openxlsx::addWorksheet(wb, "report_notes", gridLines = FALSE, tabColour = "#D2A65A")
  openxlsx::writeData(wb, "report_notes", notes, headerStyle = st$header,
    keepNA = TRUE, na.string = "NA", withFilter = nrow(notes) > 0L)
  openxlsx::setColWidths(wb, "report_notes", cols = seq_len(ncol(notes)), widths = c(14, 24, 24, 24, 100))
  openxlsx::freezePane(wb, "report_notes", firstActiveRow = 2, firstActiveCol = 3)
  apply_style("report_notes", st$body, seq_len(nrow(notes)) + 1L, seq_len(ncol(notes)))
  if (nrow(notes)) openxlsx::setRowHeights(wb, "report_notes", rows = seq_len(nrow(notes)) + 1L, heights = 32)
  # Raw sheets keep their original row order, values, column names and types.
  for (name in names(sheets)) {
    sheet <- name
    if (tolower(sheet) %in% tolower(names(wb))) sheet <- paste0("raw_", sheet)
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)
    x <- sheets[[name]]
    openxlsx::writeData(wb, sheet, x, headerStyle = st$header, keepNA = TRUE,
      na.string = "NA", withFilter = nrow(x) > 0L)
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(x)), widths = 20)
    openxlsx::setRowHeights(wb, sheet, rows = 1, heights = 32)
  }
  # This writer creates tables only. Some openxlsx versions retain relationships
  # to drawing files that were never created; discard only empty placeholders.
  for (i in seq_along(wb$worksheets_rels)) {
    if (!length(wb$drawings[[i]]) && !length(wb$vml[[i]]) && !length(wb$comments[[i]]))
      wb$worksheets_rels[[i]] <- wb$worksheets_rels[[i]][
        !grepl('Type="[^"]*/relationships/(drawing|vmlDrawing)"', wb$worksheets_rels[[i]])]
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  invisible(file)
}

# 8. Execute and save compact outputs ----------------------------------------
main <- function() {
  needed <- c("data.table", "R.utils", "TwoSampleMR", "openxlsx")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages first: ", paste(missing, collapse = ", "), ". See README.")
  suppressPackageStartupMessages(library(data.table))
  if (!is.finite(p_cut) || p_cut <= 0 || p_cut >= 1 || !is.finite(clump_r2) || clump_r2 < 0 || clump_r2 > 1 ||
      !is.finite(clump_kb) || clump_kb <= 0 || !is.finite(min_f) || min_f < 0) stop("Invalid instrument-selection settings.")
  if (!fdr_scope %in% c("family", "direction")) stop("Invalid fdr_scope.")
  if (!harmonise_action %in% 1:3) stop("harmonise_action must be 1, 2 or 3.")
  if (!file.exists(plink)) stop("PLINK not found: ", plink)
  ref_files <- paste0(ld_prefix, c(".bed", ".bim", ".fam"))
  if (any(!file.exists(ref_files))) stop("Missing LD reference: ", paste(ref_files[!file.exists(ref_files)], collapse = "; "))
  if (run_presso && !requireNamespace("MRPRESSO", quietly = TRUE)) stop("run_presso=TRUE requires MRPRESSO.")
  registry <- discover_traits(data_dirs)
  a <- select_traits(trait1, registry)
  b <- select_traits(trait2, registry)
  pairs <- make_pairs(a, b)
  traits <- rbind(a, b)
  info <- load_trait_info(traits)
  base <- file.path(outdir, paste(paste(sort(unique(a$group)), collapse = "+"),
    paste(sort(unique(b$group)), collapse = "+"), sep = "__"))
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  dir.create(instrument_cache_dir, recursive = TRUE, showWarnings = FALSE)
  outdir <- tempfile(paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = base)
  if (!dir.create(outdir)) stop("Cannot create output directory: ", outdir)
  fwrite(pairs, file.path(outdir, "planned_pairs.csv"))
  fwrite(info, file.path(outdir, "steiger_trait_info_used.csv"))
  bwmr <- load_bwmr()
  set.seed(seed)
  result <- run_batch(pairs, traits, info, bwmr$fun)
  report <- make_report_tables(result, pairs, info, bwmr$status)
  csv <- list(uvmr_results = result$mr, uvmr_hits = result$hits, uvmr_bidirectional = result$comparison,
    uvmr_snp_counts = result$counts, uvmr_failures = result$failures, uvmr_steiger = result$steiger,
    uvmr_main_table = report$main, uvmr_sensitivity_table = report$sensitivity)
  for (name in names(csv)) {
    tab <- csv[[name]]
    if (!ncol(tab)) tab <- data.table(status = character(), message = character())
    fwrite(tab, file.path(outdir, paste0(name, ".csv")))
  }
  settings <- data.table(setting = c("trait1", "trait2", "direction", "fdr_scope", "p_cut", "clump_kb", "clump_r2",
    "min_f", "harmonise_action", "LD_prefix", "BWMR", "MR_PRESSO", "leave_one_out", "seed", "TwoSampleMR_version"),
    value = as.character(c(paste(trait1, collapse = ";"), paste(trait2, collapse = ";"), direction, fdr_scope,
      p_cut, clump_kb, clump_r2, min_f, harmonise_action, ld_prefix, bwmr$status, run_presso, run_loo, seed,
      as.character(utils::packageVersion("TwoSampleMR")))))
  sheets <- list(results = result$mr, primary = result$mr[primary == TRUE], hits = result$hits,
    bidirectional = result$comparison, snp_counts = result$counts, heterogeneity = result$heterogeneity,
    egger_intercept = result$egger, steiger = result$steiger, failures = result$failures,
    optional_method_issues = result$issues, planned_pairs = pairs, trait_info = info, settings = settings)
  if (run_loo) sheets$leave_one_out <- result$leave_one_out
  sheets <- lapply(sheets, function(x) if (ncol(x)) as.data.frame(x) else data.frame(status = "No rows / not applicable"))
  # Excel has a finite row limit; complete CSV/RDS outputs remain available.
  max_excel_rows <- 1048575L
  too_long <- names(sheets)[vapply(sheets, nrow, integer(1)) > max_excel_rows]
  if (length(too_long)) for (name in too_long) {
    saveRDS(sheets[[name]], file.path(outdir, paste0("uvmr_", name, ".rds")))
    sheets[[name]] <- data.frame(status = paste("Exceeds Excel row limit; see uvmr_", name, ".rds", sep = ""))
  }
  write_results_workbook(sheets, report, file.path(outdir, "uvmr_results.xlsx"))
  if (save_harmonised) saveRDS(result$harmonised, file.path(outdir, "uvmr_harmonised.rds"))
  if (run_presso) saveRDS(result$presso, file.path(outdir, "uvmr_presso.rds"))
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  message("Completed: ", sum(result$counts$status == "ok"), "/", nrow(pairs), " pairs; ", nrow(result$hits), " IVW FDR hits.")
  message("Results saved: ", outdir)
  print(report$main[, .(direction, exposure, outcome, method, nsnp, F_mean, effect_ci, pval, FDR_method, status)])
  invisible(list(outdir = outdir, results = result, report = report))
}

# Tests may set options(uvmr.skip_main=TRUE) before sourcing the functions.
if (!isTRUE(getOption("uvmr.skip_main", FALSE))) main()
