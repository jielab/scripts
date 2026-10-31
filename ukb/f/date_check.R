#!/usr/bin/env Rscript
# Read-only audit of saved UKB dates and their underlying date evidence.
# This script never rewrites Rdata. Run through phe.sh date-check.
suppressPackageStartupMessages(library(data.table))
argv <- commandArgs(trailingOnly = TRUE)
arg <- list()
if (length(argv) %% 2L) stop("Arguments must be --name value pairs")
if (length(argv)) for (i in seq(1L, length(argv), 2L)) arg[[sub("^--", "", argv[i])]] <- argv[i + 1L]
need <- function(nm) { x <- arg[[nm]]; if (is.null(x) || !nzchar(x)) stop("Missing --", nm); x }
trait <- need("Y"); indir <- need("phenotype-dir"); work <- need("output-dir")
if (!grepl("^[A-Za-z0-9_]+$", trait)) stop("Y must be a trait name from common/*.lst")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script_file)), "date_utils.R"))
asof <- if (!is.null(arg[["as-of"]])) ukb_date_parse(arg[["as-of"]]) else Sys.Date()
if (length(asof) != 1L || is.na(asof)) stop("--as-of must be YYYY-MM-DD")
dir.create(work, showWarnings = FALSE, recursive = TRUE)
rdir <- file.path(indir, "Rdata")
defs <- fread(file.path(work, "definitions.tsv"), colClasses = "character", na.strings = c("", "NA"))
getdef <- function(src, type) unique(defs[source == src & kind == type & !is.na(value), value])
initial_meta <- if (file.exists(file.path(work, "metadata.tsv")))
    fread(file.path(work, "metadata.tsv"), colClasses = "character") else data.table(key = character(), value = character())
input_available <- function(key) any(initial_meta$key == key & tolower(initial_meta$value) == "true")
# Keep malformed/special raw strings visible. Conversion is only for validated
# calendar-date columns, never for the sole preserved copy of the evidence.
date_text <- function(x) {
    if (inherits(x, "Date") || inherits(x, "POSIXt") || is.numeric(x) || is.logical(x))
        return(as.character(ukb_date_parse(x)))
    z <- trimws(as.character(x)); z[z %in% c("", "NA", "NaN")] <- NA_character_; z
}
read_tab <- function(path) {
    if (grepl("\\.rds$", path, ignore.case = TRUE)) return(as.data.table(readRDS(path)))
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        # fread's older gzip support requires R.utils; keep this audit's dependency
        # limited to data.table. Linked extracts are uncompressed staging files.
        return(as.data.table(read.delim(gzfile(path), check.names = FALSE,
            colClasses = "character", na.strings = c("", "NA"), quote = "", comment.char = "")))
    }
    fread(path, colClasses = "character", na.strings = c("", "NA"))
}
inputs <- list(); limitations <- character(); roster <- character()
register <- function(path, role, state = "used") {
    info <- file.info(path)
    inputs[[length(inputs) + 1L]] <<- list(role = role, path = normalizePath(path, winslash = "/", mustWork = FALSE),
        status = state, bytes = if (nrow(info) && !is.na(info$size)) unname(info$size) else NULL,
        modified = if (nrow(info) && !is.na(info$mtime)) format(info$mtime, tz = "UTC", usetz = TRUE) else NULL)
}
check_id <- function(x, label) {
    if (!"eid" %in% names(x)) stop(label, ": missing eid")
    x[, eid := as.character(eid)]
    if (anyNA(x$eid) || any(!nzchar(x$eid)) || anyDuplicated(x$eid)) stop(label, ": missing or duplicate eid")
    roster <<- union(roster, x$eid)
    x
}
json <- function(x) {
    esc <- function(s) {
        s <- gsub("\\", "\\\\", s, fixed = TRUE)
        s <- gsub('"', '\\"', s, fixed = TRUE)
        s <- gsub("\n", "\\n", s, fixed = TRUE); s <- gsub("\r", "\\r", s, fixed = TRUE)
        s <- gsub("\t", "\\t", s, fixed = TRUE)
        paste0('"', s, '"')
    }
    if (is.null(x)) return("null")
    if (is.list(x)) {
        y <- vapply(x, json, character(1))
        if (!is.null(names(x))) return(paste0("{", paste(paste0(esc(names(x)), ":", y), collapse = ","), "}"))
        return(paste0("[", paste(y, collapse = ","), "]"))
    }
    if (length(x) != 1L) return(paste0("[", paste(vapply(as.list(x), json, character(1)), collapse = ","), "]"))
    if (is.na(x)) return("null")
    if (is.logical(x)) return(if (x) "true" else "false")
    if (is.numeric(x)) return(if (is.finite(x)) as.character(x) else "null")
    esc(as.character(x))
}
meta_add <- function(key, value) {
    cat(key, "\t", json(value), "\n", file = file.path(work, "metadata.tsv"), append = TRUE, sep = "")
}
clean_code <- function(code, src) {
    z <- trimws(as.character(code))
    z[z %in% c("", "NA", "NaN", ".")] <- NA_character_
    if (src %in% c("icd10", "icd9", "opcs4", "death")) z <- toupper(gsub(".", "", z, fixed = TRUE))
    z
}
matches <- function(code, src, official = FALSE) {
    source <- if (src == "death") "icd10" else src
    code <- clean_code(code, src)
    defined <- getdef(source, if (official) "official_defined" else "local_defined")
    if (length(defined) && !any(tolower(defined) %in% c("1", "true"))) return(rep(NA, length(code)))
    if (official) {
        keys <- getdef(source, "official_code")
        if (!length(keys) && !length(defined)) return(rep(NA, length(code)))
        if (source == "icd10") return(!is.na(code) & substr(code, 1L, 3L) %in% keys)
        return(!is.na(code) & code %in% keys)
    }
    if (source %in% c("gp2", "gp3")) return(!is.na(code) & code %in% getdef(source, "local_code"))
    pats <- getdef(source, "local_pattern")
    if (!length(pats) && !length(defined)) return(rep(NA, length(code)))
    if (source %in% c("icd10", "icd9", "opcs4")) pats <- toupper(gsub(".", "", pats, fixed = TRUE))
    out <- rep(FALSE, length(code))
    for (pat in pats) out <- out | (!is.na(code) & grepl(paste0("^(?:", pat, ")", if (source == "srd") "$" else ""), code, perl = TRUE))
    out
}
empty_ev <- function() data.table(eid = character(), source = character(), code = character(), date_raw = character(),
    origin = character(), granularity = character(), match_local = logical(), match_official = logical())
parts <- list()
add_events <- function(x) {
    if (!nrow(x)) return(invisible(NULL))
    for (nm in names(empty_ev())) if (!nm %in% names(x)) x[, (nm) := NA_character_]
    x <- x[, names(empty_ev()), with = FALSE]
    x[, eid := as.character(eid)]
    x[, date_raw := as.character(date_raw)]
    if (anyNA(x$eid) || any(!nzchar(x$eid))) stop("Event evidence contains missing eid")
    parts[[length(parts) + 1L]] <<- x
}
cache <- list(); source_state <- list()
sources <- c("icd10", "icd9", "opcs4", "srd", "gp2", "gp3", "gp", "ref")
message("Reading saved first-occurrence dates for ", trait)
for (src in sources) {
    path <- file.path(rdir, paste0("fod.", src, ".rds"))
    if (!file.exists(path) && src %in% c("gp2", "gp3")) path <- file.path(rdir, "fod.gp.rds")
    if (!file.exists(path)) { source_state[[src]] <- list(source = src, status = "file_missing", path = path); next }
    x <- check_id(read_tab(path), path); register(path, paste0("saved_", src))
    col <- paste0("fod_", src, "_", trait)
    if (!col %in% names(x)) { source_state[[src]] <- list(source = src, status = "trait_column_missing", path = path); next }
    cc <- paste0("cnt_", src, "_", trait); uc <- paste0("cnt_undated_", src, "_", trait)
    raw <- date_text(x[[col]])
    count <- if (cc %in% names(x)) suppressWarnings(as.numeric(x[[cc]])) else as.numeric(!is.na(raw))
    undated <- if (uc %in% names(x)) suppressWarnings(as.numeric(x[[uc]])) else rep(NA_real_, nrow(x))
    cache[[src]] <- data.table(eid = x$eid, raw = raw, count = count, undated = undated)
    ix <- which(!is.na(raw) | (!is.na(count) & count > 0))
    if (length(ix)) add_events(data.table(eid = x$eid[ix], source = src, code = col,
        date_raw = raw[ix], origin = basename(path),
        granularity = if (src == "ref") "official_ref_cache" else "fod_only", match_local = if (src == "ref") NA else TRUE, match_official = NA))
    # AOD can show dates omitted by FOD. Retain this evidence even when raw
    # input is absent, but never describe these lists as complete event histories.
    ac <- paste0("aodtxt_", src, "_", trait); al <- paste0("aod_", src, "_", trait)
    if (ac %in% names(x)) {
        ix <- which(!is.na(x[[ac]]) & nzchar(x[[ac]]))
        if (length(ix)) {
            vals <- strsplit(as.character(x[[ac]][ix]), "[;|,[:space:]]+")
            add_events(data.table(eid = rep(x$eid[ix], lengths(vals)), source = src, code = ac,
                date_raw = unlist(vals, use.names = FALSE), origin = basename(path),
                granularity = "cached_aod", match_local = if (src == "ref") NA else TRUE, match_official = NA))
        }
    } else if (al %in% names(x) && is.list(x[[al]])) {
        vals <- lapply(x[[al]], function(z) as.character(ukb_date_parse(z)))
        add_events(data.table(eid = rep(x$eid, lengths(vals)), source = src, code = al,
            date_raw = unlist(vals, use.names = FALSE), origin = basename(path),
            granularity = "cached_aod", match_local = if (src == "ref") NA else TRUE, match_official = NA))
    }
    source_state[[src]] <- list(source = src, status = "saved_fod", path = path, date_column = col)
    rm(x); invisible(gc(FALSE))
}

message("Tracing code/date slots and reference provenance")
raw_available <- setNames(rep(FALSE, 7L), c("icd10", "icd9", "opcs4", "srd", "gp2", "gp3", "ref"))
for (src in c("icd10", "icd9", "opcs4", "srd")) {
    cp <- file.path(rdir, paste0(src, ".code0.rds")); dp <- file.path(rdir, paste0(src, ".date0.rds"))
    if (!file.exists(cp) || !file.exists(dp)) {
        limitations <- c(limitations, paste(src, "raw code/date RDS pair unavailable; saved AOD/FOD is used if present.")); next
    }
    c0 <- check_id(read_tab(cp), cp); d0 <- check_id(read_tab(dp), dp)
    if (!setequal(c0$eid, d0$eid)) stop(src, ": code/date participant sets differ")
    d0 <- d0[match(c0$eid, d0$eid)]
    cn <- grep("^i[0-9]+$", names(c0), value = TRUE)
    dn <- paste0("d", sub("^i", "", cn))
    if (!length(cn) || !setequal(dn, setdiff(names(d0), "eid"))) stop(src, ": mismatched labelled iN/dN slots; rebuild icd stage")
    register(cp, paste0(src, "_codes")); register(dp, paste0(src, "_dates")); raw_available[src] <- TRUE
    gran <- if (src == "srd") "self_report_approximate" else "summary_first_per_code"
    mp <- file.path(rdir, "srd.date.meta.rds")
    if (src == "srd" && file.exists(mp)) {
        m <- readRDS(mp); register(mp, "srd_date_method"); meta_add("srd_date_method", m)
        if (any(grepl("20008|interpol", unlist(m), ignore.case = TRUE))) gran <- "self_report_interpolated"
    }
    for (j in seq_along(cn)) {
        code <- clean_code(c0[[cn[j]]], src)
        lm <- matches(code, src); om <- matches(code, src, TRUE); ix <- which((lm %in% TRUE) | (om %in% TRUE))
        if (!length(ix)) next
        add_events(data.table(eid = c0$eid[ix], source = src, code = code[ix],
            date_raw = date_text(d0[[dn[j]]][ix]),
            origin = paste0(basename(cp), ":", cn[j], "/", dn[j]),
            granularity = gran, match_local = lm[ix], match_official = om[ix]))
    }
    rm(c0, d0); invisible(gc(FALSE))
}

# Add complete, filtered GP and HES records from the streaming Python reader.
lp <- file.path(work, "linked-events.tsv")
if (file.exists(lp)) {
    x <- fread(lp, colClasses = "character", na.strings = c("", "NA"))
    if (nrow(x)) {
        bool_flag <- function(v) ifelse(is.na(v) | !nzchar(v), NA, tolower(v) %in% c("true", "1"))
        x[, match_local := bool_flag(match_local)]
        x[, match_official := bool_flag(match_official)]
        raw_available[intersect(unique(x$source), c("gp2", "gp3"))] <- TRUE
        add_events(x)
    }
    rm(x)
}
for (src in c("gp2", "gp3"))
    raw_available[src] <- raw_available[src] || input_available(paste0("raw_", src, "_available"))
# A diagnosis table without episode dates is evidence of disease, but cannot
# establish that a dated summary cache is wrong merely because its date is absent.
if (input_available("raw_hospital_episode_dates_available"))
    for (src in c("icd10", "icd9")) raw_available[src] <- raw_available[src] || input_available(paste0("raw_", src, "_available"))

# Raw reference preserves coding819 values and the companion source coding2171.
refs <- getdef("ref", "ref_field")
if (!length(refs)) refs <- getdef("ref", "field")
refs <- sub("^p?", "p", refs)
raw_ref <- NULL; raw_src <- NULL
for (path in c(file.path(rdir, "fod.date0.rds"), file.path(indir, "rap/fod.date.tab.gz"))) {
    if (file.exists(path)) { raw_ref <- check_id(read_tab(path), path); register(path, "raw_reference"); break }
}
for (path in c(file.path(rdir, "fod.ref.source0.rds"), file.path(indir, "rap/fod.source.tab.gz"))) {
    if (file.exists(path)) {
        candidate <- check_id(read_tab(path), path)
        if (ncol(candidate) > 1L) { raw_src <- candidate; register(path, "reference_source"); break }
        register(path, "reference_source", "no_source_columns")
    }
}
ref_rows <- list()
if (!is.null(raw_ref) && length(refs)) {
    fids <- sub("_.*$", "", names(raw_ref))
    ref_cols <- names(raw_ref)[fids %in% refs]
    raw_available["ref"] <- length(ref_cols) > 0L
    missing_ref_fields <- setdiff(refs, fids)
    if (length(missing_ref_fields)) limitations <- c(limitations,
        paste("Some expected reference fields were not exported:", paste(missing_ref_fields, collapse = ", "), "; grouped reference dates use available fields only."))
    for (cc in ref_cols) {
        r <- date_text(raw_ref[[cc]]); ix <- which(!is.na(r))
        if (!length(ix)) next
        srcfield <- paste0("p", as.integer(sub("^p", "", sub("_.*$", "", cc))) + 1L)
        sc <- rep(NA_character_, length(ix))
        scols <- if (!is.null(raw_src)) names(raw_src)[sub("_.*$", "", names(raw_src)) == srcfield] else character()
        exact_source <- sub("^p[0-9]+", srcfield, cc)
        if (exact_source %in% scols) scols <- exact_source
        if (length(scols) == 1L) sc <- as.character(raw_src[[scols]][match(raw_ref$eid[ix], raw_src$eid)])
        else if (length(scols) > 1L) limitations <- c(limitations, paste("Ambiguous source instance for", cc, "; source attribution left unavailable."))
        ref_rows[[length(ref_rows) + 1L]] <- data.table(eid = raw_ref$eid[ix], raw = r[ix], field = cc, source_code = sc)
        add_events(data.table(eid = raw_ref$eid[ix], source = "ref", code = cc, date_raw = r[ix],
            origin = paste0("UKB ", cc, "; ", srcfield, "=", ifelse(is.na(sc), "unavailable", sc)),
            granularity = "official_ref", match_local = NA, match_official = NA))
    }
}
if (is.null(raw_src)) limitations <- c(limitations, "Reference source fields are unavailable; date coincidences alone do not prove source attribution. Export the odd companion fields and rebuild fod_ref.")

# Disease-coded death records are separate from all-cause date_death. Both
# underlying and contributory causes are relevant to UKB first occurrence.
death_path <- file.path(indir, "rap/death.tab.gz")
death_context <- NULL
if (file.exists(death_path)) {
    dt <- check_id(read_tab(death_path), death_path); register(death_path, "death_registrations")
    dcols <- grep("^p40000(?:_i[0-9]+)?$", names(dt), value = TRUE, perl = TRUE)
    if (length(dcols)) {
        death_context <- data.table(eid = dt$eid, date_death = do.call(ukb_date_min, lapply(dt[, ..dcols], ukb_date_parse)))
        for (dc in dcols) {
            suffix <- sub("^p40000", "", dc)
            cc <- grep(paste0("^p4000[12]", suffix, "(?:_a[0-9]+)?$"), names(dt), value = TRUE, perl = TRUE)
            for (cx in cc) {
                vals <- strsplit(ifelse(is.na(dt[[cx]]), "", dt[[cx]]), "[|;[:space:]]+")
                code <- clean_code(unlist(vals, use.names = FALSE), "death")
                idx <- rep(seq_len(nrow(dt)), lengths(vals))
                lm <- matches(code, "death"); om <- matches(code, "death", TRUE); keep <- which((lm %in% TRUE) | (om %in% TRUE))
                if (length(keep)) add_events(data.table(eid = dt$eid[idx[keep]], source = "death", code = code[keep],
                    date_raw = date_text(dt[[dc]][idx[keep]]), origin = paste0("death.tab.gz:", cx, "/", dc),
                    granularity = "registration", match_local = lm[keep], match_official = om[keep]))
            }
        }
    }
    rm(dt)
} else {
    limitations <- c(limitations, "Raw death registration fields unavailable; disease-specific death dates cannot be reconstructed from all-cause death dates.")
    dp <- file.path(rdir, "death.rds")
    if (file.exists(dp)) {
        x <- check_id(read_tab(dp), dp); register(dp, "saved_death_context")
        if ("date_death" %in% names(x)) death_context <- data.table(eid = x$eid, date_death = ukb_date_valid(x$date_death))
    }
}

ev <- if (length(parts)) rbindlist(parts, use.names = TRUE, fill = TRUE) else empty_ev()
rm(parts); invisible(gc(FALSE))
if (!nrow(ev)) stop("No disease evidence found for ", trait, ". Check common/*.lst and the phenotype input paths.")
ev[, date_raw := fifelse(is.na(date_raw), "", date_raw)]
ev[, date := ukb_date_parse(date_raw)]
ev[, date_status := fifelse(!nzchar(date_raw), "unknown", fifelse(is.na(date), "unparseable", "valid"))]
special_names <- c("unknown_date", "before_birth_placeholder", "birth_date_placeholder", "birth_year_placeholder", "future_placeholder", "future_placeholder")
for (j in seq_along(ukb_special_dates)) ev[date == ukb_special_dates[j] & !is.na(date), date_status := special_names[j]]
ev[date_status != "valid", date := as.Date(NA)]
ev[!is.na(date) & date > asof, `:=`(date_status = "after_as_of", date = as.Date(NA))]
ev[, match_local := as.logical(match_local)]
ev[, match_official := as.logical(match_official)]

# No drop_na(ref): include the union of local, reference and official-map-only
# cases, even if every matched record has an unknown or protected date.
patients <- data.table(eid = sort(unique(ev$eid)))
context <- file.path(rdir, "phe.rds")
patients[, `:=`(date_birth = as.Date(NA), date_assess = as.Date(NA), date_death = as.Date(NA))]
if (file.exists(context)) {
    x <- check_id(read_tab(context), context); register(context, "baseline_context")
    for (pair in list(c("date_birth", "birth_date", "date_birth"), c("date_assess", "date_attend", "date_assess"))) {
        col <- intersect(pair[-1L], names(x))
        if (length(col)) set(patients, j = pair[1L], value = ukb_date_valid(x[[col[1L]]][match(patients$eid, x$eid)]))
    }
    rm(x)
}
if (!is.null(death_context)) patients[, date_death := death_context$date_death[match(eid, death_context$eid)]]
birth <- patients$date_birth[match(ev$eid, patients$eid)]
badbirth <- which(!is.na(ev$date) & !is.na(birth) & ev$date < birth)
if (length(badbirth)) ev[badbirth, `:=`(date_status = "before_birth", date = as.Date(NA))]

# Recomputed raw minima are kept separate from saved FOD to expose stale
# caches, code-list changes, date alignment errors and extra HES episodes.
raw_gran <- !ev$granularity %in% c("fod_only", "cached_aod", "official_ref", "official_ref_cache")
recompute <- function(src, official = FALSE) {
    keep <- ev$source == src & raw_gran & !is.na(ev$date) & if (official) (ev$match_official %in% TRUE) else (ev$match_local %in% TRUE)
    if (official && src == "death") keep <- keep & (is.na(ev$date) | ev$date >= as.Date("2006-01-01"))
    if (official && src == "srd") keep <- keep & (is.na(ev$date) | ev$date >= as.Date("1930-01-01"))
    ix <- which(keep)
    if (!length(ix)) return(rep(as.Date(NA), nrow(patients)))
    xx <- ev[ix, .(first = min(date)), by = eid]
    xx$first[match(patients$eid, xx$eid)]
}
for (src in c("icd10", "icd9", "opcs4", "srd", "gp2", "gp3", "gp", "death", "ref")) {
    observed <- recompute(src)
    set(patients, j = paste0("computed_", src), value = observed)
    chosen <- observed
    if (!is.null(cache[[src]])) {
        cc <- cache[[src]]; ix <- match(patients$eid, cc$eid)
        chosen <- ukb_date_valid(cc$raw[ix])
        bad <- !is.na(chosen) & (chosen > asof | (!is.na(patients$date_birth) & chosen < patients$date_birth))
        chosen[bad] <- as.Date(NA)
        set(patients, j = paste0("cnt_", src), value = cc$count[ix])
        set(patients, j = paste0("cnt_undated_", src), value = cc$undated[ix])
    }
    set(patients, j = src, value = chosen)
    if (src %in% c("icd10", "icd9", "srd", "gp2", "gp3", "death"))
        set(patients, j = paste0("official_", src), value = recompute(src, TRUE))
}
patients[, computed_gp := ukb_date_min(computed_gp2, computed_gp3)]
if (is.null(cache$gp)) patients[, gp := ukb_date_min(gp2, gp3)]
patients[, `:=`(me = ukb_date_min(icd10, icd9, opcs4, srd, gp),
               local_min = ukb_date_min(icd10, icd9, srd, gp, death),
               reconstructed_ref = ukb_date_min(official_icd10, official_icd9, official_srd, official_gp2, official_gp3, official_death))]
patients[, reconstructed_source_code := NA_character_]
source_dates <- list(death = patients$official_death,
    gp = ukb_date_min(patients$official_gp2, patients$official_gp3),
    hospital = ukb_date_min(patients$official_icd10, patients$official_icd9), srd = patients$official_srd)
n_sources <- rowSums(do.call(cbind, lapply(source_dates, function(x) !is.na(x))))
for (j in seq_along(source_dates)) {
    take <- which(is.na(patients$reconstructed_source_code) & !is.na(source_dates[[j]]) & source_dates[[j]] == patients$reconstructed_ref)
    patients$reconstructed_source_code[take] <- as.character(c(20L, 30L, 40L, 50L)[j] + as.integer(n_sources[take] > 1L))
}

# Collapse reference groups after preserving every component field as an event.
# saved_ref is retained separately: raw reference fields can be newer than the
# saved aggregate, or can expose an invalid date that the old aggregate accepted.
patients[, saved_ref := if (!is.null(cache$ref)) ref else as.Date(NA)]
patients[, `:=`(ref_raw = "", ref_source = "", ref_status = "missing", ref_expected_source = "", ref_source_match = "source_unavailable")]
if (isTRUE(raw_available["ref"])) {
    patients[, `:=`(ref = as.Date(NA), computed_ref = as.Date(NA), ref_status = "missing_in_raw_reference")]
}
if (length(ref_rows)) {
    rr <- rbindlist(ref_rows); rr[, clean := ukb_date_valid(raw)]
    rb <- patients$date_birth[match(rr$eid, patients$eid)]
    invalid <- which(!is.na(rr$clean) & (rr$clean > asof | (!is.na(rb) & rr$clean < rb)))
    if (length(invalid)) rr[invalid, clean := as.Date(NA)]
    collapsed <- rr[, {
        dated <- which(!is.na(clean))
        md <- if (length(dated)) min(clean[dated]) else as.Date(NA)
        min_sources <- source_code[which(!is.na(clean) & clean == md)]
        min_sources <- min_sources[!is.na(min_sources) & nzchar(min_sources)]
        list(ref = md, ref_raw = paste(unique(paste0(field, "=", raw)), collapse = ";"),
             ref_source = paste(sort(unique(min_sources)), collapse = ";"),
             ref_status = if (!length(dated)) "special_or_invalid_date" else
                 if (any(raw %in% as.character(ukb_special_dates))) "valid_with_special_component" else "valid")
    }, by = eid]
    ix <- match(collapsed$eid, patients$eid)
    for (cc in setdiff(names(collapsed), "eid")) set(patients, i = ix, j = cc, value = collapsed[[cc]])
    set(patients, i = ix, j = "computed_ref", value = collapsed$ref)
} else if (!isTRUE(raw_available["ref"]) && !is.null(cache$ref)) {
    cr <- cache$ref$raw[match(patients$eid, cache$ref$eid)]
    patients[, ref_raw := fifelse(is.na(cr), "", cr)]
    patients[!is.na(ref), ref_status := "valid_cache_only"]
    patients[is.na(ref) & nzchar(ref_raw), ref_status := "special_or_invalid_date"]
    if ("cnt_ref" %in% names(patients)) patients[is.na(ref) & !nzchar(ref_raw) & !is.na(cnt_ref) & cnt_ref > 0, ref_status := "undated_cache_evidence"]
}
src_code <- c(`20` = "death", `21` = "death", `30` = "gp", `31` = "gp", `40` = "hospital", `41` = "hospital", `50` = "srd", `51` = "srd")
for (i in which(nzchar(patients$ref_source))) {
    cs <- strsplit(patients$ref_source[i], ";", fixed = TRUE)[[1L]]
    cats <- unique(unname(src_code[cs])); cats <- cats[!is.na(cats)]
    patients$ref_expected_source[i] <- paste(cats, collapse = ";")
    # Multiple grouped reference fields can legitimately have different sources.
    actual <- patients$reconstructed_source_code[i]
    patients$ref_source_match[i] <- if (!length(cats)) "unrecognised_source_code" else if (is.na(actual)) "no_reconstructed_source" else
        if (length(cs) > 1L) "grouped_multiple_sources" else if (cs == actual) "same_code" else if (src_code[cs] == src_code[actual]) "same_first_source_other_sources_differ" else "different_first_source"
}
cmp_cols <- c("icd10", "icd9", "opcs4", "srd", "gp", "death")
patients[, ref_match_sources := ""]
for (s in cmp_cols) {
    ix <- which(!is.na(patients$ref) & !is.na(patients[[s]]) & patients[[s]] == patients$ref)
    patients$ref_match_sources[ix] <- paste0(patients$ref_match_sources[ix], ifelse(nzchar(patients$ref_match_sources[ix]), ";", ""), s)
}
patients[, ref_explanation := fifelse(ref_status == "special_or_invalid_date", "reference_special_date",
    fifelse(is.na(ref) & !is.na(me), "local_only",
    fifelse(!is.na(ref) & is.na(me), "reference_only_or_only_official_mapping",
    fifelse(!is.na(ref) & !is.na(me) & ref == me, "same_local_minimum",
    fifelse(!is.na(ref) & !is.na(reconstructed_ref) & ref == reconstructed_ref, "matched_by_official_mapping_candidate",
    fifelse(!is.na(ref) & !is.na(me) & ref < me, "reference_earlier_unresolved",
    fifelse(!is.na(ref) & !is.na(me) & ref > me, "local_earlier_unresolved", "undated_evidence")))))))]
patients[, cache_mismatch_sources := ""]
ss <- intersect(names(cache), c("icd10", "icd9", "opcs4", "srd", "gp2", "gp3", "ref"))
ss <- ss[vapply(ss, function(s) isTRUE(raw_available[s]), logical(1))]
for (s in ss) {
    a <- patients[[if (s == "ref") "saved_ref" else s]]; b <- patients[[paste0("computed_", s)]]
    hit <- which(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b))
    patients$cache_mismatch_sources[hit] <- paste0(patients$cache_mismatch_sources[hit], ifelse(nzchar(patients$cache_mismatch_sources[hit]), ";", ""), s)
}
ec <- ev[, .(n_evidence_rows = .N, n_unique_valid_dates = uniqueN(date[!is.na(date)]),
             n_undated_rows = sum(is.na(date)), n_record_rows = sum(!granularity %in% c("fod_only", "cached_aod", "official_ref", "official_ref_cache"))), by = eid]
patients <- merge(patients, ec, by = "eid", all.x = TRUE, sort = FALSE)
patients[, me_minus_ref_days := as.numeric(me - ref)]
patients[, reconstructed_minus_ref_days := as.numeric(reconstructed_ref - ref)]
setorder(patients, eid)
setcolorder(ev, c("eid", "source", "code", "date", "date_raw", "date_status", "origin", "granularity", "match_local", "match_official"))
setorder(ev, eid, source, date, code, origin, na.last = TRUE)
fwrite(patients, file.path(work, "patients.tsv.gz"), sep = "\t", na = "", quote = "auto", compress = "gzip")
fwrite(ev, file.path(work, "events.tsv.gz"), sep = "\t", na = "", quote = "auto", compress = "gzip")
meta_add("as_of", as.character(asof)); meta_add("cohort_ids_in_available_inputs", length(roster))
meta_add("case_union_n", nrow(patients)); meta_add("r_inputs", inputs); meta_add("saved_sources", unname(source_state))
meta_add("date_rules", list(me = "min(saved icd10, icd9, opcs4, srd, gp); reference excluded",
    local_min = "min(local icd10, icd9, srd, gp, disease-coded death)",
    reconstructed_ref = "min(dates selected with official code maps from available records); candidate, not guaranteed full reproduction",
    first_source_tie_order = "death > primary care > hospital > self-report",
    special_dates = "Retained as raw evidence, not plotted as calendar dates or silently converted to disease-free controls",
    cached_dates = "Source columns show saved FOD when available; computed_* columns independently use raw evidence and current definitions",
    reference_date = "ref uses validated raw UKB reference fields when any requested reference field is available; otherwise saved reference fallback. saved_ref preserves the cleaned saved aggregate and computed_ref the raw reference minimum."))
meta_add("audit_limitations", unique(c(limitations,
    "UKB first-occurrence fields and local records may come from different releases; file timestamps do not establish matching clinical data cutoffs.",
    "Official GP/SRD maps bundled here are the published January 2022 maps. Date equality does not prove complete reconstruction or true disease onset.",
    "Hospital summary AOD contains first date per diagnosis code, not all admissions. Only linked HES episodes and raw GP rows expose repeated event records.",
    "Missing disease evidence in a source is not proof of disease absence or availability of follow-up in that source.",
    "Cached FOD/AOD and raw record evidence can describe the same event; evidence-row counts are not independent episode counts.",
    "Self-reported dates and birth day 15 are approximate; coding819 birth-related placeholders remain flagged rather than assigned exact days.")))
message("Extracted ", nrow(patients), " people and ", nrow(ev), " evidence rows (all retained; no sampling).")
