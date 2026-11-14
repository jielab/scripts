# CKM phenotype builder for jielab/scripts, repository contract 971a61b.
# Reads UKB outputs; never writes or modifies UKB/Rdata. One baseline per person.
# Units: HbA1c mmol/mol; glucose, TG, TC, HDL mmol/L; creatinine umol/L;
# urine albumin mg/L; urine creatinine umol/L; UACR mg/g; eGFR mL/min/1.73m2.

ckm_ph_default <- function(x, fallback) if (is.null(x) || !length(x)) fallback else x
ckm_ph_num <- function(x) {
    z <- suppressWarnings(as.numeric(as.character(x)))
    z[!is.finite(z)] <- NA_real_; z
}
ckm_ph_date <- function(x) {
    ckm_shared()$.phenotype_date(x, special_to_na = TRUE)
}
ckm_ph_flag <- function(x) {
    v <- tolower(trimws(as.character(x)))
    out <- rep(NA, length(v)); out[v %in% c("1", "true", "yes")] <- TRUE
    out[v %in% c("0", "false", "no")] <- FALSE; out
}
ckm_ph_any <- function(...) {
    xs <- list(...)
    if (length(xs) == 1L && (is.data.frame(xs[[1L]]) || is.matrix(xs[[1L]]))) xs <- as.list(as.data.frame(xs[[1L]]))
    m <- do.call(cbind, xs)
    if (is.null(dim(m))) m <- matrix(m, ncol = 1L)
    positive <- rowSums(m == TRUE, na.rm = TRUE) > 0
    ifelse(positive, TRUE, ifelse(rowSums(is.na(m)) > 0, NA, FALSE))
}
ckm_ph_min_date <- function(...) {
    xs <- list(...)
    if (length(xs) == 1L && is.data.frame(xs[[1L]])) xs <- as.list(xs[[1L]])
    xs <- lapply(xs, function(x) as.numeric(ckm_ph_date(x)))
    z <- do.call(pmin, c(xs, list(na.rm = TRUE)))
    z[!is.finite(z)] <- NA_real_; as.Date(z, origin = "1970-01-01")
}
ckm_ph_read <- function(path) {
    if (!file.exists(path)) stop("Missing CKM input: ", path)
    x <- if (grepl("[.]rds$", path, ignore.case = TRUE)) readRDS(path) else
        data.table::fread(path, colClasses = c(eid = "character"), showProgress = FALSE)
    if (!is.data.frame(x) || !"eid" %in% names(x)) stop("Input requires a participant table with eid: ", path)
    x <- data.table::as.data.table(x); x[, eid := as.character(eid)]
    if (anyNA(x$eid) || any(!nzchar(x$eid)) || anyDuplicated(x$eid))
        stop("eid must be unique and nonmissing: ", path)
    x
}
ckm_ph_overlay <- function(x, y, columns = setdiff(names(y), "eid"), clear_absent = FALSE) {
    j <- match(x$eid, y$eid)
    for (nm in columns) {
        if (!nm %in% names(y)) next
        value <- y[[nm]][j]
        if (!isTRUE(clear_absent) && nm %in% names(x)) value[is.na(j)] <- x[[nm]][is.na(j)]
        data.table::set(x, j = nm, value = value)
    }
    x
}

ckm_ph_medications <- function(d) {
    n <- nrow(d)
    cols <- grep("^([.]baseline_med_|[.]ckm_raw_med_|(p)?(6153|6177)_i0(_a[0-9]+)?$)", names(d), value = TRUE)
    # phe0.rds preserves the original -7 none response before negative-code QC.
    if (!length(cols)) cols <- grep("^drug[.]big3", names(d), value = TRUE)
    cols <- cols[!grepl("[_.]i[1-9]", cols)]
    out <- matrix(NA, n, 3L, dimnames = list(NULL, c("statin", "bp_tx", "dm_tx")))
    if (!length(cols)) return(as.data.frame(out))
    patterns <- do.call(paste, c(lapply(d[, ..cols], as.character), sep = " "))
    keys <- unique(patterns)
    lookup <- t(vapply(keys, function(s) {
        tokens <- suppressWarnings(as.integer(regmatches(s, gregexpr("-?[0-9]+", s))[[1L]]))
        tokens <- tokens[!is.na(tokens)]
        conflict <- (-7L %in% tokens) && any(tokens %in% 1:5)
        if (conflict) return(rep(NA, 3L))
        known <- any(tokens %in% c(-7L, 1:5)) && !any(tokens %in% c(-1L, -3L))
        vapply(1:3, function(code) if (code %in% tokens) TRUE else if (known) FALSE else NA, logical(1))
    }, logical(3)))
    out[,] <- lookup[match(patterns, keys), , drop = FALSE]
    as.data.frame(out)
}

ckm_ph_diagnosis <- function(d, traits, sources) {
    n <- nrow(d)
    candidates <- unique(unlist(lapply(sources, function(s) paste0("fod_", s, "_", traits))))
    cols <- intersect(candidates, names(d))
    if (!length(cols)) return(list(date = rep(as.Date(NA), n),
        baseline = rep(NA, n), unknown = rep(TRUE, n), columns = character()))
    dates <- lapply(cols, function(nm) {
        z <- ckm_ph_date(d[[nm]])
        z[!is.na(z) & !is.na(d$birth_date) & z < d$birth_date] <- as.Date(NA)
        z[!is.na(z) & !is.na(d$death_date) & z > d$death_date] <- as.Date(NA)
        z
    })
    first <- do.call(ckm_ph_min_date, dates)
    known_prior <- !is.na(first) & !is.na(d$baseline_date) & first <= d$baseline_date
    unknown <- rep(FALSE, n)
    for (nm in cols) {
        suffix <- sub("^fod_", "", nm)
        count_name <- paste0("cnt_", suffix)
        undated_name <- paste0("cnt_undated_", suffix)
        raw <- d[[nm]]
        raw_present <- !is.na(raw) & nzchar(trimws(as.character(raw)))
        valid_date <- dates[[match(nm, cols)]]
        unknown <- unknown | (raw_present & is.na(valid_date))
        if (count_name %in% names(d)) unknown <- unknown |
            ((!is.na(ckm_ph_num(d[[count_name]])) & ckm_ph_num(d[[count_name]]) > 0) & is.na(valid_date))
        if (undated_name %in% names(d)) unknown <- unknown |
            (!is.na(ckm_ph_num(d[[undated_name]])) & ckm_ph_num(d[[undated_name]]) > 0)
    }
    baseline <- ifelse(known_prior, TRUE, ifelse(unknown | is.na(d$baseline_date), NA, FALSE))
    list(date = first, baseline = baseline, unknown = unknown, columns = cols)
}

ckm_ph_censor <- function(d, cfg) {
    p <- ckm_ph_default(cfg$phenotype, list())
    path <- ckm_ph_default(cfg$paths$censor, "")
    shared_end <- ckm_ph_date(ckm_shared()$date_follow_end)
    if (length(shared_end) != 1L || is.na(shared_end))
        stop("Shared 0f/phenotype.R::date_follow_end must be one valid date")
    end <- rep(shared_end, nrow(d))
    source <- rep("0f/phenotype.R::date_follow_end", nrow(d))
    if (nzchar(path)) {
        z <- if (grepl("[.]rds$", path, ignore.case = TRUE)) as.data.frame(readRDS(path)) else
            as.data.frame(data.table::fread(path, showProgress = FALSE))
        if (!"censor_date" %in% names(z)) stop("Censor file requires censor_date")
        keys <- intersect(c("eid", "center", "region"), names(z))
        if (length(keys) != 1L) stop("Censor table requires exactly one key: eid, center, or region")
        key <- keys[1L]
        if (!key %in% names(d) || anyNA(z[[key]]) || anyDuplicated(as.character(z[[key]])))
            stop("Invalid or absent censor matching key: ", key)
        j <- match(as.character(d[[key]]), as.character(z[[key]]))
        # Coverage may shorten the shared endpoint, never extend it. An absent
        # coverage match/date remains unobserved, rather than falling back to
        # the shared endpoint or a date_lost value.
        covered_end <- ckm_ph_date(z$censor_date[j])
        end <- pmin(end, covered_end, na.rm = FALSE)
        source[] <- paste0("0f/phenotype.R::date_follow_end+coverage_by_", key)
        source[is.na(end)] <- "missing_coverage_match_or_end"
        if ("coverage_start" %in% names(z)) {
            start <- ckm_ph_date(z$coverage_start[j])
            bad <- is.na(start) | is.na(d$baseline_date) | start > d$baseline_date
            end[bad] <- as.Date(NA); source[bad] <- "coverage_starts_after_baseline_or_unknown"
        }
    } else if (isTRUE(p$trust_input_censor)) {
        if (!"censor_date" %in% names(d))
            stop("trust_input_censor requires a censor_date column in the phenotype input")
        end <- pmin(end, ckm_ph_date(d$censor_date), na.rm = FALSE)
        source[] <- "0f/phenotype.R::date_follow_end+explicit_input_censor"
        source[is.na(end)] <- "missing_explicit_input_censor"
    }
    list(end = end, source = source, shared_end = shared_end)
}

ckm_build_phenotype <- function(cfg) {
    p <- ckm_ph_default(cfg$phenotype, list())
    paths <- ckm_ph_default(cfg$paths, list())
    root <- ckm_ph_default(paths$ukb_dir, Sys.getenv("DIR_UKB", ""))
    phe_path <- ckm_ph_default(paths$phe, file.path(root, "Rdata", "all.rds"))
    d <- ckm_ph_read(phe_path)
    # AOD text/list columns are large and not needed for a baseline phenotype.
    drop <- grep("^(aod|lod)_", names(d), value = TRUE)
    if (length(drop)) d[, (drop) := NULL]
    n <- nrow(d)
    source_files <- phe_path
    map <- list(baseline_date = c("baseline_date", "date_attend"),
        birth_date = "birth_date", death_date = c("death_date", "date_death"),
        date_lost = "date_lost", age = "age", sex = "sex", bmi = "bmi", waist = "waist",
        sbp = "sbp", dbp = "dbp", hba1c = c("hba1c", "bb_HBA1C"),
        glucose = c("glucose", "bb_GLU"), tg = c("tg", "bb_TG"),
        hdl = c("hdl", "bb_HDL"), total_chol = c("total_chol", "bb_TC"),
        creatinine = c("creatinine", "bb_CRE"), ualb = c("ualb", "bb_UALB"),
        ucre = c("ucre", "bb_UCRE"), uacr = "uacr", egfr = "egfr",
        fasting_hours = c("fasting_hours", "p74_i0"))
    if (!is.null(p$columns)) for (nm in names(p$columns)) map[[nm]] <- unlist(p$columns[[nm]], use.names = FALSE)
    for (nm in c("baseline_date", "birth_date", "death_date", "date_lost")) {
        hit <- intersect(map[[nm]], names(d))
        data.table::set(d, j = nm, value = if (length(hit)) ckm_ph_date(d[[hit[1L]]]) else rep(as.Date(NA), n))
    }
    if (all(is.na(d$baseline_date))) stop("No valid baseline date (date_attend) in phenotype input")
    # Only raw, non-imputed assays may determine clinical threshold criteria.
    biom_path <- ckm_ph_default(paths$biom, file.path(root, "Rdata", "bbc.rds"))
    assay_names <- unique(unlist(map[c("hba1c", "glucose", "tg", "hdl", "total_chol", "creatinine", "ualb", "ucre", "uacr", "egfr")]))
    assay_source <- "explicit_standardized_baseline_input"
    if (nzchar(biom_path) && file.exists(biom_path)) {
        raw <- ckm_ph_read(biom_path)
        # Clear same-named clinical fields in the merged, potentially imputed data.
        # Canonical aliases can otherwise shadow the raw bb_* field after overlay.
        for (nm in intersect(assay_names, names(d))) data.table::set(d, j = nm, value = rep(NA_real_, n))
        d <- ckm_ph_overlay(d, raw, intersect(assay_names, names(raw)), clear_absent = TRUE)
        source_files <- c(source_files, biom_path); assay_source <- if (grepl("[.]raw[.]rds$", biom_path)) "raw_unimputed_bbc" else "user_selected_clinical_assays_preprocessing_unverified"
        rm(raw)
    } else if (!isTRUE(p$allow_standardized_input)) {
        stop("Missing raw clinical assays: ", biom_path,
            ". Run ukb/phe.sh phe --step biom, or explicitly configure standardized baseline input.")
    }
    # Baseline medication responses before global negative-code cleaning.
    base_path <- ckm_ph_default(paths$baseline, file.path(root, "Rdata", "phe0.rds"))
    if (nzchar(base_path) && file.exists(base_path) && normalizePath(base_path) != normalizePath(phe_path)) {
        base <- ckm_ph_read(base_path)
        mc <- grep("^([.]baseline_med_|drug[.]big3|(p)?(6153|6177)_i0(_a[0-9]+)?$)", names(base), value = TRUE)
        mc <- mc[!grepl("[_.]i[1-9]", mc)]
        if (length(mc)) d <- ckm_ph_overlay(d, base, mc)
        source_files <- c(source_files, base_path); rm(base)
    }
    # Optional explicit baseline component file: same canonical units as above.
    comp_path <- ckm_ph_default(paths$components, "")
    comp_names <- character()
    if (nzchar(comp_path)) {
        comp <- ckm_ph_read(comp_path)
        comp_names <- names(comp)
        if ("component_date" %in% names(comp)) {
            jj <- match(comp$eid, d$eid); stamp <- ckm_ph_date(comp$component_date)
            bad <- !is.na(jj) & (is.na(stamp) | is.na(d$baseline_date[jj]) | stamp > d$baseline_date[jj])
            if (any(bad)) stop("Component measurements must be observed on/before baseline")
        }
        prohibited <- intersect(c("baseline_date", "birth_date", "death_date", "date_lost", "date_attend"), names(comp))
        if (length(prohibited)) stop("components is for baseline risk measurements, not date overrides: ", paste(prohibited, collapse = ", "))
        d <- ckm_ph_overlay(d, comp); source_files <- c(source_files, comp_path); rm(comp)
    }
    get_numeric <- function(nm) {
        hits <- intersect(map[[nm]], names(d))
        if (!length(hits)) return(rep(NA_real_, n))
        # Explicit canonical component columns take precedence when present;
        # otherwise coalesce raw aliases only when the earlier alias is missing.
        z <- ckm_ph_num(d[[hits[1L]]])
        if (length(hits) > 1L) for (hh in hits[-1L]) {
            zz <- ckm_ph_num(d[[hh]]); z[is.na(z)] <- zz[is.na(z)]
        }
        z[z < 0] <- NA_real_; z
    }
    for (nm in setdiff(names(map), c("baseline_date", "birth_date", "death_date", "date_lost", "sex")))
        data.table::set(d, j = nm, value = get_numeric(nm))
    sx <- intersect(map$sex, names(d))
    s <- if (length(sx)) tolower(trimws(as.character(d[[sx[1L]]]))) else rep(NA_character_, n)
    data.table::set(d, j = "sex", value = factor(ifelse(s %in% c("f", "female", "0"), "F",
        ifelse(s %in% c("m", "male", "1"), "M", NA_character_)), levels = c("F", "M")))
    # Explicit baseline readings have priority over previously derived summaries.
    for (nm in c("sbp", "dbp")) {
        cols <- grep(paste0("^", nm, "_.*_i0([_.]|$)"), names(d), value = TRUE)
        if (length(cols) && !nm %in% comp_names) {
            mat <- as.matrix(d[, ..cols]); storage.mode(mat) <- "double"; mat[mat <= 0] <- NA_real_
            z <- rowMeans(mat, na.rm = TRUE); z[!is.finite(z)] <- NA_real_
            data.table::set(d, j = nm, value = z)
        }
    }
    med <- ckm_ph_medications(d)
    for (nm in names(med)) {
        # Canonical flags are accepted only from a declared standardized input or
        # component override, never from derived drug.htn/drug.lipid ever-use.
        if (nm %in% names(d) && (isTRUE(p$allow_standardized_input) || nm %in% comp_names)) data.table::set(d, j = nm, value = ckm_ph_flag(d[[nm]])) else
            data.table::set(d, j = nm, value = med[[nm]])
    }
    if ("smoking" %in% names(d)) d[, smoking := ckm_ph_flag(smoking)] else {
        ss <- if ("smoke_status" %in% names(d)) ckm_ph_num(d$smoke_status) else rep(NA_real_, n)
        d[, smoking := ifelse(ss %in% c(0, 1), FALSE, ifelse(ss == 2, TRUE, NA))]
    }
    # If needed, read raw p6153/p6177 i0 responses directly; this preserves -7.
    if (all(is.na(d$bp_tx)) || all(is.na(d$statin))) {
        vip <- ckm_ph_default(paths$raw_vip, file.path(root, "rap", "vip.tab.gz"))
        if (nzchar(vip) && file.exists(vip)) {
            hdr <- names(data.table::fread(vip, nrows = 0L, showProgress = FALSE))
            mc <- grep("^(p)?(6153|6177)(_i0(_a[0-9]+)?|[.]0[.][0-9]+|-0[.][0-9]+)$", hdr, value = TRUE)
            if ("eid" %in% hdr && length(mc)) {
                z <- data.table::fread(vip, select = c("eid", mc), colClasses = c(eid = "character"), showProgress = FALSE)
                if (anyDuplicated(z$eid) || anyNA(z$eid)) stop("Duplicate or missing eid in raw medication input")
                data.table::setnames(z, mc, paste0(".ckm_raw_med_", seq_along(mc)))
                d <- ckm_ph_overlay(d, z); med <- ckm_ph_medications(d)
                for (nm in names(med)) { v <- d[[nm]]; fill <- is.na(v); v[fill] <- med[[nm]][fill]; data.table::set(d, j = nm, value = v) }
                source_files <- c(source_files, vip)
            }
        }
    }
    # CKD-EPI 2021 race-free creatinine equation; no floor/cap before staging.
    scr <- d$creatinine / 88.4; k <- ifelse(d$sex == "F", .7, .9)
    alpha <- ifelse(d$sex == "F", -.241, -.302)
    estimated <- 142 * pmin(scr/k, 1)^alpha * pmax(scr/k, 1)^(-1.2) *
        .9938^d$age * ifelse(d$sex == "F", 1.012, 1)
    estimated[!is.finite(estimated) | scr <= 0] <- NA_real_
    fill <- which(is.na(d$egfr))
    if (length(fill)) data.table::set(d, i = fill, j = "egfr", value = estimated[fill])
    ratio <- d$ualb * 8840 / d$ucre
    ratio[!is.finite(ratio) | d$ucre <= 0 | d$ualb < 0] <- NA_real_
    fill <- which(is.na(d$uacr))
    if (length(fill)) data.table::set(d, i = fill, j = "uacr", value = ratio[fill])
    d[, hba1c_pct := .09148 * hba1c + 2.152]
    d[, tg_mgdl := tg * 88.57]
    d[, nonhdl := total_chol - hdl]
    d[!is.finite(nonhdl) | nonhdl < 0, nonhdl := NA_real_]
    d[, uacr_category := as.character(cut(uacr, c(-Inf, 30, 300, Inf), right = FALSE, labels = c("A1", "A2", "A3")))]
    d[, egfr_category := as.character(cut(egfr, c(-Inf, 15, 30, 45, 60, 90, Inf), right = FALSE,
        labels = c("G5", "G4", "G3b", "G3a", "G2", "G1")))]
    # Clinical CVD domains and records known at recruitment.
    domains <- list(chd = c("cvd_cad", "mi", "ihd"), hf = "cvd_hfail",
        stroke = c("stroke", "cvd_stroke", "cvd_stroke_i", "cvd_stroke_h", "cvd_stroke_ih", "cvd_stroke_sh"),
        pad = "cvd_pad", af = "cvd_afib")
    if (!is.null(p$cvd_domains)) domains <- lapply(p$cvd_domains, unlist, use.names = FALSE)
    history_sources <- ckm_ph_default(p$history_sources, c("icd10", "icd9", "srd", "gp", "ref"))
    event_sources <- ckm_ph_default(p$event_sources, c("icd10", "icd9"))
    history <- lapply(domains, function(tr) ckm_ph_diagnosis(d, tr, history_sources))
    follow <- lapply(domains, function(tr) ckm_ph_diagnosis(d, tr, event_sources))
    missing_domains <- names(follow)[!vapply(follow, function(x) length(x$columns) > 0L, logical(1))]
    # With the standard UKB map, all three stroke classes must be observable
    # unless an explicitly aggregated stroke date is available. An ischemic-only
    # column must not silently stand in for all stroke.
    if (is.null(p$cvd_domains)) {
        has_trait <- function(tr) any(unlist(lapply(event_sources, function(s) paste0("fod_",s,"_",tr))) %in% names(d))
        if (!has_trait(c("stroke","cvd_stroke"))) {
            stroke_missing <- c(if (!has_trait("cvd_stroke_i")) "stroke_ischemic",
                if (!has_trait(c("cvd_stroke_ih","cvd_stroke_h"))) "stroke_intracerebral",
                if (!has_trait(c("cvd_stroke_sh","cvd_stroke_h"))) "stroke_subarachnoid")
            missing_domains <- unique(c(missing_domains,stroke_missing))
        }
        if (!has_trait(c("cvd_cad","ihd"))) missing_domains <- unique(c(missing_domains,"chd_beyond_mi"))
    }
    standardized_cvd <- isTRUE(p$allow_standardized_input) && "cvd_date_input" %in% names(d)
    if (length(missing_domains) && !isTRUE(p$allow_partial_cvd) && !standardized_cvd)
        stop("Missing incident CVD domain(s): ", paste(missing_domains, collapse = ", "),
            ". Check phenotype.cvd_domains and UKB fod columns; do not silently change the composite.")
    for (nm in names(domains)) {
        data.table::set(d, j = paste0("date_", nm, "_history"), value = history[[nm]]$date)
        data.table::set(d, j = paste0("date_", nm, "_incident_source"), value = follow[[nm]]$date)
    }
    prev <- do.call(ckm_ph_any, lapply(history, `[[`, "baseline"))
    d[, prevalent_cvd := prev]
    d[, cvd_history_date := do.call(ckm_ph_min_date, lapply(history, `[[`, "date"))]
    d[, cvd_date := do.call(ckm_ph_min_date, lapply(follow, `[[`, "date"))]
    # A standardized supplied CVD flag/date is useful for a curated phenotype and
    # demo; do not consume stale ckm.stage or prior pipeline t2e columns.
    if (standardized_cvd) {
        dd <- ckm_ph_date(d$cvd_date_input)
        d[, cvd_date := dd]
        d[, cvd_history_date := dd]
        d[, prevalent_cvd := ifelse(is.na(baseline_date), NA, !is.na(dd) & dd <= baseline_date)]
        if ("cvd_present_input" %in% names(d)) {
            present <- ckm_ph_flag(d$cvd_present_input)
            d[is.na(dd) & present %in% TRUE, prevalent_cvd := NA]
        }
    }
    dm_dx <- ckm_ph_diagnosis(d, c("t2dm", "t1dm", "diabetes"), history_sources)
    ht_dx <- ckm_ph_diagnosis(d, c("cvd_htn", "htn"), history_sources)
    # If no diagnosis-date columns exist, negative lab/medication evidence still
    # cannot establish the absence of known treated diabetes or hypertension.
    dm_lab <- d$hba1c >= 48
    fasting <- !is.na(d$fasting_hours) & d$fasting_hours >= 8
    dm_glucose <- ifelse(fasting, d$glucose >= 7, FALSE)
    diabetes <- ckm_ph_any(dm_dx$baseline, d$dm_tx, dm_lab, dm_glucose)
    if ("diabetes" %in% names(d) && (isTRUE(p$allow_standardized_input) || "diabetes" %in% comp_names)) diabetes <- ckm_ph_flag(d$diabetes)
    d[, diabetes := diabetes]
    hypertension <- ckm_ph_any(ht_dx$baseline, d$bp_tx, d$sbp >= 130, d$dbp >= 80)
    if ("hypertension" %in% names(d) && (isTRUE(p$allow_standardized_input) || "hypertension" %in% comp_names)) hypertension <- ckm_ph_flag(d$hypertension)
    d[, hypertension := hypertension]
    # Baseline adiposity and glycaemia. Random glucose is not fasting glucose.
    asian <- if ("ethnic.c" %in% names(d)) as.character(d$ethnic.c) %in% c("Asian", "Chinese") else rep(FALSE, n)
    bmi_cut <- ifelse(asian, 23, 25)
    waist_cut <- ifelse(d$sex == "M", ifelse(asian, 90, 102), ifelse(asian, 80, 88))
    central <- d$waist >= waist_cut
    adiposity <- ckm_ph_any(d$bmi >= bmi_cut, central)
    prediabetes <- ckm_ph_any(d$hba1c >= 39 & d$hba1c < 48,
        ifelse(fasting, d$glucose >= 5.6 & d$glucose < 7, FALSE))
    definition <- tolower(ckm_ph_default(p$definition, "aha2026"))
    if (!definition %in% c("aha2023", "aha2026")) stop("definition must be aha2023 or aha2026")
    tg_cut <- if (definition == "aha2023") 135 else 150
    # AHA2026 Table 4 uses 150 mg/dL; retain one transparent definition. An
    # explicitly requested nonfasting sensitivity uses 175 when fasting <8h.
    if (definition == "aha2026" && isTRUE(p$nonfasting_tg175)) tg_cut <- ifelse(fasting, 150, 175)
    hyper_tg <- d$tg_mgdl >= tg_cut
    low_hdl <- d$hdl * 38.67 < ifelse(d$sex == "M", 40, 50)
    mets_glucose <- ifelse(fasting, d$glucose >= 5.6, NA)
    mets_glucose <- ckm_ph_any(mets_glucose, d$diabetes)
    mets_modified <- FALSE
    if (isTRUE(p$mets_hba1c_proxy)) { mets_glucose <- ckm_ph_any(mets_glucose, d$hba1c >= 39); mets_modified <- TRUE }
    mets_components <- cbind(central, d$tg_mgdl >= 150, low_hdl,
        ckm_ph_any(d$sbp >= 130, d$dbp >= 80, d$bp_tx), mets_glucose)
    mets_min <- rowSums(mets_components == TRUE, na.rm = TRUE)
    mets_max <- mets_min + rowSums(is.na(mets_components))
    mets <- ifelse(mets_min >= 3, TRUE, ifelse(mets_max < 3, FALSE, NA))
    # CKD is an operational single-visit lab phenotype. Chronicity is unproven.
    ckd_any <- ckm_ph_any(d$egfr < 60, d$uacr >= 30)
    ckd_vhigh <- ckm_ph_any(d$egfr < 30,
        d$egfr >= 30 & d$egfr < 45 & d$uacr >= 30,
        d$egfr >= 45 & d$egfr < 60 & d$uacr >= 300)
    ckd_vhigh[is.finite(d$egfr) & d$egfr >= 60] <- FALSE
    kidney_failure <- d$egfr < 15
    if ("kidney_failure" %in% names(d) && (isTRUE(p$allow_standardized_input) || "kidney_failure" %in% comp_names)) kidney_failure <- ckm_ph_any(kidney_failure, ckm_ph_flag(d$kidney_failure))
    d[, `:=`(adiposity = adiposity, prediabetes = prediabetes, hypertriglyceridemia = hyper_tg,
        metabolic_syndrome = mets, mets_count_min = mets_min, mets_count_max = mets_max,
        ckd_lab = ckd_any, ckd_very_high = ckd_vhigh, kidney_failure = kidney_failure)]
    d[, criterion_stage1 := ckm_ph_any(adiposity, prediabetes)]
    d[, criterion_stage2 := ckm_ph_any(diabetes, hypertension, hypertriglyceridemia, metabolic_syndrome, ckd_lab)]
    # Optional subclinical evidence is only usable at/before recruitment.
    subclinical_parts <- list(subclinical_ascvd = rep(NA, n), pre_hf = rep(NA, n))
    for (flag in c("subclinical_ascvd", "pre_hf")) if (flag %in% names(d)) {
        stamp <- paste0(flag, "_date")
        if (!stamp %in% names(d)) stop(flag, " requires ", stamp, " to prevent future imaging leakage")
        v <- ckm_ph_flag(d[[flag]]); dt <- ckm_ph_date(d[[stamp]])
        v[is.na(dt) | is.na(d$baseline_date) | dt > d$baseline_date] <- NA
        subclinical_parts[[flag]] <- v
        data.table::set(d, j = flag, value = v)
    }
    if ("cac" %in% names(d)) {
        if (!"cac_date" %in% names(d)) stop("cac requires cac_date")
        dt <- ckm_ph_date(d$cac_date); v <- ckm_ph_num(d$cac) >= if (definition == "aha2026") 100 else .Machine$double.eps
        v[is.na(dt) | is.na(d$baseline_date) | dt > d$baseline_date] <- NA
        a <- subclinical_parts$subclinical_ascvd
        a[is.na(a)] <- v[is.na(a)]
        a[v %in% TRUE] <- TRUE
        subclinical_parts$subclinical_ascvd <- a
    }
    subclinical <- do.call(ckm_ph_any, subclinical_parts)
    d[, subclinical_cvd := subclinical]
    if ("risk_index_date" %in% names(d)) d[, risk_index_date := ckm_ph_date(risk_index_date)]
    if ("prevent_10y" %in% names(d)) d[, prevent_10y := ckm_ph_num(prevent_10y)]
    if ("score2_10y" %in% names(d)) d[, score2_10y := ckm_ph_num(score2_10y)]
    stage_method <- toupper(ckm_ph_default(p$stage3_method, "PREVENT"))
    if (!stage_method %in% c("PREVENT", "SCORE2")) stop("stage3_method must be PREVENT or SCORE2")
    risk_p <- ckm_prevent_risk(as.data.frame(d))
    risk_s <- ckm_score2_risk(as.data.frame(d))
    for (nm in names(risk_p)) data.table::set(d, j = nm, value = risk_p[[nm]])
    for (nm in names(risk_s)) data.table::set(d, j = nm, value = risk_s[[nm]])
    d[, risk_high_prevent := prevent_10y >= .20]
    d[, risk_high_score2 := score2_10y >= ifelse(age < 50, .025, ifelse(age < 70, .05, NA_real_))]
    # Externally computed SCORE2-DM/CKD models must provide their own validated
    # high-risk flag; the original SCORE2 cutoffs must not be silently reused.
    if ("score2_high_risk" %in% names(d)) {
        q <- ckm_ph_flag(d$score2_high_risk)
        if (any(!is.na(q))) {
            if (!all(c("risk_provenance", "risk_index_date") %in% names(d)))
                stop("score2_high_risk requires risk_provenance and risk_index_date")
            ok <- !is.na(d$risk_provenance) & nzchar(trimws(d$risk_provenance)) &
                !is.na(d$risk_index_date) & d$risk_index_date == d$baseline_date
            if (any(!is.na(q) & !(ok %in% TRUE))) stop("score2_high_risk has invalid baseline provenance")
        }
        d[!is.na(q), risk_high_score2 := q[!is.na(q)]]
    }
    d[, criterion_stage3_prevent := ckm_ph_any(ckd_very_high, subclinical_cvd, risk_high_prevent)]
    d[, criterion_stage3_score2 := ckm_ph_any(ckd_very_high, subclinical_cvd, risk_high_score2)]
    assign_stage <- function(c3) {
        st <- rep(0L, n)
        st[d$criterion_stage1 %in% TRUE] <- 1L
        st[d$criterion_stage2 %in% TRUE] <- 2L
        st[c3 %in% TRUE] <- 3L
        st[d$prevalent_cvd %in% TRUE] <- 4L
        no_info <- is.na(d$criterion_stage1) & is.na(d$criterion_stage2) & is.na(c3) & !(d$prevalent_cvd %in% TRUE)
        st[no_info] <- NA_integer_; st
    }
    d[, stage_prevent := assign_stage(criterion_stage3_prevent)]
    d[, stage_score2 := assign_stage(criterion_stage3_score2)]
    d[, stage := if (stage_method == "PREVENT") stage_prevent else stage_score2]
    risk_high <- if (stage_method == "PREVENT") d$risk_high_prevent else d$risk_high_score2
    core_c3 <- ckm_ph_any(d$ckd_very_high, risk_high)
    core_complete <- !is.na(d$criterion_stage1) & !is.na(d$criterion_stage2) & !is.na(core_c3) & !is.na(d$prevalent_cvd)
    core_complete[d$stage %in% c(3L, 4L)] <- !is.na(d$prevalent_cvd[d$stage %in% c(3L, 4L)])
    d[, stage_core_complete := core_complete]
    d[, stage_lower := stage]
    upper <- ifelse(is.na(d$stage), 0L, d$stage)
    upper[is.na(d$criterion_stage1)] <- pmax(upper[is.na(d$criterion_stage1)], 1L)
    upper[is.na(d$criterion_stage2)] <- pmax(upper[is.na(d$criterion_stage2)], 2L)
    missing_c3 <- is.na(if (stage_method == "PREVENT") d$criterion_stage3_prevent else d$criterion_stage3_score2)
    upper[missing_c3] <- pmax(upper[missing_c3], 3L)
    upper[is.na(d$prevalent_cvd)] <- 4L
    d[, stage_upper := as.integer(upper)]
    d[, stage_definite := !is.na(stage) & stage == stage_upper]
    d[, stage4_substage := ifelse(stage == 4L, ifelse(kidney_failure %in% TRUE, "4b",
        ifelse(kidney_failure %in% FALSE, "4a", "4_unknown_kidney_failure")), NA_character_)]
    risk_any <- ckm_ph_any(d$criterion_stage1, d$criterion_stage2, core_c3)
    d[, clinical_cvd_without_current_risk := prevalent_cvd %in% TRUE & risk_any %in% FALSE]
    d[, stage_aha := stage]
    # Clinical CVD without measured concurrent CKM risks is retained in the
    # operational Yang group, and left unassigned in the stricter AHA sensitivity.
    d[prevalent_cvd %in% TRUE & !(risk_any %in% TRUE), stage_aha := NA_integer_]
    censor <- ckm_ph_censor(d, cfg)
    d[, administrative_end := censor$end]
    d[, censor_source := censor$source]
    d[, censor_date := ckm_ph_min_date(administrative_end, date_lost)]
    valid <- !is.na(d$baseline_date) & !is.na(d$administrative_end) &
        !is.na(d$censor_date) & d$censor_date > d$baseline_date &
        (is.na(d$birth_date) | d$baseline_date > d$birth_date) &
        (is.na(d$death_date) | d$death_date > d$baseline_date)
    eligible <- valid & d$prevalent_cvd %in% FALSE & !is.na(d$stage)
    # Reuse the repository's t2e convention after applying the CKM-specific
    # baseline/unknown-evidence/coverage exclusions. Pass only date columns, not
    # the wide phenotype, and only rows with valid administrative end dates.
    shared <- ckm_shared()
    date_columns <- c("baseline_date", "birth_date", "date_lost", "death_date", "cvd_date")
    incident_time <- rep(NA_real_, n)
    incident_event <- rep(NA_integer_, n)
    ii <- which(eligible)
    if (length(ii)) {
        td <- as.data.frame(d[ii, ..date_columns])
        outcome <- shared$t2e(dat = td, domain = NA_character_, Y_date = "cvd_date",
            birth_date = "birth_date", date_attend = "baseline_date", date_lost = "date_lost",
            date_death = "death_date", date_end = d$administrative_end[ii],
            prefix = "ckm", time_unit = "year")
        incident_time[ii] <- outcome$ckm.t2e
        incident_event[ii] <- as.integer(outcome$ckm.Yt2e)
    }
    eligible <- eligible & is.finite(incident_time) & incident_time > 0 & !is.na(incident_event)
    # Shared t2e counts a CVD event on the death date. Death is a competing
    # event only when the shared CVD event indicator is zero.
    competing <- eligible & incident_event %in% 0L & !is.na(d$death_date) &
        d$death_date > d$baseline_date & d$death_date <= d$censor_date
    status <- ifelse(eligible, ifelse(incident_event == 1L, 1L,
        ifelse(competing, 2L, 0L)), NA_integer_)
    d[, status := as.integer(status)]
    d[, event := as.integer(status == 1L)]
    d[, time := ifelse(eligible, incident_time, NA_real_)]
    d[, incident_eligible := eligible]
    # Mortality uses the same helper on all valid baselines, including stage 4.
    # Prior CVD must not remove participants from the Yang mortality analysis.
    mortality_time <- rep(NA_real_, n)
    mortality_event <- rep(NA_integer_, n)
    ii <- which(valid)
    if (length(ii)) {
        td <- as.data.frame(d[ii, ..date_columns])
        outcome <- shared$t2e(dat = td, domain = NA_character_, Y_date = "death_date",
            birth_date = "birth_date", date_attend = "baseline_date", date_lost = "date_lost",
            date_death = "death_date", date_end = d$administrative_end[ii],
            prefix = "death", time_unit = "year")
        mortality_time[ii] <- outcome$death.t2e
        mortality_event[ii] <- as.integer(outcome$death.Yt2e)
    }
    d[, death_event := mortality_event]
    d[, event_death := death_event]
    d[, time_death := mortality_time]
    d[, followup_reason := ifelse(!valid, "invalid_or_uncovered_followup",
        ifelse(prevalent_cvd %in% TRUE, "prevalent_cvd_Yang",
        ifelse(is.na(prevalent_cvd), "uncertain_baseline_cvd_timing",
        ifelse(is.na(stage), "insufficient_stage_information",
        ifelse(!eligible, "shared_t2e_invalid", "ok")))))]
    if (!"alcohol_g_day" %in% names(d) && all(c("alcohol_days", "alcohol_units", "alcohol_status") %in% names(d))) {
        a <- ckm_ph_num(d$alcohol_days) * ckm_ph_num(d$alcohol_units) * 8 / 7
        ss <- ckm_ph_num(d$alcohol_status)
        a[ss %in% c(0, 1)] <- 0
        a[!is.finite(a) | a < 0 | is.na(ss)] <- NA_real_
        d[, alcohol_g_day := a]
    }
    # Expose source-specific endpoint dates for inspection, retain the phenotype
    # covariates, and drop bulky cached record counts after timing checks.
    drop <- grep("^(fod_|cnt_|[.]baseline_med_|[.]ckm_raw_med_)", names(d), value = TRUE)
    if (length(drop)) d[, (drop) := NULL]
    attr(d, "ckm_metadata") <- list(schema = "ckm-baseline-v2", source_commit = "971a61b44b84500ab46e2f5592478505a0eae50e",
        date_follow_end = as.character(censor$shared_end),
        followup_end_source = "0f/phenotype.R::date_follow_end", shared_helper = ckm_shared_info(),
        source_files = unique(source_files), definition = definition, stage3_method = stage_method,
        clinical_assay_source = assay_source, stage2_TG_mgdl = if (length(tg_cut) == 1L) tg_cut else "150 fasting / 175 nonfasting",
        statin_source = if ("statin" %in% comp_names || isTRUE(p$allow_standardized_input)) "explicit_baseline_flag" else "p6153_p6177_cholesterol_lowering_medication_proxy",
        mets_hba1c_proxy = mets_modified, ckd_chronicity = "single_visit_lab_proxy_not_confirmed_chronic_CKD",
        stage4 = "UKB operational: any recorded CHD, HF, stroke, PAD, AF; stage_aha requires concurrent CKM risk",
        stage4_tie_rule = "CVD event wins same-day CVD/death tie; same-day recruitment CVD is prevalent",
        history_sources = history_sources, event_sources = event_sources,
        cvd_domains = domains, missing_incident_domains = missing_domains,
        field_columns = map, risk_model_scope = "PREVENT-CVD ASCVD+HF, not the broad stage-4 endpoint",
        stage_uncertainty = "stage is highest observed stage; stage_upper accounts for unmeasured higher-stage criteria; missing imaging is not normal",
        score2_scope = "Adapted UK sensitivity: low-risk European SCORE2, age40-69, no diabetes, no lab CKD; no SCORE2-DM/CKD age substitution",
        source_date_columns = lapply(history, `[[`, "columns"),
        temporal_scope = "one baseline assessment; no inferred stage0-3 longitudinal transitions")
    d[]
}
