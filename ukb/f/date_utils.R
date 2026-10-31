# Date helpers for UKB phenotype construction. Raw inputs are retained separately.
ukb_special_dates <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02",
                              "1903-03-03", "1909-09-09", "2037-07-07"))

ukb_date_parse <- function(x) {
    if (inherits(x, "Date")) {
        y <- as.numeric(x)
    } else if (inherits(x, "POSIXt")) {
        return(as.Date(x))
    } else if (is.factor(x) || is.character(x)) {
        x <- trimws(as.character(x))
        y <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
        slash <- is.na(y) & !is.na(x) & grepl("^[0-9]{2}/[0-9]{2}/[0-9]{4}$", x)
        if (any(slash)) y[slash] <- suppressWarnings(as.Date(x[slash], format = "%d/%m/%Y"))
        return(y)
    } else {
        y <- suppressWarnings(as.numeric(x))
    }
    y[!is.finite(y)] <- NA_real_
    as.Date(y, origin = "1970-01-01")
}

ukb_date_valid <- function(x) {
    x <- ukb_date_parse(x)
    x[x %in% ukb_special_dates] <- as.Date(NA)
    x
}

ukb_date_min <- function(...) {
    xs <- list(...)
    if (length(xs) == 1L && is.data.frame(xs[[1L]])) xs <- as.list(xs[[1L]])
    xs <- Filter(Negate(is.null), xs)
    if (!length(xs)) return(as.Date(NA))
    y <- do.call(pmin, c(lapply(xs, function(x) as.numeric(ukb_date_valid(x))), list(na.rm = TRUE)))
    y[!is.finite(y)] <- NA_real_
    as.Date(y, origin = "1970-01-01")
}

# A UKB array slot is defined within an assessment instance. Flattening all
# instances and padding only at row end moves repeat-visit codes onto i0 dates.
ukb_slot_key <- function(nms) {
    inst <- ifelse(grepl("_i[0-9]+", nms), sub(".*_i([0-9]+).*", "\\1", nms), "0")
    arr <- ifelse(grepl("_a[0-9]+", nms), sub(".*_a([0-9]+).*", "\\1", nms), NA_character_)
    data.frame(instance = as.integer(inst), array = as.integer(arr), stringsAsFactors = FALSE)
}

ukb_expand_codes <- function(code, date, source) {
    code <- as.data.frame(code, stringsAsFactors = FALSE)
    date <- as.data.frame(date, stringsAsFactors = FALSE)
    if (!all(c("eid") %in% names(code)) || !"eid" %in% names(date)) stop(source, ": missing eid")
    code$eid <- as.character(code$eid); date$eid <- as.character(date$eid)
    if (anyNA(code$eid) || anyNA(date$eid) || anyDuplicated(code$eid) || anyDuplicated(date$eid))
        stop(source, ": duplicate/missing eid in code/date export")
    if (!setequal(code$eid, date$eid)) stop(source, ": code/date exports have different participant sets")
    code <- code[match(date$eid, code$eid), , drop = FALSE]
    cn <- setdiff(names(code), "eid"); dn <- setdiff(names(date), "eid")
    if (!length(cn) || !length(dn)) stop(source, ": code/date export has no event columns")
    ck <- ukb_slot_key(cn); dk <- ukb_slot_key(dn)
    dk$array[is.na(dk$array)] <- 0L
    if (anyDuplicated(paste(dk$instance, dk$array))) stop(source, ": duplicate date instance/array slot")
    out <- data.frame(eid = date$eid, stringsAsFactors = FALSE)
    for (j in seq_along(dn)) out[[paste0("i", j)]] <- rep(NA_character_, nrow(date))
    for (ins in unique(c(ck$instance, dk$instance))) {
        ci <- which(ck$instance == ins); di <- which(dk$instance == ins)
        if (!length(ci)) stop(source, ": no code field for date instance ", ins)
        packed <- ci[is.na(ck$array[ci])]
        if (length(packed)) {
            if (length(ci) != 1L) stop(source, ": mixed packed/expanded code columns in instance ", ins)
            values <- as.character(code[[cn[ci]]])
            values[is.na(values)] <- ""
            pieces <- strsplit(values, "|", fixed = TRUE)
            if (!length(di)) {
                if (any(nzchar(values))) stop(source, ": coded instance ", ins, " has no date columns")
                next
            }
            maxslot <- max(c(lengths(pieces), 0L))
            unmatched <- setdiff(seq_len(maxslot), dk$array[di] + 1L)
            for (k in unmatched) {
                extra <- vapply(pieces, function(z) if (length(z) >= k) z[k] else NA_character_, character(1))
                present <- !is.na(extra) & nzchar(trimws(extra)) & !toupper(trimws(extra)) %in% c("NA", "NAN")
                if (any(present)) stop(source, ": a coded array slot has no corresponding date column in instance ", ins)
            }
            for (j in di) {
                k <- dk$array[j] + 1L
                out[[paste0("i", j)]] <- vapply(pieces, function(z) if (length(z) >= k) z[k] else NA_character_, character(1))
            }
        } else {
            if (anyDuplicated(ck$array[ci])) stop(source, ": duplicate expanded code slot in instance ", ins)
            missing_date <- ci[!ck$array[ci] %in% dk$array[di]]
            if (length(missing_date) && any(vapply(code[cn[missing_date]], function(z) any(!is.na(z) & nzchar(as.character(z))), logical(1))))
                stop(source, ": a coded array slot has no corresponding date column")
            for (j in di) {
                hit <- ci[ck$array[ci] == dk$array[j]]
                if (length(hit)) out[[paste0("i", j)]] <- as.character(code[[cn[hit]]])
            }
        }
    }
    out
}

# Death remains an observed date. Administrative censoring belongs in t2e().
# Primary CVD death uses the code and date from the same registration instance.
ukb_death_dates <- function(dat) {
    n <- nrow(dat)
    dates <- grep("^p40000(?:_i[0-9]+)?$", names(dat), value = TRUE, perl = TRUE)
    if (!length(dates)) stop("Death export has no p40000 date fields")
    dd <- lapply(dat[dates], ukb_date_valid)
    cvd <- lapply(seq_along(dates), function(j) {
        cause <- sub("^p40000", "p40001", dates[j])
        z <- dd[[j]]
        codes <- if (cause %in% names(dat)) as.character(dat[[cause]]) else rep(NA_character_, n)
        is_cvd <- !is.na(codes) & grepl("(^|[[:space:]|])I[0-9]", toupper(codes), perl = TRUE)
        z[!is_cvd] <- as.Date(NA)
        z
    })
    data.frame(date_death = do.call(ukb_date_min, dd), date_cvd_death = do.call(ukb_date_min, cvd))
}

# The caller supplies the dataset's year limit; it must not drift with today's date.
ukb_srd_date <- function(x, birth_year, interpolated = FALSE, year_cap) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    birth_year <- suppressWarnings(as.numeric(as.character(birth_year)))
    if (interpolated) {
        yr <- floor(x)
        ok <- is.finite(x) & yr >= 1800 & yr <= year_cap & (is.na(birth_year) | yr >= birth_year)
        out <- rep(as.Date(NA), length(x))
        if (any(ok)) {
            jan <- as.Date(paste0(yr[ok], "-01-01"))
            nextjan <- as.Date(paste0(yr[ok] + 1L, "-01-01"))
            out[ok] <- jan + floor((x[ok] - yr[ok]) * as.numeric(nextjan - jan))
        }
        return(out)
    }
    # Legacy p87 has either a calendar year or an age. Both have limited
    # precision; retain the established July-1 convention, never an exact day.
    yr <- ifelse(is.finite(x) & x >= 1800, floor(x),
                 ifelse(is.finite(x) & x >= 0 & x <= 120, floor(birth_year + x), NA_real_))
    ok <- is.finite(yr) & yr <= year_cap & (is.na(birth_year) | yr >= birth_year)
    out <- rep(as.Date(NA), length(x))
    if (any(ok)) out[ok] <- as.Date(paste0(yr[ok], "-07-01"))
    out
}

# CAD includes angina/chronic IHD and cannot stand in for myocardial infarction
# when building a conventional MI/stroke/CVD-death composite.
ukb_mi_date <- function(dat) {
    cc <- intersect(c("fod_icd10_mi", "fod_icd9_mi", "fod_srd_mi"), names(dat))
    if (!length(cc)) {
        warning("No MI-specific date columns found. Rebuild ICD/SRD with the mi mappings; date_mi remains NA (CAD is retained as date_cad).", call. = FALSE)
        return(rep(as.Date(NA), nrow(dat)))
    }
    ukb_date_min(as.data.frame(dat[, cc, drop = FALSE]))
}
