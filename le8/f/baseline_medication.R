.le8_med_cache <- new.env(parent = emptyenv())
le8_attach_baseline_medication <- function(x, root = get0("indir", ifnotfound = Sys.getenv("YGB_PHE", os_path("DIR_YGB")))) {
	# Raw i0 responses retain -7 (none of the above), which phe.R may have
	# converted to NA. Never infer a negative response from that processed NA.
	path <- Sys.getenv("LE8_BASELINE_MED_FILE", file.path(root, "rap/vip.tab.gz"))
	if (!file.exists(path)) return(x)
	explicit <- nzchar(Sys.getenv("LE8_BASELINE_MED_FILE", ""))
	key <- paste(path, file.info(path)$size, as.numeric(file.info(path)$mtime), sep = "|")
	if (!identical(.le8_med_cache$key, key)) {
		hdr <- if (grepl("[.]rds$", path, ignore.case = TRUE)) names(readRDS(path)) else names(data.table::fread(path, nrows = 0, showProgress = FALSE))
		cols <- grep("^(p)?(6153|6177)(_i0(_a[0-9]+)?|[.]0[.][0-9]+|-0[.][0-9]+)$", hdr, value = TRUE)
		ids <- intersect(c("eid", "f.eid", "IID"), hdr)
		if (!length(cols) || !length(ids)) {
			if (explicit) stop("LE8_BASELINE_MED_FILE needs eid and explicit baseline p6153/p6177 fields")
			return(x)
		}
		med <- if (grepl("[.]rds$", path, ignore.case = TRUE)) as.data.frame(readRDS(path))[, c(ids[1], cols), drop = FALSE] else
			as.data.frame(data.table::fread(path, select = c(ids[1], cols), showProgress = FALSE))
		names(med)[1] <- "eid" ; med$eid <- as.character(med$eid)
		if (anyNA(med$eid) || anyDuplicated(med$eid)) stop("Raw baseline medication input has missing/duplicate eid")
		names(med)[ - 1] <- paste0(".le8_raw_med_i0_", seq_along(cols))
		.le8_med_cache$data <- med ; .le8_med_cache$key <- key ; .le8_med_cache$source <- paste(path, paste(cols, collapse = ","), sep = " : ")
	}
	med <- .le8_med_cache$data ; i <- match(as.character(x$eid), med$eid)
	for (v in setdiff(names(med), "eid")) x[[v]] <- med[[v]][i]
	attr(x, "baseline_med_source") <- .le8_med_cache$source ; x
}
le8_baseline_med_values <- function(x) {
	n <- nrow(x)
	# Only explicitly tagged first-visit responses, or explicitly supplied columns.
	medcols <- trimws(strsplit(Sys.getenv("LE8_BASELINE_MED_COLUMNS", ""), ",", fixed = TRUE)[[1]])
	medcols <- medcols[nzchar(medcols)]
	if (!length(medcols)) medcols <- grep("^([.]le8_raw_med_i0_|drug[.]big3.*[_.]i0([_.]|$)|(p)?(6153|6177)_i0(_a[0-9]+)?$)", names(x), value = TRUE)
	if (any(!medcols %in% names(x))) stop("Unknown baseline medication column")
	meds <- matrix(NA_integer_, n, 3, dimnames = list(NULL, c("drug.lipid", "drug.htn", "drug.dm")))
	if (length(medcols)) {
		patterns <- do.call(paste, c(lapply(x[medcols], as.character), sep = " "))
		keys <- unique(patterns)
		lookup <- t(vapply(keys, function(s) {
			tokens <- suppressWarnings(as.integer(regmatches(s, gregexpr("-?[0-9]+", s))[[1]]))
			known <- any(tokens %in% c( - 7L, 1 : 5)) && !any(tokens %in% c( - 1L, - 3L)) &&
				!(any(tokens ==  - 7L, na.rm = TRUE) && any(tokens %in% 1 : 5))
			if (known) as.integer(1 : 3 %in% tokens) else rep(NA_integer_, 3)
		}, integer(3)))
		meds[, ] <- lookup[match(patterns, keys), , drop = FALSE]
	}
	list(values = meds, columns = medcols)
}
