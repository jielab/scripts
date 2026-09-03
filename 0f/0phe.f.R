date_follow_end <- as.Date("2023-04-01")

`%||%` <- function(a, b) if (!is.null(a)) a else b

log1 <- function(x) ifelse(is.finite(x) & x > 0, log(x), NA_real_)

set.seed(12345)


# 🚩 Runtime and file-output helpers
setwd2 <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  setwd(path)
}

step_header <- function(x) {
  cat("\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("# ", x, "\n", sep = "")
  cat("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  invisible(NULL)
}

.is_file_path <- function(x) is.character(x) && length(x) == 1 && !is.na(x)

.as_xlsx_sheets <- function(x) {
  if (is.data.frame(x)) x <- list(data = x)
  if (!is.list(x)) stop("write_xlsx() can only write a data frame or a list.", call. = FALSE)
  x <- x[!vapply(x, is.null, logical(1))]
  if (!length(x)) x <- list(data = data.frame())
  x <- lapply(x, function(z) {
    if (inherits(z, "sf")) z <- sf::st_drop_geometry(z)
    as.data.frame(z)
  })
  if (is.null(names(x))) names(x) <- paste0("sheet", seq_along(x))
  bad_names <- is.na(names(x)) | !nzchar(names(x))
  names(x)[bad_names] <- paste0("sheet", which(bad_names))
  names(x) <- make.unique(names(x), sep = "_")
  x
}

write_xlsx <- function(x, file, ...) {
  if (.is_file_path(x) && !.is_file_path(file)) {
    tmp <- x
    x <- file
    file <- tmp
  }
  if (!.is_file_path(file)) stop("write_xlsx() needs one character file path.", call. = FALSE)
  sheets <- .as_xlsx_sheets(x)
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(sheets, path = file, ...)
  } else if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(sheets, file = file, overwrite = TRUE, ...)
  } else {
    stop("write_xlsx() requires package 'writexl' or 'openxlsx'.", call. = FALSE)
  }
  invisible(file)
}

write_csv <- function(x, file, ..., row.names = FALSE) {
  if (.is_file_path(x) && !.is_file_path(file)) {
    tmp <- x
    x <- file
    file <- tmp
  }
  if (!.is_file_path(file)) stop("write_csv() needs one character file path.", call. = FALSE)
  utils::write.csv(x, file = file, row.names = row.names, ...)
  invisible(file)
}

write_raw_csv <- function(x, file, rawdir = get0("rawdir", ifnotfound = getwd()), ...) {
  write_csv(x, file.path(rawdir, file), ...)
}

save_rds <- function(x, file, ..., compress = TRUE) {
  if (.is_file_path(x) && !.is_file_path(file)) {
    tmp <- x
    x <- file
    file <- tmp
  }
  if (!.is_file_path(file)) stop("save_rds() needs one character file path.", call. = FALSE)
  saveRDS(x, file = file, ..., compress = compress)
  invisible(file)
}

save_plot <- function(p, file, width = 8, height = 6, w = width, h = height,
                      dpi = 300, bg = "white", units = "in",
                      limitsize = FALSE, ...) {
  ggplot2::ggsave(
    filename = file, plot = p, width = w, height = h,
    units = units, dpi = dpi, bg = bg, limitsize = limitsize, ...
  )
  invisible(file)
}

save_editable_plot <- function(p, file_stem, width, height, bg = "white", ...) {
  file <- if (grepl("\\.(png|jpe?g|tiff?|pdf|svg)$", file_stem, ignore.case = TRUE)) file_stem else paste0(file_stem, ".png")
  save_plot(p, file, width = width, height = height, bg = bg, ...)
}


# 🚩 Shared phenotype constants
center.lst <- c(
	"c10003" = "Stockport", "c11001" = "Manchester", "c11002" = "Oxford", "c11003" = "Cardiff", "c11004" = "Glasgow",
	"c11005" = "Edinburgh", "c11006" = "Stoke", "c11007" = "Reading", "c11008" = "Bury", "c11009" = "Newcastle",
	"c11010" = "Leeds", "c11011" = "Bristol", "c11012" = "Barts", "c11013" = "Nottingham", "c11014" = "Sheffield",
	"c11016" = "Liverpool", "c11017" = "Middlesbrough", "c11018" = "Hounslow", "c11020" = "Croydon", "c11021" = "Birmingham",
	"c11022" = "Swansea", "c11023" = "Wrexham", "c11024" = "Cheadle", "c11025" = "Cheadle", "c11026" = "Reading",
	"c11027" = "Newcastle", "c11028" = "Bristol"
)

dx.lst <- c(
	cvd_htn = "Hypertension", cvd_cad = "Coronary artery disease", mi = "Myocardial infarction",
	cvd_stroke_i = "Ischemic stroke", cvd_stroke_ih = "Intracerebral hemorrhage", cvd_stroke_sh = "Subarachnoid hemorrhage",
	cvd_hfail = "Heart failure", cvd_afib = "Atrial fibrillation",
	cvd_pad = "Peripheral arterial disease", cvd_dvt = "Deep vein thrombosis", cvd_pe = "Pulmonary embolism",
	t2dm = "Type 2 diabetes",
	acd = "All-cause dementia", ad = "Alzheimer's disease", depress = "Depression",
	copd = "Chronic obstructive pulmonary disease", asthma = "Asthma",
	ckd = "Chronic kidney disease", ra = "Rheumatoid arthritis",
	death = "All-cause mortality"
)

cvd.lst <- c(
	cvd_cad = "Coronary artery disease", cvd_stroke_i = "Ischemic stroke",
	cvd_hfail = "Heart failure", cvd_afib = "Atrial fibrillation",
	cvd_aosten = "Aortic valve stenosis", cvd_thaneu = "Thoracic aneurysm", cvd_abaneu = "Abdominal aneurysm",
	cvd_dvt = "Deep vein thrombosis", cvd_pe = "Pulmonary embolism",
	cvd_htn = "Hypertension",
	cvd_pad = "Peripheral arterial disease", cvd_cmyo = "Cardiomyopathy",
	cvd_stroke_ih = "Intracerebral hemorrhage", cvd_stroke_sh = "Subarachnoid hemorrhage", cvd_tia = "Transient ischemic attack",
	carrest = "Cardiac arrest", death = "All-cause mortality"
)

# surv.objs <- paste0("Surv(", Ys, ".t2e, ", Ys, ".Yt2e)")
vars.basic <- c("age", "sex", "tdi", "PC1", "PC2")
vars.basic1 <- c(vars.basic, "center")
vars.le8 <- c("diet.pts", "pa.pts", "smoke.pts", "bmi.pts", "nonhdl.pts", "hba1c.pts", "bp.pts", "sleep.pts")
names.le8 <- gsub("\\.pts$", "", vars.le8)


# 🚩 Data cleaning and variant matching
expo <- function(x) 1.1^x
rb <- function(x) if (is.null(x)) NA else round(x, 3)
rp <- function(x) if (is.null(x)) NA else signif(x, 2)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))

SNP_match <- function(datA, datA.CHR, datA.POS, datA.REF, datA.ALT = NULL,
                      datB, datB.CHR, datB.POS, datB.REF, datB.ALT = NULL,
                      by = c("CHR", "POS", "REF")){
	library(data.table)
	by <- toupper(by)
	if(!all(c("CHR","POS") %in% by)) stop("by must at least contain CHR and POS")
	if(!all(by %in% c("CHR","POS","REF","ALT"))) stop("by can only use CHR, POS, REF, ALT")

	getv <- function(d, nm) if(is.null(nm) || is.na(nm) || !nm %in% names(d)) rep(NA_character_, nrow(d)) else d[[nm]]
	nchr <- function(x){ x <- toupper(sub("^CHR", "", as.character(x))); x[x == "X"] <- "23"; x[x == "Y"] <- "24"; x }
	nal  <- function(x){ x <- toupper(trimws(as.character(x))); x[x %in% c("",".","NA","N")] <- NA_character_; x }

	A0 <- as.data.table(datA); B0 <- as.data.table(datB)
	A <- data.table(.A_id = seq_len(nrow(A0)), .chr = nchr(getv(A0, datA.CHR)), .pos = as.integer(getv(A0, datA.POS)),
	                .refA = nal(getv(A0, datA.REF)), .altA = nal(getv(A0, datA.ALT)))
	B <- data.table(.B_id = seq_len(nrow(B0)), .chr = nchr(getv(B0, datB.CHR)), .pos = as.integer(getv(B0, datB.POS)),
	                .refB = nal(getv(B0, datB.REF)), .altB = nal(getv(B0, datB.ALT)))

	z <- merge(A, B, by = c(".chr",".pos"), allow.cartesian = TRUE)
	if(!nrow(z)) return(data.table())

	in_set <- function(x, a, b) !is.na(x) & ((!is.na(a) & x == a) | (!is.na(b) & x == b))

	if(setequal(by, c("CHR","POS"))){
		z[, `:=`(keep = TRUE, match_type = "CHR_POS")]
	} else if("ALT" %in% by){
		z[, keep := !is.na(.refA) & !is.na(.altA) & !is.na(.refB) & !is.na(.altB) &
			((.refA == .refB & .altA == .altB) | (.refA == .altB & .altA == .refB))]
		z[, match_type := fifelse(.refA == .refB & .altA == .altB, "REF_ALT_exact",
			fifelse(.refA == .altB & .altA == .refB, "REF_ALT_swapped", NA_character_))]
	} else {
		z[, keep := in_set(.refB, .refA, .altA) | in_set(.refA, .refB, .altB)]
		z[, match_type := fifelse(.refB == .refA, "REF_REF",
			fifelse(.refB == .altA, "B_REF_matches_A_ALT",
			fifelse(!is.na(.altB) & .refA == .altB, "A_REF_matches_B_ALT", "allele_set")))]
	}

	z <- z[keep == TRUE]
	if(!nrow(z)) return(data.table())

	AA <- copy(A0); BB <- copy(B0)
	setnames(AA, paste0("A.", names(AA))); setnames(BB, paste0("B.", names(BB)))
	cbind(z[, .(match_chr = .chr, match_pos = .pos, match_type, A_id = .A_id, B_id = .B_id)], AA[z$.A_id], BB[z$.B_id])
}

neg2NA <- function(dat1, neg_list) {
	check_cols <- names(dat1)[vapply(dat1, function(x) is.numeric(x) || inherits(x, "integer64") || inherits(x, "Date") || is.character(x), logical(1))]
	neg_cols <- check_cols[vapply(dat1[check_cols], function(x) {
		if (is.character(x)) {
			v <- suppressWarnings(as.numeric(x))
			nv <- v[!is.na(v) & v < 0]
			length(nv) > 0 && all(nv %in% neg_list)
		} else {
			nv <- x[!is.na(x) & x < 0]
			length(nv) > 0 && all(nv %in% neg_list)
		}
	}, logical(1))]
	for (col in neg_cols) {
		if (is.character(dat1[[col]])) {
			v <- suppressWarnings(as.numeric(dat1[[col]]))
			dat1[[col]][!is.na(v) & v < 0] <- NA
		} else {
			dat1[[col]][!is.na(dat1[[col]]) & dat1[[col]] < 0] <- NA
		}
	}
	dat1
}

miss_cnt <- function(dat, eid_col = "eid", cuts = seq(0.1, 0.9, 0.1), right = TRUE) {
	bio_cols <- setdiff(names(dat), eid_col)
	miss <- colMeans(is.na(dat[bio_cols]))
	brks <- c(0, cuts)
	labs <- paste0("(", head(brks, -1) * 100, "%,", tail(brks, -1) * 100, "%]")
	bin <- rep(NA_character_, length(miss))
	bin[miss == 0] <- "==0%"
	bin[miss > 0] <- as.character(cut(miss[miss > 0], breaks = brks, include.lowest = FALSE, right = right, labels = labs))
	lvls <- c("==0%", labs)
	bin <- factor(bin, levels = lvls)
	out <- as.data.frame(table(bin), stringsAsFactors = FALSE)
	names(out) <- c("thr", "n_var")
	out$vars <- I(lapply(lvls, function(lv) names(miss)[bin == lv]))
	out
}

drop_10pct <- function(dat1, vars, imp = FALSE, mis_rate = 0.1, rowmax = mis_rate, colmax = mis_rate, k = 10, seed = 1) {
	mis <- colMeans(is.na(dat1[, vars]))
	Ydrop <- names(mis)[mis > mis_rate]
	dat1 <- dat1 %>% dplyr::select(-all_of(Ydrop))
	new_vars <- setdiff(vars, Ydrop)
	if (imp && length(new_vars) > 0) {
		message("Performing KNN imputation...")
		imp_res <- impute::impute.knn(as.matrix(dat1[, new_vars]), rowmax = rowmax, colmax = colmax, k = k, rng.seed = seed)
		dat1[, new_vars] <- as.data.frame(imp_res$data)
	}
	message(paste("Dropped:", length(Ydrop), "vars. Remaining:", length(new_vars)))
	return(list(dat = dat1, vars = new_vars))
}

clean_biom <- function(dat, eid_col = "eid", miss_col = 0.1, miss_row = 0.1, var_thr = 1e-6, imp = TRUE) {
	dat <- as.data.frame(dat)
	stopifnot(eid_col %in% names(dat))
	bio_cols <- setdiff(names(dat), eid_col)
	dat[bio_cols] <- lapply(dat[bio_cols], function(x) {
		x <- suppressWarnings(as.numeric(x))
		x[!is.finite(x)] <- NA_real_
		x
	})
	dat[[eid_col]] <- as.character(dat[[eid_col]])
	bio_cols <- setdiff(names(dat), eid_col)
	missv <- colMeans(is.na(dat[bio_cols]))
	keep_cols <- names(missv)[missv <= miss_col]
	dat <- dat[, c(eid_col, keep_cols), drop = FALSE]
	message("Cols kept after miss_col: ", length(keep_cols), " / ", length(missv))
	bio_cols <- setdiff(names(dat), eid_col)
	keep_row <- rowMeans(is.na(dat[bio_cols])) <= miss_row
	dat <- dat[keep_row, , drop = FALSE]
	message("Rows kept after miss_row: ", nrow(dat))
	bio_cols <- setdiff(names(dat), eid_col)
	vv <- sapply(dat[bio_cols], var, na.rm = TRUE)
	vdrop <- names(vv)[is.na(vv) | vv < var_thr]
	if (length(vdrop) > 0)
		dat <- dat[, setdiff(names(dat), vdrop), drop = FALSE]
	message("Cols dropped by var_thr: ", length(vdrop), "; remaining: ", ncol(dat) - 1)
	if (isTRUE(imp)) {
		bio_cols <- setdiff(names(dat), eid_col)
		datDrop <- drop_10pct(dat, vars = bio_cols, imp = TRUE, mis_rate = 0.1)
		dat <- datDrop$dat
	}
	dat
}

# Missing-phenotype imputation
impPhe <- function(dat1, Y, Xs, postfix = NULL) {
	  if (!Y %in% names(dat1) || all(is.na(dat1[[Y]]))) return(dat1)
	  dat1[[paste0(Y, "_mis")]] <- ifelse(is.na(dat1[[Y]]), 1, 0)
	  fit  <- lm(reformulate(Xs, Y), data = dat1)
	  pred <- predict(fit, newdata = dat1[, Xs, drop = FALSE])
	  new_name <- if (is.null(postfix)) Y else paste0(Y, postfix)
	  missing_idx <- which(is.na(dat1[[Y]]))
	  if (length(missing_idx) > 0) {  dat1[missing_idx, new_name] <- pred[missing_idx] }
	  return(dat1)
}

f3c <- function(x, keepNA = FALSE) { # Split into low, middle, and high groups.
	qunt <- ntile(x, 4)
	res <- ifelse(qunt == 1, "low",  ifelse(qunt == 4, "high", "middle"))
	if (keepNA) {res[is.na(x)] <- "TBD"; return(factor(res, levels = c("TBD", "low", "middle", "high")))
	} else {return(factor(res, levels = c("low", "middle", "high")))}
}

is.Date1 <- function(x, format = "%Y-%m-%d") {
	if (!is.character(x) || length(x) != 1) return(FALSE)
	out <- tryCatch(as.Date(x, format = format), error = function(e) NA)
	!is.na(out) && format(out, format) == x
}

range_sd <- function(dat, n = 3) {
	X_mean <- mean(dat, na.rm = TRUE); X_sd <- sd(dat, na.rm = TRUE); return(c(X_mean - n * X_sd, X_mean + n * X_sd))
}
del_sd <- function(x, n=3) (
	ifelse((x > (mean(x, na.rm = TRUE) + n * sd(x,na.rm = TRUE)) | x < (mean(x, na.rm = TRUE) - n * sd(x, na.rm = TRUE))), NA, x)
)
cap_sd <- function(x, n = 3) {
	mu <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE); function(v) pmin(pmax(v, mu - n * s), mu + n * s)
}
cap_pval <- function(x, p = 0.05) {
	qs <- quantile(x, c(p, 1 - p), na.rm = TRUE); function(v) pmin(pmax(v, qs[1]), qs[2])
}

check_dup <- function(x) {table(x[x %in% x[duplicated(x)]])}

split_multiple <- function(x, sep = "|", fun = "max") {
	parts <- strsplit(as.character(x), sep, fixed = TRUE)
	sapply(parts, function(v) {
		v <- v[nzchar(v)] # Drop empty string entries.
		if (length(v) == 0) return(NA_real_)
		nums <- suppressWarnings(as.numeric(v))
		if (fun %in% c("max","min","mean")) { return(do.call(fun, list(nums))) }
		if (grepl("^-?\\d+$", fun)) { # If fun is numeric, return the Nth element.
			N <- as.integer(fun)
			if (N > 0 && N <= length(nums)) return(nums[N])
			if (N < 0 && abs(N) <= length(nums)) return(nums[length(nums)+N+1])
			return(NA_real_)
		}
		stop("Invalid 'fun' argument.")
	})
}

score_wts <- function(dat1, mat.wts, thresh) {
	mat <- as.matrix(dat1[, names(mat.wts)])
	pw_sum <- (!is.na(mat)) %*% as.matrix(mat.wts)
	m_zero <- mat; m_zero[is.na(m_zero)] <- 0
	f_score <- (m_zero %*% as.matrix(mat.wts)) / pw_sum
	f_score[rowMeans(!is.na(mat)) < thresh] <- NA
	return(as.vector(f_score))
}

# P-value thinning for large association scans.
thinP0 <- function(dat1, P = 1e-3, p_col = "P", base_keep = 0.10, decade_factor = 0.5, min_keep = 0.01, seed = 1234) {
	stopifnot(is.data.table(dat1) || is.data.frame(dat1))
	set.seed(seed)
	p <- as.numeric(dat1[[p_col]])
	ok <- !is.na(p) & p > 0 & p <= 1; keep <- rep(TRUE, length(p))
	idx_thin <- which(ok & p > P)
	if (length(idx_thin) > 0) {
		p2 <- p[idx_thin]
		decade <- floor(log10(p2 / P))
		prob <- base_keep * (decade_factor ^ decade); prob <- pmax(prob, min_keep)
		keep[idx_thin] <- (runif(length(idx_thin)) < prob)
	}
	dat1[keep, ]
}

# Compress very small p-values for plotting.
thinP1 <- function(p, min_p = 1e-300, max_log10 = 50, compress_threshold = 10) {
	  p <- as.numeric(p)
	  p[!is.na(p) & (p < min_p | p == 0)] <- min_p
	  log10p <- -log10(p)
	  current_max <- max(log10p[is.finite(log10p)], na.rm = TRUE)
	  if (current_max <= max_log10) { return(p) }
	  scaling_factor <- (max_log10 - compress_threshold) / (current_max - compress_threshold)
	  log10p_new <- ifelse(log10p <= compress_threshold, log10p, compress_threshold + (log10p - compress_threshold) * scaling_factor)
	  return(10^(-log10p_new))
}

strsplit1 <- function(x, sep = "\\s+") {
	strsplit(x, split = sep)[[1]]
}
rowMeans1 <- function(x, per = 0.5) {
	x <- data.matrix(as.data.frame(x))
	NA.pt <- rowMeans(is.na(x))
	out <- rowMeans(x, na.rm = TRUE)
	out[NA.pt > per] <- NA_real_
	out[is.nan(out)] <- NA_real_
	return(out)
}
rowSums1 <- function(x, per = 0.5) {
	x <- data.matrix(as.data.frame(x))
	NA.pt <- rowMeans(is.na(x))
	out <- rowSums(x, na.rm = TRUE)
	out[NA.pt > per] <- NA_real_
	out[is.nan(out)] <- NA_real_
	return(out)
}
min1 <- function(..., na.rm = FALSE) {
	m <- suppressWarnings(min(..., na.rm = na.rm))
	if (!is.finite(m)) NA_real_ else m
}
bulk_rowMeans <- function(dat, fields, per = 0.8) {
	for (f in fields) {dat[[paste0(f, ".mi")]] <- rowMeans1(dplyr::select(dat, dplyr::starts_with(f)), per = per)}
	return(dat)
}

# Column-name harmonization for GWAS files.
replacement = c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')
pattern = c('^snp$|^rsid$|variant_id', '^chr$|^chrom', '^bp$|^pos$|^position|^base_pair', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele', '^eaf$|^freq$|a1freq|effect_allele_freq', '^n$|Neff', '^beta$|^effect$', '^se$|standard_error', '^p$|^pval$|^p_bolt_lmm')
replace_code <- function(code, sep, mapping_df) {
	base_part <- sub(paste0("\\", sep, ".*"), "", code)
	word <- mapping_df$V2[mapping_df$V1 == base_part]
	suffix <- sub(paste0("^[^", sep, "]*"), "", code)
	if (length(word) == 1 && !is.na(word)) return(paste0(word, suffix)) else return(code)
}

# Time-to-event helper.
#
# In addition to the baseline-relative variables, `bi2e` records attained age
# at diagnosis/censoring (birth -> event/censor).  It is intended for a Cox
# model with delayed entry at the participant's baseline age, for example
# Surv(age_at_baseline, bi2e, Yt2e).  It does *not* imply that an omic value
# measured in adulthood was observed, or can be projected, at birth.
t2e <- function(dat, domain, Y_date, birth_date, date_attend, date_lost, date_death,
	              date_end = date_follow_end, prefix = NA, time_unit = "year") {
	if (!is.na(prefix) && all(paste0(prefix, c(".Yt2e", ".Yr2e", ".t2e", ".r2e", ".b2e", ".bi2e")) %in% names(dat))) return(dat)
	dt_att <- if (inherits(date_attend, "Date")) date_attend else as.Date(dat[[date_attend]])
	for (v in c(birth_date, Y_date, date_lost, date_death, domain)) if (!is.null(v) && !is.na(v) && v %in% names(dat)) dat[[v]] <- as.Date(dat[[v]])
	date_end0 <- if (inherits(date_end, "Date")) date_end else as.Date(date_end)
	# Use the earliest observed censoring date.  The administrative end date is
	# always present, so pmin() cannot return Inf when loss/death are both NA.
	censor_dt <- pmin(dat[[date_lost]], dat[[date_death]], date_end0, na.rm = TRUE)
	d_var <- if (!is.na(domain) && domain %in% names(dat)) domain else if (!is.na(domain) && paste0("icd10Ct_", domain) %in% names(dat)) paste0("icd10Ct_", domain) else NA_character_
	d_is_date <- !is.na(d_var) && d_var == domain
	dat <- dat %>% mutate(
	  Yd = if_else(.data[[Y_date]] == dt_att, as.Date(NA), .data[[Y_date]]),
	  is_aft = !is.na(Yd) & Yd > dt_att, is_bef = !is.na(Yd) & Yd < dt_att, is_nev = is.na(Yd),
	  d_happ = if (!is.na(d_var)) { if (d_is_date) (!is.na(.data[[d_var]]) & .data[[d_var]] > dt_att) else (.data[[d_var]] > 0) } else FALSE,
	  is_dty = is_nev & d_happ,
	  Yt2e = case_when(is_aft ~ 1, is_dty ~ NA_real_, is_nev ~ 0, TRUE ~ NA_real_),
	  Yr2e = case_when(is_bef ~ 1, is_dty ~ NA_real_, is_nev ~ 0, TRUE ~ NA_real_),
	  f_dt = if_else(is_aft, Yd, censor_dt),
	  r_dt = if_else(is_bef, Yd, as.Date(.data[[birth_date]])),
	  bi_dt = case_when(is_bef ~ Yd, !is.na(Yt2e) ~ f_dt, TRUE ~ as.Date(NA)),
	  t2e = if_else(!is.na(Yt2e), as.numeric(f_dt - dt_att), NA_real_),
	  r2e = if_else(!is.na(Yr2e), as.numeric(dt_att - r_dt), NA_real_),
	  b2e = case_when(Yr2e == 1 ~ -r2e, Yt2e == 1 ~  t2e, TRUE ~ NA_real_),
	  bi2e = if_else(!is.na(bi_dt) & !is.na(.data[[birth_date]]),
	                   as.numeric(bi_dt - .data[[birth_date]]), NA_real_)
	) %>% { if (time_unit == "year") mutate(., t2e = t2e/365.25, r2e = r2e/365.25,
	                                                b2e = b2e/365.25, bi2e = bi2e/365.25) else . }
	if (!is.na(prefix)) dat <- dat %>% rename_with(~ paste0(prefix, ".", .), c(Yt2e, Yr2e, t2e, r2e, b2e, bi2e))
	dat %>% dplyr::select(-Yd, -is_aft, -is_bef, -is_nev, -d_happ, -is_dty, -f_dt, -r_dt, -bi_dt)
}

# ICD-10 helper.
get_icd <- function(ver, dat.code = NA, dat.date = NA, dat.this, code.this, label.this,
                    save_aod_list = TRUE, save_aod_text = TRUE) {
	fod_col <- paste0("fod_", ver, "_", label.this)
	lod_col <- paste0("lod_", ver, "_", label.this)
	ct_col  <- paste0("cnt_", ver, "_", label.this)
	aod_col <- paste0("aod_", ver, "_", label.this)
	aodtxt_col <- paste0("aodtxt_", ver, "_", label.this)

	need_cols <- c(fod_col, lod_col, ct_col)
	if (isTRUE(save_aod_list)) need_cols <- c(need_cols, aod_col)
	if (isTRUE(save_aod_text)) need_cols <- c(need_cols, aodtxt_col)
	if (is.na(code.this) || all(need_cols %in% names(dat.this))) return(dat.this)

	if (!is.data.frame(dat.code)) dat.code <- readRDS(paste0(indir, "/Rdata/", ver, ".code0.rds"))
	if (!is.data.frame(dat.date)) dat.date <- readRDS(paste0(indir, "/Rdata/", ver, ".date0.rds"))
	code_mat <- as.matrix(dat.code[, -1, drop = FALSE])
	date_mat <- as.matrix(dat.date[, -1, drop = FALSE])
	match_mask <- matrix(grepl(code.this, code_mat, perl = TRUE), nrow = nrow(code_mat))
	date_mat[!match_mask] <- NA
	code_mat[!match_mask] <- NA

	out <- data.frame(eid = as.character(dat.code$eid), stringsAsFactors = FALSE)
	out[[ct_col]] <- rowSums(!is.na(code_mat))

	temp_dates <- as.data.frame(date_mat)
	temp_dates[] <- lapply(temp_dates, function(x) {
		if (inherits(x, "Date")) return(x)
		if (inherits(x, "IDate")) return(as.Date(x))
		if (is.character(x)) return(as.Date(x))
		y <- suppressWarnings(as.numeric(unclass(x)))
		y[!is.finite(y)] <- NA_real_
		as.Date(y, origin = "1970-01-01")
	})

	min_dates <- do.call(pmin, c(temp_dates, list(na.rm = TRUE)))
	max_dates <- do.call(pmax, c(temp_dates, list(na.rm = TRUE)))
	min_dates <- as.Date(min_dates, origin = "1970-01-01")
	max_dates <- as.Date(max_dates, origin = "1970-01-01")
	min_dates[!is.finite(as.numeric(min_dates))] <- as.Date(NA)
	max_dates[!is.finite(as.numeric(max_dates))] <- as.Date(NA)

	out[[fod_col]] <- min_dates
	out[[lod_col]] <- max_dates

	if (isTRUE(save_aod_list) || isTRUE(save_aod_text)) {
		aod_list <- lapply(seq_len(nrow(temp_dates)), function(i) {
			z <- unlist(temp_dates[i, ], use.names = FALSE)
			z <- as.Date(z, origin = "1970-01-01")
			z <- sort(unique(z[!is.na(z)]))
			z
		})
		if (isTRUE(save_aod_list)) out[[aod_col]] <- I(aod_list)
		if (isTRUE(save_aod_text)) {
			out[[aodtxt_col]] <- vapply(aod_list, function(z) {
				if (!length(z)) NA_character_ else paste(format(z, "%Y-%m-%d"), collapse = "|")
			}, character(1))
		}
	}

	# Remove overlapping columns before merge, otherwise base::merge may make .x/.y columns.
	replace_cols <- intersect(setdiff(names(out), "eid"), names(dat.this))
	if (length(replace_cols)) dat.this <- dat.this[, setdiff(names(dat.this), replace_cols), drop = FALSE]
	merge(dat.this, out, by = "eid", all.x = TRUE, sort = FALSE)
}


# Build all-trait ICD outputs from one sparse event table per source.
get_icd_fast <- function(ver, dat.code = NA, dat.date = NA, fn,
                         save_aod_list = FALSE, save_aod_text = TRUE,
                         normalize_dot = TRUE, verbose = TRUE) {
	if (!requireNamespace("data.table", quietly = TRUE)) stop("get_icd_fast() requires data.table.", call. = FALSE)
	DT <- data.table::as.data.table
	as_date_fast <- function(x) {
		if (inherits(x, "Date")) return(x)
		if (inherits(x, "IDate")) return(as.Date(x))
		if (inherits(x, "POSIXt")) return(as.Date(x))
		if (is.character(x)) return(suppressWarnings(as.Date(x)))
		y <- suppressWarnings(as.numeric(unclass(x)))
		y[!is.finite(y)] <- NA_real_
		as.Date(y, origin = "1970-01-01")
	}
	clean_code_fast <- function(x) {
		x <- toupper(trimws(as.character(x)))
		x[x %in% c("", "NA", "NAN")] <- NA_character_
		if (isTRUE(normalize_dot)) x <- gsub("\\.", "", x)
		x
	}
	clean_pattern_fast <- function(x) {
		x <- toupper(trimws(as.character(x)))
		if (isTRUE(normalize_dot)) x <- gsub("\\.", "", x)
		x
	}
	wide_one <- function(ans, value, prefix, all_eids, all_traits, default = NA) {
		out <- data.table::data.table(eid = all_eids)
		if (nrow(ans) > 0 && value %in% names(ans)) {
			w <- data.table::dcast(ans, eid ~ trait, value.var = value)
			old <- setdiff(names(w), "eid")
			data.table::setnames(w, old, paste0(prefix, old))
			out <- merge(out, w, by = "eid", all.x = TRUE, sort = FALSE)
		}
		need <- paste0(prefix, all_traits)
		miss <- setdiff(need, names(out))
		for (cc in miss) out[[cc]] <- default
		out <- out[, c("eid", need), with = FALSE]
		out
	}

	if (!is.data.frame(dat.code)) dat.code <- readRDS(paste0(indir, "/Rdata/", ver, ".code0.rds"))
	if (!is.data.frame(dat.date)) dat.date <- readRDS(paste0(indir, "/Rdata/", ver, ".date0.rds"))
	code0 <- DT(dat.code); date0 <- DT(dat.date)
	code0[, eid := as.character(eid)]; date0[, eid := as.character(eid)]
	stopifnot(identical(code0$eid, date0$eid))

	fn0 <- DT(fn)[, 1:2, with = FALSE]
	data.table::setnames(fn0, c("pattern", "trait"))
	fn0 <- fn0[!is.na(pattern) & nzchar(as.character(pattern)) & !is.na(trait) & nzchar(as.character(trait))]
	fn0[, `:=`(pattern = as.character(pattern), trait = as.character(trait))]
	fn0[, pattern_clean := clean_pattern_fast(pattern)]
	all_traits <- unique(fn0$trait)
	all_eids <- code0$eid
	make_aodtxt <- isTRUE(save_aod_text) || isTRUE(save_aod_list)

	if (isTRUE(verbose)) message("get_icd_fast(): ", ver, " | traits=", nrow(fn0), " | participants=", nrow(code0))
	if (!nrow(fn0)) return(as.data.frame(data.table::data.table(eid = all_eids)))

	code_cols <- setdiff(names(code0), "eid")
	date_cols <- setdiff(names(date0), "eid")
	nslot <- min(length(code_cols), length(date_cols))
	code_cols <- code_cols[seq_len(nslot)]
	date_cols <- date_cols[seq_len(nslot)]

	code_mat <- as.matrix(code0[, ..code_cols])
	idx <- which(!is.na(code_mat) & code_mat != "", arr.ind = TRUE)
	if (!nrow(idx)) {
		ans <- data.table::data.table(eid = character(), trait = character(), cnt = integer(), fod = as.Date(character()), lod = as.Date(character()), aodtxt = character())
	} else {
		long <- data.table::data.table(eid = code0$eid[idx[, 1]], row_idx = idx[, 1], slot = idx[, 2], code = clean_code_fast(code_mat[idx]))
		rm(code_mat); gc(verbose = FALSE)
		long[, date := as.Date(NA)]
		for (j in seq_len(nslot)) {
			ii <- which(long$slot == j)
			if (length(ii)) data.table::set(long, i = ii, j = "date", value = as_date_fast(date0[[date_cols[j]]])[long$row_idx[ii]])
		}
		long[, row_idx := NULL]
		long <- long[!is.na(code) & nzchar(code)]
		if (isTRUE(verbose)) message("get_icd_fast(): ", ver, " | coded events=", format(nrow(long), big.mark = ","), " | dated events=", format(sum(!is.na(long$date)), big.mark = ","), " | unique codes=", length(unique(long$code)))

		if (!nrow(long)) {
			ans <- data.table::data.table(eid = character(), trait = character(), cnt = integer(), fod = as.Date(character()), lod = as.Date(character()), aodtxt = character())
		} else {
			data.table::setkey(long, code)
			codes <- unique(long$code)
			hits_cnt <- vector("list", nrow(fn0))
			hits_date <- vector("list", nrow(fn0))
			for (i in seq_len(nrow(fn0))) {
				if (isTRUE(verbose) && (i == 1 || i == nrow(fn0) || i %% 25 == 0)) message("  ", ver, " trait ", i, "/", nrow(fn0), ": ", fn0$trait[i])
				hit_codes <- codes[grepl(fn0$pattern_clean[i], codes, perl = TRUE)]
				if (length(hit_codes)) {
					z <- long[data.table::J(hit_codes), .(eid, slot, code, date), nomatch = 0L]
					if (nrow(z)) {
						z[, trait := fn0$trait[i]]
						hits_cnt[[i]] <- z[, .(eid, trait, slot, code)]
						zd <- z[!is.na(date)]
						if (nrow(zd)) hits_date[[i]] <- zd[, .(eid, trait, slot, code, date)]
					}
				}
			}
			hit_cnt <- data.table::rbindlist(hits_cnt, use.names = TRUE, fill = TRUE)
			hit_date <- data.table::rbindlist(hits_date, use.names = TRUE, fill = TRUE)

			if (nrow(hit_cnt)) {
				hit_cnt <- unique(hit_cnt, by = c("eid", "trait", "slot", "code"))
				ans_cnt <- hit_cnt[, .(cnt = .N), by = .(eid, trait)]
			} else {
				ans_cnt <- data.table::data.table(eid = character(), trait = character(), cnt = integer())
			}

			if (nrow(hit_date)) {
				hit_date <- unique(hit_date, by = c("eid", "trait", "slot", "code", "date"))
				ans_date <- hit_date[, .(
					fod = min(date, na.rm = TRUE),
					lod = max(date, na.rm = TRUE),
					aodtxt = if (isTRUE(make_aodtxt)) paste(format(sort(unique(date)), "%Y-%m-%d"), collapse = "|") else NA_character_
				), by = .(eid, trait)]
			} else {
				ans_date <- data.table::data.table(eid = character(), trait = character(), fod = as.Date(character()), lod = as.Date(character()), aodtxt = character())
			}

			ans <- merge(ans_cnt, ans_date, by = c("eid", "trait"), all = TRUE, sort = FALSE)
			ans[is.na(cnt), cnt := 0L]
		}
	}

	out <- data.table::data.table(eid = all_eids)
	out <- merge(out, wide_one(ans, "fod", paste0("fod_", ver, "_"), all_eids, all_traits, default = as.Date(NA)), by = "eid", all.x = TRUE, sort = FALSE)
	out <- merge(out, wide_one(ans, "lod", paste0("lod_", ver, "_"), all_eids, all_traits, default = as.Date(NA)), by = "eid", all.x = TRUE, sort = FALSE)
	out <- merge(out, wide_one(ans, "cnt", paste0("cnt_", ver, "_"), all_eids, all_traits, default = 0L), by = "eid", all.x = TRUE, sort = FALSE)
	cnt_cols <- grep(paste0("^cnt_", ver, "_"), names(out), value = TRUE)
	for (cc in cnt_cols) out[is.na(get(cc)), (cc) := 0L]

	if (isTRUE(make_aodtxt)) {
		aodtxt_wide <- wide_one(ans, "aodtxt", paste0("aodtxt_", ver, "_"), all_eids, all_traits, default = NA_character_)
		if (isTRUE(save_aod_text)) out <- merge(out, aodtxt_wide, by = "eid", all.x = TRUE, sort = FALSE)
		if (isTRUE(save_aod_list)) {
			# Optional list columns. Keep disabled by default because they greatly increase all.rds size.
			aodtxt_cols <- setdiff(names(aodtxt_wide), "eid")
			for (cc in aodtxt_cols) {
				newcc <- sub("^aodtxt_", "aod_", cc)
				out[[newcc]] <- I(lapply(aodtxt_wide[[cc]], function(x) {
					if (is.na(x) || !nzchar(x)) as.Date(character()) else as.Date(strsplit(x, "\\|", perl = TRUE)[[1]])
				}))
			}
		}
	}
	as.data.frame(out)
}

# FOD helper.
get_fod <- function(dat.fod, dat.this, code.this, label.this) {
	new_col_name <- paste0("fod_ref_", label.this)
	if (new_col_name %in% names(dat.this)) return(dat.this)
	if (!is.data.frame(dat.fod)) dat.fod <- readRDS(paste0("D:/data/ukb/phe/Rdata/fod.date0.rds"))
	  target_cols <- grep(code.this, names(dat.fod), value = TRUE)
	  if (length(target_cols) == 0) { warning(paste("No columns found for pattern:", code.this)); return(dat.this) }
	  out <- dat.fod[, "eid", drop = FALSE]
	  subset_dates <- dat.fod[, target_cols, drop = FALSE]
	  subset_dates[] <- lapply(subset_dates, as.Date)
	  out[[new_col_name]] <- do.call(pmin, c(subset_dates, na.rm = TRUE))
	  out <- merge(dat.this, out, by = "eid", all.x = TRUE, sort = FALSE)
	  return(out)
}

exclude_top_mis <- function(dat, th = 0.25) {
	mis <- dat %>% summarise(across(everything(), ~ mean(is.na(.)))) %>% gather(variable, mis_rate) %>% arrange(desc(mis_rate))
	exclude <- mis[mis$mis_rate > th, "variable"]
	out <- dat %>% dplyr::select(-all_of(exclude))
	return(out)
}

keep_top_sd <- function(dat, th = 0.25) {
	means <- apply(dat, 2, mean, na.rm = TRUE)
	sds   <- apply(dat, 2, sd, na.rm = TRUE)
	cv <- sds / abs(means); cv[is.nan(cv) | is.infinite(cv)] <- 0
	keep <- order(cv, decreasing = TRUE)[1:ceiling(th * sum(cv > 0))]
	dat <- dat[, keep, drop = FALSE]
	return(dat)
}

rel_seek <- function(eid, cutoff = 0.0882, id1_col = "ID1", id2_col = "ID2", kin_col = "Kinship") {
	rel <- read.table(paste0(dir0, '/data/ukb/gen/typ/ukb_rel.dat'), header = TRUE)
	rel_sub <- rel[rel[[kin_col]] > cutoff & rel[[id1_col]] %in% eid & rel[[id2_col]] %in% eid, c(id1_col, id2_col, kin_col)]
	rel_sub$rowid <- seq_len(nrow(rel_sub))
	work <- rel_sub[, c(id1_col, id2_col, "rowid")]; names(work) <- c("ID1", "ID2", "rowid")
	ids_all <- sort(unique(c(work$ID1, work$ID2))); n_ids   <- length(ids_all)
	idx_all <- match(c(work$ID1, work$ID2), ids_all)
	deg     <- tabulate(idx_all, nbins = n_ids)
	remove_ids    <- vector("integer", 0L)
	removed_pairs <- rel_sub[0, c(id1_col, id2_col, kin_col)]
	removed_pairs$removed_id <- removed_pairs[[id1_col]][0]

	while (nrow(work) > 0L) { # Greedy loop: each step removes the ID with the most relatives.
		i_max <- which.max(deg)
		id_remove <- ids_all[i_max]
		rows <- work$ID1 == id_remove | work$ID2 == id_remove # Find all edges containing this ID.
		if (!any(rows)) {deg[i_max] <- 0L; next}
		rm_rowid <- work$rowid[rows]
		tmp <- rel_sub[match(rm_rowid, rel_sub$rowid), c(id1_col, id2_col, kin_col)]
		tmp$removed_id <- id_remove
		removed_pairs <- rbind(removed_pairs, tmp)
		partners <- ifelse(work$ID1[rows] == id_remove, work$ID2[rows], work$ID1[rows])
		tab_partners <- table(partners)
		idx_partners <- match(as.integer(names(tab_partners)), ids_all)
		deg[idx_partners] <- deg[idx_partners] - as.integer(tab_partners)
		deg[i_max] <- 0L # Zero out the degree for the removed ID.
		work <- work[!rows, , drop = FALSE] # Remove these edges from the working table.
		remove_ids <- c(remove_ids, id_remove) # Record the removed ID.
	}
	remove_ids <- unique(remove_ids)
	list(remove_ids = remove_ids, pairs_removed = removed_pairs)
}



merge_check <- function(..., by = "eid", by.x = NULL, by.y = NULL,
                        update_parent = TRUE, verbose = TRUE) {
	exprs <- as.list(substitute(list(...)))[-1]
	dots <- list(...)
	if (!length(dots)) stop("merge_check() needs at least one data frame or one list of data frames.", call. = FALSE)

	is_df_list <- length(dots) == 1L && is.list(dots[[1]]) && !is.data.frame(dots[[1]]) &&
		all(vapply(dots[[1]], is.data.frame, logical(1)))
	if (is_df_list) {
		dfs <- dots[[1]]
		nm <- names(dfs)
		if (is.null(nm) || any(!nzchar(nm))) nm <- paste0("dat", seq_along(dfs))
		names(dfs) <- nm
	} else {
		if (!all(vapply(dots, is.data.frame, logical(1)))) stop("merge_check() inputs must be data frames, or one list of data frames.", call. = FALSE)
		dfs <- dots
		nm <- names(dfs)
		ex_nm <- vapply(exprs, function(z) paste(deparse(z), collapse = ""), character(1))
		if (is.null(nm)) nm <- rep("", length(dfs))
		nm[!nzchar(nm)] <- ex_nm[!nzchar(nm)]
		nm[!nzchar(nm)] <- paste0("dat", which(!nzchar(nm)))
		names(dfs) <- nm
	}

	if (!is.null(by.x) || !is.null(by.y)) {
		if (length(dfs) != 2L) stop("merge_check(): by.x/by.y are supported only for two data frames.", call. = FALSE)
		if (is.null(by.x) || is.null(by.y)) stop("merge_check(): use both by.x and by.y, or use by only.", call. = FALSE)
		key_list <- list(as.character(by.x), as.character(by.y))
		key_label <- paste0("by.x=", paste(by.x, collapse = "+"), ", by.y=", paste(by.y, collapse = "+"))
	} else {
		key_list <- rep(list(as.character(by)), length(dfs))
		key_label <- paste(by, collapse = "+")
	}

	fmt <- function(x) {
		x <- unique(as.character(x))
		x <- x[nzchar(x)]
		if (!length(x)) return("")
		if (length(x) == 1L) return(x)
		if (length(x) == 2L) return(paste(x, collapse = " and "))
		paste0(paste(head(x, -1), collapse = ", "), ", and ", tail(x, 1))
	}

	first_from <- setNames(character(0), character(0))
	seen_order <- character(0)
	drop_log <- list()
	log_i <- 0L

	for (i in seq_along(dfs)) {
		df <- as.data.frame(dfs[[i]], stringsAsFactors = FALSE, check.names = FALSE)
		keys <- key_list[[i]]
		missing_keys <- setdiff(keys, names(df))
		if (length(missing_keys)) stop("merge_check(): key column(s) missing in ", nm[i], ": ", paste(missing_keys, collapse = ", "), call. = FALSE)

		# Remove duplicated column names within the same data frame, keeping the first occurrence.
		internal_dup <- duplicated(names(df))
		if (any(internal_dup)) {
			dropped <- names(df)[internal_dup]
			log_i <- log_i + 1L
			drop_log[[log_i]] <- data.frame(
				variable = dropped,
				dropped_from = nm[i],
				kept_from = nm[i],
				reason = "duplicated_within_data_frame",
				stringsAsFactors = FALSE
			)
			df <- df[, !internal_dup, drop = FALSE]
		}

		# Remove variables that appeared in earlier data frames, except merge keys needed by base::merge().
		non_key <- setdiff(names(df), keys)
		dup_vars <- intersect(non_key, names(first_from))
		if (length(dup_vars)) {
			log_i <- log_i + 1L
			drop_log[[log_i]] <- data.frame(
				variable = dup_vars,
				dropped_from = nm[i],
				kept_from = unname(first_from[dup_vars]),
				reason = "duplicated_across_data_frames",
				stringsAsFactors = FALSE
			)
			df <- df[, !(names(df) %in% dup_vars), drop = FALSE]
		}

		new_vars <- setdiff(names(df), keys)
		new_vars <- setdiff(new_vars, names(first_from))
		if (length(new_vars)) {
			first_from[new_vars] <- nm[i]
			seen_order <- c(seen_order, new_vars)
		}
		dfs[[i]] <- df
	}

	report <- if (length(drop_log)) do.call(rbind, drop_log) else data.frame(variable = character(), dropped_from = character(), kept_from = character(), reason = character())
	attr(dfs, "merge_check_report") <- report

	if (verbose) {
		if (!nrow(report)) {
			message("merge_check(): all columns except the by column(s) are unique [", key_label, "].")
		} else {
			message("merge_check(): duplicated non-key columns detected [", key_label, "]; keeping the first occurrence.")

			# Across-data-frame duplicated variables.
			rep2 <- report[report$reason == "duplicated_across_data_frames", , drop = FALSE]
			if (nrow(rep2)) {
				for (v in unique(rep2$variable)) {
					rv <- rep2[rep2$variable == v, , drop = FALSE]
					kept <- rv$kept_from[1]
					dropped <- unique(rv$dropped_from)
					appeared <- unique(c(kept, dropped))
					message("  - ", v, " appeared in ", fmt(appeared), "; keeping ", v, " in ", kept, "; dropping ", v, " from ", fmt(dropped), ".")
				}
			}

			# Internal duplicated variables.
			rep1 <- report[report$reason == "duplicated_within_data_frame", , drop = FALSE]
			if (nrow(rep1)) {
				for (ii in seq_len(nrow(rep1))) {
					message("  - ", rep1$variable[ii], " appeared more than once in ", rep1$dropped_from[ii], "; keeping the first ", rep1$variable[ii], " in ", rep1$kept_from[ii], ".")
				}
			}
		}
	}

	if (isTRUE(update_parent)) {
		penv <- parent.frame()
		for (i in seq_along(dfs)) {
			obj <- nm[i]
			if (nzchar(obj) && make.names(obj) == obj) assign(obj, dfs[[i]], envir = penv)
		}
	}

	return(dfs)
}
