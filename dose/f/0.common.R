source(file.path(dirname(dirname(dirname(normalizePath(tail(Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile)), 1)[[1]])))), "0f", "png_output.R"), local = TRUE)

# Shared paths for DOSE. DOSE_PROJECT is retained as an alias for the results root.
# Batch graphics helpers may open the default device while measuring grobs.
# Discard implicit pages instead of creating Rplots.pdf; explicit exports still work.
if (!interactive()) options(device=function(...) grDevices::pdf(file=NULL, ...))

dose_runtime_dir <- function() Sys.getenv("DOSE_RUNTIME_DIR", file.path(path.expand("~"), ".local/share/dose"))
if (dir.exists(file.path(dose_runtime_dir(), "R-library")))
    .libPaths(unique(c(file.path(dose_runtime_dir(), "R-library"), .libPaths())))
dose_drive <- function() if (.Platform$OS.type == "windows") "D:" else "/mnt/d"
dose_output_dir <- function() Sys.getenv("DOSE_OUTDIR", Sys.getenv("DOSE_PROJECT", file.path(dose_drive(), "analysis/dose")))
dose_script_dir <- function() Sys.getenv("DOSE_SCRIPT_DIR", file.path(dose_drive(), "scripts/dose/f"))
dose_ukb_dir <- function() Sys.getenv("DOSE_UKB", file.path(dose_drive(), "data/ukb"))
dose_files_dir <- function() Sys.getenv("DOSE_FILES_DIR", file.path(dose_drive(), "files/diet"))
# Study-specific dictionaries never fall back to the generic UKB lists.
dose_definition_file <- function(name) {
    stopifnot(length(name) == 1L, name %in% c("diet", "vip", "icd10_other"))
    file.path(dose_ukb_dir(), "common", paste0("DOSE.", name, ".lst"))
}
dose_helper_dir <- function() Sys.getenv("DOSE_HELPER_DIR", file.path(dose_drive(), "scripts/0f"))
dose_mri_dictionary <- function() Sys.getenv("DOSE_MRI_DICTIONARY", file.path(dose_files_dir(), "MRI_analysis_dictionary.csv"))
dose_python <- function() Sys.getenv("DOSE_PYTHON", Sys.getenv("RETICULATE_PYTHON", if (.Platform$OS.type == "windows") "py" else "python3"))

dose_replace <- function() identical(toupper(Sys.getenv("DOSE_REPLACE", "FALSE")), "TRUE")

# Commit an RDS only after serialization succeeds. An interrupt leaves the last
# complete file intact; checkpoint manifests are committed in the same way.
dose_save_rds <- function(object, file, ...) {
    if (inherits(file, "connection")) return(base::saveRDS(object, file, ...))
    temp <- tempfile(paste0(".", basename(file), "."), tmpdir=dirname(file))
    on.exit(unlink(temp), add=TRUE)
    base::saveRDS(object, temp, ...)
    if (!file.rename(temp, file)) stop("Could not commit RDS: ", file)
    invisible(NULL)
}

dose_cached_phase <- function(name, inputs, outputs, code) {
    expr <- substitute(code)
    stamp <- file.path(dose_output_dir(), ".checkpoints/foodwas/preprocess", paste0(name, ".rds"))
    info <- function(paths) {
        z <- file.info(paths)
        data.frame(path=as.character(paths), size=z$size, mtime=as.numeric(z$mtime), row.names=NULL)
    }
    if (any(!file.exists(inputs))) stop("Missing phase input: ", paste(inputs[!file.exists(inputs)], collapse=", "))
    identity <- list(code=expr, inputs=info(inputs),
        helpers=tools::md5sum(c(file.path(dose_script_dir(), "0.common.R"),
            file.path(dose_helper_dir(), c("phenotype.R", "association.R", "plotting.R")))))
    old <- if (file.exists(stamp)) tryCatch(readRDS(stamp), error=function(e) NULL) else NULL
    if (!dose_replace() && identical(old$identity, identity) && all(file.exists(outputs))) {
        current <- info(outputs)
        if (identical(old$output_info, current) || identical(old$md5, tools::md5sum(outputs))) {
            message("SKIPPED preprocess phase: ", name, " (verified)")
            return(invisible(NULL))
        }
    }
    message("Running preprocess phase: ", name)
    eval(expr, envir=new.env(parent=parent.frame()))
    invisible(gc())
    if (!identical(identity$inputs, info(inputs))) stop("Phase inputs changed during analysis: ", name)
    if (any(!file.exists(outputs))) stop("Phase did not produce all outputs: ", name)
    dir.create(dirname(stamp), recursive=TRUE, showWarnings=FALSE)
    dose_save_rds(list(identity=identity, output_info=info(outputs), md5=tools::md5sum(outputs)), stamp)
    message("Saved preprocess checkpoint: ", name)
    invisible(NULL)
}

# Cache small results after each expensive fit/block. Hash prepared inputs once,
# and keep fitted data/model frames out of the checkpoint itself.
dose_result_cache <- function(directory, identity, packages=character()) {
    scripts <- file.path(dose_script_dir(), c("0.common.R", "mechanism.R"))
    versions <- vapply(packages,
        function(p) as.character(utils::packageVersion(p)), character(1))
    token <- digest::digest(list(identity, R.version.string, versions,
                                tools::md5sum(scripts)), algo="sha256")
    rm(identity)
    directory <- file.path(directory, token)
    replace <- dose_replace()
    function(key, code) {
        path <- file.path(directory, paste0(digest::digest(key, algo="sha256"), ".rds"))
        old <- if (!replace && file.exists(path)) tryCatch(readRDS(path), error=function(e) NULL) else NULL
        if (is.list(old) && identical(old$key, key) && "value" %in% names(old)) {
            message("Reusing completed fit/block: ", key)
            return(old$value)
        }
        value <- eval(substitute(code), envir=parent.frame())
        dir.create(directory, recursive=TRUE, showWarnings=FALSE)
        dose_save_rds(list(key=key, value=value), path)
        value
    }
}

dose_add_family_history <- function(dat) {
    fields <- c("fhist_father", "fhist_mother", "fhist_siblings")
    missing <- setdiff(fields, names(dat))
    if (!length(missing)) return(dat)
    path <- file.path(dose_output_dir(), "ukb/00data/inputs/family_history.tsv")
    message("Extracting baseline family history using ", dose_definition_file("vip"))
    status <- system2(dose_python(), shQuote(file.path(dose_script_dir(), "0.common.py")),
                      env="DOSE_STEP=family_history")
    if (status != 0L) stop("Family-history extraction failed (exit ", status, ")")
    extra <- utils::read.delim(path, colClasses="character", check.names=FALSE,
                              quote="", na.strings=c("", "NA"))
    if (anyDuplicated(extra$eid)) stop("Duplicate IDs in family-history input")
    index <- match(as.character(dat$eid), extra$eid)
    if (anyNA(index)) stop("Family-history source does not cover the cohort IDs")
    for (field in missing) {
        value <- gsub('[\\[\\]" ]', "", extra[[field]][index], perl=TRUE)
        dat[[field]] <- gsub(",", "|", value, fixed=TRUE)
    }
    dat
}


# Module source registry: loading a module defines its steps without running analyses.
# Each selected step is executed in its own R process by dose.sh.
if (!exists(".dose_steps", inherits=FALSE)) .dose_steps <- new.env(parent=emptyenv())
dose_register_step <- function(module, step, code) {
    stopifnot(is.expression(code), length(code) == 1L)
    assign(paste(module, step, sep="/"), code, envir=.dose_steps)
    invisible(NULL)
}
dose_step_code <- function(module, step) {
    key <- paste(module, step, sep="/")
    if (!exists(key, envir=.dose_steps, inherits=FALSE)) {
        if (!module %in% c("foodwas","ml","score","mechanism","signature","nhanes"))
            stop("Unknown module: ", module)
        sys.source(file.path(dose_script_dir(), paste0(module, ".R")), envir=globalenv())
    }
    if (!exists(key, envir=.dose_steps, inherits=FALSE)) stop("Unknown step: ", key)
    block <- get(key, envir=.dose_steps, inherits=FALSE)[[1L]]
    as.expression(as.list(block)[-1L])
}
dose_source_step <- function(module, step, envir=parent.frame()) {
    eval(dose_step_code(module, step), envir=envir)
}
dose_run_module <- function(module) {
    step <- Sys.getenv("DOSE_STEP")
    if (!nzchar(step)) stop("Run ./dose.sh ", module, " (see ./dose.sh --help)")
    dose_source_step(module, step, envir=globalenv())
}

# Read one participant table at a time: retaining every input alongside the
# growing full join can exceed the preprocess scope's memory budget.
dose_merge_participant_rds <- function(paths) {
    stopifnot(length(paths) > 0L)
    out <- NULL
    for (path in paths) {
        message("Merging participant data: ", basename(path))
        next_table <- readRDS(path)
        if (!is.data.frame(next_table) || !"eid" %in% names(next_table))
            stop("Participant table lacks eid: ", path)
        next_table$eid <- as.character(next_table$eid)
        if (anyNA(next_table$eid) || anyDuplicated(next_table$eid))
            stop("Participant table has missing/duplicate eid: ", path)
        out <- if (is.null(out)) next_table else
            merge(out, next_table, by="eid", all=TRUE)
        rm(next_table)
        invisible(gc())
        message("  Merged: ", nrow(out), " participants, ", ncol(out), " columns")
    }
    out
}

# ===========================================================================
# Food association helpers
# ===========================================================================
# Shared helpers for DOSE s1 analyses

std_10_90 <- function(x){
	x <- as.numeric(x); q <- quantile(x, c(.1,.9), na.rm=TRUE, names=FALSE)
	if(!all(is.finite(q)) || q[1]==q[2]) return(rep(NA_real_, length(x)))
	(x-q[1])/(q[2]-q[1])
}

std_z <- function(x){
	x <- as.numeric(x); s <- sd(x, na.rm=TRUE)
	if(!is.finite(s) || s==0) return(rep(NA_real_, length(x)))
	(x-mean(x, na.rm=TRUE))/s
}

qcat <- function(x, n=4){ # MYTH-like quartiles; identical scores are not forcibly split
	x <- as.numeric(x); out <- rep(NA_integer_, length(x)); ok <- is.finite(x)
	if(sum(ok)<2) return(factor(out))
	br <- unique(quantile(x[ok], probs=seq(0,1,length.out=n+1), na.rm=TRUE, names=FALSE))
	if(length(br)<2) return(factor(out))
	out[ok] <- as.integer(cut(x[ok], breaks=br, include.lowest=TRUE, labels=FALSE))
	k <- max(out, na.rm=TRUE)
	factor(out, levels=seq_len(k), labels=paste0("Q",seq_len(k)), ordered=TRUE)
}

adaptive_food_cat <- function(x){
	x <- as.numeric(x); out <- rep(NA_integer_, length(x)); ok <- is.finite(x)
	if(!any(ok)) return(list(cat=factor(out), zero_pct=NA_real_, n_cat=0L))
	z <- mean(x[ok]==0)
	if(z>=.8){
		out[ok] <- ifelse(x[ok]==0,1L,2L)
	} else if(z<.2){
		br <- unique(quantile(x[ok], seq(0,1,length.out=6), na.rm=TRUE, names=FALSE))
		if(length(br)>=3) out[ok] <- as.integer(cut(x[ok], br, include.lowest=TRUE, labels=FALSE))
	} else {
		npos <- if(z<.4) 4L else if(z<.6) 3L else 2L
		pos <- ok & x>0; out[ok & x==0] <- 1L
		br <- unique(quantile(x[pos], seq(0,1,length.out=npos+1), na.rm=TRUE, names=FALSE))
		if(length(br)>=2) out[pos] <- 1L + as.integer(cut(x[pos], br, include.lowest=TRUE, labels=FALSE))
	}
	# robust fallback when quantile cut-points collapse
	if(length(unique(out[!is.na(out)]))<2 && any(ok & x==0) && any(ok & x>0)) out[ok] <- ifelse(x[ok]==0,1L,2L)
	lv <- sort(unique(out[!is.na(out)])); out <- match(out, lv)
	k <- length(lv); list(cat=factor(out, levels=seq_len(k), labels=paste0("C",seq_len(k)), ordered=TRUE), zero_pct=z, n_cat=k)
}


cox_cont <- function(dat, Y, X, covs){
	y <- paste0(Y,".Yt2e"); t <- paste0(Y,".t2e")
	f <- reformulate(c(X,covs), response=paste0("Surv(",t,",",y,")"))
	fit <- survival::coxph(f, data=dat)
	e <- broom::tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>% dplyr::filter(term==X)
	tibble::tibble(total=fit$n, case=fit$nevent, HR=e$estimate, LCI=e$conf.low, HCI=e$conf.high, P=e$p.value)
}

cox_cat <- function(dat, Y, Xcat, Xnum, covs, all_refs=FALSE){
	y <- paste0(Y,".Yt2e"); t <- paste0(Y,".t2e")
	d <- dat %>% dplyr::select(dplyr::all_of(c(t,y,Xcat,Xnum,covs))) %>% tidyr::drop_na()
	d$pheno <- droplevels(factor(d[[Xcat]], ordered=FALSE)); d$xnum <- as.numeric(d[[Xnum]])
	lv <- levels(d$pheno); if(length(lv)<2) return(NULL)
	f0 <- reformulate(covs, response=paste0("Surv(",t,",",y,")"))
	f1 <- reformulate(c("pheno",covs), response=paste0("Surv(",t,",",y,")"))
	fit0 <- survival::coxph(f0, data=d); fit1 <- survival::coxph(f1, data=d)
	lrt <- anova(fit0, fit1, test="LRT"); p.global <- as.numeric(lrt[2,ncol(lrt)])
	ph <- tryCatch(survival::cox.zph(fit1)$table["pheno","p"], error=function(e) NA_real_)
	med <- tapply(d$xnum, d$pheno, median, na.rm=TRUE); d$trend <- unname(med[as.character(d$pheno)])
	fit.tr <- survival::coxph(reformulate(c("trend",covs), response=paste0("Surv(",t,",",y,")")), data=d)
	p.trend <- broom::tidy(fit.tr) %>% dplyr::filter(term=="trend") %>% dplyr::pull(p.value)

	b <- c(setNames(0,lv[1]), setNames(stats::coef(fit1)[grep("^pheno",names(stats::coef(fit1)))], lv[-1]))
	optimal <- names(which.min(b))
	global <- tibble::tibble(total=fit1$n, case=fit1$nevent, n_cat=length(lv), P_global=p.global, P_trend=p.trend, PH_P=ph, optimal_cat=optimal)

	stats <- d %>% dplyr::group_by(pheno) %>% dplyr::summarise(n=dplyr::n(), case=sum(.data[[y]]==1), median=median(xnum), min=min(xnum), max=max(xnum), .groups="drop") %>% dplyr::rename(category=pheno)
	refs <- if(all_refs) lv else lv[1]
	contr <- lapply(refs, function(ref){
		d$pheno <- stats::relevel(d$pheno, ref=ref)
		fit <- survival::coxph(f1, data=d)
		td <- broom::tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>% dplyr::filter(grepl("^pheno",term))
		lev2 <- levels(d$pheno)
		res <- tibble::tibble(ref=ref, category=lev2, HR=1, LCI=1, HCI=1, P=1)
		if(nrow(td)){
			for(i in seq_len(nrow(td))){
				cat_i <- sub("^pheno","",td$term[i]); j <- match(cat_i,res$category)
				if(!is.na(j)) res[j,c("HR","LCI","HCI","P")] <- list(td$estimate[i],td$conf.low[i],td$conf.high[i],td$p.value[i])
			}
		}
		res %>% dplyr::left_join(stats, by="category")
	}) %>% dplyr::bind_rows()
	list(global=global, contrast=contr, stats=stats)
}

# ===========================================================================
# Mediation helpers
# ===========================================================================
# Shared functions for DOSE multi-omics mediation analyses.
# This file performs no data import and no omics-specific preprocessing.

# Derive the manuscript mediation rows from the upstream analysis, without Excel.
# b is the Cox mediator coefficient adjusted for DOSE, as in the product method.
s4_mediation_table <- function(associations, product, domain) {
    domain <- match.arg(domain, c("metabolite", "inflammation"))
    normalize_key <- function(x) {
        if (!domain %in% names(x) && "mediator" %in% names(x)) names(x)[names(x) == "mediator"] <- domain
        x
    }
    associations <- normalize_key(as.data.frame(associations))
    product <- normalize_key(as.data.frame(product))
    assoc_cols <- c(domain, "candidate", "beta_XM", "P_XM", "FDR_XM", "beta_MY", "P_MY", "FDR_MY")
    cox_cols <- c(domain, "a_estimate", "a_P", "b_adjusted_estimate", "b_adjusted_P",
                  "ACME", "P_ACME", "FDR_ACME", "Proportion_mediated", "P_proportion", "FDR_proportion")
    check_input <- function(x, required, label) {
        # A mediation run without candidates may export an empty note-only table.
        if (!nrow(x)) for (v in setdiff(required, names(x))) x[[v]] <- numeric()
        missing <- setdiff(required, names(x))
        if (length(missing)) stop(label, " lacks columns: ", paste(missing, collapse = ", "))
        if (anyNA(x[[domain]]) || any(!nzchar(x[[domain]])) || anyDuplicated(x[[domain]]))
            stop(label, " has missing/duplicate marker IDs")
        x
    }
    associations <- check_input(associations, assoc_cols, "Association results")
    product <- check_input(product, cox_cols, "Cox mediation results")
    keep <- is.finite(product$Proportion_mediated) & product$Proportion_mediated > 0 &
        is.finite(product$P_proportion) & product$P_proportion >= 0 & product$P_proportion < .05
    x <- product[keep, , drop = FALSE]
    ia <- match(x[[domain]], associations[[domain]])
    if (anyNA(ia)) stop("Selected mediator missing from association results")
    a <- associations[ia, , drop = FALSE]
    if (anyNA(a$candidate) || any(!a$candidate)) stop("Selected mediator is not an upstream candidate")
    coefficients <- c("a_estimate", "a_P", "b_adjusted_estimate", "b_adjusted_P")
    if (any(!vapply(x[coefficients], function(z) all(is.finite(z)), logical(1))))
        stop("Nonfinite coefficients in selected Cox mediation results")
    marker <- as.character(x[[domain]])
    display <- if (domain == "metabolite") sub("^met_", "", marker) else marker
    if (domain == "inflammation") {
        labels <- c(bc_NEUT = "Neutrophil count", bc_LYMPH = "Lymphocyte count",
                    bc_MONO = "Monocyte count", bc_PLT = "Platelet count", bb_CRP = "CRP")
        matched <- marker %in% names(labels)
        display[matched] <- unname(labels[marker[matched]])
    }
    out <- data.frame(marker = marker, display_name = display,
                      domain = rep(domain, length(marker)),
                      proportion_percent = 100 * x$Proportion_mediated,
                      P_proportion = x$P_proportion,
                      beta_a = x$a_estimate, P_a = x$a_P,
                      beta_b = x$b_adjusted_estimate, P_b = x$b_adjusted_P,
                      stringsAsFactors = FALSE)
    for (v in setdiff(assoc_cols, domain)) out[[v]] <- a[[v]]
    for (v in c("ACME", "P_ACME", "FDR_ACME", "Proportion_mediated", "FDR_proportion")) out[[v]] <- x[[v]]
    out$source_row <- which(keep) + 1L
    out <- out[order(-out$proportion_percent, out$marker), , drop = FALSE]
    rownames(out) <- NULL
    out
}

s4_read_mediation_table <- function(dir_out, domain) {
    paths <- c(associations = file.path(dir_out, paste0("01_all_", domain, "_associations.csv")),
               mediation = file.path(dir_out, "07_mediation_Cox_product.csv"))
    if (any(!file.exists(paths))) stop("Run the upstream ", domain, " analysis first. Missing: ",
                                      paste(paths[!file.exists(paths)], collapse = "; "))
    out <- s4_mediation_table(read.csv(paths[[1]], check.names = FALSE),
                              read.csv(paths[[2]], check.names = FALSE), domain)
    out$source_file <- rep(paths[[2]], nrow(out))
    out
}

s4_required_packages <- c(
    "dplyr", "purrr", "survival", "broom", "ggplot2", "ggrepel",
    "openxlsx", "mediation"
)

s4_load_packages <- function(extra = character()) {
    pkgs <- unique(c(s4_required_packages, extra))
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
        stop("Install required packages first: ", paste(missing, collapse = ", "))
    }
    suppressPackageStartupMessages({
        library(dplyr)
        library(purrr)
        library(survival)
        library(broom)
        library(ggplot2)
    })
    invisible(pkgs)
}

s4_safe_scale <- function(x, label = deparse(substitute(x))) {
    x <- suppressWarnings(as.numeric(x))
    if (any(!is.finite(x))) stop("Non-finite values remain in ", label, ".")
    sx <- stats::sd(x)
    if (!is.finite(sx) || sx <= 0) stop("Cannot standardize ", label, ".")
    as.numeric((x - mean(x)) / sx)
}

s4_transform_direct_z <- function(x, label) {
    s4_safe_scale(x, label)
}

s4_transform_log1p_z <- function(x, label) {
    x <- suppressWarnings(as.numeric(x))
    if (any(!is.finite(x))) stop("Non-finite values remain in ", label, ".")
    if (any(x < 0)) stop("Negative values found in ", label,
                         "; log(x+1) is not valid.")
    s4_safe_scale(log1p(x), paste0("log1p(", label, ")"))
}

s4_valid_covariate <- function(x) {
    if (is.factor(x) || is.character(x) || is.logical(x)) {
        return(dplyr::n_distinct(x, na.rm = TRUE) > 1L)
    }
    sx <- suppressWarnings(stats::sd(as.numeric(x), na.rm = TRUE))
    is.finite(sx) && sx > 0
}

s4_formula <- function(response, terms) {
    stats::reformulate(termlabels = terms, response = response)
}

s4_extract_term <- function(fit, term, exponentiate = FALSE) {
    z <- broom::tidy(fit, conf.int = TRUE, exponentiate = exponentiate)
    z <- z[z$term == term, , drop = FALSE]
    if (nrow(z) != 1L) stop("Could not uniquely extract term: ", term)
    z
}

s4_med_component <- function(object, name, default = NA_real_) {
    x <- object[[name]]
    if (is.null(x) || !length(x)) return(default)
    as.numeric(x[1])
}

s4_med_ci <- function(object, name) {
    x <- object[[name]]
    if (is.null(x) || length(x) < 2L) return(c(NA_real_, NA_real_))
    as.numeric(x[1:2])
}

s4_write_csv <- function(x, path) {
    if (!ncol(x)) x <- data.frame(note = "No results")
    utils::write.csv(x, path, row.names = FALSE, na = "")
}

s4_sheet_safe <- function(x) {
    if (!ncol(x)) data.frame(note = "No results") else x
}

s4_significance_row <- function(result, p, fdr) {
    data.frame(
        result = result,
        N_tested = sum(!is.na(p)),
        N_P05 = sum(p < 0.05, na.rm = TRUE),
        N_FDR05 = sum(fdr < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}

s4_significance_from_data <- function(result, data, p_col, fdr_col) {
    p <- if (p_col %in% names(data)) data[[p_col]] else numeric()
    fdr <- if (fdr_col %in% names(data)) data[[fdr_col]] else numeric()
    s4_significance_row(result, p, fdr)
}

s4_validate_ids <- function(x, name) {
    if (!"eid" %in% names(x)) stop(name, " has no eid column.")
    x$eid <- as.character(x$eid)
    if (anyNA(x$eid) || anyDuplicated(x$eid)) {
        stop(name, " has missing or duplicated eid values.")
    }
    x
}

s4_validate_dose <- function(dose_object) {
    component_cols <- attr(dose_object, "component_columns", exact = TRUE)
    meta <- list(
        version = attr(dose_object, "score_version", exact = TRUE),
        method = attr(dose_object, "primary_scoring_method", exact = TRUE),
        unit = attr(dose_object, "score_intake_unit", exact = TRUE),
        range = attr(dose_object, "score_theoretical_range", exact = TRUE),
        step = attr(dose_object, "score_step", exact = TRUE),
        component_steps = attr(dose_object, "component_steps", exact = TRUE)
    )
    dose <- s4_validate_ids(as.data.frame(dose_object), "dose")
    valid_meta <- identical(meta$version, "dose") &&
        identical(meta$unit, "servings/day") &&
        length(meta$method) == 1L && !is.na(meta$method) &&
        length(meta$range) == 2L && all(is.finite(meta$range)) &&
        diff(meta$range) > 0 && length(meta$step) == 1L &&
        is.finite(meta$step) && meta$step > 0 && length(component_cols) > 0L &&
        !anyDuplicated(component_cols) &&
        all(c("dose", component_cols) %in% names(dose)) &&
        identical(names(meta$component_steps), component_cols) &&
        all(is.finite(meta$component_steps) & meta$component_steps > 0)
    if (!valid_meta) {
        stop("Invalid primary DOSE metadata. Run the current score.R first.")
    }
    valid_components <- vapply(component_cols, function(v) {
        is.numeric(dose[[v]]) && all(
            is.na(dose[[v]]) |
                dose[[v]] %in% seq(0, 1, by = meta$component_steps[[v]])
        )
    }, logical(1))
    valid_dose <- is.numeric(dose$dose) && all(
        is.na(dose$dose) |
            dose$dose %in% seq(meta$range[1], meta$range[2], by = meta$step)
    ) && isTRUE(all.equal(
        dose$dose, rowSums(dose[component_cols]), check.attributes = FALSE
    ))
    if (!all(valid_components) || !valid_dose) {
        stop("DOSE must match its component sum, range and step metadata.")
    }
    if ("diet.dose.sum" %in% names(dose) &&
        !identical(dose$dose, dose$diet.dose.sum)) {
        stop("The compatibility alias diet.dose.sum differs from dose.")
    }
    list(data = dose, meta = meta, component_cols = component_cols)
}

s4_get_m2_covariates <- function(cfg) {
    alcohol_covariates <- unique(c(cfg$covs.alcohol, "alcohol.c"))
    covs <- setdiff(unique(c(cfg$covs.m1, cfg$covs.m2.add)),
                    alcohol_covariates)
    if (!length(covs) || anyNA(covs) || any(!nzchar(covs)) ||
        any(alcohol_covariates %in% covs)) {
        stop("Invalid M2 covariates: separate alcohol adjustment is not permitted.")
    }
    covs
}

s4_prepare_base <- function(dat, dose, mediator_data, mediator_vars,
                            time_var, event_var, covariates, cohort_label) {
    dat <- s4_validate_ids(as.data.frame(dat), "dat")
    dose <- s4_validate_ids(as.data.frame(dose), "dose")
    mediator_data <- s4_validate_ids(as.data.frame(mediator_data), cohort_label)
    if (!setequal(dat$eid, dose$eid)) stop("Score/cohort ID mismatch.")
    dat$dose <- dose$dose[match(dat$eid, dose$eid)]
    if (anyNA(dat$dose) || any(!is.finite(dat$dose))) {
        stop("The primary dose score contains missing or non-finite values.")
    }
    required <- unique(c("eid", time_var, event_var, covariates))
    if (!all(required %in% names(dat))) {
        stop("Missing cohort variables: ",
             paste(setdiff(required, names(dat)), collapse = ", "))
    }
    if (!all(mediator_vars %in% names(mediator_data))) {
        stop("Missing mediator variables: ",
             paste(setdiff(mediator_vars, names(mediator_data)), collapse = ", "))
    }
    base_vars <- unique(c("eid", "dose", time_var, event_var, covariates))
    dat <- dat[, base_vars, drop = FALSE]
    datp <- dplyr::inner_join(
        dat, mediator_data[c("eid", mediator_vars)], by = "eid"
    )
    ok <- complete.cases(datp[base_vars]) & is.finite(datp[[time_var]]) &
        datp[[time_var]] > 0 & datp[[event_var]] %in% 0:1
    for (v in covariates) {
        if (is.numeric(datp[[v]]) || is.integer(datp[[v]])) {
            ok <- ok & is.finite(datp[[v]])
        }
    }
    if (!any(ok)) stop("No eligible ", cohort_label, " M2 participants.")
    datp0 <- droplevels(datp[ok, , drop = FALSE])
    datp0$dose_z <- s4_safe_scale(datp0$dose, "dose")
    bad_covs <- covariates[
        !vapply(datp0[covariates], s4_valid_covariate, logical(1))
    ]
    if (length(bad_covs)) {
        stop("M2 covariates without variation in the ", cohort_label,
             " cohort: ", paste(bad_covs, collapse = ", "))
    }
    list(data = datp0, overlap_N = nrow(datp), covariates = covariates)
}

s4_empty_result <- function() data.frame()

s4_capture_error <- function(stage, mediator, e) {
    data.frame(stage = stage, mediator = mediator,
               error = conditionMessage(e), stringsAsFactors = FALSE)
}

s4_run_xm <- function(dat, mediator_vars, transform_mediator, covariates,
                      min_n, fdr_n, label) {
    errors <- list()
    out <- purrr::map(seq_along(mediator_vars), function(i) {
        m <- mediator_vars[i]
        if (i == 1L || i %% 100L == 0L || i == length(mediator_vars)) {
            message(sprintf("Flow 1, DOSE -> %s: %d/%d | %s",
                            label, i, length(mediator_vars), m))
        }
        tryCatch({
            vars <- unique(c("eid", "dose_z", covariates, m))
            d <- droplevels(dat[complete.cases(dat[vars]), vars, drop = FALSE])
            if (nrow(d) < min_n) stop("Insufficient N")
            d$mediator_z <- transform_mediator(d[[m]], m)
            fit <- stats::lm(s4_formula("mediator_z", c("dose_z", covariates)),
                             data = d)
            a <- s4_extract_term(fit, "dose_z")
            data.frame(mediator = m, N_XM = nrow(d), beta_XM = a$estimate,
                       SE_XM = a$std.error, LCI_XM = a$conf.low,
                       HCI_XM = a$conf.high, P_XM = a$p.value)
        }, error = function(e) {
            errors[[length(errors) + 1L]] <<- s4_capture_error("flow1_XM", m, e)
            NULL
        })
    }) %>% dplyr::bind_rows()
    if (!nrow(out)) stop("All flow-1 models failed.")
    out <- out %>%
        dplyr::mutate(FDR_XM = p.adjust(P_XM, method = "BH", n = fdr_n)) %>%
        dplyr::arrange(FDR_XM, P_XM)
    list(results = out, errors = dplyr::bind_rows(errors))
}

s4_run_my <- function(dat, mediator_vars, transform_mediator, covariates,
                      time_var, event_var, min_n, min_events, fdr_n, label) {
    errors <- list()
    out <- purrr::map(seq_along(mediator_vars), function(i) {
        m <- mediator_vars[i]
        if (i == 1L || i %% 100L == 0L || i == length(mediator_vars)) {
            message(sprintf("Flow 2, %s -> IS Cox: %d/%d | %s",
                            label, i, length(mediator_vars), m))
        }
        tryCatch({
            vars <- unique(c("eid", time_var, event_var, covariates, m))
            d <- droplevels(dat[complete.cases(dat[vars]), vars, drop = FALSE])
            if (nrow(d) < min_n || sum(d[[event_var]]) < min_events) {
                stop("Insufficient N/events")
            }
            d$mediator_z <- transform_mediator(d[[m]], m)
            d$time <- d[[time_var]]
            d$event <- as.integer(d[[event_var]])
            # Screening model deliberately excludes DOSE.
            fit <- survival::coxph(
                s4_formula("Surv(time, event)", c("mediator_z", covariates)),
                data = d, ties = "efron", x = TRUE
            )
            b <- s4_extract_term(fit, "mediator_z")
            hr <- s4_extract_term(fit, "mediator_z", exponentiate = TRUE)
            data.frame(mediator = m, N_MY = nrow(d),
                       Events_MY = sum(d$event), beta_MY = b$estimate,
                       SE_MY = b$std.error, HR_MY = hr$estimate,
                       LCI_HR_MY = hr$conf.low, HCI_HR_MY = hr$conf.high,
                       P_MY = b$p.value)
        }, error = function(e) {
            errors[[length(errors) + 1L]] <<- s4_capture_error("flow2_MY_Cox", m, e)
            NULL
        })
    }) %>% dplyr::bind_rows()
    if (!nrow(out)) stop("All flow-2 models failed.")
    out <- out %>%
        dplyr::mutate(FDR_MY = p.adjust(P_MY, method = "BH", n = fdr_n)) %>%
        dplyr::arrange(FDR_MY, P_MY)
    list(results = out, errors = dplyr::bind_rows(errors))
}

s4_mediation_one <- function(dat, mediator, transform_mediator, covariates,
                             time_var, event_var, method,
                             min_n, min_events, sims, seed) {
    vars <- unique(c("dose_z", time_var, event_var, covariates, mediator))
    d <- droplevels(dat[complete.cases(dat[vars]), vars, drop = FALSE])
    if (nrow(d) < min_n || sum(d[[event_var]]) < min_events) {
        stop("Insufficient N/events")
    }
    d$mediator_z <- transform_mediator(d[[mediator]], mediator)
    d$time <- d[[time_var]]
    d$event <- as.integer(d[[event_var]])
    fit.X2M <- stats::lm(
        s4_formula("mediator_z", c("dose_z", covariates)),
        data = d, x = TRUE, y = TRUE, model = TRUE
    )
    a <- s4_extract_term(fit.X2M, "dose_z")

    if (method == "lhr_aft") {
        # No dist argument: this matches the referenced LHR code.
        fit.X2Y <- survival::survreg(
            s4_formula("Surv(time, event)", c("dose_z", covariates)),
            data = d, x = TRUE, y = TRUE, model = TRUE
        )
        fit.M2Y <- survival::survreg(
            s4_formula("Surv(time, event)",
                       c("dose_z", "mediator_z", covariates)),
            data = d, x = TRUE, y = TRUE, model = TRUE
        )
        b_reported <- b_adjusted <- s4_extract_term(fit.M2Y, "mediator_z")
        effect_scale <- "survreg AFT model-based effect (default distribution)"
    } else if (method == "modern_logistic") {
        fit_b <- stats::glm(
            s4_formula("event", c("mediator_z", covariates)), data = d,
            family = stats::binomial(link = "logit"),
            x = TRUE, y = TRUE, model = TRUE
        )
        fit.X2Y <- stats::glm(
            s4_formula("event", c("dose_z", covariates)), data = d,
            family = stats::binomial(link = "logit"),
            x = TRUE, y = TRUE, model = TRUE
        )
        fit.M2Y <- stats::glm(
            s4_formula("event", c("dose_z", "mediator_z", covariates)),
            data = d, family = stats::binomial(link = "logit"),
            x = TRUE, y = TRUE, model = TRUE
        )
        b_reported <- s4_extract_term(fit_b, "mediator_z")
        b_adjusted <- s4_extract_term(fit.M2Y, "mediator_z")
        effect_scale <- "logistic model-based effect"
    } else if (method == "cox_product") {
        fit_b <- survival::coxph(
            s4_formula("Surv(time, event)", c("mediator_z", covariates)),
            data = d, ties = "efron", x = TRUE
        )
        fit.X2Y <- survival::coxph(
            s4_formula("Surv(time, event)", c("dose_z", covariates)),
            data = d, ties = "efron", x = TRUE
        )
        fit.M2Y <- survival::coxph(
            s4_formula("Surv(time, event)",
                       c("dose_z", "mediator_z", covariates)),
            data = d, ties = "efron", x = TRUE
        )
        b_reported <- s4_extract_term(fit_b, "mediator_z")
        b_adjusted <- s4_extract_term(fit.M2Y, "mediator_z")
    } else stop("Unknown mediation method: ", method)

    c_total <- s4_extract_term(fit.X2Y, "dose_z")
    c_direct <- s4_extract_term(fit.M2Y, "dose_z")

    if (method == "cox_product") {
        indirect <- a$estimate * b_adjusted$estimate
        se_indirect <- sqrt(
            b_adjusted$estimate^2 * a$std.error^2 +
                a$estimate^2 * b_adjusted$std.error^2
        )
        z <- indirect / se_indirect

        # Approximate proportion mediated on the Cox log-hazard scale:
        #     proportion = (a * b) / c
        # where c is the total DOSE coefficient from the Cox X -> Y model.
        # The first-order delta variance below includes uncertainty in a, b,
        # and c, but assumes their cross-model covariances are zero.  This is
        # an exploratory approximation; it can be unstable when c is near 0.
        c_est <- c_total$estimate
        c_se <- c_total$std.error
        proportion <- proportion_se <- proportion_z <- proportion_p <- NA_real_
        proportion_lci <- proportion_hci <- NA_real_
        denominator_ok <- is.finite(c_est) &&
            abs(c_est) > sqrt(.Machine$double.eps)
        inputs_ok <- all(is.finite(c(
            a$estimate, a$std.error,
            b_adjusted$estimate, b_adjusted$std.error,
            c_est, c_se
        )))
        if (denominator_ok && inputs_ok) {
            proportion <- indirect / c_est
            proportion_var <-
                (b_adjusted$estimate / c_est)^2 * a$std.error^2 +
                (a$estimate / c_est)^2 * b_adjusted$std.error^2 +
                (indirect / c_est^2)^2 * c_se^2
            if (is.finite(proportion_var) && proportion_var > 0) {
                proportion_se <- sqrt(proportion_var)
                proportion_z <- proportion / proportion_se
                proportion_p <- 2 * stats::pnorm(-abs(proportion_z))
                proportion_lci <- proportion -
                    stats::qnorm(0.975) * proportion_se
                proportion_hci <- proportion +
                    stats::qnorm(0.975) * proportion_se
            }
        }
        return(data.frame(
            mediator = mediator, method = method,
            effect_scale = "approximate log-hazard",
            effect_contrast = "DOSE +1 SD (dose_z 0 to 1)",
            N = nrow(d), Events = sum(d$event),
            a_estimate = a$estimate, a_SE = a$std.error, a_P = a$p.value,
            b_reported_estimate = b_reported$estimate,
            b_reported_P = b_reported$p.value,
            b_adjusted_estimate = b_adjusted$estimate,
            b_adjusted_P = b_adjusted$p.value,
            c_total_model_estimate = c_total$estimate,
            c_total_model_P = c_total$p.value,
            c_direct_model_estimate = c_direct$estimate,
            c_direct_model_P = c_direct$p.value,
            ACME = indirect,
            ACME_LCI = indirect - qnorm(0.975) * se_indirect,
            ACME_HCI = indirect + qnorm(0.975) * se_indirect,
            P_ACME = 2 * pnorm(-abs(z)),
            ADE = c_direct$estimate, ADE_LCI = c_direct$conf.low,
            ADE_HCI = c_direct$conf.high, P_ADE = c_direct$p.value,
            Total_effect = c_total$estimate, Total_LCI = c_total$conf.low,
            Total_HCI = c_total$conf.high, P_total = c_total$p.value,
            Proportion_mediated = proportion,
            Proportion_SE = proportion_se,
            Proportion_z = proportion_z,
            Proportion_LCI = proportion_lci,
            Proportion_HCI = proportion_hci,
            P_proportion = proportion_p,
            proportion_inference = paste0(
                "first-order delta method for (a*b)/c; ",
                "cross-model covariances assumed zero"
            ),
            simulations = 0L,
            bootstrap = FALSE, seed = seed
        ))
    }

    set.seed(seed)
    med <- mediation::mediate(
        model.m = fit.X2M, model.y = fit.M2Y,
        treat = "dose_z", mediator = "mediator_z",
        control.value = 0, treat.value = 1,
        sims = sims, boot = FALSE, long = FALSE
    )
    acme_ci <- s4_med_ci(med, "d.avg.ci")
    ade_ci <- s4_med_ci(med, "z.avg.ci")
    total_ci <- s4_med_ci(med, "tau.ci")
    prop_ci <- s4_med_ci(med, "n.avg.ci")
    data.frame(
        mediator = mediator, method = method, effect_scale = effect_scale,
        effect_contrast = "DOSE +1 SD (dose_z 0 to 1)",
        N = nrow(d), Events = sum(d$event),
        a_estimate = a$estimate, a_SE = a$std.error, a_P = a$p.value,
        b_reported_estimate = b_reported$estimate,
        b_reported_P = b_reported$p.value,
        b_adjusted_estimate = b_adjusted$estimate,
        b_adjusted_P = b_adjusted$p.value,
        c_total_model_estimate = c_total$estimate,
        c_total_model_P = c_total$p.value,
        c_direct_model_estimate = c_direct$estimate,
        c_direct_model_P = c_direct$p.value,
        ACME = s4_med_component(med, "d.avg"),
        ACME_LCI = acme_ci[1], ACME_HCI = acme_ci[2],
        P_ACME = s4_med_component(med, "d.avg.p"),
        ADE = s4_med_component(med, "z.avg"),
        ADE_LCI = ade_ci[1], ADE_HCI = ade_ci[2],
        P_ADE = s4_med_component(med, "z.avg.p"),
        Total_effect = s4_med_component(med, "tau.coef"),
        Total_LCI = total_ci[1], Total_HCI = total_ci[2],
        P_total = s4_med_component(med, "tau.p"),
        Proportion_mediated = s4_med_component(med, "n.avg"),
        Proportion_LCI = prop_ci[1], Proportion_HCI = prop_ci[2],
        P_proportion = s4_med_component(med, "n.avg.p"),
        simulations = sims, bootstrap = FALSE, seed = seed
    )
}

# Three named wrappers are intentionally retained so each analysis method is
# explicit and independently callable.
s4_mediation_lhr_aft <- function(...) {
    s4_mediation_one(..., method = "lhr_aft")
}

s4_mediation_modern_logistic <- function(...) {
    s4_mediation_one(..., method = "modern_logistic")
}

s4_mediation_cox_product <- function(...) {
    s4_mediation_one(..., method = "cox_product")
}

s4_run_mediation_methods <- function(dat, candidates, transform_mediator,
                                     covariates, time_var, event_var,
                                     min_n, min_events, sims, seed,
                                     cache = function(key, code) force(code)) {
    methods <- list(
        lhr_aft = s4_mediation_lhr_aft,
        modern_logistic = s4_mediation_modern_logistic,
        cox_product = s4_mediation_cox_product
    )
    errors <- list()
    results <- lapply(names(methods), function(method) {
        message("Starting mediation method: ", method,
                "; candidates=", length(candidates))
        if (!length(candidates)) return(data.frame())
        out <- purrr::map(seq_along(candidates), function(i) {
            m <- candidates[i]
            if (i == 1L || i %% 10L == 0L || i == length(candidates)) {
                message(sprintf("%s: %d/%d | %s",
                                method, i, length(candidates), m))
            }
            tryCatch(
                cache(paste(method, m, sep="/"), methods[[method]](
                    dat = dat, mediator = m,
                    transform_mediator = transform_mediator,
                    covariates = covariates, time_var = time_var,
                    event_var = event_var, min_n = min_n,
                    min_events = min_events, sims = sims, seed = seed
                )),
                error = function(e) {
                    errors[[length(errors) + 1L]] <<-
                        s4_capture_error(paste0("mediation_", method), m, e)
                    NULL
                }
            )
        }) %>% dplyr::bind_rows()
        if (nrow(out)) {
            n_fdr <- length(candidates)
            out <- out %>% dplyr::mutate(
                FDR_total = p.adjust(P_total, method = "BH", n = n_fdr),
                FDR_ACME = p.adjust(P_ACME, method = "BH", n = n_fdr),
                FDR_ADE = p.adjust(P_ADE, method = "BH", n = n_fdr),
                FDR_proportion = if (all(is.na(P_proportion))) NA_real_ else
                    p.adjust(P_proportion, method = "BH", n = n_fdr)
            ) %>% dplyr::arrange(P_ACME)
        }
        method_counts <- dplyr::bind_rows(
            s4_significance_from_data(
                paste0(method, ": ACME/indirect effect"),
                out, "P_ACME", "FDR_ACME"
            ),
            s4_significance_from_data(
                paste0(method, ": proportion mediated"),
                out, "P_proportion", "FDR_proportion"
            )
        )
        message("Flow 3 method completed: ", method)
        print(method_counts, row.names = FALSE)
        out
    })
    names(results) <- names(methods)
    list(results = results, errors = dplyr::bind_rows(errors))
}

s4_compact_method <- function(x, prefix) {
    if (!nrow(x)) {
        z <- data.frame(
            mediator = character(), ACME = numeric(), P_ACME = numeric(),
            FDR_ACME = numeric(), Proportion_mediated = numeric(),
            P_proportion = numeric(), FDR_proportion = numeric()
        )
    } else {
        z <- x %>% dplyr::select(
            mediator, ACME, P_ACME, FDR_ACME,
            Proportion_mediated, P_proportion, FDR_proportion
        )
    }
    names(z)[-1] <- paste0(prefix, "_", names(z)[-1])
    z
}

s4_method_to_long <- function(x, n_fdr) {
    if (!nrow(x)) return(data.frame())
    dplyr::bind_rows(
        x %>% dplyr::transmute(method, mediator, effect = "TE",
            estimate = Total_effect, LCI = Total_LCI,
            HCI = Total_HCI, p.value = P_total),
        x %>% dplyr::transmute(method, mediator, effect = "ACME",
            estimate = ACME, LCI = ACME_LCI,
            HCI = ACME_HCI, p.value = P_ACME),
        x %>% dplyr::transmute(method, mediator, effect = "ADE",
            estimate = ADE, LCI = ADE_LCI,
            HCI = ADE_HCI, p.value = P_ADE),
        x %>% dplyr::transmute(method, mediator, effect = "prop",
            estimate = Proportion_mediated, LCI = Proportion_LCI,
            HCI = Proportion_HCI, p.value = P_proportion)
    ) %>%
        dplyr::group_by(method, effect) %>%
        dplyr::mutate(FDR = if (all(is.na(p.value))) NA_real_ else
            p.adjust(p.value, method = "BH", n = n_fdr)) %>%
        dplyr::ungroup()
}

s4_run_analysis <- function(dat, mediator_vars, transform_mediator,
                            mediator_label, time_var, event_var, covariates,
                            min_n = 1000L, min_events = 20L,
                            xm_fdr = 0.05, my_p = 0.05, my_fdr = 0.05,
                            sims = 100L, seed = 12345L,
                            fdr_family_n = length(mediator_vars)) {
    stopifnot(is.function(transform_mediator), length(mediator_vars) > 0L)
    fdr_family_n <- max(as.integer(fdr_family_n), length(mediator_vars))
    cache <- dose_result_cache(file.path(dose_output_dir(), ".checkpoints/mechanism", mediator_label, "fits"),
        list(dat, mediator_vars, body(transform_mediator), formals(transform_mediator),
             mediator_label, time_var, event_var, covariates, min_n, min_events,
             xm_fdr, my_p, my_fdr, sims, seed, fdr_family_n),
        packages=c("survival", "mediation", "dplyr", "broom"))

    total_data <- dat %>% dplyr::transmute(
        time = .data[[time_var]], event = .data[[event_var]],
        dose_z = dose_z, dplyr::across(dplyr::all_of(covariates))
    )
    total_fit <- survival::coxph(
        s4_formula("Surv(time, event)", c("dose_z", covariates)),
        data = total_data, ties = "efron"
    )
    tt <- s4_extract_term(total_fit, "dose_z", exponentiate = TRUE)
    total_effect <- data.frame(
        N = nrow(total_data), Events = sum(total_data$event),
        HR_DOSE_per_SD = tt$estimate, LCI = tt$conf.low,
        HCI = tt$conf.high, P = tt$p.value
    )

    xm <- cache("xm", s4_run_xm(dat, mediator_vars, transform_mediator, covariates,
                    min_n, fdr_family_n, mediator_label))
    message("Flow 1 completed:")
    flow1_counts <- s4_significance_row(
        paste0("Flow 1: DOSE -> ", mediator_label, " linear association"),
        xm$results$P_XM, xm$results$FDR_XM
    )
    print(flow1_counts, row.names = FALSE)

    my <- cache("my", s4_run_my(dat, mediator_vars, transform_mediator, covariates,
                    time_var, event_var, min_n, min_events,
                    fdr_family_n, mediator_label))
    message("Flow 2 completed:")
    flow2_counts <- s4_significance_row(
        paste0("Flow 2: ", mediator_label,
               " -> IS Cox association without DOSE"),
        my$results$P_MY, my$results$FDR_MY
    )
    print(flow2_counts, row.names = FALSE)

    candidates <- dplyr::inner_join(
        xm$results %>% dplyr::filter(FDR_XM < xm_fdr),
        my$results %>% dplyr::filter(P_MY < my_p),
        by = "mediator"
    ) %>% dplyr::arrange(FDR_XM, P_MY)
    associations <- dplyr::full_join(
        xm$results, my$results, by = "mediator"
    ) %>% dplyr::mutate(candidate = mediator %in% candidates$mediator) %>%
        dplyr::arrange(FDR_XM, P_MY)
    message("Intersection retained for pathway analysis and mediation: ",
            nrow(candidates), " ", mediator_label,
            "(s) (X-M FDR<0.05 and M-Y P<0.05).")

    med <- s4_run_mediation_methods(
        dat, candidates$mediator, transform_mediator, covariates,
        time_var, event_var, min_n, min_events, sims, seed, cache=cache
    )
    lhr <- med$results$lhr_aft
    modern <- med$results$modern_logistic
    product <- med$results$cox_product
    combined <- candidates %>%
        dplyr::left_join(s4_compact_method(lhr, "LHR"), by = "mediator") %>%
        dplyr::left_join(s4_compact_method(modern, "MODERN"), by = "mediator") %>%
        dplyr::left_join(s4_compact_method(product, "PRODUCT"), by = "mediator")
    long <- dplyr::bind_rows(
        s4_method_to_long(lhr, nrow(candidates)),
        s4_method_to_long(modern, nrow(candidates)),
        s4_method_to_long(product, nrow(candidates))
    )
    if (!ncol(long)) long <- data.frame(
        method = character(), mediator = character(), effect = character(),
        estimate = numeric(), LCI = numeric(), HCI = numeric(),
        p.value = numeric(), FDR = numeric()
    )
    positive_prop <- long %>% dplyr::filter(
        method %in% c("lhr_aft", "modern_logistic"), effect == "prop",
        p.value < 0.05, estimate > 0
    ) %>% dplyr::arrange(method, p.value)
    product_positive <- if (nrow(product)) product %>% dplyr::filter(
        P_ACME < 0.05, is.finite(Proportion_mediated),
        Proportion_mediated > 0
    ) %>% dplyr::arrange(P_ACME) else data.frame(mediator = character())

    association_counts <- dplyr::bind_rows(flow1_counts, flow2_counts)
    mediation_counts <- dplyr::bind_rows(
        s4_significance_from_data("LHR AFT: total effect", lhr,
                                  "P_total", "FDR_total"),
        s4_significance_from_data("LHR AFT: ACME", lhr,
                                  "P_ACME", "FDR_ACME"),
        s4_significance_from_data("LHR AFT: ADE", lhr,
                                  "P_ADE", "FDR_ADE"),
        s4_significance_from_data("LHR AFT: proportion mediated", lhr,
                                  "P_proportion", "FDR_proportion"),
        s4_significance_from_data("MODERN logistic: total effect", modern,
                                  "P_total", "FDR_total"),
        s4_significance_from_data("MODERN logistic: ACME", modern,
                                  "P_ACME", "FDR_ACME"),
        s4_significance_from_data("MODERN logistic: ADE", modern,
                                  "P_ADE", "FDR_ADE"),
        s4_significance_from_data("MODERN logistic: proportion mediated", modern,
                                  "P_proportion", "FDR_proportion"),
        s4_significance_from_data("Cox product: total effect", product,
                                  "P_total", "FDR_total"),
        s4_significance_from_data("Cox product: indirect effect", product,
                                  "P_ACME", "FDR_ACME"),
        s4_significance_from_data("Cox product: direct effect", product,
                                  "P_ADE", "FDR_ADE"),
        s4_significance_from_data("Cox product: proportion mediated", product,
                                  "P_proportion", "FDR_proportion")
    )
    message("Mediation-stage significance counts:")
    print(mediation_counts, row.names = FALSE)
    summary <- data.frame(
        N_base_complete_M2 = nrow(dat), IS_events = sum(dat[[event_var]]),
        HR_DOSE_IS_overall = total_effect$HR_DOSE_per_SD,
        P_DOSE_IS_overall = total_effect$P,
        N_mediators_tested = nrow(associations),
        N_XM_P05 = sum(associations$P_XM < 0.05, na.rm = TRUE),
        N_XM_FDR05 = sum(associations$FDR_XM < 0.05, na.rm = TRUE),
        N_MY_P05 = sum(associations$P_MY < 0.05, na.rm = TRUE),
        N_MY_FDR05 = sum(associations$FDR_MY < my_fdr, na.rm = TRUE),
        N_candidates_XM_FDR05_MY_P05 = nrow(candidates),
        N_LHR_completed = nrow(lhr),
        N_LHR_ACME_P05 = sum(lhr$P_ACME < 0.05, na.rm = TRUE),
        N_LHR_ACME_FDR05 = sum(lhr$FDR_ACME < 0.05, na.rm = TRUE),
        N_LHR_positive_prop_P05 = sum(
            lhr$P_proportion < 0.05 & lhr$Proportion_mediated > 0,
            na.rm = TRUE),
        N_MODERN_completed = nrow(modern),
        N_MODERN_ACME_P05 = sum(modern$P_ACME < 0.05, na.rm = TRUE),
        N_MODERN_ACME_FDR05 = sum(modern$FDR_ACME < 0.05, na.rm = TRUE),
        N_MODERN_positive_prop_P05 = sum(
            modern$P_proportion < 0.05 & modern$Proportion_mediated > 0,
            na.rm = TRUE),
        N_PRODUCT_completed = nrow(product),
        N_PRODUCT_indirect_P05 = sum(product$P_ACME < 0.05, na.rm = TRUE),
        N_PRODUCT_indirect_FDR05 = sum(product$FDR_ACME < 0.05, na.rm = TRUE),
        N_PRODUCT_proportion_P05 = sum(
            product$P_proportion < 0.05, na.rm = TRUE),
        N_PRODUCT_proportion_FDR05 = sum(
            product$FDR_proportion < 0.05, na.rm = TRUE),
        N_PRODUCT_positive_proportion_P05 = sum(
            product$P_proportion < 0.05 &
                product$Proportion_mediated > 0,
            na.rm = TRUE)
    )
    list(
        total_effect = total_effect, xm = xm$results, my = my$results,
        candidates = candidates, associations = associations,
        lhr = lhr, modern = modern, product = product,
        combined = combined, long = long,
        positive_prop = positive_prop, product_positive = product_positive,
        significance_counts = dplyr::bind_rows(association_counts,
                                                mediation_counts),
        summary = summary,
        errors = dplyr::bind_rows(xm$errors, my$errors, med$errors)
    )
}

s4_export_results <- function(res, dir_out, mediator_label,
                              workbook_name, candidate_prefix = mediator_label,
                              extra_sheets = list(), notes = character()) {
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    p <- function(name) file.path(dir_out, name)
    label <- mediator_label
    domain_columns <- function(x) {
        if ("mediator" %in% names(x)) {
            names(x)[names(x) == "mediator"] <- label
        }
        x
    }
    xm_p <- res$xm %>% dplyr::filter(P_XM < 0.05)
    xm_fdr <- res$xm %>% dplyr::filter(FDR_XM < 0.05)
    my_p <- res$my %>% dplyr::filter(P_MY < 0.05)
    my_fdr <- res$my %>% dplyr::filter(FDR_MY < 0.05)
    s4_write_csv(res$total_effect, p("00_overall_DOSE_IS_total_effect.csv"))
    s4_write_csv(domain_columns(res$associations),
                 p(paste0("01_all_", label, "_associations.csv")))
    s4_write_csv(domain_columns(res$xm),
                 p(paste0("01a_DOSE_", label, "_XM_all.csv")))
    s4_write_csv(domain_columns(res$my),
                 p(paste0("01b_", label, "_IS_MY_all.csv")))
    s4_write_csv(domain_columns(xm_p),
                 p(paste0("02a_DOSE_", label, "_XM_P05.csv")))
    s4_write_csv(domain_columns(xm_fdr),
                 p(paste0("02b_DOSE_", label, "_XM_FDR05.csv")))
    s4_write_csv(domain_columns(my_p),
                 p(paste0("03a_", label, "_IS_MY_P05.csv")))
    s4_write_csv(domain_columns(my_fdr),
                 p(paste0("03b_", label, "_IS_MY_FDR05.csv")))
    s4_write_csv(domain_columns(res$candidates),
                 p(paste0("04_candidate_", label,
                          "_XM_FDR05_MY_P05.csv")))
    utils::write.table(
        res$candidates$mediator,
        p(paste0(candidate_prefix, "_candidates.txt")),
        row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    dose_save_rds(res$candidates$mediator,
            p(paste0(candidate_prefix, "_candidates.rds")))
    s4_write_csv(domain_columns(res$lhr), p("05_mediation_LHR_AFT.csv"))
    s4_write_csv(domain_columns(res$modern),
                 p("06_mediation_MODERN_logistic.csv"))
    s4_write_csv(domain_columns(res$product),
                 p("07_mediation_Cox_product.csv"))
    if (label %in% c("metabolite", "inflammation")) {
        table_name <- if (label == "metabolite") "Table S17" else "Table S19"
        mediation_table <- s4_mediation_table(res$associations, res$product, label)
        s4_write_csv(mediation_table, p(paste0(gsub(" ", "_", table_name), ".csv")))
        extra_sheets[[table_name]] <- mediation_table
    }
    s4_write_csv(res$significance_counts,
                 p("07a_significance_counts_P05_FDR05.csv"))
    s4_write_csv(domain_columns(res$combined),
                 p("08_candidate_mediation_combined.csv"))
    s4_write_csv(domain_columns(res$long),
                 p("09_all_mediation_effects_long.csv"))
    s4_write_csv(domain_columns(res$positive_prop),
                 p("10_positive_proportion_P05_LHR_MODERN.csv"))
    s4_write_csv(domain_columns(res$product_positive),
                 p("11_product_indirect_P05_same_direction.csv"))
    s4_write_csv(res$summary, p("12_summary_counts.csv"))
    s4_write_csv(res$errors, p("13_model_errors.csv"))
    mediator_any <- sort(unique(c(
        res$positive_prop$mediator, res$product_positive$mediator
    )))
    for (method in c("lhr_aft", "modern_logistic")) {
        vals <- res$positive_prop %>% dplyr::filter(.data$method == method) %>%
            dplyr::pull(mediator)
        utils::write.table(vals,
            p(paste0(candidate_prefix, "_mediators_",
                     ifelse(method == "lhr_aft", "LHR", "MODERN"), ".txt")),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    utils::write.table(res$product_positive$mediator,
        p(paste0(candidate_prefix, "_mediators_PRODUCT.txt")),
        row.names = FALSE, col.names = FALSE, quote = FALSE)
    utils::write.table(mediator_any,
        p(paste0(candidate_prefix, "_mediators_any_method.txt")),
        row.names = FALSE, col.names = FALSE, quote = FALSE)

    volcano <- res$associations %>% dplyr::mutate(
        neglog10_P_XM = -log10(pmax(P_XM, .Machine$double.xmin)),
        group = dplyr::case_when(
            candidate ~ "Candidate: X-M FDR<0.05; M-Y P<0.05",
            FDR_XM < 0.05 ~ "X-M FDR<0.05", TRUE ~ "Other"
        )
    )
    labels <- volcano %>% dplyr::filter(candidate) %>%
        dplyr::arrange(P_MY, FDR_XM) %>% dplyr::slice_head(n = 20)
    fig1 <- ggplot2::ggplot(
        volcano, ggplot2::aes(beta_XM, neglog10_P_XM, color = group)
    ) + ggplot2::geom_point(alpha = 0.78, size = 1.7) +
        ggplot2::geom_vline(xintercept = 0, linetype = 3, color = "grey45") +
        ggrepel::geom_text_repel(
            data = labels, ggplot2::aes(label = mediator), size = 3,
            max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
        ) + ggplot2::scale_color_manual(values = c(
            "Candidate: X-M FDR<0.05; M-Y P<0.05" = "#D73027",
            "X-M FDR<0.05" = "#4575B4", "Other" = "grey78"
        )) + ggplot2::labs(
            x = paste0("Beta: DOSE (per SD) to ", label, " (per SD)"),
            y = expression(-log[10](P)), color = NULL,
            title = paste0("DOSE-", label, " associations"),
            subtitle = "Candidates: X-M BH-FDR<0.05 and M-Y Cox P<0.05"
        ) + ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(legend.position = "top")
    bounded_ggsave(p(paste0("Fig1_DOSE_", label, "_volcano.pdf")),
                    fig1, width = 7.5, height = 5.8, device = cairo_pdf)

    if (nrow(res$combined)) {
        bubble <- res$combined %>% dplyr::mutate(
            abs_beta_XM = abs(beta_XM), abs_beta_MY = abs(beta_MY),
            proportion_plot = pmin(pmax(PRODUCT_Proportion_mediated, 0), 1),
            ACME_significant = ifelse(
                !is.na(LHR_FDR_ACME) & LHR_FDR_ACME < 0.05,
                "LHR AFT ACME FDR<0.05", "Not LHR-FDR significant"
            )
        )
        bubble_labels <- bubble %>% dplyr::arrange(
            LHR_P_ACME, PRODUCT_P_ACME
        ) %>% dplyr::slice_head(n = 12)
        fig2 <- ggplot2::ggplot(
            bubble, ggplot2::aes(abs_beta_XM, abs_beta_MY)
        ) + ggplot2::geom_point(
            ggplot2::aes(size = proportion_plot, color = ACME_significant),
            alpha = 0.85
        ) + ggrepel::geom_text_repel(
            data = bubble_labels, ggplot2::aes(label = mediator), size = 3,
            max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE
        ) + ggplot2::scale_color_manual(values = c(
            "LHR AFT ACME FDR<0.05" = "#D73027",
            "Not LHR-FDR significant" = "#4575B4"
        )) + ggplot2::scale_size_continuous(range = c(2.5, 8),
                                             limits = c(0, 1)) +
            ggplot2::labs(
                x = paste0("|beta a|: DOSE-", label),
                y = paste0("|beta|: screening ", label, "-IS Cox"),
                size = "Cox approximate\nproportion", color = NULL,
                title = paste0("Candidate ", label, " mediation")
            ) + ggplot2::theme_classic(base_size = 12) +
            ggplot2::theme(legend.position = "right")
        bounded_ggsave(
            p("Fig2_candidate_mediation_bubble.pdf"), fig2,
            width = 7.2, height = 5.8, device = cairo_pdf
        )
    }

    sheets <- c(list(
        "Summary" = s4_sheet_safe(res$summary),
        "P-FDR counts" = s4_sheet_safe(res$significance_counts),
        "Overall X-Y" = s4_sheet_safe(res$total_effect)
    ), extra_sheets, list(
        "All associations" = s4_sheet_safe(domain_columns(res$associations)),
        "X-M all" = s4_sheet_safe(domain_columns(res$xm)),
        "M-Y all" = s4_sheet_safe(domain_columns(res$my)),
        "X-M P05" = s4_sheet_safe(domain_columns(xm_p)),
        "X-M FDR05" = s4_sheet_safe(domain_columns(xm_fdr)),
        "M-Y P05" = s4_sheet_safe(domain_columns(my_p)),
        "M-Y FDR05" = s4_sheet_safe(domain_columns(my_fdr)),
        "Candidates" = s4_sheet_safe(domain_columns(res$candidates)),
        "LHR AFT" = s4_sheet_safe(domain_columns(res$lhr)),
        "MODERN logistic" = s4_sheet_safe(domain_columns(res$modern)),
        "Cox product" = s4_sheet_safe(domain_columns(res$product)),
        "All effects long" = s4_sheet_safe(domain_columns(res$long)),
        "Combined" = s4_sheet_safe(domain_columns(res$combined)),
        "Model errors" = s4_sheet_safe(domain_columns(res$errors))
    ))
    openxlsx::write.xlsx(sheets, p(workbook_name), overwrite = TRUE,
                         asTable = TRUE)
    common_notes <- c(
        "Flow 1: DOSE -> M linear regression; BH-FDR<0.05.",
        "Flow 2: M -> incident IS Cox regression adjusted for M2 but not DOSE; raw P<0.05.",
        "Candidates are the intersection of flow-1 FDR<0.05 and flow-2 P<0.05.",
        "Formal LHR and MODERN outcome models contain both dose_z and mediator_z plus M2.",
        "LHR uses lm + survreg + mediation::mediate; no dist argument is specified.",
        "MODERN uses lm + logistic regression + mediation::mediate.",
        "The Cox product method uses a*b with the delta/Sobel standard error.",
        paste0(
            "For the Cox product method, proportion=(a*b)/c. Its SE, 95% CI, ",
            "P value and BH-FDR use a first-order delta approximation including ",
            "SE(a), SE(b) and SE(c), while assuming cross-model covariances are zero."
        ),
        paste0(
            "Cox product proportions are exploratory and may be unstable when ",
            "the total-effect coefficient c is near zero or imprecisely estimated."
        ),
        "sims=100 and seed=12345 unless settings in the calling script are changed.",
        "DOSE contains Wine; separate alcohol covariates are excluded from M2."
    )
    writeLines(c(common_notes, notes), p("README_three_method_analysis.txt"))
    invisible(res)
}

s4_write_manifest <- function(input_files, dir_out) {
    paths <- unname(input_files)
    exists <- file.exists(paths)
    out <- data.frame(
        name = names(input_files), path = paths, exists = exists,
        md5 = NA_character_, stringsAsFactors = FALSE
    )
    out$md5[exists] <- unname(tools::md5sum(paths[exists]))
    utils::write.csv(out, file.path(dir_out, "input_manifest.csv"),
                     row.names = FALSE)
    capture.output(sessionInfo(), file = file.path(dir_out, "sessionInfo.txt"))
    invisible(out)
}

# ===========================================================================
# Signature helpers
# ===========================================================================
# Shared functions for DOSE-related proteomic/metabolomic signatures.
# Omics-specific LASSO fitting is performed by the two Python scripts.
# This file provides score QC, quartile/continuous Cox models and three
# score-level mediation methods through the validated s4 functions.

s5_required_packages <- c("survival", "broom", "dplyr", "openxlsx", "mediation")

s5_load_packages <- function() {
    missing <- s5_required_packages[
        !vapply(s5_required_packages, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing)) {
        stop("Install required R packages first: ", paste(missing, collapse = ", "))
    }
    suppressPackageStartupMessages({
        library(survival)
        library(dplyr)
    })
    invisible(TRUE)
}

s5_validate_ids <- function(x, label) {
    x <- as.data.frame(x)
    if (!"eid" %in% names(x)) stop(label, " has no eid column.")
    x$eid <- as.character(x$eid)
    if (anyNA(x$eid) || anyDuplicated(x$eid)) {
        stop(label, " has missing or duplicated eid values.")
    }
    x
}

s5_prepare_lasso_input <- function(domain, dir_project = NULL, dir_ukb = NULL) {
    domain <- match.arg(domain, c("protein", "metabolite"))
    if (is.null(dir_project)) {
        dir_project <- dose_output_dir()
    }
    if (is.null(dir_ukb)) {
        dir_ukb <- dose_ukb_dir()
    }
    is_protein <- identical(domain, "protein")
    association_file <- file.path(
        dir_project, "ukb/05mechanism", domain, "M2",
        paste0("02b_DOSE_", domain, "_XM_FDR05.csv")
    )
    omics_file <- file.path(dir_ukb, "Rdata",
                            if (is_protein) "prot.rds" else "met.rds")
    dose_file <- file.path(dir_project, "ukb/04score/dose_score.rds")
    inputs <- c(association = association_file, omics = omics_file, dose = dose_file)
    if (any(!file.exists(inputs))) {
        stop("Missing LASSO preparation input(s): ",
             paste(inputs[!file.exists(inputs)], collapse = "; "))
    }
    assoc <- utils::read.csv(association_file, check.names = FALSE)
    feature_col <- domain
    if (!all(c(feature_col, "FDR_XM") %in% names(assoc))) {
        stop("Invalid s4 association file: ", association_file)
    }
    features <- unique(as.character(assoc[[feature_col]][
        !is.na(assoc$FDR_XM) & assoc$FDR_XM < 0.05
    ]))
    features <- features[!is.na(features) & nzchar(features)]
    if (!length(features)) stop("No ", domain, " has X-M BH-FDR<0.05.")
    omics <- s5_validate_ids(readRDS(omics_file), paste0(domain, " RDS"))
    dose <- s5_validate_ids(readRDS(dose_file), "dose score")
    if (!"dose" %in% names(dose)) stop("dose_score.rds has no dose column.")
    missing <- setdiff(features, names(omics))
    if (length(missing)) {
        stop("FDR-selected ", domain, " variables missing from RDS: ",
             paste(head(missing, 20), collapse = ", "))
    }
    out <- dplyr::inner_join(
        omics[c("eid", features)], dose[c("eid", "dose")], by = "eid"
    )
    if (anyDuplicated(out$eid)) stop("Duplicated eid after LASSO input join.")
    out <- out[is.finite(out$dose), , drop = FALSE]
    rm(omics, dose)
    invisible(gc())
    dir_out <- file.path(dir_project, "ukb/06signature", domain)
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    output_file <- file.path(dir_out, paste0(domain, "_LASSO_input.csv.gz"))
    con <- gzfile(output_file, open = "wt")
    on.exit(close(con), add = TRUE)
    utils::write.csv(out, con, row.names = FALSE, na = "")
    close(con)
    on.exit(NULL, add = FALSE)
    writeLines(features, file.path(dir_out, paste0(domain, "_FDR05_features.txt")))
    utils::write.csv(
        data.frame(domain = domain, N = nrow(out), N_features = length(features),
                   association_file = association_file, omics_file = omics_file,
                   dose_file = dose_file),
        file.path(dir_out, paste0(domain, "_LASSO_input_QC.csv")), row.names = FALSE
    )
    message("Prepared ", domain, " LASSO input: N=", nrow(out),
            "; features=", length(features), "; ", output_file)
    invisible(output_file)
}

s5_formula <- function(response, terms) {
    stats::reformulate(termlabels = terms, response = response)
}

s5_valid_covariate <- function(x) {
    if (is.factor(x) || is.character(x) || is.logical(x)) {
        return(dplyr::n_distinct(x, na.rm = TRUE) > 1L)
    }
    sx <- suppressWarnings(stats::sd(as.numeric(x), na.rm = TRUE))
    is.finite(sx) && sx > 0
}

s5_extract_cox <- function(fit, term, increment = 1) {
    b <- unname(stats::coef(fit)[term])
    se <- sqrt(stats::vcov(fit)[term, term])
    if (!is.finite(b) || !is.finite(se) || se <= 0) {
        stop("Invalid Cox coefficient for ", term)
    }
    data.frame(
        beta = b * increment,
        SE = se * abs(increment),
        HR = exp(b * increment),
        LCI = exp((b - stats::qnorm(0.975) * se) * increment),
        HCI = exp((b + stats::qnorm(0.975) * se) * increment),
        P = 2 * stats::pnorm(-abs(b / se)),
        stringsAsFactors = FALSE
    )
}

s5_lr_p <- function(full, reduced, df) {
    stats::pchisq(
        max(0, 2 * (as.numeric(stats::logLik(full)) -
                    as.numeric(stats::logLik(reduced)))),
        df = df, lower.tail = FALSE
    )
}

s5_ph_p <- function(fit, term) {
    tryCatch({
        z <- survival::cox.zph(fit, terms = TRUE)$table
        if (!term %in% rownames(z)) return(NA_real_)
        unname(z[term, "p"])
    }, error = function(e) NA_real_)
}

s5_model_covariates <- function(cfg, model) {
    # The signature is trained to reproduce DOSE, which includes Wine.
    # Match the s4/DOSE convention by excluding separate alcohol adjustment.
    alcohol <- unique(c(cfg$covs.alcohol, "alcohol.c"))
    covs <- switch(
        model,
        M1 = unique(cfg$covs.m1),
        M2 = unique(c(cfg$covs.m1, cfg$covs.m2.add)),
        M3 = unique(c(cfg$covs.m1, cfg$covs.m2.add, cfg$covs.m3.add)),
        stop("Unknown model: ", model)
    )
    covs <- setdiff(covs, alcohol)
    if (!length(covs) || anyNA(covs) || any(!nzchar(covs))) {
        stop("Invalid covariate definition for ", model)
    }
    covs
}

s5_make_quartiles <- function(score) {
    if (!is.numeric(score) || any(!is.finite(score))) {
        stop("Signature score must be finite numeric values.")
    }
    breaks <- unname(stats::quantile(
        score, probs = c(0, 0.25, 0.50, 0.75, 1), type = 7, na.rm = TRUE
    ))
    if (anyDuplicated(breaks)) {
        stop("Fewer than four distinct signature-score quantile intervals.")
    }
    group <- cut(
        score, breaks = breaks, include.lowest = TRUE, right = TRUE,
        labels = paste0("Q", 1:4)
    )
    group <- stats::relevel(group, ref = "Q1")
    list(group = group, breaks = breaks)
}

s5_prepare_analysis_data <- function(dat, dose, score, score_col,
                                     time_var, event_var, all_covariates) {
    dat <- s5_validate_ids(dat, "cohort")
    dose <- s5_validate_ids(dose, "dose")
    score <- s5_validate_ids(score, "signature score")
    if (!score_col %in% names(score)) stop("Missing score column: ", score_col)
    if (!is.numeric(score[[score_col]])) stop(score_col, " must be numeric.")
    required <- unique(c("eid", time_var, event_var, all_covariates))
    if (!all(required %in% names(dat))) {
        stop("Missing cohort fields: ", paste(setdiff(required, names(dat)), collapse = ", "))
    }
    if (!"dose" %in% names(dose)) stop("dose_score.rds has no dose column.")
    out <- dat[, required, drop = FALSE]
    out$dose <- dose$dose[match(out$eid, dose$eid)]
    names(score)[names(score) == score_col] <- "signature_score"
    out <- dplyr::inner_join(out, score[c("eid", "signature_score")], by = "eid")
    out
}

s5_run_score_cox <- function(dat, time_var, event_var, cfg,
                             signature_label, models = c("M1", "M2", "M3")) {
    complete_covs <- unique(unlist(lapply(models, function(m) s5_model_covariates(cfg, m))))
    required <- unique(c("signature_score", time_var, event_var, complete_covs))
    ok <- complete.cases(dat[required]) & is.finite(dat[[time_var]]) &
        dat[[time_var]] > 0 & dat[[event_var]] %in% 0:1 &
        is.finite(dat$signature_score)
    d0 <- droplevels(dat[ok, , drop = FALSE])
    if (!nrow(d0) || sum(d0[[event_var]]) < 20L) stop("Insufficient Cox N/events.")
    d0$time <- d0[[time_var]]
    d0$event <- as.integer(d0[[event_var]])
    d0$score_z <- as.numeric(scale(d0$signature_score))
    q <- s5_make_quartiles(d0$signature_score)
    d0$score_q <- q$group
    medians <- tapply(d0$signature_score, d0$score_q, stats::median)
    d0$score_trend <- unname(medians[as.character(d0$score_q)])

    quartile_definitions <- do.call(rbind, lapply(1:4, function(j) {
        idx <- d0$score_q == paste0("Q", j)
        data.frame(
            category = paste0("Q", j), cut_lower = q$breaks[j],
            cut_upper = q$breaks[j + 1L], lower_inclusive = j == 1L,
            upper_inclusive = TRUE, N = sum(idx), Events = sum(d0$event[idx]),
            score_min = min(d0$signature_score[idx]),
            score_max = max(d0$signature_score[idx]),
            score_median = stats::median(d0$signature_score[idx]),
            stringsAsFactors = FALSE
        )
    }))

    continuous <- list()
    categories <- list()
    overall <- list()
    model_spec <- list()
    for (model in models) {
        covs <- s5_model_covariates(cfg, model)
        bad <- covs[!vapply(d0[covs], s5_valid_covariate, logical(1))]
        if (length(bad)) stop(model, " covariates without variation: ", paste(bad, collapse = ", "))
        null <- survival::coxph(
            s5_formula("Surv(time, event)", covs), data = d0,
            ties = "efron", singular.ok = FALSE
        )
        linear <- survival::coxph(
            s5_formula("Surv(time, event)", c("score_z", covs)), data = d0,
            ties = "efron", singular.ok = FALSE
        )
        categorical <- survival::coxph(
            s5_formula("Surv(time, event)", c("score_q", covs)), data = d0,
            ties = "efron", singular.ok = FALSE
        )
        trend <- survival::coxph(
            s5_formula("Surv(time, event)", c("score_trend", covs)), data = d0,
            ties = "efron", singular.ok = FALSE
        )
        continuous[[model]] <- cbind(
            data.frame(signature = signature_label, outcome = event_var,
                       model = model, scale = "Per 1-SD", N = linear$n,
                       Events = linear$nevent, PH_P = s5_ph_p(linear, "score_z")),
            s5_extract_cox(linear, "score_z")
        )
        categories[[model]] <- do.call(rbind, lapply(1:4, function(j) {
            level <- paste0("Q", j)
            idx <- d0$score_q == level
            effect <- if (j == 1L) data.frame(
                beta = 0, SE = NA_real_, HR = 1, LCI = 1, HCI = 1, P = NA_real_
            ) else s5_extract_cox(categorical, paste0("score_q", level))
            cbind(data.frame(
                signature = signature_label, outcome = event_var, model = model,
                category = level, reference = "Q1", N = categorical$n,
                Events = categorical$nevent, Group_N = sum(idx),
                Group_Events = sum(d0$event[idx]),
                Person_years = sum(d0$time[idx]),
                Incidence_percent = 100 * mean(d0$event[idx]),
                Rate_per_100000 = 1e5 * sum(d0$event[idx]) / sum(d0$time[idx]),
                score_min = min(d0$signature_score[idx]),
                score_max = max(d0$signature_score[idx]),
                score_median = stats::median(d0$signature_score[idx])
            ), effect)
        }))
        overall[[model]] <- data.frame(
            signature = signature_label, outcome = event_var, model = model,
            N = categorical$n, Events = categorical$nevent,
            P_global = s5_lr_p(categorical, null, 3L),
            P_trend = s5_extract_cox(trend, "score_trend")$P,
            PH_P = s5_ph_p(categorical, "score_q"),
            AIC_linear = stats::AIC(linear),
            AIC_categorical = stats::AIC(categorical)
        )
        model_spec[[model]] <- data.frame(
            model = model, covariates = paste(covs, collapse = ";"),
            separate_alcohol_adjustment = FALSE
        )
    }
    list(
        data = d0,
        score_summary = data.frame(
            signature = signature_label, N = nrow(d0), Events = sum(d0$event),
            mean = mean(d0$signature_score), SD = stats::sd(d0$signature_score),
            min = min(d0$signature_score), Q25 = q$breaks[2],
            median = q$breaks[3], Q75 = q$breaks[4], max = max(d0$signature_score)
        ),
        quartile_definitions = quartile_definitions,
        continuous = dplyr::bind_rows(continuous),
        categories = dplyr::bind_rows(categories),
        overall = dplyr::bind_rows(overall),
        model_specification = dplyr::bind_rows(model_spec)
    )
}

s5_run_score_mediation <- function(dat, time_var, event_var, cfg,
                                    sims = 100L, seed = 12345L) {
    required_functions <- c(
        "s4_mediation_lhr_aft", "s4_mediation_modern_logistic",
        "s4_mediation_cox_product", "s4_transform_direct_z"
    )
    if (!all(vapply(required_functions, exists, logical(1), mode = "function"))) {
        stop("Source 0.common.R before running score mediation.")
    }
    covs <- s5_model_covariates(cfg, "M2")
    required <- unique(c("dose", "signature_score", time_var, event_var, covs))
    ok <- complete.cases(dat[required]) & is.finite(dat$dose) &
        is.finite(dat$signature_score) & is.finite(dat[[time_var]]) &
        dat[[time_var]] > 0 & dat[[event_var]] %in% 0:1
    d <- droplevels(dat[ok, , drop = FALSE])
    if (nrow(d) < 1000L || sum(d[[event_var]]) < 20L) {
        stop("Insufficient M2 N/events for score mediation.")
    }
    d$dose_z <- as.numeric(scale(d$dose))
    args <- list(
        dat = d, mediator = "signature_score",
        transform_mediator = s4_transform_direct_z,
        covariates = covs, time_var = time_var, event_var = event_var,
        min_n = 1000L, min_events = 20L, sims = sims, seed = seed
    )
    list(
        lhr_aft = do.call(s4_mediation_lhr_aft, args),
        modern_logistic = do.call(s4_mediation_modern_logistic, args),
        cox_product = do.call(s4_mediation_cox_product, args)
    )
}

s5_write_analysis <- function(cox, mediation = NULL, dir_out, prefix, input_files) {
    dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
    has_mediation <- !is.null(mediation)
    med_all <- if (has_mediation) {
        dplyr::bind_rows(mediation)
    } else {
        data.frame(note = "Mediation not run for this Cox-only sensitivity analysis")
    }
    utils::write.csv(cox$score_summary,
                     file.path(dir_out, paste0(prefix, "_score_summary.csv")), row.names = FALSE)
    utils::write.csv(cox$quartile_definitions,
                     file.path(dir_out, paste0(prefix, "_quartile_definitions.csv")), row.names = FALSE)
    utils::write.csv(cox$continuous,
                     file.path(dir_out, paste0(prefix, "_IS_continuous.csv")), row.names = FALSE)
    utils::write.csv(cox$categories,
                     file.path(dir_out, paste0(prefix, "_IS_quartiles.csv")), row.names = FALSE)
    utils::write.csv(cox$overall,
                     file.path(dir_out, paste0(prefix, "_IS_overall_tests.csv")), row.names = FALSE)
    utils::write.csv(cox$model_specification,
                     file.path(dir_out, paste0(prefix, "_model_specification.csv")), row.names = FALSE)
    utils::write.csv(med_all,
                     file.path(dir_out, paste0(prefix, "_three_method_mediation.csv")), row.names = FALSE)
    sheets <- list(
        "Score summary" = cox$score_summary,
        "Quartile definitions" = cox$quartile_definitions,
        "Continuous Cox" = cox$continuous,
        "Quartile Cox" = cox$categories,
        "Overall tests" = cox$overall,
        "Model specification" = cox$model_specification
    )
    if (has_mediation) {
        sheets <- c(sheets, list(
            "LHR AFT mediation" = mediation$lhr_aft,
            "MODERN mediation" = mediation$modern_logistic,
            "Cox product mediation" = mediation$cox_product
        ))
    } else {
        sheets[["Mediation status"]] <- med_all
    }
    openxlsx::write.xlsx(
        sheets, file.path(dir_out, paste0(prefix, "_IS_Cox_and_mediation.xlsx")),
    overwrite = TRUE, asTable = TRUE)
    manifest <- data.frame(
        name = names(input_files), path = unname(input_files),
        exists = file.exists(input_files), md5 = NA_character_
    )
    manifest$md5[manifest$exists] <- unname(tools::md5sum(input_files[manifest$exists]))
    utils::write.csv(
        manifest, file.path(dir_out, paste0(prefix, "_input_manifest.csv")),
        row.names = FALSE
    )
    capture.output(
        sessionInfo(),
        file = file.path(dir_out, paste0(prefix, "_sessionInfo.txt"))
    )
    invisible(list(cox = cox, mediation = mediation))
}

s5_run_signature_analysis <- function(score_file, score_col, signature_label,
                                      output_prefix, dir_project,
                                      sims = 100L, seed = 12345L,
                                      run_mediation = TRUE) {
    s5_load_packages()
    dir_data <- file.path(dir_project, "ukb/00data")
    dir_score <- file.path(dir_project, "ukb/04score")
    file_s4 <- file.path(dose_script_dir(), "0.common.R")
    files <- c(
        s4_functions = file_s4,
        cohort = file.path(dir_data, "dat.base.imp.rds"),
        config = file.path(dir_data, "analysis_config.rds"),
        dose = file.path(dir_score, "dose_score.rds"),
        signature = score_file
    )
    if (any(!file.exists(files))) {
        stop("Missing input(s): ", paste(files[!file.exists(files)], collapse = "; "))
    }

    cfg <- readRDS(files[["config"]])
    outcome <- "cvd_stroke_i"
    time_var <- paste0(outcome, ".t2e")
    event_var <- paste0(outcome, ".Yt2e")
    if (!outcome %in% cfg$Ys ||
        !identical(cfg$primary_development_outcome, outcome) ||
        !identical(cfg$primary_model, "M2")) {
        stop("analysis_config.rds does not define cvd_stroke_i/M2 as primary.")
    }
    dat <- readRDS(files[["cohort"]])
    dose <- readRDS(files[["dose"]])
    score <- utils::read.csv(score_file, check.names = FALSE)
    all_covs <- unique(unlist(lapply(c("M1", "M2", "M3"), function(m) {
        s5_model_covariates(cfg, m)
    })))
    analysis_data <- s5_prepare_analysis_data(
        dat, dose, score, score_col, time_var, event_var, all_covs
    )
    rm(dat, dose, score)
    invisible(gc())
    cox <- s5_run_score_cox(
        analysis_data, time_var, event_var, cfg, signature_label
    )
    mediation <- if (isTRUE(run_mediation)) {
        s5_run_score_mediation(
            analysis_data, time_var, event_var, cfg, sims, seed
        )
    } else NULL
    dir_out <- file.path(dirname(score_file), "association_mediation")
    s5_write_analysis(cox, mediation, dir_out, output_prefix, files)
    message("Completed ", signature_label, " analysis: ", dir_out,
            if (isTRUE(run_mediation)) "" else " (Cox only)")
    invisible(list(cox = cox, mediation = mediation, output_dir = dir_out))
}
