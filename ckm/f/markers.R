# Cross-sectional molecular stage markers and optional offline gene-set ORA.
# Source associations.R before this file. Markers do not establish drivers.
ckm_run_markers <- function(d, cfg, exposures) {
    a <- ckm_a_options(cfg)
    opt <- utils::modifyList(list(reference_stage = 0L, stages = 0:4, min_n = 100L, min_stage_n = 20L,
                                 covariates = a$covariates), ckm_a_default(cfg$markers, list()))
    ref <- as.integer(opt$reference_stage); stages <- sort(unique(as.integer(unlist(opt$stages))))
    if (length(ref) != 1L || is.na(ref) || !ref %in% stages || anyNA(stages) || any(!stages %in% 0:4)) stop("Invalid marker stage/reference configuration")
    d <- as.data.frame(d); m <- ckm_a_manifest(exposures$manifest); x <- as.data.frame(exposures$data)
    required <- c("eid", a$stage_column, as.character(unlist(opt$covariates)))
    if (length(setdiff(required, names(d)))) stop("Marker adjustment/data columns missing: ", paste(setdiff(required, names(d)), collapse = ","))
    if (!"eid" %in% names(x) || anyNA(d$eid) || anyNA(x$eid) || anyDuplicated(d$eid) || anyDuplicated(x$eid)) stop("Markers require unique, nonmissing eid in both tables")
    mi <- match(as.character(d$eid), as.character(x$eid)); sg <- suppressWarnings(as.integer(as.character(d[[a$stage_column]])))
    covars <- as.character(unlist(opt$covariates)); cv <- d[, covars, drop = FALSE]
    for (v in covars) {
        if (is.character(cv[[v]]) || is.factor(cv[[v]])) cv[[v]] <- factor(cv[[v]])
        if (is.numeric(cv[[v]])) cv[[v]][!is.finite(cv[[v]])] <- NA_real_
    }
    base <- !is.na(sg) & sg %in% stages & stats::complete.cases(cv)
    if (!any(base)) stop("No complete observations for configured marker adjustment")
    cv <- droplevels(cv[base, , drop = FALSE]); ss <- sg[base]
    constants <- names(cv)[vapply(cv, function(v) length(unique(v)) < 2L, logical(1))]
    cv <- cv[, setdiff(names(cv), constants), drop = FALSE]
    if (ncol(cv)) names(cv) <- paste0("c", seq_along(cv))
    Xc <- if (ncol(cv)) stats::model.matrix(~ ., data = cv) else matrix(1, nrow = sum(base), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
    targets <- setdiff(stages, ref); Xs <- vapply(targets, function(s) as.numeric(ss == s), numeric(length(ss)))
    if (is.null(dim(Xs))) Xs <- matrix(Xs, ncol = 1L)
    colnames(Xs) <- paste0("stage", targets); X <- cbind(Xc, Xs)
    rows <- lapply(seq_len(nrow(m)), function(i) {
        r <- m[i, ]; meta <- data.frame(source = r$source, feature = r$feature, column = r$column, reference_stage = ref,
            transform = r$transform, interpretation = "cross_sectional_stage_marker", stringsAsFactors = FALSE)
        out <- cbind(meta, data.frame(n = 0L, p = NA_real_, df_stage = NA_integer_, df_residual = NA_integer_,
            tested_stages = "", scale_sd = NA_real_, status = "not_estimable", dropped_constant_covariates = paste(constants, collapse = ",")))
        effects <- cbind(meta[rep(1L, length(targets)), , drop = FALSE], data.frame(stage = targets, n_reference = 0L, n_stage = 0L,
            beta = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, status = "not_estimable"))
        fail <- function(reason) {out$status <- reason; effects$status <- reason; list(marker = out, stage = effects)}
        if (r$type != "continuous") return(fail("requires_continuous_assay"))
        if (!r$column %in% names(x)) return(fail("missing_exposure_column"))
        prep <- ckm_a_prepare_exposure(x[[r$column]][mi], r, !is.na(sg) & sg %in% stages)
        if (!prep$ok) return(fail(prep$status))
        y <- as.numeric(prep$x[base, 1L]); ok <- is.finite(y)
        counts <- table(factor(ss[ok], levels = stages)); effects$n_reference <- as.integer(counts[as.character(ref)])
        effects$n_stage <- as.integer(counts[as.character(targets)]); out$scale_sd <- prep$scale
        if (counts[as.character(ref)] < opt$min_stage_n) return(fail("reference_stage_unavailable_or_sparse"))
        tested <- stages[counts >= opt$min_stage_n]; ok <- ok & ss %in% tested
        out$n <- sum(ok); out$tested_stages <- paste(tested, collapse = ",")
        if (sum(ok) < opt$min_n || length(tested) < 2L) return(fail("insufficient_marker_sample"))
        xx <- X[ok, , drop = FALSE]; active <- colSums(abs(xx)) > 0
        xx <- xx[, active, drop = FALSE]; fit <- stats::lm.fit(xx, y[ok]); rank <- fit$rank; df <- sum(ok) - rank
        if (rank < ncol(xx)) return(fail("aliased_marker_design"))
        if (df <= 0L) return(fail("insufficient_residual_df"))
        keep <- fit$qr$pivot[seq_len(rank)]; variance <- sum(fit$residuals^2) / df
        V <- variance * chol2inv(qr.R(fit$qr)[seq_len(rank), seq_len(rank), drop = FALSE]); dimnames(V) <- list(colnames(xx)[keep], colnames(xx)[keep])
        stage_terms <- intersect(paste0("stage", setdiff(tested, ref)), names(fit$coefficients)[keep])
        if (!length(stage_terms) || !is.finite(variance) || variance <= 0) return(fail("nonestimable_stage_effect"))
        b <- fit$coefficients[stage_terms]; vv <- V[stage_terms, stage_terms, drop = FALSE]
        stat <- tryCatch(as.numeric(crossprod(b, solve(vv, b))) / length(b), error = function(e) NA_real_)
        out$p <- stats::pf(stat, length(b), df, lower.tail = FALSE); out$df_stage <- length(b); out$df_residual <- df
        out$status <- if (is.finite(out$p)) "ok" else "nonestimable_stage_test"
        for (j in seq_along(targets)) {
            term <- paste0("stage", targets[j])
            if (!targets[j] %in% tested) {effects$status[j] <- "sparse_stage"; next}
            if (!term %in% stage_terms) {effects$status[j] <- "aliased_stage_effect"; next}
            be <- fit$coefficients[term]; se <- sqrt(V[term, term]); cutoff <- stats::qt(.975, df)
            effects$beta[j] <- be; effects$ci_low[j] <- be - cutoff * se; effects$ci_high[j] <- be + cutoff * se
            effects$p[j] <- 2 * stats::pt(-abs(be / se), df); effects$status[j] <- "ok"
        }
        list(marker = out, stage = effects)
    })
    marker <- ckm_a_adjust(ckm_a_bind(lapply(rows, `[[`, "marker")), "p", "q", "source")
    effects <- ckm_a_adjust(ckm_a_bind(lapply(rows, `[[`, "stage")), "p", "q", c("source", "stage"))
    list(marker = marker, stage = effects, settings = opt)
}

ckm_run_enrichment <- function(results, cfg) {
    opt <- utils::modifyList(list(min_size = 5L, max_size = 500L, min_selected = 3L,
                                 fdr = ckm_a_options(cfg)$fdr), ckm_a_default(cfg$enrichment, list()))
    status_row <- function(status, message, source = NA_character_, set = NA_character_) data.table::data.table(
        source = source, set = set, pathway = NA_character_, description = NA_character_,
        background_genes = NA_integer_, selected_genes = NA_integer_, pathway_genes = NA_integer_, overlap = NA_integer_,
        fold_enrichment = NA_real_, p = NA_real_, q = NA_real_, genes = "", status = status, message = message)
    annotation <- cfg$paths$annotation; gmt <- cfg$paths$gmt
    if (is.null(annotation) || is.null(gmt) || !nzchar(annotation) || !nzchar(gmt) || !file.exists(annotation) || !file.exists(gmt)) {
        return(status_row("not_configured", "Configure readable annotation (source,feature,gene) and GMT files to run offline ORA."))
    }
    tab <- results$assoc
    if (is.null(tab) || !nrow(tab)) return(status_row("associations_not_available", "No prospective association results supplied."))
    if (opt$min_size < 1 || opt$max_size < opt$min_size || opt$min_selected < 1) stop("Invalid enrichment gene-set size configuration")
    ann <- if (grepl("[.]rds$", annotation, ignore.case = TRUE)) as.data.frame(readRDS(annotation)) else as.data.frame(data.table::fread(annotation))
    if (!all(c("source", "feature", "gene") %in% names(ann))) stop("Annotation requires source, feature, gene columns")
    ann <- unique(ann[, c("source", "feature", "gene")]); ann[] <- lapply(ann, as.character)
    ann <- ann[stats::complete.cases(ann) & nzchar(ann$gene) & nzchar(ann$feature), , drop = FALSE]
    connection <- if (grepl("[.]gz$", gmt)) gzfile(gmt, "rt") else file(gmt, "rt")
    lines <- readLines(connection, warn = FALSE); close(connection)
    raw <- strsplit(lines[nzchar(lines) & !grepl("^#", lines)], "\t", fixed = FALSE)
    raw <- raw[lengths(raw) >= 3L]
    if (!length(raw)) stop("GMT has no valid gene sets")
    ids <- vapply(raw, `[`, character(1), 1L)
    if (anyDuplicated(ids)) stop("GMT repeats pathway identifiers")
    sets <- stats::setNames(lapply(raw, function(v) unique(v[-c(1L, 2L)][nzchar(v[-c(1L, 2L)])])), ids)
    desc <- stats::setNames(vapply(raw, `[`, character(1), 2L), ids)
    tab <- as.data.frame(tab); out <- list()
    needed <- c("source", "feature", "status", "p_shared", "q_shared", "p_interaction", "p_equivalence", "pattern", "direction")
    if (length(setdiff(needed, names(tab)))) stop("Association table lacks enrichment fields: ", paste(setdiff(needed, names(tab)), collapse = ","))
    for (source in unique(tab$source)) {
        z <- tab[tab$source == source, , drop = FALSE]; aa <- ann[ann$source == source, , drop = FALSE]
        main <- z$status == "ok" & is.finite(z$p_shared)
        candidates <- list(
            shared_risk = list(tested = main & is.finite(z$p_equivalence), selected = z$pattern == "shared_equivalent_within_margin" & z$direction == "risk_associated"),
            shared_protective = list(tested = main & is.finite(z$p_equivalence), selected = z$pattern == "shared_equivalent_within_margin" & z$direction == "protective_associated"),
            stage_dependent = list(tested = main & is.finite(z$p_interaction), selected = grepl("^stage_dependent", z$pattern)),
            pooled_risk = list(tested = main, selected = is.finite(z$q_shared) & z$q_shared < opt$fdr & z$direction == "risk_associated"),
            pooled_protective = list(tested = main, selected = is.finite(z$q_shared) & z$q_shared < opt$fdr & z$direction == "protective_associated"))
        for (group in names(candidates)) {
            cnd <- candidates[[group]]; cnd$tested[is.na(cnd$tested)] <- FALSE; cnd$selected[is.na(cnd$selected)] <- FALSE
            background <- unique(aa$gene[aa$feature %in% z$feature[cnd$tested]])
            selected <- intersect(background, unique(aa$gene[aa$feature %in% z$feature[cnd$tested & cnd$selected]]))
            if (!length(background) || length(selected) < opt$min_selected) {
                rr <- status_row(if (!length(background)) "no_mapped_tested_genes" else "insufficient_selected_genes",
                                 "Gene identifiers must match the GMT; background is the successfully tested, mapped assay universe.", source, group)
                rr$background_genes <- length(background); rr$selected_genes <- length(selected); out[[length(out) + 1L]] <- rr; next
            }
            rows <- lapply(names(sets), function(id) {
                pathway <- intersect(sets[[id]], background); K <- length(pathway); N <- length(background); n <- length(selected)
                if (K < opt$min_size || K > opt$max_size) return(NULL)
                overlap <- intersect(selected, pathway); k <- length(overlap)
                data.table::data.table(source = source, set = group, pathway = id, description = desc[[id]], background_genes = N,
                    selected_genes = n, pathway_genes = K, overlap = k, fold_enrichment = (k / n) / (K / N),
                    p = stats::phyper(k - 1L, K, N - K, n, lower.tail = FALSE), genes = paste(sort(overlap), collapse = ";"),
                    status = "ok", message = "Observational gene-set association; zero-overlap pathways remain in the FDR family.")
            })
            rr <- ckm_a_bind(rows)
            if (!nrow(rr)) rr <- status_row("no_eligible_gene_sets", "No GMT set meets the size limits within the tested assay background.", source, group)
            out[[length(out) + 1L]] <- rr
        }
    }
    ans <- ckm_a_bind(out)
    if (nrow(ans)) ans <- ckm_a_adjust(ans, "p", "q", c("source", "set"))
    ans
}
