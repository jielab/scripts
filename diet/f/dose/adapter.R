#!/usr/bin/env Rscript
# Diet coordinator adapter. The scientific engines retain their original models.
# DOSE_ACTION is supplied by diet/1.disco.sh or diet/diet.sh.

.entry <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))[1L]
if (is.na(.entry)) stop("Run this adapter with Rscript")
.scripts <- dirname(normalizePath(.entry, winslash="/", mustWork=TRUE))
Sys.setenv(DOSE_SCRIPT_DIR=.scripts)
source(file.path(.scripts, "0.common.R"), local=globalenv())

.action <- Sys.getenv("DOSE_ACTION", "")
.diet_root <- Sys.getenv("DIET_OUT", dirname(dose_output_dir()))
.part <- function(cohort) Sys.getenv("DIET_DISCO_PART",
    file.path(.diet_root, ".cache", "disco", paste0(cohort, ".dose.rds")))
.require <- function(paths, instruction) {
    absent <- paths[!file.exists(paths)]
    if (length(absent)) stop(instruction, "\nMissing: ", paste(absent, collapse="\n"))
}
.save <- function(x, path) {
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
    dose_save_rds(x, path)
    invisible(path)
}
.with_phases <- function(phases, code) {
    old <- Sys.getenv("DOSE_PREPROCESS_PHASES", unset=NA_character_)
    on.exit(if (is.na(old)) Sys.unsetenv("DOSE_PREPROCESS_PHASES") else
        Sys.setenv(DOSE_PREPROCESS_PHASES=old), add=TRUE)
    Sys.setenv(DOSE_PREPROCESS_PHASES=paste(phases, collapse=","))
    force(code)
}
.run <- function(module, step, env=new.env(parent=globalenv())) {
    message("DOSE ", module, "/", step)
    dose_source_step(module, step, envir=env)
    invisible(env)
}
.score_functions <- function() {
    old <- options(dose.score.run=FALSE)
    on.exit(options(old), add=TRUE)
    .run("score", "build")
}
.nh_functions <- function() {
    old <- options(dose.nhanes.run=FALSE)
    on.exit(options(old), add=TRUE)
    env <- .run("nhanes", "score")
    .run("nhanes", "patterns", env)
    env
}
.canonical <- function(cohort) {
    path <- Sys.getenv("DIET_DISCO_FILE", file.path(.diet_root, "disco", paste0(cohort, ".disco.rds")))
    .require(path, paste0("Run 1.disco.sh ", cohort, " before the phenotype stage."))
    contract <- new.env(parent=globalenv())
    source(file.path(dirname(.scripts), "common.R"), local=contract)
    d <- contract$diet_profile(readRDS(path), "dose", cohort)
    if (anyNA(d$eid) || anyDuplicated(d$eid)) stop("Canonical DOSE profile has invalid eid values")
    attr(d, "canonical_file") <- normalizePath(path, winslash="/", mustWork=TRUE)
    attr(d, "canonical_md5") <- unname(tools::md5sum(path))
    d
}
.overlay <- function(d, profile, columns=setdiff(intersect(names(d), names(profile)), c("id", "eid", "SEQN", "cycle"))) {
    # columns is native-name -> canonical-profile-name; default names are identical.
    if (is.null(names(columns))) names(columns) <- columns
    if (!all(unname(columns) %in% names(profile)))
        stop("Canonical DOSE profile lacks score fields; rerun 1.disco.sh: ",
             paste(setdiff(unname(columns), names(profile)), collapse=", "))
    ids <- if ("eid" %in% names(d)) as.character(d$eid) else paste(d$cycle, d$SEQN, sep="_")
    j <- match(ids, as.character(profile$eid))
    if (anyNA(j)) stop("Canonical DOSE profile does not cover the native score cache; rerun 1.disco.sh")
    for (nm in names(columns)) d[[nm]] <- profile[[columns[[nm]]]][j]
    attr(d, "canonical_md5") <- attr(profile, "canonical_md5")
    d
}
.require_current <- function(cohort) {
    profile <- .canonical(cohort)
    path <- file.path(.diet_root, "phe", paste0(cohort, ".dose.rds"))
    .require(path, paste0("Run diet.sh ", cohort, " before association."))
    old <- readRDS(path)$data
    if (!identical(attr(old, "canonical_md5"), attr(profile, "canonical_md5")))
        stop("Canonical scores changed after phenotype preparation. Rerun diet.sh ", cohort,
             " before association so that the same score artifact drives every model.")
    invisible(profile)
}
.score_attrs <- function(d) attributes(d)[setdiff(names(attributes(d)), c("names", "row.names", "class"))]
.save_part <- function(d, cohort, definitions, metadata=NULL) {
    if (anyNA(d$eid) || anyDuplicated(d$eid)) stop("Score part has invalid participant IDs")
    attr(d, "cohort") <- cohort
    attr(d, "profile") <- "dose"
    attr(d, "score_definitions") <- definitions
    if (!is.null(metadata)) attr(d, "score_metadata") <- metadata
    attr(d, "score_columns") <- grep("^diet[.].*[.]sum$|^dose(11|9|9half)?$|^(dose(11|9|9half)?|dash|mind|amed)_(day1|two)$", names(d), value=TRUE)
    attr(d, "outcome_filter_applied") <- FALSE
    attr(d, "source_revision") <- "971a61b44b84500ab46e2f5592478505a0eae50e"
    .save(d, .part(cohort))
    message("Saved outcome-free scores: ", .part(cohort), " (", nrow(d), " participants)")
}

disco_ukb <- function() {
    .with_phases(c("qc", "medito", "webq", "patterns", "diet_merge", "upf"),
                 .run("foodwas", "preprocess"))
    source(file.path(dose_helper_dir(), "phenotype.R"), local=globalenv())
    food <- readRDS(dose_rdata_file("diet.fudan"))
    food$eid <- as.character(food$eid)
    field <- as.data.frame(readxl::read_excel(file.path(dose_files_dir(), "foods.xlsx"), sheet="Sheet4"))
    groups <- dose_group_intakes(food, field)
    .save(groups, dose_rdata_file("dose.groups"))
    scoring <- .score_functions()
    score_dir <- file.path(dose_rdata_dir(), "scores")
    result <- scoring$run_dose_scoring(dose_rdata_file("dose.groups"), score_dir)
    patterns <- readRDS(dose_rdata_file("diet"))
    patterns$eid <- as.character(patterns$eid)
    totals <- patterns[, c("eid", grep("\\.(sum|pts|q5)$", names(patterns), value=TRUE)), drop=FALSE]
    # Include all records with a score, including touchscreen-only MEDI records.
    part <- merge(totals, result$primary, by="eid", all=TRUE, sort=FALSE)
    for (version in names(result$scores)) {
        d <- result$scores[[version]]
        columns <- c("eid", setdiff(names(d), names(part)))
        part <- merge(part, d[, columns, drop=FALSE], by="eid", all.x=TRUE, sort=FALSE)
    }
    metadata <- food[, intersect(c("eid", "energy.kJ", "enroll", "start_date"), names(food)), drop=FALSE]
    part <- merge(part, metadata, by="eid", all.x=TRUE, sort=FALSE)
    defs <- list(primary="dose9half", range=c(0,9), step=.5,
        components=attr(result$primary, "food_scoring_rules"),
        comparator_profile="DOSE manuscript WebQ preprocessing and scoring",
        exposure_date="Last valid WebQ recall; touchscreen MEDI is separately named medito")
    metadata <- c(lapply(result$scores, .score_attrs), list(dose=.score_attrs(result$primary)))
    .save_part(part, "ukb", defs, metadata)
}

cohort_ukb <- function() {
    .require(c(dose_rdata_file("diet"), dose_rdata_file("dose.groups"),
               file.path(dose_rdata_dir(), "scores", "dose_score.rds")),
             "Run 1.disco.sh ukb before diet.sh ukb.")
    canonical <- .canonical("ukb")
    cached_patterns <- readRDS(dose_rdata_file("diet"))
    cached_patterns <- .overlay(cached_patterns, canonical)
    .save(cached_patterns, dose_rdata_file("diet"))
    .with_phases("clinical_merge", .run("foodwas", "preprocess"))
    .run("foodwas", "cohort")
    data_dir <- file.path(dose_output_dir(), "ukb/00data")
    base <- readRDS(file.path(data_dir, "dat.base.imp.rds"))
    attr(base, "canonical_md5") <- attr(canonical, "canonical_md5")
    attr(base, "canonical_file") <- attr(canonical, "canonical_file")
    .save(base, file.path(data_dir, "dat.base.imp.rds"))
    cfg <- readRDS(file.path(data_dir, "analysis_config.rds"))
    flow_path <- file.path(dose_rdata_dir(), "flow.csv")
    if (file.exists(flow_path)) {
        flow <- utils::read.csv(flow_path, stringsAsFactors=FALSE)
        flow <- rbind(flow, data.frame(stage="No prevalent stroke and valid follow-up", n=nrow(base)))
        utils::write.csv(flow, file.path(data_dir, "flow.csv"), row.names=FALSE)
    }
    groups <- readRDS(dose_rdata_file("dose.groups"))
    idx <- match(as.character(base$eid), as.character(groups$eid))
    if (anyNA(idx)) stop("DOSE intake cache does not cover the prepared cohort")
    keep <- unique(c("eid", "region", "region_code", "cvd_stroke_i.t2e", "cvd_stroke_i.Yt2e",
                     cfg$covs.m1, cfg$covs.m2.add))
    dat.group <- base[, keep, drop=FALSE]
    dat.group[setdiff(names(groups), "eid")] <- groups[idx, setdiff(names(groups), "eid"), drop=FALSE]
    attr(dat.group, "intake_unit") <- "servings/day"
    group_dir <- file.path(dose_output_dir(), "ukb/02foodwas")
    .save(dat.group, file.path(group_dir, "dat.group.rds"))
    field <- as.data.frame(readxl::read_excel(file.path(dose_files_dir(), "foods.xlsx"), sheet="Sheet4"))
    utils::write.csv(field, file.path(group_dir, "food_group_map.csv"), row.names=FALSE)

    # Subset precomputed fixed scores by ID. No disease-dependent recalibration.
    score_from <- file.path(dose_rdata_dir(), "scores")
    score_to <- file.path(dose_output_dir(), "ukb/04score")
    dir.create(score_to, recursive=TRUE, showWarnings=FALSE)
    subsets <- list()
    for (version in c("dose11", "dose9", "dose9half", "dose")) {
        d <- readRDS(file.path(score_from, paste0(version, "_score.rds")))
        d <- .overlay(d, canonical, setdiff(names(d), "eid"))
        meta <- attr(canonical, "score_metadata")[[version]]
        if (!is.null(meta)) for (nm in names(meta)) attr(d, nm) <- meta[[nm]]
        components <- attr(d, "component_columns")
        total <- if (version == "dose") "dose" else paste0("diet.", version, ".sum")
        if (!isTRUE(all.equal(as.numeric(d[[total]]), as.numeric(rowSums(d[components])), check.attributes=FALSE)))
            stop("Canonical DOSE totals and component scores disagree: ", version)
        j <- match(as.character(base$eid), as.character(d$eid))
        if (anyNA(j)) stop("DOSE score cache is missing cohort participants")
        attrs <- attributes(d)[setdiff(names(attributes(d)), c("names", "row.names", "class"))]
        z <- d[j, , drop=FALSE]
        rownames(z) <- NULL
        for (nm in names(attrs)) attr(z, nm) <- attrs[[nm]]
        .save(z, file.path(score_to, paste0(version, "_score.rds")))
        subsets[[version]] <- z
    }
    for (name in c("DOSE_scoring_rules.csv", "DOSE_food_scoring_rules.csv", "DOSE_version_definitions.csv")) {
        if (!file.copy(file.path(score_from, name), file.path(score_to, name), overwrite=TRUE))
            stop("Cannot copy score definition: ", name)
    }
    scoring <- .score_functions()
    distribution <- do.call(rbind, lapply(subsets[c("dose11", "dose9", "dose9half")], scoring$dose_distribution))
    utils::write.csv(distribution, file.path(score_to, "DOSE_score_distribution.csv"), row.names=FALSE)
    dir.create(file.path(score_to, "figures"), showWarnings=FALSE)
    for (version in c("dose11", "dose9", "dose9half")) {
        scoring$plot_dose_rules(scoring$dose_food_rules(version),
            file.path(score_to, "figures", paste0("Fig_", toupper(version), ".pdf")))
    }
    .save(base, file.path(dose_output_dir(), "ukb", "phenotypes.rds"))
    .save(list(data=base, profile="dose", id="eid"), file.path(.diet_root, "phe", "ukb.dose.rds"))
    message("Prepared DOSE cohort and selected existing scores: ", nrow(base), " participants")
}

disco_nhanes <- function() {
    nh <- .run("nhanes", "score")
    .run("nhanes", "patterns", nh)
    patterns <- nh$nh_patterns()$scores
    cache <- file.path(dose_output_dir(), ".cache/nhanes")
    one <- readRDS(file.path(cache, "nhanes_day1_score.rds"))
    two <- readRDS(file.path(cache, "nhanes_two_day_score.rds"))
    dash <- readRDS(file.path(cache, "nhanes_dash_score.rds"))
    score_fields <- setdiff(names(nh$nh_scoring$calculate_dose(
        structure(one[1, c("eid", nh$nh_scoring$dose_food_rules()$Food)], intake_unit="servings/day"))), "eid")
    fields <- c(score_fields, "energy_kJ", "valid")
    part <- one[, c("eid", "SEQN", "cycle", fields), drop=FALSE]
    names(part)[match(fields, names(part))] <- paste0(fields, "_day1")
    j <- match(part$eid, two$eid)
    for (nm in fields) part[[paste0(nm, "_two")]] <- two[[nm]][j]
    k <- match(paste(part$cycle, part$SEQN), paste(dash$cycle, dash$SEQN))
    for (nm in c("dash_day1", "dash_two")) part[[nm]] <- dash[[nm]][k]
    k <- match(part$eid, patterns$eid)
    for (nm in c("mind_day1", "mind_two", "amed_day1", "amed_two")) part[[nm]] <- patterns[[nm]][k]
    # The pre-specified main recall mode is two-day mean in 2003-2018.
    # Day-one values remain explicit and are never substituted for missing day two.
    part$dose <- part$dose_two
    part$dose9half <- part$dose
    part$diet.dose.sum <- part$dose
    part$diet.dose9half.sum <- part$dose9half
    part$diet.dash.sum <- part$dash_two
    part$diet.mind.sum <- part$mind_two
    part$diet.amed.sum <- part$amed_two
    defs <- list(primary="dose9half", range=c(0,9), step=.5,
        recall_mode="two_day_2003_2018", day1_available="1999-2018", calculated_fields=score_fields,
        dose_aggregation="score(mean of reliable day-one and day-two intakes)",
        dash_aggregation="dietaryindex native mean of two daily scores",
        med_label="aMED (0-9), distinct from UKB MedDiet medi24",
        mind_label="MIND-adapted")
    .save_part(part, "nhanes", defs)
}

cohort_nhanes <- function() {
    cache <- file.path(dose_output_dir(), ".cache/nhanes")
    .require(file.path(cache, c("nhanes_day1_score.rds", "nhanes_two_day_score.rds", "nhanes_patterns.rds")),
             "Run 1.disco.sh nhanes before diet.sh nhanes.")
    canonical <- .canonical("nhanes")
    score_fields <- attr(canonical, "score_definitions")$calculated_fields
    if (is.null(score_fields)) stop("Rerun 1.disco.sh nhanes: canonical calculated score fields are missing")
    for (recall in c("day1", "two")) {
        filename <- if (recall == "day1") "nhanes_day1_score.rds" else "nhanes_two_day_score.rds"
        d <- readRDS(file.path(cache, filename))
        fields <- c(score_fields, "energy_kJ", "valid")
        d <- .overlay(d, canonical, setNames(paste0(fields, "_", recall), fields))
        d$diet.dose.sum <- d$dose
        if (!identical(d$dose, d$dose9half)) stop("Canonical NHANES primary score aliases disagree")
        .save(d, file.path(cache, filename))
    }
    dash <- readRDS(file.path(cache, "nhanes_dash_score.rds"))
    dash <- .overlay(dash, canonical, c("dash_day1", "dash_two"))
    .save(dash, file.path(cache, "nhanes_dash_score.rds"))
    patterns <- readRDS(file.path(cache, "nhanes_patterns.rds"))
    patterns$scores <- .overlay(patterns$scores, canonical, c("mind_day1", "mind_two", "amed_day1", "amed_two"))
    .save(patterns, file.path(cache, "nhanes_patterns.rds"))
    nh <- .nh_functions()
    cov <- dplyr::bind_rows(lapply(nh$NH_CYCLES$cycle, nh$nh_covariates))
    if (anyDuplicated(cov[c("SEQN", "cycle")])) stop("Duplicate NHANES clinical IDs")
    .save(cov, file.path(cache, "nhanes_covariates.rds"))
    one <- readRDS(file.path(cache, "nhanes_day1_score.rds"))
    phe <- dplyr::left_join(one, cov, by=c("SEQN", "cycle"))
    attr(phe, "canonical_md5") <- attr(canonical, "canonical_md5")
    attr(phe, "canonical_file") <- attr(canonical, "canonical_file")
    .save(phe, file.path(dose_output_dir(), "nhanes", "phenotypes.rds"))
    .save(list(data=phe, profile="dose", id="eid"), file.path(.diet_root, "phe", "nhanes.dose.rds"))
    message("Prepared NHANES clinical variables separately from dietary scores")
}

export_effects <- function(cohort) {
    base <- dose_output_dir()
    tables <- list()
    if (cohort == "ukb") {
        manuscript <- "ukb/04score/pattern/pattern_continuous_10patterns.csv"
        main <- if (file.exists(file.path(base, manuscript)))
            c(stroke_manuscript_10=manuscript) else
            c(stroke_all_patterns_14="ukb/04score/pattern/pattern_continuous.csv")
        sources <- c(main,
                     outcome_wide="ukb/04score/otherOutcomes/DOSE_other_continuous.csv")
        for (family in names(sources)) {
            path <- file.path(base, sources[[family]])
            if (!file.exists(path)) next
            d <- utils::read.csv(path, check.names=FALSE)
            if (!all(c("beta", "se") %in% names(d))) stop("Rerun DOSE association to export actual coefficient SE: ", path)
            tables[[family]] <- data.frame(project="dose", cohort="ukb", family=family,
                score=d$diet, outcome=d$outcome, model=d$model, contrast=d$scale,
                beta=d$beta, se=d$se, effect="HR", estimate=d$HR, lower=d$LCI,
                upper=d$HCI, p=d$P, fdr=if ("FDR" %in% names(d)) d$FDR else NA_real_,
                n=d$N, events=d$Events, increment=d$increment, df=Inf,
                ci_method="normal", source=sources[[family]])
        }
    } else {
        path <- file.path(base, "nhanes/NHANES_M2_results.csv")
        if (file.exists(path)) {
            d <- utils::read.csv(path, check.names=FALSE)
            d <- d[d$Contrast %in% c("Per point", "Per SD"), , drop=FALSE]
            if (!all(c("beta", "se") %in% names(d))) stop("Rerun DOSE NHANES association to export actual survey SE")
            tables$stroke <- data.frame(project="dose", cohort="nhanes", family=d$Analysis,
                score=d$Score, outcome="prevalent_stroke", model=d$Model, contrast=d$Contrast,
                beta=d$beta, se=d$se, effect="OR", estimate=d$OR, lower=d$LCI,
                upper=d$UCI, p=d$p, fdr=NA_real_, n=d$N, events=d$Cases,
                increment=d$Increment, df=d$df, ci_method="survey t",
                source="nhanes/NHANES_M2_results.csv")
        }
    }
    if (!length(tables)) stop("No DOSE association estimates are available for ", cohort)
    out <- do.call(rbind, tables)
    rownames(out) <- NULL
    path <- file.path(base, cohort, "effects.tsv")
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
    utils::write.table(out, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
    message("Saved actual coefficients and SEs: ", path)
    invisible(out)
}

assoc_ukb <- function() {
    .require_current("ukb")
    .require(file.path(dose_output_dir(), c("ukb/00data/dat.base.imp.rds", "ukb/04score/dose_score.rds")),
             "Run diet.sh ukb before diet.sh assoc ukb.")
    steps <- Sys.getenv("DOSE_ASSOC_STEPS", "distribution,stroke,sensitivity,plot,outcomes")
    for (step in strsplit(steps, ",", fixed=TRUE)[[1L]]) .run("score", trimws(step))
    export_effects("ukb")
}
assoc_nhanes <- function() {
    .require_current("nhanes")
    .require(file.path(dose_output_dir(), ".cache/nhanes/nhanes_covariates.rds"),
             "Run diet.sh nhanes before diet.sh assoc nhanes.")
    for (step in c("stroke", "plot", "flow")) .run("nhanes", step)
    export_effects("nhanes")
}

if (.action %in% c("disco_ukb", "cohort_ukb", "disco_nhanes", "cohort_nhanes", "assoc_ukb", "assoc_nhanes")) {
    get(.action, mode="function")()
} else if (.action == "effects") {
    cohort <- Sys.getenv("DIET_COHORT", "ukb")
    if (!cohort %in% c("ukb", "nhanes")) stop("DOSE effects are available for ukb or nhanes")
    export_effects(cohort)
} else if (.action == "step") {
    module <- Sys.getenv("DOSE_MODULE", "")
    step <- Sys.getenv("DOSE_STEP", "")
    if (!nzchar(module) || !nzchar(step)) stop("DOSE_MODULE and DOSE_STEP are required")
    if (module == "nhanes" && step == "patterns") {
        nh <- .nh_functions(); nh$nh_patterns()
    } else .run(module, step)
} else stop("Unknown DOSE_ACTION: ", .action)
