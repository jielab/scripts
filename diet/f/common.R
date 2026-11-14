# Shared data contracts. Scientific scoring/model definitions live in their engines.
diet_env <- function(key, default = "", required = FALSE) {
    value <- Sys.getenv(key, default)
    if (required && !nzchar(value)) stop("Missing environment setting: ", key)
    value
}
diet_dir <- function(path) { dir.create(path, recursive=TRUE, showWarnings=FALSE); path }
diet_save <- function(x, path) {
    diet_dir(dirname(path)); tmp <- tempfile(".diet-", tmpdir=dirname(path))
    on.exit(unlink(tmp), add=TRUE); saveRDS(x, tmp, compress=TRUE)
    if (!file.rename(tmp, path)) stop("Cannot publish ", path)
    invisible(path)
}
diet_tsv <- function(x, path) {
    diet_dir(dirname(path)); write.table(x, path, sep="\t", quote=TRUE,
        row.names=FALSE, na="NA", fileEncoding="UTF-8")
    invisible(path)
}
diet_read <- function(path) {
    if (!file.exists(path)) stop("Missing input: ", path)
    if (grepl("\\.rds$", path, ignore.case=TRUE)) return(readRDS(path))
    reader <- if(grepl("\\.csv(\\.gz)?$", path, ignore.case=TRUE)) read.csv else read.delim
    header <- names(reader(path, nrows=0L, check.names=FALSE))
    classes <- setNames(rep(NA_character_,length(header)),header)
    classes[intersect(header,c("id","ID","eid","iid","SEQN","IDIND","participant_id"))] <- "character"
    reader(path, check.names=FALSE, stringsAsFactors=FALSE, colClasses=classes)
}
diet_id <- function(x, cohort) {
    candidates <- switch(cohort, ukb=c("eid","id"), nhanes=c("SEQN","id"),
                         chns=c("IDIND","id"), charls=c("id","ID"))
    found <- intersect(candidates, names(x))
    if (!length(found)) stop("Missing participant ID for ", cohort)
    ids <- trimws(as.character(x[[found[1L]]]))
    if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
        stop(cohort, ": IDs must be nonempty and unique; resolve repeated waves explicitly")
    ids
}
diet_score_columns <- function(x) {
    declared <- attr(x, "score_columns", exact=TRUE)
    if (!is.null(declared)) return(intersect(declared, names(x)))
    grep("(^|\\.)diet\\..*\\.sum$|^(DOSE9half|DOSE11|DOSE9|DOSE|dose9half|dose11|dose9|dose|MAHA|DASH|MIND|aMED|MEDI)$",
         names(x), value=TRUE)
}
diet_publish_disco <- function(cohort, output, parts, selection=diet_env("DIET_SCORE_SELECTION", "all")) {
    present <- parts[file.exists(parts)]
    if (!length(present)) stop("No completed score parts for ", cohort, "; run 1.disco.sh first")
    out <- NULL; catalog <- list(); profiles <- list()
    for (profile in names(present)) {
        x <- readRDS(present[[profile]])
        if (!is.data.frame(x)) stop("Score part must be a data.frame: ", present[[profile]])
        ids <- diet_id(x, cohort); scores <- diet_score_columns(x)
        source_attrs <- attributes(x)[setdiff(names(attributes(x)), c("names","row.names","class"))]
        x$id <- ids
        native <- switch(cohort, ukb="eid", nhanes="SEQN", chns="IDIND", charls="ID")
        x[native] <- NULL
        original <- names(x)
        collide <- if (is.null(out)) character() else intersect(setdiff(names(x), "id"), names(out))
        names(x)[names(x) %in% collide] <- paste0(profile, ".", collide)
        mapped <- setNames(names(x), original)
        if (is.null(out)) out <- x else out <- merge(out, x, by="id", all=TRUE, sort=FALSE)
        for (score in scores) {
            v <- suppressWarnings(as.numeric(x[[mapped[[score]]]])); ok <- is.finite(v)
            q <- if (any(ok)) quantile(v[ok], c(.1,.5,.9), names=FALSE) else rep(NA_real_,3)
            catalog[[length(catalog)+1L]] <- data.frame(cohort=cohort, profile=profile,
                variable=mapped[[score]], native_variable=score, n=nrow(x), n_scored=sum(ok),
                missing=sum(!ok), min=if(any(ok)) min(v[ok]) else NA_real_,
                q10=q[1], median=q[2], q90=q[3], max=if(any(ok)) max(v[ok]) else NA_real_,
                definition=paste0(cohort,"/",profile,"/",score), stringsAsFactors=FALSE)
        }
        profiles[[profile]] <- list(source=present[[profile]], attributes=source_attrs,
            column_map=mapped, participant_ids=ids)
    }
    # Native aliases support each cohort's established analysis code.
    native <- switch(cohort, ukb="eid", nhanes="SEQN", chns="IDIND", charls="ID")
    out[[native]] <- out$id
    out <- out[c("id", native, setdiff(names(out), c("id",native)))]
    cat <- if (length(catalog)) do.call(rbind, catalog) else data.frame()
    attr(out,"schema_version") <- 1L
    attr(out,"cohort") <- cohort
    attr(out,"profiles") <- profiles
    attr(out,"score_columns") <- if (nrow(cat)) unique(cat$variable) else character()
    attr(out,"score_catalog") <- cat
    attr(out,"source_commit") <- "971a61b44b84500ab46e2f5592478505a0eae50e"
    foods <- diet_env("DIET_FOODS")
    attr(out,"foods_md5") <- if(file.exists(foods)) unname(tools::md5sum(foods)) else NA_character_
    attr(out,"created_at") <- format(Sys.time(), tz="UTC", usetz=TRUE)
    # Keep original names, components, intake/recall metadata and profile maps in
    # the raw artifact. Equal score names across engines need not mean equal rules.
    raw <- file.path(output,"disco",paste0(cohort,".raw.rds"))
    diet_save(out,raw)
    totals <- grep("(^|\\.)diet\\.[^.]+\\.sum$", names(out), value=TRUE)
    score_names <- sub("\\.sum$", "", sub("^.*diet\\.", "", totals))
    requested <- if(identical(selection,"all")) unique(score_names) else unique(strsplit(selection,",",fixed=TRUE)[[1L]])
    selected <- totals[score_names %in% requested]
    usable <- vapply(selected, function(nm) any(is.finite(suppressWarnings(as.numeric(out[[nm]])))), logical(1))
    for(nm in selected[!usable]) message("SKIP ",cohort,": ",nm,": no finite score values (missing inputs or components)")
    selected <- selected[usable]
    available <- sub("\\.sum$", "", sub("^.*diet\\.", "", selected))
    for(nm in setdiff(requested,unique(available)))
        message("SKIP ",cohort,": ",nm,": unavailable in the selected cohort/profile outputs")
    compact <- out[c("id",native,selected)]
    # Full reconstruction metadata belongs to raw, not the slim score table.
    attr(compact,"profiles") <- NULL
    attr(compact,"schema_version") <- 2L
    attr(compact,"raw_file") <- normalizePath(raw,winslash="/",mustWork=TRUE)
    attr(compact,"raw_md5") <- unname(tools::md5sum(raw))
    attr(compact,"requested_scores") <- requested
    attr(compact,"score_columns") <- selected
    selected_catalog <- if(nrow(cat)) cat[cat$variable %in% selected,,drop=FALSE] else cat
    attr(compact,"score_catalog") <- selected_catalog
    dest <- file.path(output,"disco",paste0(cohort,".disco.rds"))
    diet_save(compact,dest)
    diet_tsv(selected_catalog,file.path(output,"disco",paste0(cohort,".scores.tsv")))
    message("Saved ",raw," (full scoring data)")
    message("Saved ",dest," (",nrow(compact)," participants; ",length(selected)," selected score totals)")
    invisible(compact)
}
diet_profile <- function(x, profile, cohort) {
    # Native analyses need components and their original source populations.
    # Follow the paired raw artifact and check it, rather than silently using a
    # stale engine cache after the public score table has been regenerated.
    raw <- attr(x,"raw_file",exact=TRUE)
    if(!is.null(raw)) {
        if(!file.exists(raw) || !identical(unname(tools::md5sum(raw)),attr(x,"raw_md5",exact=TRUE)))
            stop("Missing or changed raw score artifact: ",raw,"; rerun 1.disco.sh ",cohort)
        compact <- x
        x <- readRDS(raw)
        if(!identical(x$id,compact$id)) stop("Raw and compact score participant IDs disagree")
        for(nm in diet_score_columns(compact)) x[[nm]] <- compact[[nm]]
    }
    meta <- attr(x,"profiles",exact=TRUE)[[profile]]
    if(is.null(meta) || is.null(meta$column_map) || is.null(meta$participant_ids))
        stop("Canonical ",cohort," scores do not contain profile ",profile,
             "; rerun 1.disco.sh ",cohort," --score all")
    map <- meta$column_map
    if(!all(unname(map) %in% names(x))) stop("Canonical score profile has missing columns")
    idx <- match(as.character(meta$participant_ids),as.character(x$id))
    if(anyNA(idx)) stop("Canonical score profile is missing participant rows")
    y <- x[idx,unname(map),drop=FALSE]; names(y) <- names(map)
    native <- switch(cohort,ukb="eid",nhanes="SEQN",chns="IDIND",charls="ID")
    y[[native]] <- y$id
    for(nm in names(meta$attributes)) attr(y,nm) <- meta$attributes[[nm]]
    y
}
diet_bundle_phe <- function(cohort, output) {
    files <- file.path(output,"phe",paste0(cohort,c(".maha.rds",".dose.rds",".charls.rds")))
    names(files) <- c("maha","dose","charls"); files <- files[file.exists(files)]
    if(!length(files)) stop("No prepared phenotype parts for ",cohort)
    x <- lapply(files, readRDS)
    attr(x,"cohort") <- cohort
    attr(x,"schema_version") <- 1L
    attr(x,"note") <- "Profiles retain their own eligibility and time origins. Do not pool rows across profiles."
    diet_save(x,file.path(output,"phe",paste0(cohort,".phe.rds")))
}
if (sys.nframe()==0L) {
    args <- commandArgs(trailingOnly=TRUE)
    if(length(args)!=2L) stop("Usage: Rscript common.R disco|phe COHORT")
    output <- diet_env("DIET_OUTDIR",required=TRUE); cohort <- args[2]
    if(args[1]=="disco") {
        parts <- file.path(output,".cache","disco",paste0(cohort,c(".maha.rds",".dose.rds",".charls.rds")))
        names(parts)<-c("maha","dose","charls")
        selected_profiles <- diet_env("DIET_SCORE_PROFILES")
        if(nzchar(selected_profiles)) {
            profiles <- strsplit(selected_profiles,",",fixed=TRUE)[[1L]]
            parts <- parts[profiles]
            if(anyNA(parts) || !all(file.exists(parts))) stop("Missing selected scoring engine output")
        }
        diet_publish_disco(cohort,output,parts)
    } else if(args[1]=="phe") diet_bundle_phe(cohort,output) else stop("Unknown action: ",args[1])
}
