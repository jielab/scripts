# Shared provenance and supplemental-data contracts. No participant data are shipped.
BMONTH_CODE_VERSION <- "2026-10-20-exposure-focus"
bm_truth <- function(x) toupper(x) %in% c("1", "TRUE", "YES")

phewas_input_path <- function(validate=TRUE) {
  explicit <- Sys.getenv("BMONTH_PHEWAS_INPUT", "")
  candidates <- c(file.path(outdir,"phewas","PheWAS.rds"),file.path(indir,"Rdata","PheWAS.rds"))
  path <- if(nzchar(explicit)) explicit else {
    found <- candidates[file.exists(candidates)]
    if(!length(found)) stop("No existing phenotype panel. Missing: ",paste(candidates,collapse="; "),
      "; optional rebuild: ./bmonth.sh rebuild --dir-in ",shQuote(indir)," --dir-out ",shQuote(outdir))
    found[[1]]
  }
  if (!file.exists(path)) stop("Specified phenotype panel does not exist: ", path)
  selected <- Sys.getenv("BMONTH_PHEWAS_SELECTED_INPUT", "")
  if(nzchar(selected) && !identical(normalizePath(path),normalizePath(selected,mustWork=FALSE)))
    stop("Shell/R phenotype selection changed; restart the stage: ",selected," -> ",path)
  if (validate) invisible(phewas_provenance(path))
  path
}

phewas_provenance <- local({
  cached <- new.env(parent=emptyenv())
  function(input=phewas_input_path(validate=FALSE)) {
    manifest_path <- Sys.getenv("BMONTH_PHEWAS_MANIFEST", "")
    if (!nzchar(manifest_path)) manifest_path <- sub("\\.rds$", ".manifest.rds", input)
    if (identical(manifest_path, input)) manifest_path <- paste0(input, ".manifest.rds")
    meta <- if (file.exists(manifest_path)) tryCatch(readRDS(manifest_path), error=function(e) NULL) else NULL
    sidecar_names <- c("PheWAS.dates.rds", "PheWAS.roster.rds", "PheWAS.uncertainty.rds")
    sidecar_paths <- file.path(dirname(input), sidecar_names)
    # Independent validation hashes real contents; pass this provenance explicitly
    # to immediately following consumers to avoid redundant scans.
    # mtime/size alone cannot detect a same-size replacement with preserved dates.
    sidecar_hashes <- rep(NA_character_,length(sidecar_paths))
    present <- file.exists(sidecar_paths)
    sidecar_hashes[present] <- unname(tools::md5sum(sidecar_paths[present]))
    panel_hash <- unname(tools::md5sum(input))
    manifest_hash <- if(file.exists(manifest_path)) unname(tools::md5sum(manifest_path)) else "missing"
    key <- paste(normalizePath(input, mustWork=TRUE), panel_hash, manifest_hash,
                 paste(sidecar_hashes, collapse="|"), sep="|")
    if (exists(key, envir=cached, inherits=FALSE)) return(get(key, envir=cached))
    # The existing isolated packer validates the wide table and publishes counts.
    # The parent reads only metadata; it never deserializes the original panel.
    if(!exists("bm_panel_isolated",mode="function"))
      source(file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R"))
    panel_cache <- file.path(rawdir,"phewas_panel_cache",bm_hash(panel_hash))
    packed <- bm_panel_isolated(input,panel_cache,
      file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R"),input_md5=panel_hash)
    counts <- packed$counts; n_participants <- packed$n
    errors <- character()
    if (!is.list(meta)) errors <- "missing or unreadable PheWAS.manifest.rds"
    else {
      if (!identical(meta$format_version, 2L)) errors <- c(errors, "unsupported manifest version")
      if (!identical(meta$vocabulary_id, "ICD10")) errors <- c(errors, "vocabulary is not UK ICD10")
      if (!identical(meta$mapping_name, "phecode_map_icd10")) errors <- c(errors, "wrong mapping table")
      for (field in c("roster_validated", "control_exclusions_applied", "sex_restrictions_applied"))
        if (!isTRUE(meta[[field]])) errors <- c(errors, paste("unverified", field))
      if (!identical(meta$panel_md5, panel_hash)) errors <- c(errors, "panel hash does not match manifest")
      if (!is.character(meta$coverage_definition) || length(meta$coverage_definition)!=1L ||
          is.na(meta$coverage_definition) || !nzchar(trimws(meta$coverage_definition)))
        errors <- c(errors, "missing coverage definition")
      expected <- meta$sidecar_md5
      # The builder records original absolute paths; resolve by basename so the
      # whole verified panel bundle can be moved without accepting unrelated files.
      if (!is.character(expected) || is.null(names(expected)) ||
          anyDuplicated(basename(names(expected))) ||
          !all(sidecar_names %in% basename(names(expected)))) {
        errors <- c(errors, "missing or invalid sidecar hashes")
      } else {
        expected <- unname(expected[match(sidecar_names, basename(names(expected)))])
        bad <- is.na(sidecar_hashes) | is.na(expected) | sidecar_hashes!=expected
        if (any(bad)) errors <- c(errors, paste("missing or mismatched sidecar:", paste(sidecar_names[bad],collapse=", ")))
      }
    }
    # Missing/old/mismatched external provenance is descriptive, not a gate.
    # Do not pass untrusted metadata to consumers as validated source evidence.
    out <- list(validation_status=if(length(errors)) "existing_panel" else "validated_ICD10_v2",
                input_path=normalizePath(input,mustWork=TRUE),
                manifest_path=manifest_path, panel_md5=panel_hash, panel_cache=panel_cache,
                issues=errors, manifest=if(length(errors)) NULL else meta,
                n_participants=n_participants,phenotype_counts=counts,
                source_note=if(length(errors)) "Using existing phenotype file; optional provenance not verified"
                  else "Panel and optional ICD10 provenance bundle agree")
    message("[PheWAS] panel: ",out$input_path,"; MD5: ",panel_hash,"; provenance: ",out$validation_status)
    assign(key, out, envir=cached); out
  }
})

# Only use with provenance just validated in this process. Its packer already
# checked the binary hash; do not scan it again for each exposure/plot.
phewas_packed_panel <- function(provenance) {
  panel <- bm_cached_panel(provenance$panel_cache,provenance$panel_md5,verify_binary=FALSE)
  if(is.null(panel)) stop("Packed phenotype metadata changed; restart this stage")
  panel
}

merge_supplemental_covariates <- function(dat) {
  path <- Sys.getenv("BMONTH_SUPPLEMENTAL_COVARIATES", "")
  if (!nzchar(path)) return(dat)
  if (!file.exists(path)) stop("Missing supplemental covariates RDS: ", path)
  extra <- validate_ids(readRDS(path), "eid", "supplemental covariates")
  reserved<-unique(c(BASE_EXPOSURES,EXPLORATORY_ALL,"birth_year_bm","age_fractional",
    "assessment_year","assessment_month","adopted_num","sibling_n","sibship_size",
    "sibship_size_cap","sibship_size_cat","birth_rank","impossible_rank","eligible_border",
    "move","move_km","tdi.birth","tdi.birth.nn","lon.birth","lat.birth","lon.home","lat.home",
    "bmonth_cos_x_north","bmonth_sin_x_north","bmonth_cos_x_east","bmonth_sin_x_east",
    "bmonth_cos_x_move","bmonth_sin_x_move"))
  collision<-intersect(setdiff(names(extra),"eid"),reserved)
  if(length(collision))stop("Supplemental columns use reserved derived-field names: ",paste(collision,collapse=", "),
    "; supply explicitly named raw covariates instead")
  overlap <- intersect(setdiff(names(extra), "eid"), names(dat))
  if (length(overlap)) stop("Supplemental columns must not silently replace cohort columns: ",
                             paste(overlap, collapse=", "))
  dplyr::left_join(dat, extra, by="eid")
}

# This is age calculated from the existing month-15 birth-date proxy, not exact DOB.
add_retained_covariates <- function(dat) {
  if (all(c("birth_date", "date_attend") %in% names(dat))) {
    dat$age_fractional <- as.numeric(as_date(dat$date_attend)-as_date(dat$birth_date))/365.25
    dat$age_fractional[!is.finite(dat$age_fractional) | dat$age_fractional<0 | dat$age_fractional>110] <- NA_real_
  }
  if ("date_attend" %in% names(dat)) {
    dat$assessment_year <- as_num(format(as_date(dat$date_attend), "%Y"))
    dat$assessment_month <- as_num(format(as_date(dat$date_attend), "%m"))
  }
  dat
}

bm_planned_adjust <- function(p, eligible=rep(TRUE, length(p)), method="BH") {
  if (length(p)!=length(eligible)) stop("Eligibility and P values must have the same length")
  eligible <- !is.na(eligible) & eligible
  answer <- rep(NA_real_, length(p))
  pp <- p[eligible]; valid <- is.finite(pp) & pp>=0 & pp<=1
  pp[!valid] <- 1
  adjusted <- p.adjust(pp, method)
  answer[which(eligible)[valid]] <- adjusted[valid]
  answer
}

bm_write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile(".csv-", tmpdir=dirname(path)); on.exit(unlink(tmp))
  write.csv(x, tmp, row.names=FALSE, na="")
  if(!file.rename(tmp, path)) stop("Cannot publish CSV: ", path)
  invisible(path)
}
