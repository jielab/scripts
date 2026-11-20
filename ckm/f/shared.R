# Reuse the repository's current phenotype helpers in a private namespace.
# The shared source owns date_follow_end; CKM never carries a second date literal.
.ckm_repository <- new.env(parent = emptyenv())

ckm_load_shared <- function(file) {
  if (!file.exists(file)) stop("Missing shared phenotype helper: ", file,
      ". Place ckm/ beside the existing scripts/0f/ directory.")
  file <- normalizePath(file, winslash = "/", mustWork = TRUE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  # The shared helper calls set.seed(); preserve the caller's random stream.
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  shared <- new.env(parent = globalenv())
  base::source(file, local = shared, chdir = FALSE)
  needed <- c("date_follow_end", "vars.basic", ".phenotype_date", "t2e")
  missing <- needed[!vapply(needed, exists, logical(1), envir = shared, inherits = FALSE)]
  if (length(missing)) stop("Shared phenotype helper lacks: ", paste(missing, collapse = ", "))
  cutoff <- shared$date_follow_end
  if (!inherits(cutoff, "Date") || length(cutoff) != 1L || is.na(cutoff) || !is.finite(as.numeric(cutoff)))
    stop("0f/phenotype.R date_follow_end must be one valid Date")
  if (!is.character(shared$vars.basic) || !length(shared$vars.basic) ||
      anyNA(shared$vars.basic) || any(!nzchar(shared$vars.basic)) || anyDuplicated(shared$vars.basic))
    stop("0f/phenotype.R vars.basic must contain unique covariate names")
  .ckm_repository$shared <- shared
  .ckm_repository$file <- file
  .ckm_repository$md5 <- unname(tools::md5sum(file))
  invisible(shared)
}

ckm_shared <- function() {
  if (is.null(.ckm_repository$shared))
    stop("Load the sibling 0f/phenotype.R with ckm_load_shared() before CKM analysis")
  .ckm_repository$shared
}

ckm_shared_info <- function() {
  shared <- ckm_shared()
  list(file = .ckm_repository$file, md5 = .ckm_repository$md5,
    date_follow_end = as.character(shared$date_follow_end), vars.basic = shared$vars.basic,
    reused = c("date_follow_end", "vars.basic", ".phenotype_date", "t2e"))
}
