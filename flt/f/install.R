# Explicit setup for flt.sh run. No scientific analyses or bulk package updates.
source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(args, "--check"))) stop("Unknown installer argument")
check <- "--check" %in% args
cran <- c("tidyverse", "data.table", "survival", "readxl", "ggrepel", "lubridate",
           "broom", "R.utils", "openxlsx", "MendelianRandomization", "robustbase", "ieugwasr",
          "patchwork", "png", "deTS")
bioc <- c("limma", "impute")
universe <- c("TwoSampleMR", "MVMR", "MRPRESSO")
github <- c(HDL = "zhenin/HDL/HDL", GenomicSEM = "GenomicSEM/GenomicSEM")
if (Sys.getenv("FLT_RUN_HDL", "TRUE") == "FALSE") github <- github[names(github) != "HDL"]
required <- unique(c(cran, bioc, universe, names(github)))
missing <- function(packages) packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
cat("R: ", R.home(), "\nLibraries: ", paste(.libPaths(), collapse = ":"), "\n", sep = "")
absent <- missing(required)
if (length(absent)) cat("Missing/unloadable: ", paste(absent, collapse = ", "), "\n", sep = "")
if (!check && length(absent)) {
    lib <- Sys.getenv("FLT_R_LIBRARY")
    if (!nzchar(lib)) stop("Use ./install.sh run to set FLT_R_LIBRARY")
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(lib) || file.access(lib, 2L) != 0L) stop("R library is not writable: ", lib)
    .libPaths(c(lib, .libPaths()))
    system_prefix <- Sys.getenv("DIET_SYSTEM_PREFIX")
    if (nzchar(system_prefix) && dir.exists(file.path(system_prefix, "include"))) {
        # Configure scripts such as gmp also need non-pkg-config header paths.
        for (key in c("CPATH", "LIBRARY_PATH")) {
            path <- file.path(system_prefix, if (key == "CPATH") "include" else "lib")
            if (key == "CPATH") path <- c(path, file.path(system_prefix, "include/freetype2"))
            old <- Sys.getenv(key)
            do.call(Sys.setenv, setNames(list(paste(c(path, old[nzchar(old)]), collapse=":")), key))
        }
    }
    offline <- file.path(Sys.getenv("FLT_RESOURCE_DIR"), "R-repo")
    if (Sys.getenv("submit") == "hpc" && file.exists(file.path(offline, "src/contrib/PACKAGES"))) {
        # Package startup notices may try the internet during R CMD INSTALL's
        # load tests; computation nodes are offline, so bound those requests.
        Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT="10")
        install.packages(absent, lib=lib, repos=paste0("file://", offline), type="source")
        absent <- missing(required)
        if (length(absent)) stop("Offline installation failed: ", paste(absent, collapse=", "))
    }
    cran_url <- Sys.getenv("FLT_CRAN", if (Sys.getenv("submit") == "hpc")
        "https://mirrors.tuna.tsinghua.edu.cn/CRAN" else "https://cloud.r-project.org")
    options(repos = c(CRAN = cran_url), timeout = 120)
    if (Sys.getenv("submit") == "hpc") options(BioC_mirror =
        Sys.getenv("FLT_BIOC_MIRROR", "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"))
    need <- missing(cran)
    if (length(need)) install.packages(need, lib = lib)
    need <- missing(bioc)
    if (length(need)) {
        if (length(missing("BiocManager"))) install.packages("BiocManager", lib = lib)
        BiocManager::install(need, lib = lib, ask = FALSE, update = FALSE)
    }
    need <- missing(universe)
    if (length(need)) install.packages(need, lib = lib,
        repos = c("https://mrcieu.r-universe.dev", getOption("repos")))
    if (length(missing(names(github)))) {
        if (length(missing("remotes"))) install.packages("remotes", lib = lib)
        for (pkg in names(github)) if (length(missing(pkg)))
            remotes::install_github(github[[pkg]], lib = lib, upgrade = "never",
                                    dependencies = NA, build_vignettes = FALSE)
    }
    absent <- missing(required)
}
if (length(absent)) stop("Missing/unloadable: ", paste(absent, collapse = ", "),
                        ". Run ./install.sh run with the same runtime/library.", call. = FALSE)
bwmr <- Sys.getenv("FLT_BWMR")
if (!nzchar(bwmr)) stop("FLT_BWMR is not set")
if (!file.exists(bwmr) && !check) {
    dir.create(dirname(bwmr), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(tmpdir = dirname(bwmr))
    url <- "https://raw.githubusercontent.com/jiazhao97/BWMR/master/updates/BWMR_updated.R"
    download.file(url, tmp, mode = "wb")
    invisible(parse(tmp))
    if (!file.rename(tmp, bwmr)) stop("Cannot save BWMR: ", bwmr)
    writeLines(c(url, paste("MD5", unname(tools::md5sum(bwmr)))), paste0(bwmr, ".source.txt"))
}
if (!file.exists(bwmr)) stop("Missing BWMR source; run ./install.sh run")
env <- new.env(parent = globalenv())
sys.source(bwmr, envir = env)
if (!is.function(env$BWMR)) stop("Invalid BWMR source: missing BWMR()")
credential <- file.path(Sys.getenv("USERPROFILE", unset = path.expand("~")), ".config/opengwas/.Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(credential)) readRenviron(credential)
# Fresh runs have no OpenGWAS cache. Check presence without exposing the secret.
if (tolower(Sys.getenv("FLT_PHEWAS_ONLINE", "true")) == "true" &&
    !nzchar(ieugwasr::get_opengwas_jwt()) && !flt_phewas_cached())
    stop("Fresh PheWAS needs OPENGWAS_JWT or ~/.config/opengwas/.Renviron. No token was printed.")
versions <- data.frame(package = required,
    version = vapply(required, function(x) as.character(packageVersion(x)), character(1)))
if (!check) {
    dir.create(Sys.getenv("FLT_R_LIBRARY"), recursive = TRUE, showWarnings = FALSE)
    write.csv(versions, file.path(Sys.getenv("FLT_R_LIBRARY"), "flt-packages.csv"), row.names = FALSE)
}
cat("Analysis R dependencies and BWMR passed. Online token validity is checked by OpenGWAS during queries.\n")
