# Download source packages on the HPC login node; compilation remains in LSF.
source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
repo <- file.path(Sys.getenv("FLT_RESOURCE_DIR"), "R-repo")
dest <- file.path(repo, "src/contrib")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
options(timeout = 120, download.file.method="curl",
        download.file.extra="--fail --location --connect-timeout 15 --max-time 120")
# Use the system TLS stack for downloads, outside the conda R library path.
Sys.unsetenv("LD_LIBRARY_PATH")
Sys.setenv(PATH=paste("/usr/bin:/bin", Sys.getenv("PATH"), sep=":"))
bwmr <- Sys.getenv("FLT_BWMR")
if (nzchar(bwmr) && !file.exists(bwmr)) {
    dir.create(dirname(bwmr), recursive=TRUE, showWarnings=FALSE)
    staged <- tempfile(tmpdir=dirname(bwmr))
    download.file("https://raw.githubusercontent.com/jiazhao97/BWMR/master/updates/BWMR_updated.R",
                  staged, mode="wb")
    invisible(parse(staged))
    if (!file.rename(staged, bwmr)) stop("Cannot save BWMR source")
}
cran <- Sys.getenv("FLT_CRAN", "https://mirrors.tuna.tsinghua.edu.cn/CRAN")
upstream <- c(TwoSampleMR="MRCIEU/TwoSampleMR", MVMR="WSpiller/MVMR",
    MRPRESSO="rondolab/MR-PRESSO", HDL="zhenin/HDL", GenomicSEM="GenomicSEM/GenomicSEM",
    limma="bioc/limma", impute="bioc/impute", MRMix="gqi/MRMix", RadialMR="WSpiller/RadialMR")
for (pkg in names(upstream)) {
    if (requireNamespace(pkg, quietly=TRUE) || length(Sys.glob(file.path(dest, paste0(pkg, "_*.tar.gz"))))) next
    tmp <- tempfile(); dir.create(tmp)
    archive <- file.path(tmp, "source.tar.gz")
    ref <- if (pkg %in% c("limma", "impute")) "RELEASE_3_22" else "HEAD"
    url <- paste0("https://api.github.com/repos/", upstream[[pkg]], "/tarball/", ref)
    cat("Fetch source: ", pkg, "\n", sep="")
    download.file(url, archive, mode="wb")
    untar(archive, exdir=tmp)
    roots <- list.dirs(tmp, recursive=FALSE)
    root <- roots[1L]
    if (pkg == "HDL") root <- file.path(root, "HDL")
    desc <- read.dcf(file.path(root, "DESCRIPTION"))
    stopifnot(desc[1L,"Package"] == pkg)
    file.rename(root, file.path(tmp, pkg))
    old <- setwd(tmp)
    tar(file.path(dest, paste0(pkg, "_", desc[1L,"Version"], ".tar.gz")), files=pkg,
        compression="gzip", tar="internal")
    setwd(old)
    unlink(tmp, recursive=TRUE)
}
tools::write_PACKAGES(dest, type="source")
local <- available.packages(contriburl=paste0("file://", dest), filters=list())
remote <- available.packages(repos=c(CRAN=cran))
db <- rbind(local, remote[!rownames(remote) %in% rownames(local), , drop=FALSE])
required <- c("tidyverse", "data.table", "survival", "readxl", "ggrepel", "lubridate",
    "broom", "R.utils", "openxlsx", "MendelianRandomization", "robustbase", "ieugwasr", "patchwork", "png", "deTS", names(upstream))
deps <- unique(c(required, unlist(tools::package_dependencies(required, db=db,
    which=c("Depends", "Imports", "LinkingTo"), recursive=TRUE))))
need <- deps[!vapply(deps, requireNamespace, logical(1), quietly=TRUE)]
# Installed namespaces can still be older than a new source package requires.
specs <- unlist(strsplit(na.omit(as.vector(db[intersect(deps, rownames(db)),
    c("Depends", "Imports", "LinkingTo"), drop=FALSE])), ","))
for (spec in specs) {
    spec <- trimws(gsub("[[:space:]]+", " ", spec))
    hit <- regmatches(spec, regexec("^([^ (]+) *\\((>=|<=|==|>|<) *([^ )]+)\\)", spec))[[1L]]
    if (length(hit) != 4L || hit[2L] == "R") next
    pkg <- hit[2L]
    if (!requireNamespace(pkg, quietly=TRUE) ||
        !do.call(hit[3L], list(packageVersion(pkg), package_version(hit[4L])))) need <- c(need, pkg)
}
need <- unique(need)
need <- setdiff(need, rownames(local))
unknown <- setdiff(need, rownames(remote))
if (length(unknown)) stop("Source repository missing: ", paste(unknown, collapse=", "))
if (length(need)) download.packages(need, destdir=dest, available=remote, repos=c(CRAN=cran), type="source")
tools::write_PACKAGES(dest, type="source")
cached <- available.packages(contriburl=paste0("file://", dest), filters=list())
if (length(setdiff(need, rownames(cached)))) stop("Incomplete source downloads; rerun install.sh run")
cat("HPC source cache ready: ", repo, "\n", sep="")
