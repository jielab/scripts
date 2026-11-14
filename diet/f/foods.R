# The authoritative sheets remain untouched. Supplements supply only columns
# missing from those sheets; even blank authoritative cells must not be replaced.
diet_read_foods <- function(path, sheet, n_max=Inf) {
    read_sheet <- function(name) {
        x <- as.data.frame(readxl::read_excel(path,sheet=name),stringsAsFactors=FALSE)
        if(nrow(x)) x <- x[rowSums(!is.na(x))>0,,drop=FALSE]
        if(!"name" %in% names(x) || anyNA(x$name) || anyDuplicated(x$name))
            stop("Invalid dictionary names in ",path,"/",name)
        x
    }
    base <- read_sheet(sheet)
    extra <- switch(sheet,Sheet2="Sheet5",Sheet3="Sheet6",NULL)
    if(!is.null(extra) && extra %in% readxl::excel_sheets(path)) {
        supplement <- read_sheet(extra)
        idx <- match(base$name,supplement$name)
        if(anyNA(idx)) stop(extra," does not cover all authoritative ",sheet," names")
        for(nm in setdiff(names(supplement),names(base))) base[[nm]] <- supplement[[nm]][idx]
        # Extra bread subtype rows are retained in the workbook for explicit
        # use. Appending them here would double-count the existing .b/.w rows.
    }
    if(is.finite(n_max)) base <- head(base,n_max)
    base
}
