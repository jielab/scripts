# Resolve the YGB dietary dictionary consistently for scoring and figures.
ygb_diet_dictionary <- function(raw_ygb) {
    override <- Sys.getenv("DIET_YGB_DIET_LST", "")
    candidates <- if (nzchar(override)) override else
        file.path(raw_ygb, "common", "diet.lst")
    available <- candidates[file.exists(candidates)]
    if (!length(available)) stop("Missing YGB dietary dictionary: ",
        paste(candidates, collapse=", "), call.=FALSE)
    available[[1L]]
}
