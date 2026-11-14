# Resolve the UKB dietary dictionary consistently for scoring and figures.
ukb_diet_dictionary <- function(raw_ukb) {
    override <- Sys.getenv("DIET_UKB_DIET_LST", "")
    candidates <- if (nzchar(override)) override else
        file.path(raw_ukb, "common", "diet.lst")
    available <- candidates[file.exists(candidates)]
    if (!length(available)) stop("Missing UKB dietary dictionary: ",
        paste(candidates, collapse=", "), call.=FALSE)
    available[[1L]]
}
