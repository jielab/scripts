# Aggregate input-date sanity checks. Never censor or alter the source values.
ukb_input_date_range <- function(x, variable, source) {
    flatten <- function(z) {
        if (is.list(z) && !inherits(z, "POSIXlt"))
            return(unlist(lapply(z, flatten), use.names = FALSE))
        if (inherits(z, "Date") || inherits(z, "POSIXt") || is.numeric(z))
            return(as.character(ukb_date_parse(z)))
        unlist(strsplit(as.character(z), "[;|,]+"), use.names = FALSE)
    }
    if (inherits(x, "Date") || inherits(x, "POSIXt") || is.numeric(x)) {
        parsed <- ukb_date_parse(x)
        present <- !is.na(x)
    } else {
        raw <- trimws(flatten(x))
        present <- !is.na(raw) & !raw %in% c("", "NA", "NaN", "NULL", ".")
        # Parse each distinct date once; large cached columns repeat dates often.
        unique_raw <- unique(raw)
        parsed <- ukb_date_parse(unique_raw)[match(raw, unique_raw)]
    }
    valid <- !is.na(parsed)
    special <- parsed %in% as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"))
    list(variable = variable, source = basename(source), n_values = length(parsed),
         n_dates = sum(valid), n_missing = sum(!present), n_unparseable = sum(present & !valid),
         n_special = sum(special),
         min_date = if (any(valid)) as.character(min(parsed[valid])) else NULL,
         max_date = if (any(valid)) as.character(max(parsed[valid])) else NULL)
}
