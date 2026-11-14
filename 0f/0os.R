# R adapter for the single path table in 0os.sh. No data-root defaults live here.
.os_table <- file.path(dirname(normalizePath(tail(Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile)), 1L)[[1L]], winslash = "/")), "0os.sh")
os_to_windows <- function(path) {
    converted <- sub("^/mnt/([A-Za-z])(/|$)", "\\1:/", path)
    hit <- grepl("^/mnt/[A-Za-z](/|$)", path)
    substr(converted[hit], 1L, 1L) <- toupper(substr(converted[hit], 1L, 1L))
    converted
}
os_to_wsl <- function(path) {
    path <- chartr("\\", "/", path)
    ifelse(grepl("^[A-Za-z]:/", path),
           paste0("/mnt/", tolower(substr(path, 1L, 1L)), substring(path, 3L)), path)
}
os_paths_init <- function(default_mode = "auto") {
    mode <- Sys.getenv("submit", Sys.getenv("SUBMIT", default_mode))
    table <- .os_table
    native_windows <- .Platform$OS.type == "windows"
    to_wsl <- os_to_wsl
    if (native_windows) table <- to_wsl(table)
    command <- paste0("source ", shQuote(table, type = "sh"), " || exit; ",
        "for name in submit SUBMIT \"${OS_PATH_KEYS[@]}\"; do ",
        "printf '%s\\t%s\\n' \"$name\" \"${!name}\"; marker=_OS_DEFAULT_$name; if [[ -v $marker ]]; then printf '%s\\t%s\\n' \"$marker\" \"${!marker}\"; fi; done")
    # Pass explicit overrides to WSL too when called by native Windows R.
    keys <- grep("^(DIR_|OS_|_OS_DEFAULT_)|^ANALYSIS_ROOT$", names(Sys.getenv()), value = TRUE)
    command <- paste0("export submit=", shQuote(mode, type = "sh"), "; ", command)
    if (native_windows) {
        overrides <- Sys.getenv(keys); overrides <- overrides[nzchar(overrides)]
        overrides <- vapply(overrides, to_wsl, "")
        if (length(overrides)) command <- paste0("export ", paste(paste0(names(overrides), "=", shQuote(overrides, type = "sh")), collapse = " "), "; ", command)
        out <- system2("wsl.exe", c("--exec", "bash", "-c", shQuote(command, type = "cmd")), stdout = TRUE)
    } else out <- system2("bash", c("-c", shQuote(command)), stdout = TRUE)
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0L) stop("Could not load paths from ", .os_table)
    fields <- strsplit(out, "\t", fixed = TRUE)
    fields <- Filter(function(x) length(x) == 2L, fields)
    values <- setNames(vapply(fields, `[[`, "", 2L), vapply(fields, `[[`, "", 1L))
    if (!"DIR_UKB" %in% names(values)) stop("0os.sh did not return DIR_UKB")
    if (native_windows) values <- os_to_windows(values)
    do.call(Sys.setenv, as.list(values))
    invisible(values)
}
os_path <- function(name) {
    value <- Sys.getenv(name, unset = "")
    if (!nzchar(value)) stop("Missing ", name, "; initialize paths with os_paths_init() from 0f/0os.R")
    value
}
