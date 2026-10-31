# Shared PNG export: a 220-DPI render ceiling plus a verified 1,000,000-byte cap.
# PDF/SVG exports keep their original rendering settings.
.png_output_python <- file.path(dirname(normalizePath(tail(Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile)), 1)[[1]], mustWork = TRUE)), "png_output.py")

limit_png_file <- function(path) {
    if (!file.exists(path)) stop("Missing exported PNG: ", path)
    if (file.size(path) <= 1000000) return(invisible(path))
    python <- Sys.getenv("PNG_PYTHON", Sys.getenv("PYTHON_BIN", "python3"))
    status <- system2(python, c(shQuote(.png_output_python), shQuote(path)))
    if (status != 0L || !file.exists(path) || file.size(path) > 1000000)
        stop("PNG export could not meet the 1 MB limit: ", path)
    invisible(path)
}

bounded_ggsave <- function(filename, plot = ggplot2::last_plot(), ..., dpi = 300, path = NULL) {
    target <- if (is.null(path)) filename else file.path(path, filename)
    is_png <- grepl("[.]png$", target, ignore.case = TRUE)
    if (is_png) {
        if (is.character(dpi)) dpi <- c(screen = 72, print = 300, retina = 320)[[dpi]]
        dpi <- min(dpi, 220)
    }
    ggplot2::ggsave(filename, plot = plot, ..., dpi = dpi, path = path)
    if (is_png) limit_png_file(target)
    invisible(target)
}

# Explicit wrappers for base/grid graphics; compress only after closing a PNG.
if (!exists(".bounded_png_devices", inherits = FALSE)) .bounded_png_devices <- new.env(parent = emptyenv())
bounded_png <- function(filename, ...) {
    if (grepl("%", filename, fixed = TRUE)) stop("Use one explicit PNG filename per export")
    grDevices::png(filename, ...)
    .bounded_png_devices[[as.character(grDevices::dev.cur())]] <- normalizePath(filename, mustWork = FALSE)
    invisible(NULL)
}
bounded_dev_off <- function(which = grDevices::dev.cur()) {
    key <- as.character(which)
    path <- .bounded_png_devices[[key]]
    value <- grDevices::dev.off(which)
    if (!is.null(path)) {
        rm(list = key, envir = .bounded_png_devices)
        limit_png_file(path)
    }
    invisible(value)
}
