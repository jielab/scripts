# Run through ems120.sh --preflight, inside an allocated worker.
# Only workbook samples and synthetic model data are used; no production caches.
script_dir <- Sys.getenv("EMS120_SCRIPT_DIR")
main <- file.path(script_dir, "f", "ems120.R")
lines <- readLines(main, warn = FALSE, encoding = "UTF-8")
invisible(parse(text = lines))
boundary <- grep("^# .*Data harmonization, ml_phone, ml_geo, and ml_dx$", lines)
stopifnot(length(boundary) == 1L)

# Execute the real startup, helpers and constants, stopping before data processing.
# Give this startup its own temporary analysis root even if the caller selected
# an existing production analysis directory.
production_analysis <- Sys.getenv("EMS120_ANALYSIS_DIR")
probe_root <- tempfile("ems120-preflight-", tmpdir = tempdir())
dir.create(probe_root)
Sys.setenv(EMS120_ANALYSIS_DIR = probe_root, EMS120_MODEL_DIR = file.path(probe_root, "ml"))
eval(parse(text = lines[seq_len(boundary - 1L)]), envir = globalenv())
cat("PASS: complete R syntax and actual analysis startup/shared helpers\n")

for (yy in years) {
    raw <- file.path(dir.dat, "clean", paste0(yy, ".xlsx"))
    geo <- file.path(dir.dat, "szu_xia", paste0(yy, ".geo.xlsx"))
    stopifnot(file.exists(raw), file.info(raw)$size > 0, file.exists(geo), file.info(geo)$size > 0)
    h <- trimws(names(readxl::read_excel(raw, n_max = 0)))
    essential <- match(c("电话", "地址", "开始受理时刻", "性别", "年龄"), vars)
    missing <- vars[essential][!vapply(vars.alias[essential], function(p) any(grepl(p, h)), logical(1))]
    if (length(missing)) stop(yy, " missing EMS columns: ", paste(missing, collapse = ", "))
    gh <- clean_xia_colnames(names(readxl::read_excel(geo, n_max = 0)))
    missing <- setdiff(fig2_required_vars, union(h, gh))
    if (length(missing)) stop(yy, " missing geo/house-price columns: ", paste(missing, collapse = ", "))
    cat(sprintf("PASS: %d EMS and geo/house-price workbooks, required column names\n", yy))
    invisible(gc())
}

# Validate the same first 10,000 annotated 2019 rows used by a fresh full run.
train_path <- file.path(production_analysis, "dat", "2019.train_dx.xlsx")
if (file.exists(train_path)) {
    train <- readxl::read_excel(train_path)
    cat("Using existing training workbook for validation\n")
} else {
    train <- readxl::read_excel(file.path(dir.dat, "clean", "2019.xlsx"), n_max = 10000)
    original <- trimws(names(train)); renamed <- original
    for (i in seq_along(vars)) {
        matches <- grep(vars.alias[i], original, value = TRUE)
        if (!length(matches)) next
        best <- matches[1]
        if (length(matches) > 1) for (p in strsplit(vars.alias[i], "\\|")[[1]]) {
            m <- grep(p, matches, value = TRUE)
            if (length(m)) { best <- if (vars[i] == "电话" && vars.phone[["2019"]] %in% matches) vars.phone[["2019"]] else m[1]; break }
        }
        renamed[match(best, original)] <- vars[i]
    }
    names(train) <- renamed
    for (v in setdiff(vars, names(train))) train[[v]] <- NA
}
stopifnot(all(c("疾病类型", "地址", "地址类型") %in% names(train)))
geo_labels <- trimws(as.character(train$地址类型))
stopifnot(length(unique(geo_labels[!is.na(geo_labels) & nzchar(geo_labels)])) >= 2)

ensure_ml_python()
preflight_py <- reticulate::import_main(convert = FALSE)
preflight_py$ems120_preflight_training <- train
reticulate::source_python(file.path(script_dir, "f", "preflight.py"))
rm(train)

# Exercise the actual export helper on this headless node, with synthetic data.
p <- ggplot2::ggplot(data.frame(x = 1:3, y = c(1, 3, 2)), ggplot2::aes(x, y)) + ggplot2::geom_point()
png_file <- file.path(probe_root, "preflight.png")
save_plot(p, png_file, width = 3, height = 2, dpi = 150)
save_plot(p, file.path(probe_root, "preflight.pdf"), width = 3, height = 2)
stopifnot(file.info(png_file)$size > 0, file.info(png_file)$size <= 1000000)
probe_xlsx <- file.path(probe_root, "preflight.xlsx")
writexl::write_xlsx(data.frame(result = "PASS"), probe_xlsx)
stopifnot(identical(readxl::read_excel(probe_xlsx)$result, "PASS"))
asset <- file.path(dir.files, "YinYang.png")
if (!file.exists(asset)) asset <- file.path(dir.files, "YinYang_clean.png")
stopifnot(file.exists(asset))
invisible(png::readPNG(asset))
cat("PASS: headless PNG/PDF export, Excel round trip and YinYang image\n")
cat("EMS120 PREFLIGHT PASS: all 24 input workbooks, training labels, dependencies and synthetic GPU model checks\n")
