# 🚩 Packages, paths, and shared inputs
# Download NHANES diet and mortality files for score validation.
pacman::p_load(tidyverse, data.table, stringr, purrr, readr, curl)

options(timeout = 1000, width = 200)

outdir <- "D:/data/nhanes/raw"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cycles <- tribble(
	~cycle,      ~path_year, ~prefix, ~suffix, ~use_in_maha,
	"1999-2000", 1999,       "",      "",      TRUE,
	"2001-2002", 2001,       "",      "_B",    TRUE,
	"2003-2004", 2003,       "",      "_C",    TRUE,
	"2005-2006", 2005,       "",      "_D",    TRUE,
	"2007-2008", 2007,       "",      "_E",    TRUE,
	"2009-2010", 2009,       "",      "_F",    TRUE,
	"2011-2012", 2011,       "",      "_G",    TRUE,
	"2013-2014", 2013,       "",      "_H",    TRUE,
	"2015-2016", 2015,       "",      "_I",    TRUE,
	# You can keep 2017-2018 on disk, but MAHA_NHANES.R should not combine it with 2017-2020.
	"2017-2018", 2017,       "",      "_J",    FALSE,
	"2017-2020", 2017,       "P_",    "",      TRUE,
	"2021-2023", 2021,       "",      "_L",    TRUE
)

tables <- tribble(
	~root,     ~component,       ~need_for_maha,
	"DEMO",   "Demographics",   TRUE,
	"BMX",    "Examination",    TRUE,
	"BPX",    "Examination",    TRUE,
	"OHX",    "Examination",    FALSE,
	"AUX",    "Examination",    FALSE,
	"VIX",    "Examination",    FALSE,
	"SPX",    "Examination",    FALSE,
	"PAX",    "Examination",    FALSE,
	"BPQ",    "Questionnaire",  TRUE,
	"MCQ",    "Questionnaire",  TRUE,
	"DIQ",    "Questionnaire",  TRUE,
	"SMQ",    "Questionnaire",  TRUE,
	"ALQ",    "Questionnaire",  TRUE,
	"PAQ",    "Questionnaire",  TRUE,
	"DPQ",    "Questionnaire",  TRUE,
	"SLQ",    "Questionnaire",  TRUE,
	"DBQ",    "Questionnaire",  TRUE,
	"WHQ",    "Questionnaire",  TRUE,
	"HIQ",    "Questionnaire",  FALSE,
	"HUQ",    "Questionnaire",  FALSE,
	"HSQ",    "Questionnaire",  FALSE,
	"OCQ",    "Questionnaire",  FALSE,
	"INQ",    "Questionnaire",  FALSE,
	"RXQ_RX", "Questionnaire",  TRUE,
	"RXQASA", "Questionnaire",  FALSE,
	"PFQ",    "Questionnaire",  FALSE,
	"KIQ_U",  "Questionnaire",  TRUE,
	"KIQ",    "Questionnaire",  TRUE,
	"HEQ",    "Questionnaire",  FALSE,
	"RHQ",    "Questionnaire",  FALSE,
	"SXQ",    "Questionnaire",  FALSE,
	"DR1TOT", "Dietary",        TRUE,
	"DR2TOT", "Dietary",        TRUE,
	"DR1IFF", "Dietary",        TRUE,
	"DR2IFF", "Dietary",        TRUE,
	"DRXFCD", "Dietary",        TRUE,
	"CBC",    "Laboratory",     FALSE,
	"BIOPRO", "Laboratory",     TRUE,
	"TCHOL",  "Laboratory",     TRUE,
	"HDL",    "Laboratory",     TRUE,
	"TRIGLY", "Laboratory",     TRUE,
	"GLU",    "Laboratory",     TRUE,
	"GHB",    "Laboratory",     TRUE,
	"INS",    "Laboratory",     TRUE,
	"CRP",    "Laboratory",     FALSE,
	"HSCRP",  "Laboratory",     TRUE,
	"ALB_CR", "Laboratory",     TRUE,
	"UCPREG", "Laboratory",     FALSE,
	"UCOT",   "Laboratory",     FALSE,
	"COT",    "Laboratory",     FALSE,
	"VID",    "Laboratory",     FALSE,
	"VITD",   "Laboratory",     FALSE,
	"FERTIN", "Laboratory",     FALSE,
	"FOLATE", "Laboratory",     FALSE,
	"B12",    "Laboratory",     FALSE
)

mk_stem <- function(root, cycle, prefix, suffix) {
	if (cycle == "2017-2020") paste0(prefix, root) else paste0(root, suffix)
}

stem_candidates <- function(root, cycle, prefix, suffix) {
	# Early NHANES diet files used DRX* names. Second dietary recall is generally available from 2003-2004 onward.
	if (root == "DR1TOT") {
		if (cycle %in% c("1999-2000", "2001-2002")) return(paste0("DRXTOT", suffix))
		return(mk_stem(root, cycle, prefix, suffix))
	}
	if (root == "DR2TOT") {
		if (cycle %in% c("1999-2000", "2001-2002")) return(character(0))
		return(mk_stem(root, cycle, prefix, suffix))
	}
	if (root == "DR1IFF") {
		if (cycle %in% c("1999-2000", "2001-2002")) return(unique(c(paste0("DRXIFF", suffix), paste0("DR1IFF", suffix))))
		return(mk_stem(root, cycle, prefix, suffix))
	}
	if (root == "DR2IFF") {
		if (cycle %in% c("1999-2000", "2001-2002")) return(character(0))
		return(mk_stem(root, cycle, prefix, suffix))
	}
	if (root == "DRXFCD") {
		return(unique(c(mk_stem("DRXFCD", cycle, prefix, suffix), mk_stem("DRXFMT", cycle, prefix, suffix))))
	}
	# Older lab releases sometimes used LAB/L-number stems. The URL check will keep only existing files.
	if (root == "GHB")    return(unique(c(mk_stem(root, cycle, prefix, suffix), "LAB10", "L10_B", "L10_C")))
	if (root == "TCHOL")  return(unique(c(mk_stem(root, cycle, prefix, suffix), "LAB13", "L13_B", "L13_C")))
	if (root == "HDL")    return(unique(c(mk_stem(root, cycle, prefix, suffix), "LAB13", "L13_B", "L13_C")))
	if (root == "TRIGLY") return(unique(c(mk_stem(root, cycle, prefix, suffix), "LAB13AM", "L13AM_B", "L13AM_C")))
	mk_stem(root, cycle, prefix, suffix)
}

make_manifest <- function(only_maha = FALSE) {
	tab <- if (only_maha) filter(tables, need_for_maha) else tables
	crossing(cycles, tab) |>
		rowwise() |>
		mutate(stem = list(stem_candidates(root, cycle, prefix, suffix))) |>
		ungroup() |>
		unnest(stem) |>
		filter(!is.na(stem), stem != "") |>
		mutate(
			file = paste0(stem, ".xpt"),
			url = paste0("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/", path_year, "/DataFiles/", file),
			dest_dir = file.path(outdir, cycle, component),
			dest_file = file.path(dest_dir, file)
		) |>
		dplyr::select(cycle, path_year, use_in_maha, component, root, need_for_maha, stem, file, url, dest_dir, dest_file) |>
		distinct(url, dest_file, .keep_all = TRUE)
}

url_exists <- function(url) {
	h <- new_handle(nobody = TRUE, followlocation = TRUE, failonerror = FALSE)
	x <- tryCatch(curl_fetch_memory(url, handle = h), error = function(e) NULL)
	!is.null(x) && x$status_code == 200
}

download_one <- function(url, dest_file) {
	dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
	if (file.exists(dest_file) && file.info(dest_file)$size > 0) return("exists")
	for (i in 1:3) {
		message("Downloading: ", basename(dest_file), " | try ", i)
		ok <- tryCatch({
			curl_download(url, dest_file, quiet = FALSE, mode = "wb")
			file.exists(dest_file) && file.info(dest_file)$size > 0
		}, error = function(e) {
			message("Failed: ", conditionMessage(e))
			FALSE
		})
		if (ok) return("downloaded")
		if (file.exists(dest_file) && file.info(dest_file)$size == 0) unlink(dest_file)
		Sys.sleep(2 * i)
	}
	"failed"
}

download_nhanes_xpt <- function(only_maha = FALSE) {
	manifest0 <- make_manifest(only_maha = only_maha)
	write_csv(manifest0, file.path(outdir, ifelse(only_maha, "manifest_candidates_maha.csv", "manifest_candidates_all.csv")))
	message("Candidate URLs: ", nrow(manifest0))

	manifest <- manifest0 |>
		mutate(available = map_lgl(url, url_exists)) |>
		filter(available) |>
		dplyr::select(-available)

	if (nrow(manifest) == 0) {
		message("URL check found 0 files. Now trying to download all candidate URLs directly.")
		manifest <- manifest0
	}

	message("Files selected for download: ", nrow(manifest))
	print(count(manifest, cycle, component), n = 300)

	log <- manifest |>
		arrange(path_year, component, root, stem) |>
		mutate(status = map2_chr(url, dest_file, download_one)) |>
		mutate(
			file_exists = file.exists(dest_file),
			size_mb = ifelse(file_exists, round(file.info(dest_file)$size / 1024^2, 3), NA_real_)
		)

	manifest_ok <- log |>
		filter(status %in% c("downloaded", "exists"), file_exists, size_mb > 0)

	write_csv(manifest_ok, file.path(outdir, ifelse(only_maha, "manifest_maha.csv", "manifest_common.csv")))
	write_csv(log, file.path(outdir, ifelse(only_maha, "download_log_maha.csv", "download_log.csv")))
	write_csv(filter(log, status == "failed" | !file_exists), file.path(outdir, ifelse(only_maha, "failed_downloads_maha.csv", "failed_downloads.csv")))

	log |>
		count(cycle, component, status) |>
		arrange(cycle, component, status) |>
		print(n = 300)

	message("Successfully downloaded/existing XPT files: ", nrow(manifest_ok))
	invisible(log)
}

download_mortality <- function() {
	mort_dir <- file.path(outdir, "mortality_2019_public")
	dir.create(mort_dir, recursive = TRUE, showWarnings = FALSE)
	mort_cycles <- c("1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014", "2015-2016", "2017-2018")
	mort_manifest <- tibble(cycle = mort_cycles) |>
		mutate(
			file = paste0("NHANES_", str_replace_all(cycle, "-", "_"), "_MORT_2019_PUBLIC.dat"),
			url = paste0("https://ftp.cdc.gov/pub/health_statistics/nchs/datalinkage/linked_mortality/", file),
			dest_file = file.path(mort_dir, file)
		)

	log <- mort_manifest |>
		mutate(status = map2_chr(url, dest_file, download_one)) |>
		mutate(
			file_exists = file.exists(dest_file),
			size_mb = ifelse(file_exists, round(file.info(dest_file)$size / 1024^2, 3), NA_real_)
		)
	write_csv(log, file.path(mort_dir, "download_log_mortality.csv"))
	print(log, n = 50)
	invisible(log)
}

# Set only_maha = TRUE for a smaller download focused on MAHA_NHANES.R.
# Set only_maha = FALSE to download the broader common table.
log_xpt <- download_nhanes_xpt(only_maha = FALSE)
log_mort <- download_mortality()

message("Done.")
message("Output folder: ", outdir)
message("Note: for MAHA_NHANES.R, use 2017-2020 pre-pandemic instead of 2017-2018 to avoid overlapping samples.")
