source(Sys.getenv("GU_RESULTS_R", "/mnt/d/scripts/0f/results.R"))


# 🚩 Reusable method results Tables remain ordinary R data frames. Native fitted objects and exact exchange bytes are
# kept together with the permanent native result files.
gu_read_metadata <- function(path) {
	lines <- readLines(path, warn = FALSE)
	lines <- lines[nzchar(lines)]
	parts <- strsplit(gsub("\\t", "\t", lines, fixed = TRUE), "\t", fixed = TRUE)
	width <- max(c(2L, lengths(parts)))
	x <- matrix("", nrow = length(parts), ncol = width)
	for (i in seq_along(parts)) x[i, seq_along(parts[[i]])] <- parts[[i]]
	x <- as.data.frame(x, stringsAsFactors = FALSE)
	names(x) <- c("key", paste0("value", seq_len(width - 1L)))
	x$record_line <- seq_along(lines)
	x$raw_record <- lines
	stopifnot(nrow(x) == length(lines))
	x
}

gu_pack_result <- function(spec) {
 tables <- list()
 for (entry in spec$files) {
  path <- file.path(spec$source, entry$name)
  if (is.null(entry$link) && grepl('[.]tsv([.]gz)?$', path) && file.size(path) > 0) {
   x <- if (grepl('(^|/)(run|cache)[.]meta[.]tsv$', path)) gu_read_metadata(path) else
    withCallingHandlers(result_read_table(path), warning = function(w) stop('Incomplete result table: ', path, ': ', conditionMessage(w)))
   if (!is.null(x) && ncol(x)) tables[[entry$name]] <- x
  }
 }
 hashes <- vapply(tables, digest::digest, character(1), algo = 'sha256')
 tables <- tables[!duplicated(hashes)]
 if (!length(tables)) tables <- list(availability = data.frame(status = 'No completed tabular results; native inputs are retained in the raw archive'))
 result_stream_workbook(tables, spec$destination)
 invisible(spec$destination)
}


gu_restore_result <- function(source, target) {
	x <- readRDS(source)
	dir.create(target, recursive = TRUE, showWarnings = FALSE)
	writeLines(attr(x, "source_root"), file.path(target, ".result-source-root"))
	for (name in names(attr(x, "native_results"))) {
		if (startsWith(name, "/") || ".." %in% strsplit(name, "/", fixed = TRUE)[[1]])
			stop("Unsafe stored result path")
		entry <- attr(x, "native_results")[[name]]
		path <- file.path(target, name)
		dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
		# Recovery must not replace newer native results or active-run files.
		link <- Sys.readlink(path)
		if (file.exists(path) || (!is.na(link) && nzchar(link))) next
		if (!is.null(entry$link)) {
			file.symlink(entry$link, path)
		} else {
			if (!identical(digest::digest(entry$bytes, algo = "sha256", serialize = FALSE), entry$sha256))
				stop("Native result checksum mismatch")
			stage <- tempfile(paste0(".", basename(path), ".part."), tmpdir = dirname(path))
			writeBin(entry$bytes, stage)
			stopifnot(file.rename(stage, path))
			Sys.setFileTime(path, as.POSIXct(entry$modified, origin = "1970-01-01", tz = "UTC"))
		}
	}
	for (figure in attr(x, "figure_files")) {
		path <- file.path(target, figure$source)
		dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
		if (!file.exists(path)) stopifnot(file.copy(file.path(dirname(source), figure$name), path))
		workbook <- sub("[.]png$", ".xlsx", path)
		if (!file.exists(workbook)) stopifnot(file.copy(file.path(dirname(source), sub("[.]png$", ".xlsx", figure$name)), workbook))
	}
	invisible(attr(x, "figure_files"))
}


# 🚩 Relational participant results
gu_pack_database <- function(source, destination) {
	con <- DBI::dbConnect(RSQLite::SQLite(), source)
	on.exit(DBI::dbDisconnect(con), add = TRUE)
	stopifnot(DBI::dbGetQuery(con, "PRAGMA quick_check")[[1]] == "ok")
	schema <- DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master WHERE sql IS NOT NULL AND name NOT LIKE 'sqlite_%' ORDER BY type DESC, name")
	table_names <- schema$name[schema$type == "table"]
	tables <- setNames(lapply(table_names, function(name) DBI::dbReadTable(con, name)), table_names)
	attr(tables, "sqlite_schema") <- schema
	result_write_rds(tables, destination)
	message("Database tables preserved: ", length(tables), "; rows: ", sum(vapply(tables, nrow, integer(1))))
}

gu_restore_database <- function(source, destination) {
	tables <- readRDS(source)
	schema <- attr(tables, "sqlite_schema")
	temporary <- tempfile("gu-database-", tmpdir = "/tmp", fileext = ".sqlite")
	on.exit(unlink(temporary), add = TRUE)
	con <- DBI::dbConnect(RSQLite::SQLite(), temporary)
	on.exit(if (DBI::dbIsValid(con)) DBI::dbDisconnect(con), add = TRUE)
	DBI::dbExecute(con, "PRAGMA journal_mode=OFF")
	DBI::dbWithTransaction(con, {
		for (i in which(schema$type == "table")) DBI::dbExecute(con, schema$sql[i])
		for (name in names(tables)) {
			DBI::dbAppendTable(con, name, tables[[name]])
			stopifnot(identical(DBI::dbReadTable(con, name), tables[[name]]))
		}
		for (i in which(schema$type != "table")) DBI::dbExecute(con, schema$sql[i])
	})
	stopifnot(DBI::dbGetQuery(con, "PRAGMA quick_check")[[1]] == "ok")
	DBI::dbDisconnect(con)
	dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
	stage <- tempfile(paste0(".", basename(destination), ".part."), tmpdir = dirname(destination))
	on.exit(unlink(stage), add = TRUE)
	stopifnot(file.copy(temporary, stage), file.rename(stage, destination))
	Sys.chmod(destination, "0600")
}


# 🚩 Saved phylogeny workbooks
gu_result_workbooks <- function(spec) {
	if (spec$method != "phyml")
		return(invisible(NULL))
	final <- file.path(spec$source, "final")
	out <- dirname(spec$destination)
	source(Sys.getenv("GU_PHYML_R", "/mnt/d/scripts/gu/f/phyml.R"), local = FALSE)
	for (table in c("evidence_trees.tsv", "region_common_trees.tsv")) {
		d <- gu_b_read(file.path(final, table))
		if (!nrow(d))
			next
		for (i in seq_len(nrow(d))) {
			row <- d[i, , drop = FALSE]
			if (!nzchar(gu_b_value(row, "tree_newick")))
				next
			bundle <- gu_b_bundle(row, final, row$locus_id[[1]])
			for (figure in spec$figures) {
				if (!grepl(paste0(".", bundle$lineage, ".panelB"), figure$source, fixed = TRUE))
						next
				if (length(unique(d$locus_id)) > 1L && !grepl(paste0("/", bundle$locus, "/"), figure$source, fixed = TRUE))
						next
				minimum <- if (grepl(".full.png", figure$source, fixed = TRUE))
						2L else 11L
				dir.create(out, recursive = TRUE, showWarnings = FALSE)
				grDevices::png(file.path(out, figure$name), width = 2400, height = 2500, res = 240, type = "cairo")
				tryCatch(gu_draw_panel_b(bundle, minimum), finally = grDevices::dev.off())
				result_write_workbook(gu_b_tables(bundle, minimum), file.path(out, sub("[.]png$", ".xlsx", figure$name)))
			}
		}
	}
	# Association/eligibility tables are included in the main per-run workbook.
}


# 🚩 Storage commands
if (sys.nframe() == 0L) {
	args <- commandArgs(TRUE)
	if (args[1] == "pack") {
		for (spec in jsonlite::read_json(args[2])) {
			gu_pack_result(spec)
			gu_result_workbooks(spec)
			message("Packed ", spec$method, ": ", basename(spec$source))
		}
	} else if (args[1] == "plots") {
		for (spec in jsonlite::read_json(args[2])) gu_result_workbooks(spec)
	} else if (args[1] == "restore") {
		for (spec in jsonlite::read_json(args[2])) gu_restore_result(spec$source, spec$target)
	} else if (args[1] == "database-pack")
		gu_pack_database(args[2], args[3]) else if (args[1] == "database-restore")
		gu_restore_database(args[2], args[3]) else stop("Unknown GU storage command")
}
