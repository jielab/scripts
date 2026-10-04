source(Sys.getenv("GU_RESULTS_R", "/mnt/d/scripts/0f/results.R"))


# 🚩 Reusable method results Tables remain ordinary R data frames. Native fitted objects and exact exchange bytes are
# kept together so algorithm readers can resume from /tmp.
gu_pack_result <- function(spec) {
	tables <- list()
	native <- list()
	for (entry in spec$files) {
		if (!is.null(entry$link)) {
			native[[entry$name]] <- list(link = entry$link)
			next
		}
		path <- file.path(spec$source, entry$name)
		bytes <- readBin(path, "raw", n = file.size(path))
		sha <- digest::digest(bytes, algo = "sha256", serialize = FALSE)
		native[[entry$name]] <- list(bytes = bytes, sha256 = sha, modified = entry$modified)
		if (grepl("[.]tsv([.]gz)?$", path) && length(bytes)) {
			table <- tryCatch(suppressWarnings(result_read_table(path)), error = function(e) NULL)
			if (!is.null(table))
				tables[[entry$name]] <- table
		}
	}
	value <- list(tables = tables)
	attr(value, "native_results") <- native
	attr(value, "figure_files") <- spec$figures
	attr(value, "source_root") <- spec$source_root
	result_write_rds(value, spec$destination)
	# The exact bytes, links, and parsed tables have all passed readRDS equality.
	invisible(spec$destination)
}

gu_restore_result <- function(source, target) {
	if (!startsWith(normalizePath(target, mustWork = FALSE), "/tmp/"))
		stop("GU exchanges must stay in /tmp")
	x <- readRDS(source)
	dir.create(target, recursive = TRUE, showWarnings = FALSE)
	writeLines(attr(x, "source_root"), file.path(target, ".result-source-root"))
	for (name in names(attr(x, "native_results"))) {
		if (startsWith(name, "/") || ".." %in% strsplit(name, "/", fixed = TRUE)[[1]])
			stop("Unsafe stored result path")
		entry <- attr(x, "native_results")[[name]]
		path <- file.path(target, name)
		dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
		if (!is.null(entry$link)) {
			if (nzchar(Sys.readlink(path)))
				unlink(path)
			if (!file.exists(path))
				file.symlink(entry$link, path)
		} else {
			if (!identical(digest::digest(entry$bytes, algo = "sha256", serialize = FALSE), entry$sha256))
				stop("Native result checksum mismatch")
			writeBin(entry$bytes, path)
			Sys.setFileTime(path, as.POSIXct(entry$modified, origin = "1970-01-01", tz = "UTC"))
		}
	}
	for (figure in attr(x, "figure_files")) {
		path <- file.path(target, figure$source)
		dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
		stopifnot(file.copy(file.path(dirname(source), figure$name), path, overwrite = TRUE))
		stopifnot(file.copy(file.path(dirname(source), sub("[.]png$", ".xlsx", figure$name)), sub("[.]png$", ".xlsx", path),
			overwrite = TRUE))
	}
	invisible(attr(x, "figure_files"))
}


# 🚩 Relational participant results
gu_pack_database <- function(source, destination) {
	con <- DBI::dbConnect(RSQLite::SQLite(), source)
	on.exit(DBI::dbDisconnect(con), add = TRUE)
	stopifnot(DBI::dbGetQuery(con, "PRAGMA quick_check")[[1]] == "ok")
	schema <- DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master WHERE sql IS NOT NULL AND name NOT LIKE 'sqlite_%' ORDER BY type DESC, name")
	tables <- setNames(lapply(DBI::dbListTables(con), function(name) DBI::dbReadTable(con, name)), DBI::dbListTables(con))
	attr(tables, "sqlite_schema") <- schema
	result_write_rds(tables, destination)
	message("Database tables preserved: ", length(tables), "; rows: ", sum(vapply(tables, nrow, integer(1))))
}

gu_restore_database <- function(source, destination) {
	if (!startsWith(normalizePath(dirname(destination), mustWork = FALSE), "/tmp/"))
		stop("SQLite is a temporary viewing cache")
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
	stopifnot(file.copy(temporary, destination, overwrite = TRUE))
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
				result_write_workbook(gu_b_tables(bundle, minimum), file.path(out, sub("[.]png$", ".xlsx", figure$name)))
			}
		}
	}
	# Eligibility and locus association results are useful even without a tree.
	tables <- list()
	for (name in c("gwas_loci", "gwas_haplotypes", "skipped_loci")) {
		x <- gu_b_read(file.path(final, paste0(name, ".tsv")))
		if (nrow(x)) {
			x <- x[, !result_private_columns(names(x)) & !grepl("(file|path|directory)$", names(x)), drop = FALSE]
			if (ncol(x))
				tables[[sub("gwas_", "", name)]] <- x
		}
	}
	if (length(tables))
		result_write_workbook(tables, file.path(out, "phyml.association.xlsx"))
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
	} else if (args[1] == "restore") {
		for (spec in jsonlite::read_json(args[2])) gu_restore_result(spec$source, spec$target)
	} else if (args[1] == "database-pack")
		gu_pack_database(args[2], args[3]) else if (args[1] == "database-restore")
		gu_restore_database(args[2], args[3]) else stop("Unknown GU storage command")
}
