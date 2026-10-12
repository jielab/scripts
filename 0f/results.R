# 🚩 Named analysis results
result_private_columns <- function(columns) {
	names <- tolower(gsub('[^[:alnum:]]', '', columns))
	names %in% c('id', 'eid', 'iid', 'fid', 'sample', 'samples', 'sampleid', 'sampleids',
		'individualid', 'participantid', 'subjectid', 'patientid', 'copyid', 'copyids',
		'referencesampleid', 'membersamples', 'carriersamples', 'donorid',
		'memberids', 'carrierids', 'samplenames', 'individuals', 'donorids') |
		grepl('^(sample|individual|participant|subject|patient)(ids|identifiers)$', names)
}

result_read_table <- function(path) {
	if (grepl('[.]rds$', path, ignore.case = TRUE)) return(as.data.frame(readRDS(path)))
	header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
	ids <- header[result_private_columns(header)]
	as.data.frame(data.table::fread(path, colClasses = if (length(ids)) list(character = ids) else NULL,
		check.names = FALSE, showProgress = FALSE))
}

result_write_rds <- function(value, path) {
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	temporary <- tempfile('result-', tmpdir = '/tmp', fileext = '.rds')
	on.exit(unlink(temporary), add = TRUE)
	saveRDS(value, temporary, compress = 'gzip')
	if (!identical(readRDS(temporary), value)) stop('RDS verification failed: ', path)
	if (!file.copy(temporary, path, overwrite = TRUE) ||
		!identical(unname(tools::md5sum(temporary)), unname(tools::md5sum(path)))) stop('RDS publication failed: ', path)
	Sys.chmod(path, '0600')
	invisible(path)
}

result_import_table <- function(source, destination, metadata = NULL) {
	x <- result_read_table(source)
	attr(x, 'table_exchange') <- list(name = basename(source),
		bytes = readBin(source, 'raw', n = file.size(source)), sha256 = digest::digest(file = source, algo = 'sha256'))
	if (!is.null(metadata)) attr(x, 'model_info') <- metadata
	if (grepl('[.]xlsx$', destination, ignore.case = TRUE)) result_write_workbook(list(results = x), destination, source)
	else result_write_rds(x, destination)
}


# 🚩 Large review tables and migration from table-only RDS
result_stream_workbook <- function(tables, path, sources = character(), metadata = NULL) {
 directory <- tempfile('result-stream-', tmpdir = '/tmp')
 dir.create(directory)
 on.exit(unlink(directory, recursive = TRUE), add = TRUE)
 spec <- list(tables = list(), sources = unname(as.list(sources)), metadata = metadata)
 for (i in seq_along(tables)) {
  x <- as.data.frame(tables[[i]])
  source <- file.path(directory, paste0(sprintf('%04d_', i), gsub('[^[:alnum:]_.-]', '_', names(tables)[i]), '.tsv'))
  types <- vapply(names(x), function(n) {
   v <- x[[n]]
   if (result_private_columns(n) || inherits(v, c('Date','POSIXt','integer64'))) 'character'
   else if (is.numeric(v)) 'numeric' else if (is.logical(v)) 'logical' else 'character'
  }, character(1))
  for (n in names(x)) {
   if (is.list(x[[n]])) x[[n]] <- vapply(x[[n]], function(z) jsonlite::toJSON(z, auto_unbox = TRUE, null = 'null'), character(1))
   if (types[[n]] == 'character') x[[n]] <- as.character(x[[n]])
   if (types[[n]] == 'numeric') x[[n]] <- ifelse(is.na(x[[n]]), NA_character_, sprintf('%.17g', x[[n]]))
  }
  data.table::fwrite(x, source, sep = '\t', na = '')
  spec$tables[[i]] <- list(name = names(tables)[i], path = source, types = unname(as.list(types)))
 }
 job <- file.path(directory, 'workbook.json')
 jsonlite::write_json(spec, job, auto_unbox = TRUE, null = 'null')
 code <- Sys.getenv('RESULTS_PY', '/mnt/d/scripts/0f/results.py')
 status <- system2(Sys.getenv('RESULTS_PYTHON', 'python3'), c(shQuote(code), 'stream-workbook', shQuote(job), shQuote(path)))
 if (status != 0L) stop('Workbook publication failed: ', path)
 invisible(path)
}
result_rds_tables <- function(x) {
 if (is.data.frame(x)) tables <- list(results = x)
 else if (is.list(x) && is.list(x$tables)) tables <- x$tables
 else if (is.list(x)) tables <- Filter(is.data.frame, x)
 else stop('RDS contains a fitted/raw object, not review tables')
 figure <- attr(x, 'figure_data')
 if (is.data.frame(figure)) tables$figure_data <- figure
 else if (is.list(figure)) tables <- c(tables, Filter(is.data.frame, figure))
 tables <- Filter(function(t) is.data.frame(t) && ncol(t) > 0L, tables)
 hashes <- vapply(tables, digest::digest, character(1), algo = 'sha256')
 tables <- tables[!duplicated(hashes)]
 if (!length(tables)) stop('No reviewable tables in this RDS')
 tables
}
result_convert_rds <- function(source, destination) {
 x <- readRDS(source)
 tables <- result_rds_tables(x)
 result_stream_workbook(tables, destination, metadata = attr(x, 'model_info'))
 invisible(destination)
}

# 🚩 Figure workbooks
result_write_workbook <- function(tables, path, sources = character()) {
 tables <- Filter(function(x) is.data.frame(x) && ncol(x) > 0L, tables)
 if (!length(tables)) stop('No analysis results for workbook: ', path)
 result_stream_workbook(tables, path, sources)
}


# 🚩 Python and shell exchange in /tmp
if (sys.nframe() == 0L) {
	args <- commandArgs(TRUE)
	if (length(args) < 2L) stop('Usage: results.R import-table|export-table|metadata SOURCE [DESTINATION] [METADATA_JSON]')
	if (args[1] == 'metadata') {
		cat(jsonlite::toJSON(attr(readRDS(args[2]), 'model_info'), auto_unbox = TRUE, null = 'null'))
		quit(save = 'no')
	}
	if (length(args) < 3L) stop('A result destination is required')
	if (args[1] == 'convert-rds') {
		result_convert_rds(args[2], args[3])
	} else if (args[1] == 'workbook') {
		spec <- jsonlite::read_json(args[2], simplifyVector = TRUE)
		tables <- lapply(spec$tables, result_read_table)
		result_write_workbook(tables, args[3], unlist(spec$tables, use.names = FALSE))
	} else if (args[1] == 'import-object') {
		result_write_rds(jsonlite::read_json(args[2], simplifyVector = TRUE), args[3])
	} else if (args[1] == 'import-table') {
		metadata <- if (length(args) > 3L) jsonlite::read_json(args[4], simplifyVector = TRUE) else NULL
		result_import_table(args[2], args[3], metadata)
	} else if (args[1] == 'export-table') {
		x <- readRDS(args[2])
		exchange <- attr(x, 'table_exchange')
		target <- args[3]
		if (dir.exists(target)) target <- file.path(target, if (!is.null(exchange)) exchange$name else sub('[.]rds$', '.tsv', basename(args[2])))
		if (!startsWith(normalizePath(dirname(target), mustWork = TRUE), '/tmp/')) stop('Table exchanges must stay under /tmp')
		if (!is.null(exchange) && identical(basename(target), exchange$name)) {
			if (!identical(digest::digest(exchange$bytes, algo = 'sha256', serialize = FALSE), exchange$sha256)) stop('Stored source checksum mismatch')
			writeBin(exchange$bytes, target)
		} else data.table::fwrite(x, target, sep = '\t', na = 'NA')
		cat(normalizePath(target))
	} else stop('Unknown result command: ', args[1])
}
