# 🚩 Named analysis results
result_private_columns <- function(columns) {
	names <- tolower(gsub('[^[:alnum:]]', '', columns))
	names %in% c('eid', 'iid', 'fid', 'sample', 'samples', 'sampleid', 'sampleids',
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
	result_write_rds(x, destination)
}


# 🚩 Figure workbooks
result_write_workbook <- function(tables, path, sources = character()) {
	tables <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L && ncol(x) > 0L, tables)
	if (!length(tables)) stop('No analysis results for workbook: ', path)
	if (any(vapply(tables, function(x) any(result_private_columns(names(x))), logical(1))))
		stop('Individual records belong in a named RDS: ', path)
	if (any(vapply(tables, function(x) any(vapply(x, function(column) {
		if (!is.character(column)) return(FALSE)
		any(grepl('(^|[;,|[:space:]])(HG|NA)[0-9]{5}($|[;,|[:space:]])', column), na.rm = TRUE)
	}, logical(1))), logical(1)))) stop('Individual identifiers found in workbook cells: ', path)
	if (any(vapply(tables, function(x) nrow(x) > 1048575L || ncol(x) > 16384L, logical(1))))
		stop('Split this result table by analysis scope: ', path)
	sheets <- substr(gsub('[^[:alnum:]_. -]', '_', names(tables)), 1L, 31L)
	if (any(!nzchar(sheets)) || anyDuplicated(tolower(sheets))) stop('Use distinct, descriptive worksheet names')
	workbook <- openxlsx::createWorkbook()
	for (i in seq_along(tables)) {
		openxlsx::addWorksheet(workbook, sheets[i])
		openxlsx::writeData(workbook, sheets[i], tables[[i]], withFilter = TRUE)
		openxlsx::freezePane(workbook, sheets[i], firstRow = TRUE)
	}
	temporary <- tempfile('result-workbook-', tmpdir = '/tmp', fileext = '.xlsx')
	on.exit(unlink(temporary), add = TRUE)
	openxlsx::saveWorkbook(workbook, temporary, overwrite = TRUE)
	if (length(sources)) {
		# Preserve exact aggregate exports for software readers and numerical reuse.
		stage <- tempfile('result-sources-', tmpdir = '/tmp')
		dir.create(file.path(stage, 'results', 'exports'), recursive = TRUE)
		on.exit(unlink(stage, recursive = TRUE), add = TRUE)
		manifest <- list(format = 'analysis-tables-v1', files = list())
		for (source in sources) {
			if (any(result_private_columns(names(result_read_table(source))))) stop('Private workbook source: ', source)
			name <- basename(source)
			if (name %in% names(manifest$files)) stop('Duplicate source name: ', name)
			sha <- digest::digest(file = source, algo = 'sha256')
			part <- paste0('results/exports/', sha, '.bin')
			stopifnot(file.copy(source, file.path(stage, part), overwrite = TRUE))
			manifest$files[[name]] <- list(part = part, sha256 = sha)
		}
		jsonlite::write_json(manifest, file.path(stage, 'results/manifest.json'), auto_unbox = TRUE)
		utils::unzip(temporary, files = '[Content_Types].xml', exdir = stage)
		types <- file.path(stage, '[Content_Types].xml')
		xml <- paste(readLines(types, warn = FALSE), collapse = '')
		for (extension in c('json', 'bin')) if (!grepl(paste0('Extension="', extension, '"'), xml, fixed = TRUE)) {
			type <- if (extension == 'json') 'application/json' else 'application/octet-stream'
			xml <- sub('</Types>', paste0('<Default Extension="', extension, '" ContentType="', type, '"/></Types>'), xml, fixed = TRUE)
		}
		# Keep the printer-settings .bin default; give exact result exports
		# explicit types so Excel does not mistake them for printer records.
		parts <- unique(vapply(manifest$files, `[[`, character(1), 'part'))
		for (part in parts) {
			override <- paste0('<Override PartName="/', part, '" ContentType="application/octet-stream"/>')
			xml <- sub('</Types>', paste0(override, '</Types>'), xml, fixed = TRUE)
		}
		writeLines(xml, types, useBytes = TRUE)
		zip::zip_append(temporary, c('[Content_Types].xml', 'results'), root = stage,
			compression_level = 6, include_directories = FALSE)
	}
	if (any(grepl('/$', utils::unzip(temporary, list = TRUE)$Name))) stop('Unexpected directory entries in workbook: ', path)
	stopifnot(identical(openxlsx::getSheetNames(temporary), sheets))
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	if (!file.copy(temporary, path, overwrite = TRUE) ||
		!identical(unname(tools::md5sum(temporary)), unname(tools::md5sum(path)))) stop('Workbook publication failed: ', path)
	invisible(path)
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
	if (args[1] == 'workbook') {
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
