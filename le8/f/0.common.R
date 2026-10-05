# 🚩 Consolidated result tables
# Excel workbooks contain aggregate worksheets and exact source exports.
# Named RDS files contain participant data; CSV/TSV exchanges stay in /tmp.
le8_table_tmpdir <- function() if (.Platform$OS.type == 'windows') tempdir() else '/tmp'
le8_table_private_columns <- c(
	'eid', 'iid', 'fid', 'id', 'id1', 'id2', 'participantid', 'individualid',
	'subjectid', 'patientid', 'queryid', 'targetid', 'referenceid', 'donorid',
	'ukbid', 'ukbeid', 'sampleid', 'datebirth', 'birthdate'
)
le8_table_private <- function(x) {
	names <- tolower(names(x))
	# S7 aggregate metrics use reference_id for a MODEL; participant reference IDs
	# remain private in every other schema, and eid columns are always private.
	if (all(c('run_id', 'model_id', 'primary_id', 'reference_id', 'n', 'uno_c_horizon') %in% names))
		names <- names[names != 'reference_id']
	any(gsub('[^a-z0-9]', '', names) %in% le8_table_private_columns |
		grepl('(^|[_. #])(eid|iid|fid)($|[_. ])', names))
}
le8_table_blank <- function(path) {
	if (file.info(path)$size == 0) return(TRUE)
	if (file.info(path)$size > 1024) return(FALSE)
	!any(nzchar(trimws(readLines(path,warn=FALSE))))
}
le8_table_read <- function(path) {
	if (le8_table_blank(path)) return(data.frame())
	separator <- if (grepl('[.]tsv([.]gz)?$', path)) '\t' else ','
	tryCatch(as.data.frame(data.table::fread(path, sep = separator, check.names = FALSE, showProgress = FALSE)),
		error = function(e) stop('Cannot preserve table: ', path, ': ', conditionMessage(e)))
}
le8_table_entry <- function(path, table = NULL) {
	list(data = table, bytes = readBin(path, 'raw', n = file.size(path)),
		sha256 = digest::digest(file = path, algo = 'sha256'), mtime = file.info(path)$mtime)
}
le8_table_save <- function(value, path) {
	temporary <- tempfile('le8-tables-', tmpdir = le8_table_tmpdir(), fileext = '.rds')
	on.exit(unlink(temporary), add = TRUE)
	# Fast compression keeps publication responsive for large reusable tables.
	connection <- gzfile(temporary, open = 'wb', compression = 1)
	tryCatch(saveRDS(value, connection), finally = close(connection))
	check <- readRDS(temporary)
	if (!identical(check, value)) stop('RDS verification failed: ', path)
	if (!file.copy(temporary, path, overwrite = TRUE, copy.mode = FALSE)) stop('Cannot publish ', path)
	if (!identical(unname(tools::md5sum(path)), unname(tools::md5sum(temporary)))) stop('Published RDS checksum mismatch: ', path)
	invisible(path)
}
le8_table_archive_files <- function(directory) {
	paths <- list.files(directory, pattern = '[.]xlsx$', full.names = TRUE)
	paths[vapply(paths, function(path) 'le8/manifest.json' %in% utils::unzip(path, list = TRUE)$Name, logical(1))]
}
le8_table_archive_read <- function(path) {
	connection <- unz(path, 'le8/manifest.json', open = 'r')
	on.exit(close(connection), add = TRUE)
	x <- jsonlite::fromJSON(paste(readLines(connection, warn = FALSE), collapse = '\n'), simplifyVector = FALSE)
	if (!identical(x$format, 'le8-workbook-v1')) stop('Invalid result workbook: ', path)
	for (name in names(x$files)) {
		entry <- x$files[[name]]
		if (basename(name) != name || !grepl('^le8/exports/[a-f0-9]{64}[.]bin$', entry$part)) stop('Unsafe workbook source path')
		entry$workbook <- path
		entry$mtime <- as.POSIXct(entry$mtime, origin = '1970-01-01', tz = 'UTC')
		x$files[[name]] <- entry
	}
	x
}
le8_table_load <- function(directory) {
	x <- list(format = 'le8-workbook-v1', files = list(), private_files = list())
	for (path in le8_table_archive_files(directory)) {
		part <- le8_table_archive_read(path)
		if (any(names(part$files) %in% names(x$files))) stop('Duplicate sources in result workbooks: ', directory)
		x$files <- c(x$files, part$files)
		x$private_files <- c(x$private_files, part$private_files)
	}
	x
}
le8_table_bytes <- function(entry) {
	if (is.raw(entry$bytes)) bytes <- entry$bytes else {
		connection <- unz(entry$workbook, entry$part, open = 'rb')
		on.exit(close(connection), add = TRUE)
		bytes <- readBin(connection, 'raw', n = entry$size)
	}
	if (!identical(digest::digest(bytes, algo = 'sha256', serialize = FALSE), entry$sha256)) stop('Stored export checksum mismatch')
	bytes
}
le8_table_sheets <- function(path) {
	workbook <- openxlsx::loadWorkbook(path)
	setNames(lapply(names(workbook), function(sheet) withCallingHandlers(
		openxlsx::read.xlsx(workbook, sheet = sheet, check.names = FALSE),
		warning = function(w) {
			if (grepl('^No data found on worksheet[.]?$', trimws(conditionMessage(w)))) invokeRestart('muffleWarning')
		}
	)), names(workbook))
}
le8_table_materialize <- function(entry, source) {
	entry$bytes <- le8_table_bytes(entry)
	if (is.data.frame(entry$data) || length(entry$sheets) || identical(entry$kind, 'exchange')) return(entry)
	suffix <- if (grepl('[.]xlsx$', source)) '.xlsx' else {
		paste0(if (grepl('[.]tsv([.]gz)?$', source)) '.tsv' else '.csv', if (grepl('[.]gz$', source)) '.gz')
	}
	temporary <- tempfile('le8-read-table-', tmpdir = le8_table_tmpdir(), fileext = suffix)
	on.exit(unlink(temporary), add = TRUE)
	writeBin(entry$bytes, temporary)
	if (suffix == '.xlsx') entry$sheets <- le8_table_sheets(temporary) else entry$data <- le8_table_read(temporary)
	entry
}
le8_table_archive_add <- function(path, entries, private_files) {
	stage <- tempfile('le8-workbook-sources-', tmpdir = le8_table_tmpdir())
	dir.create(file.path(stage, 'le8', 'exports'), recursive = TRUE)
	on.exit(unlink(stage, recursive = TRUE), add = TRUE)
	manifest <- list(format = 'le8-workbook-v1', layout = 'topic-workbooks-v2', files = list(), private_files = private_files)
	for (source in names(entries)) {
		entry <- entries[[source]] ; bytes <- le8_table_bytes(entry)
		part <- paste0('le8/exports/', entry$sha256, '.bin')
		writeBin(bytes, file.path(stage, part))
		kind <- entry$kind
		if (is.null(kind)) kind <- if (length(entry$sheets)) 'workbook' else if (is.data.frame(entry$data)) 'table' else 'exchange'
		manifest$files[[source]] <- list(part = part, size = length(bytes), sha256 = entry$sha256, mtime = as.numeric(entry$mtime), kind = kind)
	}
	jsonlite::write_json(manifest, file.path(stage, 'le8', 'manifest.json'), auto_unbox = TRUE, null = 'null', digits = NA)
	# Standard workbook cells remain readable in Excel. Package parts retain
	# exact aggregate exports for numerical reuse and source-hash verification.
	utils::unzip(path, files = '[Content_Types].xml', exdir = stage)
	content_path <- file.path(stage, '[Content_Types].xml')
	content <- paste(readLines(content_path, warn = FALSE), collapse = '')
	for (extension in c('json', 'bin')) if (!grepl(paste0('Extension="', extension, '"'), content, fixed = TRUE)) {
		type <- if (extension == 'json') 'application/json' else 'application/octet-stream'
		content <- sub('</Types>', paste0('<Default Extension="', extension, '" ContentType="', type, '"/></Types>'), content, fixed = TRUE)
	}
	# openxlsx reserves the .bin default for printer settings. Exact data
	# exports need their own content type so Excel does not parse them as printers.
	parts <- unique(vapply(manifest$files, `[[`, character(1), 'part'))
	for (part in parts) {
		override <- paste0('<Override PartName="/', part, '" ContentType="application/octet-stream"/>')
		content <- sub('</Types>', paste0(override, '</Types>'), content, fixed = TRUE)
	}
	writeLines(content, content_path, useBytes = TRUE)
	zip::zip_append(path, files = c('[Content_Types].xml', 'le8'), root = stage,
		compression_level = 6, include_directories = FALSE)
	if (any(grepl('/$', utils::unzip(path, list = TRUE)$Name))) stop('Unexpected directory entries in workbook: ', path)
	check <- le8_table_archive_read(path)
	stopifnot(identical(as.character(names(check$files)), as.character(names(entries))))
	for (source in names(check$files)) invisible(le8_table_bytes(check$files[[source]]))
	invisible(path)
}
le8_table_private_save <- function(entry, source, directory, manifest = list()) {
	known <- names(manifest)[vapply(manifest, function(x) identical(x, source), logical(1))]
	stem <- sub('[.](csv|tsv|xlsx|jsonl)([.]gz)?$', '', source)
	filename <- if (length(known)) known[1] else paste0(stem, '.rds')
	# Do not overwrite a model or another existing data object with the same stem.
	if (!length(known) && file.exists(file.path(directory, filename))) {
		filename <- paste0(source, '.rds')
		if (file.exists(file.path(directory, filename))) stop('Private RDS name conflicts with existing data: ', filename)
	}
	data <- entry$data
	if (is.null(data)) data <- entry$sheets
	if (is.null(data)) stop('Private export has no data object: ', source)
	metadata <- entry[c('bytes', 'sha256', 'mtime')]
	metadata$format <- 'le8-private-table-v1' ; metadata$source <- source
	attr(data, 'le8_table_export') <- metadata
	le8_table_save(data, file.path(directory, filename))
	Sys.chmod(file.path(directory, filename), mode = '0600')
	manifest[[filename]] <- source
	manifest
}
le8_table_restore_entry <- function(entry, directory, name) {
	if (basename(name) != name || name %in% c('.', '..')) stop('Unsafe table name: ', name)
	target <- file.path(directory, name)
	if (file.exists(target)) return(invisible(NULL))
	writeBin(le8_table_bytes(entry), target)
	Sys.setFileTime(target, entry$mtime)
	invisible(target)
}

le8_table_sheet_names <- function(names) {
	used <- '_tables'
	vapply(names, function(name) {
		base <- sub('[.](csv|tsv)([.]gz)?$', '', basename(name))
		base <- sub('^c[1-5][.]', '', base)
		base <- sub('^out[.]', '', base)
		base <- sub('^Fig[0-9]+[.][^.]+[.]', '', base)
		base <- sub('^out[.]', '', base)
		base <- sub('^(pgs_focus|focus|questions|final[.]questions)[.]', '', base)
		base <- sub('^question_', '', base)
		base <- sub('^dandelion[._]', '', base, ignore.case = TRUE)
		base <- gsub('[^[:alnum:]_. -]', '_', base)
		if (!nzchar(base)) base <- 'table'
		candidate <- substr(base, 1, 31)
		number <- 1L
		while (tolower(candidate) %in% tolower(used)) {
			number <- number + 1L
			candidate <- paste0(substr(base, 1, 27), '_', number)
		}
		used <<- c(used, candidate)
		candidate
	}, character(1), USE.NAMES = FALSE)
}


# 🚩 Figure result workbooks
le8_table_from_data <- function(x) {
	file <- tempfile('le8-result-', tmpdir = le8_table_tmpdir(), fileext = '.csv')
	on.exit(unlink(file), add = TRUE)
	data.table::fwrite(as.data.frame(x), file, na = 'NA')
	le8_table_entry(file, as.data.frame(x))
}
le8_table_plot_exports <- function(plots, directory, stem) {
	for (i in seq_along(plots)) {
		plot <- plots[[i]]
		tables <- attr(plot, 'le8_result_tables')
		if (is.null(tables)) {
			tables <- c(list(plot$data), lapply(plot$layers, function(layer) layer$data))
			tables <- Filter(function(x) le8_table_useful(x) && !le8_table_private(x) &&
				!all(names(x) %in% c('x', 'y', 'xend', 'yend', 'label', 'colour', 'fill')), tables)
			fingerprints <- vapply(tables, function(x) digest::digest(x, algo = 'sha256'), character(1))
			tables <- tables[!duplicated(fingerprints)]
		}
		for (j in seq_along(tables)) {
			x <- as.data.frame(tables[[j]])
			if (!le8_table_useful(x)) next
			if (le8_table_private(x)) stop('Individual data cannot be exported as figure results')
			x$panel <- LETTERS[i]
			name <- paste0(stem, '.panel_', LETTERS[i], if (length(tables) > 1L) paste0('.', j), '.csv')
			data.table::fwrite(x, file.path(directory, name), na = 'NA')
		}
	}
	invisible(NULL)
}
le8_table_administrative <- function(name) {
	grepl('(^|[./_])(manifest|provenance|output_index|omission_audit|source_registry|sources|input_audit|input_features|input_feature_annotation_audit|analysis_scope|method_scope_audit|design_status|annotation_status|status|completed|panels|caption|claim_register|audit_findings|supplement_source_inventory|review_association_scope|review_interpretation)([./_]|$)|^outcome_audit[.]|^c4[.]focus[.]design[.]|^(final[.]questions[.]|question_)design[.]', name, ignore.case = TRUE)
}
le8_table_useful <- function(x) {
	is.data.frame(x) && nrow(x) > 0 && ncol(x) > 0 &&
		!all(tolower(names(x)) %in% c('note', 'status', 'detail', 'reason', 'message', 'caption', 'exit_code', 'source', 'path', 'file'))
}
le8_table_c3_results <- function(directory, public) {
	if (basename(directory) != 'c3_coloc') return(public)
	# The fitted object already retains every regional SNP and posterior. Export
	# only the loci/variants drawn in the paired figures, without duplicating it.
	cache <- file.path(directory, 'c3.res.rds')
	large <- c('c3.out.xlsx', 'c3.regional_rows.csv', 'c3.variant_posteriors.csv')
	if (any(large %in% names(public$files))) {
		if (!file.exists(cache)) stop('Cannot replace full C3 exports without the fitted result: ', directory)
		obj <- readRDS(cache)
		if (!all(c('summary', 'regional', 'variants') %in% names(obj))) stop('Incomplete C3 result: ', cache)
		x <- as.data.frame(obj$summary)
		x <- x[x$status %in% 'ok', , drop = FALSE]
		robust <- if ('PP.H4_robust_min' %in% names(x)) x$PP.H4_robust_min else rep(NA_real_,nrow(x))
		x <- x[order(-robust, -x$PP.H4, na.last = TRUE), , drop = FALSE]
		loci <- head(x, 4L)
		key <- function(z) paste(z$feature, z$locus, sep = '\r')
		regional <- as.data.frame(obj$regional)
		regional <- regional[key(regional) %in% key(loci) & is.finite(regional$p) & regional$p > 0, , drop = FALSE]
		variants <- as.data.frame(obj$variants)
		variants <- variants[key(variants) %in% key(loci), , drop = FALSE]
		if (nrow(variants)) {
			variants <- variants[order(variants$feature, variants$locus, -variants$SNP.PP.H4), , drop = FALSE]
			variants <- variants[ave(seq_len(nrow(variants)), key(variants), FUN = seq_along) <= 40L, , drop = FALSE]
		}
		tables <- list(regional_top_loci = regional, credible_set_variants = variants,
			GPU_results = obj$GPU_coloc$results, GPU_status = obj$GPU_coloc$status)
		for (name in names(tables)) if (is.data.frame(tables[[name]]) && ncol(tables[[name]]))
			public$files[[paste0('c3.', name, '.csv')]] <- le8_table_from_data(tables[[name]])
		public$files[intersect(large, names(public$files))] <- NULL
	}
	public$files[grep('^c3[.]CIGMA_|^qtl_cad_manifest[.]', names(public$files), value = TRUE)] <- NULL
	public
}
le8_table_figure_pattern <- function(file) {
	stem <- sub('[.]png$', '', basename(file))
	key <- sub('^(c[1-5][.])?Fig[0-9]+[.]', '', stem)
	module <- substr(stem, 1, 2)
	patterns <- switch(module,
		c1 = c(
			mh = '(pwas|mwas)_(pgs_|incident_adj2|prevalent_adj2|birthline_attained_age_adj2)',
			circular = 'circular_associations|mwas_(pgs_|incident_adj2|prevalent_adj2|birthline_attained_age_adj2)',
			vc = '(pwas|mwas)_(pgs_|incident_adj2|prevalent_adj2|birthline_attained_age_adj2)',
			measured_volcano = '(pwas|mwas)_(incident_adj2|prevalent_adj2)',
			temporal_profiles = 'gradient|mock_trajectories', diagnosis_timed_profiles = 'mock_(trajectories|clusters)',
			gradient_cluster = 'cluster_membership|cluster_selection|mock_(trajectories|clusters)',
			gradient_cluster_diagnostics = 'cluster_membership|cluster_selection|mock_clusters',
			quantile_top = 'quantile', enrich_sig = 'enrichment_(incident|prevalent)_sig',
			temporal_evidence = 'directionality_triage|temporal_heterogeneity',
			landmark_birthline_sensitivity = 'landmark_adj2|birthline_attained_age_adj2',
			diagnosis_window_riskset = 'diagnosis_window_riskset', paired_temporal_validation = 'paired_pgs_measured|pgs_actual_concordance|prevalent_duration_adj2',
			reverse_time_exploratory = 'prevalent_reverse_cox_adj2|prevalent_duration_adj2',
			pgs_actual_concordance = 'pgs_actual_concordance|paired_pgs_measured',
			L_VLDL_TG_pct_deep_dive = 'L_VLDL_TG_pct'
		),
		c2 = c(
			instrument_diagnostics = 'instrument|pQTL|heritability', effect_concordance = 'cis_trans_comparison|MR_all',
			mr_incident_prevalent = 'mr_incident_prevalent', mrlink2 = 'MRLink2_results|mrlink2[.]results',
			bidirectional_mr = 'bidirectional_mr|reverse_MR_all|MR_all', genetic_decomposition = 'individual_genetic_decomposition[.]csv|decomposition_calibration',
			genetic_component_leadtime = 'genetic_component_leadtime', evidence_grades = 'evidence_grades|top_candidates',
			dandelion_sensitivity = 'dandelion.*(targets|integration|sensitivity|scores)',
			dandelion_mr_integration = 'dandelion.*(integration|MR_integration|targets)',
			dandelion_network = 'dandelion.*(pairs|lead_snps|snp_gene_map|targets)',
			prots.top = 'top_candidates|MR_best|observational|cis_trans_comparison',
			sensitivity_architecture = 'heterogeneity|pleiotropy|leave_one|sensitivity|MR_all',
			directionality_causal = 'directionality_causal', observational_mr_overview = 'mr_incident_prevalent|observational|MR_all',
			qtl_variance_ranked = 'pQTL|heritability|instrument'
		),
		c3 = c(
			posterior_evidence = 'coloc_summary|selected_mr|mr_locus_selection|same_locus|susie',
			regional_top_loci = 'regional_top_loci', credible_sets = 'credible_set',
			gpu_coloc_validation = 'GPU_', pgs_coloc_triangulation = 'pgs_|triangulation[.](features|loci)'
		),
		c4 = c(
			deployed_concept_fidelity = 'deployed_concept_fidelity',
			supervised_connections = 'LE8_feature_associations|primary_pillar_assignment|supervised_module',
			mediation = 'mediation', imaging_context = 'imaging_associations|imaging_fields',
			lifestyle_omics_risk = 'lifestyle_omics_risk_display', sex_interaction = 'sex_interaction',
			LE8_component_interactions = 'LE8_pairwise_interactions', spline_patterns = 'nonlin_(curves|tests)',
			pass_fail_penalty = 'penalty_(summary|CV_by_fold|selected_gamma)',
			proxy_heatmap = 'LE8_feature_associations|proxy_membership',
			connection_bridge = 'genetic_omic_disease_bridges|matched_PGS_bridges',
			connection_evidence = 'genetic_omic_disease_bridges|matched_PGS_scan|PRS_feature_associations',
			network_globe = 'YS_edges|proxy_membership|supervised_module_membership',
			selection_mediation = 'mediation|primary_pillar_assignment|proxy_membership',
			imaging_atlas = 'imaging_associations|mock_brain_region_counts',
			state_remodeling = 'state_network_(edges|counts|hubs)',
			equal_budget = 'focus[.](metrics|contrasts|calibration|panel_members|pillar_counts|heterogeneity|proxy_accuracy)[.]'
		),
		c5 = c(cell.enrichment = 'cell[.](enrichment|coverage|panel_annotation)'),
		NULL
	)
	if (module == 'c5') key <- sub('^c5[.]', '', stem)
	if (key %in% names(patterns)) return(unname(patterns[key]))
	if (grepl('^Fig[0-9S]+$', stem)) return(paste0('^', stem, '[_./]'))
	if (grepl('^Fig6[.]question_LE8', stem)) return('final[.]questions[.](contrasts|proxy|pillars|prediction)[.]')
	if (grepl('^Fig7[.]question_ABM', stem)) return('final[.]questions[.]abm_(metrics|coverage|paired|support)[.]')
	if (grepl('^Fig9[.]question_ABM_attention', stem)) return('final[.]questions[.]abm_(metrics|coverage|paired|registry|release_audit|fit_status)[.]')
	if (grepl('^Fig8[.]question_genetic', stem)) return('final[.]questions[.](genetic|temporal|same_locus)[.]')
	if (grepl('coverage', stem, ignore.case = TRUE)) return('coverage_curve|support_error')
	if (stem == 'Fig_masked_reconstruction') return('masked_feature_metrics|masked_reconstruction')
	if (stem == 'Fig_model_comparison') return('test_metrics|approach_comparison|paired_contrasts')
	if (grepl('paired|prediction|benchmark', stem, ignore.case = TRUE)) return('paired_contrasts|test_metrics|approach_comparison')
	if (grepl('support|error', stem, ignore.case = TRUE)) return('support_error|coverage_curve')
	gsub('[.]', '[.]', key)
}
le8_table_result_groups <- function(names, directory) {
	if (!length(names)) return(character())
	module <- if (grepl('^c[1-5]_', basename(directory))) sub('_.*$', '', basename(directory)) else basename(directory)
	# Rules are ordered by scientific specificity, not input order or table count.
	rules <- switch(module,
		c1 = c(
			pgs.temporal = 'pgs_focus[.](windows|window_contrasts|landmarks|landmark_contrasts)[.]',
			pgs.decomposition = 'pgs_focus[.](bootstrap|calibration|components|contrasts|folds|lifestyle|composition)[.]',
			pgs.comparison = 'pgs_focus[.]',
			cohort = 'cohort|endpoint_',
			enrichment = 'enrich|mock_(function_terms|gene_universe|ppi_edges|tf_edges)',
			pgs = 'pgs_',
			temporal = 'birthline|landmark|window|duration|reverse_prevalent|directionality|time_heterogeneity|trajector|profiles',
			association = 'pwas_|mwas_|assoc|prevalent|incident|attenuation'
		),
		c2 = c(
			dandelion = 'dandelion|review_dan_', mrlink2 = 'mrlink2',
			instruments = 'QTL_R2|instrument|heritability',
			genetic_leadtime = 'genetic_leadtime', genetic_decomposition = 'decomp',
			mr = 'MR_|reverse_MR|effect_forest|cis_local_vs_trans', evidence = 'evidence'
		),
		c4 = c(
			validation.genetics = 'focus[.]PGS_',
			validation.models = 'focus[.](baseline_hazards|model_coefficients|preprocessing|fit_diagnostics)[.]',
			validation.panels = 'focus[.](panel_|membership|training_screen|domain_availability|inflammation_definition)',
			validation = 'focus[.]', explain = 'explain[.]', cohort = 'cohort',
			association = 'assoc|behavior_biology', connections = 'primary_assignment|genetic_omic_bridges',
			splines = 'spline|nonlin|nadir', interactions = 'interaction|prediction_surfaces',
			penalty = 'penalty', networks = 'state_network', imaging = 'imaging'
		),
		c5 = c(cellulation = 'cell'),
		abm_reference = c(validation = 'test_|subgroup|development_metrics',
			development = 'tuning|oof_fits|embedding|metric_features|mosaic_weights|token_membership',
			selective_training = 'c1[.]selective[.](training_comparison|model_tuning|crossfit_audit)',
			selective_validation = 'c1[.]selective[.](risk_stratified_gain|gate_diagnostics|audit_contrasts|decision_curve)'),
		abm_selective_attention = c(validation = 'test_|coverage|paired|development_metrics|audit_decision|audit_contrasts|censoring_sensitivity',
			training = 'model_registry|model_tuning|gate_training|crossfit_audit|fit_status',
			inputs = 'preprocessing_audit|role_audit|reference_diagnostics'),
		abm_tabicl = c(model = 'test_|tuning|learning_curve|feature_selection'),
		attention = c(interventions = 'intervention'),
		le8_annotations = c(protein_interactions = 'string_physical'),
		Yin = c(connections = 'connection'), YinYang = c(connections = 'connection'),
		NULL
	)
	if (module %in% c('final', 'shiny')) rules <- c(
		catalogue = '^(tables|figures|status|question_sources)[.]',
		abm = '(questions[.]|question_)abm_', cell = '(questions[.]|question_)cell',
		dandelion = 'dandelion', nonlinear = 'nonlinear|age_models', mediation = 'mediation',
		state_projection = 'state_projection',
		temporal = '(questions[.]|question_)temporal',
		genetics = '(^loci[.]|[._](genetic|same_locus|mr_signal_evidence|mr_scope|susie_pairs|susie_diagnostics|decomposition)[.])',
		connections = '[._](members|modules|pillars|inflammation_definition|concept_coefficients|concept_status|concept_fold_panels)[.]',
		validation = 'prediction|proxy|deployed_concept_fidelity|[._](contrasts|fit|heterogeneity|design)[.]',
		overview = 'cohort|overview|association_counts|discovery_counts',
		association = '^(candidates|effects)[.]'
	)
	group <- rep('results', length(names))
	for (key in rev(names(rules))) group[grepl(rules[[key]], names, ignore.case = TRUE)] <- key
	paste0(module, '.', group)
}

le8_table_figure_supplements <- function(figure, sources) {
	# Figure numbers can change at publication. Match the scientific topic in
	# explicit figure exports so their tables follow the corresponding image.
	stem <- sub('[.]png$', '', basename(figure))
	topic <- sub('^c[1-5][.]Fig[0-9]+[.]', '', stem)
	module <- substr(stem, 1, 2)
	if (!grepl('^c[1-5][.]Fig', stem)) return(character())
	topic_alias <- switch(topic, enrich_sig = 'enrichment', diagnosis_timed_profiles = 'profiles', topic)
	explicit <- paste0('^', module, '[.]Fig[0-9]+[.]', gsub('[.]', '[.]', topic_alias), '[.]')
	extra <- switch(paste(module, topic, sep = '.'),
		c1.enrich_sig = 'enrichment_(incident|prevalent)|mock_(function_terms|gene_universe|ppi_edges|tf_edges)',
		c1.paired_temporal_validation = 'review_paired_associations|review_time_heterogeneity',
		c1.temporal_evidence = 'out[.]directionality[.]',
		c1.reverse_time_exploratory = 'out[.]reverse_prevalent[.]',
		c2.instrument_diagnostics = 'out[.]QTL_R2[.]',
		c2.mrlink2 = 'mrlink2_(job_audit|jobs)|out[.]MRLink2_(job_audit|jobs)',
		c2.bidirectional_mr = 'out[.]MR_reverse[.]|reverse_MR_audit',
		c2.genetic_decomposition = 'genetic_decomp_summary|review_decomp_folds',
		c2.evidence_grades = 'evidence_overlap',
		c2.effect_concordance = 'effect_forest|cis_local_vs_trans_distal',
		c4.supervised_connections = 'out[.]primary_assignment[.]',
		c4.connection_bridge = 'out[.]genetic_omic_bridges[.]',
		c4.pass_fail_penalty = 'penalty_inner_CV',
		c4.state_remodeling = 'state_network_counts',
		NULL
	)
	sources[grepl(paste(c(explicit, extra), collapse = '|'), sources, ignore.case = TRUE)]
}

le8_table_workbook <- function(directory, public) {
	if (!length(public$files) && !length(public$private_files)) return(invisible(character()))
	public <- le8_table_c3_results(directory, public)
	if (basename(directory) == 'final') {
		# These are copies of module sources, not new Final analysis results.
		copied <- grepl('^[^.]+[.](prot|met|joint)[.]', names(public$files)) |
			grepl('^final[.]questions[.]panome_', names(public$files))
		public$files[copied] <- NULL
	}
	# Flatten old bundles once. Never embed a whole workbook inside another.
	for (source in names(public$files)) {
		if (basename(directory) == 'c3_coloc' && grepl('[.]xlsx$', source)) {
			public$files[[source]] <- NULL
			next
		}
		public$files[[source]] <- le8_table_materialize(public$files[[source]], source)
	}
	weights <- 'c2.genetic_score_weights.tsv'
	if (weights %in% names(public$files)) {
		public$private_files <- le8_table_private_save(public$files[[weights]], weights, directory, public$private_files)
		public$files[[weights]] <- NULL
	}
	tables <- lapply(public$files, `[[`, 'data')
	tables <- tables[vapply(tables, is.data.frame, logical(1))]
	for (source in names(public$files)) if (length(public$files[[source]]$sheets)) {
		for (name in names(public$files[[source]]$sheets)) {
			x <- public$files[[source]]$sheets[[name]]
			if (!le8_table_useful(x) || le8_table_administrative(name)) next
			same <- vapply(tables, function(y) identical(names(x), names(y)) && identical(dim(x), dim(y)) &&
				isTRUE(all.equal(x, y, check.attributes = FALSE, tolerance = 0)), logical(1))
			if (any(same)) next
			key <- paste0(sub('[.]xlsx$', '', source), '.', name, '.csv')
			if (key %in% names(public$files)) next
			public$files[[key]] <- le8_table_from_data(x)
			tables[[key]] <- x
		}
		public$files[[source]] <- NULL
	}
	if (any(vapply(tables, le8_table_private, logical(1)))) stop('Private table in aggregate workbook: ', directory)
	visible <- vapply(tables, le8_table_useful, logical(1)) & !le8_table_administrative(names(tables))
	# Shiny is a machine-readable catalogue, with no corresponding PNGs.
	if (basename(directory) == 'shiny') visible <- vapply(tables, le8_table_useful, logical(1))
	tables <- tables[visible]
	figures <- sort(list.files(directory, pattern = '[.]png$'))
	plans <- list()
	for (figure in figures) {
		selected <- names(tables)[startsWith(names(tables), paste0(sub('[.]png$', '', figure), '.panel_'))]
		if (!length(selected)) selected <- names(tables)[grepl(le8_table_figure_pattern(figure), names(tables), ignore.case = TRUE)]
		if (!length(selected)) {
			pattern <- le8_table_figure_pattern(figure)
			expected <- names(public$files)[grepl(pattern,names(public$files),ignore.case=TRUE)]
			empty <- length(expected) && all(vapply(public$files[expected],function(z) !le8_table_useful(z$data),logical(1)))
			if (!empty) stop('No analysis table mapped to figure: ', file.path(directory, figure))
			selected <- paste0(sub('[.]png$','',figure),'.availability.csv')
			tables[[selected]] <- data.frame(result_table=expected,availability='No analyzable rows; figure contains an unavailable-results panel')
			public$files[[selected]] <- le8_table_from_data(tables[[selected]])
		}
		plans[[sub('[.]png$', '.xlsx', figure)]] <- selected
	}
	for (figure in figures) {
		# Unplotted full results of a paginated topic get their own topic book;
		# never copy them into every page just because the topic matches.
		topic <- sub('^(c[1-5][.])?Fig[0-9]+[.]', '', figure)
		if (sum(sub('^(c[1-5][.])?Fig[0-9]+[.]', '', figures) == topic) != 1L) next
		remaining <- setdiff(names(tables), unique(unlist(plans)))
		file <- sub('[.]png$', '.xlsx', figure)
		plans[[file]] <- c(plans[[file]], le8_table_figure_supplements(figure, remaining))
	}
	remaining <- setdiff(names(tables), unique(unlist(plans)))
	groups <- le8_table_result_groups(remaining, directory)
	for (group in unique(groups)) {
		plans[[paste0(group, '.xlsx')]] <- remaining[groups == group]
	}
	if (!length(plans)) {
		# A non-applicable cell analysis has a scientific eligibility result.
		eligibility <- 'c5.cellulation_status.csv'
		if (eligibility %in% names(public$files)) {
			tables[[eligibility]] <- public$files[[eligibility]]$data
			plans[['c5.cellulation.xlsx']] <- eligibility
		} else stop('No analysis results to publish: ', directory)
	}
	primary <- names(plans)[1]
	owners <- setNames(rep(primary, length(public$files)), names(public$files))
	for (file in rev(names(plans))) owners[intersect(plans[[file]], names(owners))] <- file
	for (file in names(plans)) {
		selected <- tables[plans[[file]]]
		# Identical source aliases serve old internal consumers but need one sheet.
		fingerprints <- vapply(selected, function(x) digest::digest(x, algo = 'sha256'), character(1))
		selected <- selected[!duplicated(fingerprints)]
		le8_table_write_workbook(directory, selected, file, public$files[owners == file],
			if (file == primary) public$private_files else list())
	}
	invisible(names(plans))
}
le8_table_write_workbook <- function(directory, tables, filename, entries, private_files) {
	if (length(tables) && all(grepl('[.]panel_', names(tables)))) {
		# Panels with the same result schema share a worksheet, identified by
		# the panel column; six survival panels need curves + contrasts, not 12 tabs.
		groups <- vapply(tables, function(x) paste(names(x), collapse = '\r'), character(1))
		merged <- list()
		for (key in unique(groups)) {
			x <- as.data.frame(data.table::rbindlist(tables[groups == key], use.names = TRUE))
			name <- if (all(c('time', 'surv') %in% names(x))) 'curves' else
				if (all(c('contrast', 'HR') %in% names(x))) 'contrasts' else 'panel_results'
			if (name %in% names(merged)) name <- paste0(name, '_', length(merged) + 1L)
			merged[[name]] <- x
		}
		tables <- merged
	}
	if (length(tables) > 1L && !all(names(tables) %in% c('curves', 'contrasts'))) {
		groups <- vapply(tables, function(x) paste(names(x), collapse = '\r'), character(1))
		merged <- list()
		for (key in unique(groups)) {
			selected <- tables[groups == key]
			name <- names(selected)[1]
			if (length(selected) == 1L || 'result_table' %in% names(selected[[1]])) {
				merged <- c(merged, selected)
				next
			}
			x <- as.data.frame(data.table::rbindlist(selected, use.names = TRUE, idcol = 'result_table'))
			x$result_table <- sub('[.](csv|tsv)([.]gz)?$', '', x$result_table)
			name <- if ('beta_difference' %in% names(x)) 'contrasts' else
				if ('attenuation_pct' %in% names(x)) 'attenuation' else
				if ('adjusted_p' %in% names(x)) 'enrichment' else
				if (all(c('beta', 'landmark') %in% names(x))) 'time_effects' else
				if (any(c('p.value', 'pval', 'beta', 'Beta') %in% names(x))) 'associations' else
				if ('metric' %in% names(x)) 'metrics' else 'results'
			if (name %in% names(merged)) name <- paste0(name, '_', length(merged) + 1L)
			merged[[name]] <- x
		}
		tables <- merged
	}
	sheets <- le8_table_sheet_names(names(tables))
	if (!length(tables)) stop('Cannot publish an empty result workbook: ', filename)
	cells <- sum(vapply(tables, function(x) as.double(nrow(x)) * ncol(x), numeric(1)))
	if (cells > 3000000 || any(vapply(tables, function(x) nrow(x) > 1048575L || ncol(x) > 16384L, logical(1))))
		stop('Result table too large for its figure workbook: ', filename, '; split by analysis scope')
	wb <- openxlsx::createWorkbook()
	for (i in seq_along(tables)) {
		openxlsx::addWorksheet(wb, sheets[i])
		x <- tables[[i]]
		if (ncol(x)) {
			# List-valued annotations remain explicit instead of being discarded.
			for (name in names(x)) if (is.list(x[[name]])) x[[name]] <- vapply(x[[name]], function(z) paste(z, collapse = '; '), character(1))
			openxlsx::writeData(wb, sheets[i], x, withFilter = nrow(x) > 0)
			openxlsx::freezePane(wb, sheets[i], firstRow = TRUE)
		}
	}
	target <- file.path(directory, filename)
	temporary <- tempfile('le8-workbook-', tmpdir = le8_table_tmpdir(), fileext = '.xlsx')
	on.exit(unlink(temporary), add = TRUE)
	openxlsx::saveWorkbook(wb, temporary, overwrite = TRUE)
	if (!identical(openxlsx::getSheetNames(temporary), sheets)) stop('Workbook verification failed: ', target)
	le8_table_archive_add(temporary, entries, private_files)
	if (!file.copy(temporary, target, overwrite = TRUE)) stop('Cannot publish ', target)
	if (!identical(unname(tools::md5sum(temporary)), unname(tools::md5sum(target)))) stop('Workbook copy verification failed: ', target)
	invisible(target)
}
le8_table_stores <- function(root) {
	paths <- list.files(root, pattern = '[.]xlsx$', recursive = TRUE, full.names = TRUE)
	paths <- paths[!grepl('/(_history|_source_figures|_previous|le8_annotations)/', paths)]
	as.character(unlist(lapply(unique(dirname(paths)), le8_table_archive_files), use.names = FALSE))
}
le8_tables_pack <- function(root, clean = TRUE) {
	root <- normalizePath(root, winslash = '/', mustWork = TRUE)
	paths <- list.files(root, pattern = '[.](csv|tsv|jsonl)([.]gz)?$|[.]xlsx$|[.]json$', recursive = TRUE, full.names = TRUE, all.files = TRUE)
	internal <- grepl('^c5[.].*[.]json$',basename(paths)) | grepl('/abm_(reference|tabicl|selective_attention)(/attention)?/[^/]+[.]json$', paths)
	paths <- paths[!grepl('[.]json$', paths) | internal]
	paths <- paths[!grepl('/(_history|_source_figures|_previous|le8_annotations)/', paths)]
	directories <- sort(unique(c(dirname(paths), dirname(le8_table_stores(root)))))
	for (directory in directories) {
		old_workbooks <- le8_table_archive_files(directory)
		public <- le8_table_load(directory)
		changed <- any(vapply(old_workbooks, function(path) !identical(le8_table_archive_read(path)$layout, 'topic-workbooks-v2'), logical(1)))
		paired <- sub('[.]png$', '.xlsx', list.files(directory, pattern = '[.]png$'))
		changed <- changed || any(!file.exists(file.path(directory, paired)))
		files <- paths[dirname(paths) == directory]
		files <- setdiff(files, old_workbooks)
		if (Sys.getenv('LE8_TABLE_WORKSPACE') == '1') {
			removed <- setdiff(names(public$files), basename(files))
			if (length(removed)) { public$files[removed] <- NULL ; changed <- TRUE }
		}
		intermediate <- any(strsplit(directory, '/', fixed = TRUE)[[1]] %in% c('prepared', 'neural', 'quality_neural'))
		for (path in files) {
			name <- basename(path) ; stored <- public$files[[name]]
			if (!is.null(stored$sha256) && identical(stored$sha256, digest::digest(file = path, algo = 'sha256'))) next
			changed <- TRUE
			if (grepl('[.]json$', name)) {
				entry <- le8_table_entry(path) ; entry$kind <- 'exchange' ; is_private <- FALSE
			} else if (grepl('[.]jsonl([.]gz)?$', name)) {
				connection <- if (grepl('[.]gz$', name)) gzfile(path, open = 'rt') else file(path, open = 'rt')
				x <- tryCatch(jsonlite::stream_in(connection, verbose = FALSE), finally = close(connection))
				entry <- le8_table_entry(path, x) ; is_private <- TRUE
			} else if (grepl('[.]xlsx$', name)) {
				sheets <- le8_table_sheets(path)
				is_private <- any(vapply(sheets, le8_table_private, logical(1)))
				entry <- le8_table_entry(path) ; entry$sheets <- sheets
			} else {
				separator <- if (grepl('[.]tsv([.]gz)?$', path)) '\t' else ','
				header <- if (le8_table_blank(path)) data.frame() else data.table::fread(path, sep = separator, nrows = 0, showProgress = FALSE)
				is_private <- le8_table_private(header)
				x <- if (!intermediate || is_private) le8_table_read(path) else NULL
				entry <- le8_table_entry(path, x)
			}
			if (is_private) {
				public$private_files <- le8_table_private_save(entry, name, directory, public$private_files)
				public$files[[name]] <- NULL
			} else public$files[[name]] <- entry
		}
		if (changed || !length(old_workbooks)) {
			written <- le8_table_workbook(directory, public)
			unlink(setdiff(old_workbooks, file.path(directory, written)))
		}
		# A source XLSX may have the same name as its newly consolidated figure
		# workbook. Never delete the published replacement while cleaning inputs.
		if (clean && length(files)) unlink(setdiff(files,le8_table_archive_files(directory)))
		message('Packed result workbooks: ', directory)
		invisible(gc())
	}
	invisible(length(directories))
}
le8_tables_restore <- function(root) {
	paths <- le8_table_stores(root)
	paths <- paths[!duplicated(dirname(paths))]
	for (path in paths) {
		x <- le8_table_load(dirname(path))
		for (name in names(x$files)) le8_table_restore_entry(x$files[[name]], dirname(path), name)
		if (grepl('/abm_(reference|tabicl|selective_attention)/', path) && Sys.getenv('LE8_TABLE_ABM_PRIVATE', 'TRUE') != 'TRUE') next
		for (filename in names(x$private_files)) {
			if (basename(filename) != filename) stop('Unsafe private table path')
			data <- readRDS(file.path(dirname(path), filename))
			entry <- attr(data, 'le8_table_export')
			if (!identical(entry$format, 'le8-private-table-v1') || !identical(entry$source, x$private_files[[filename]])) stop('Private table metadata mismatch: ', filename)
			le8_table_restore_entry(entry, dirname(path), entry$source)
		}
	}
	invisible(length(paths))
}


# 🚩 Portable aggregate viewer
le8_tables_share <- function(root, destination) {
	root <- normalizePath(root, winslash = '/', mustWork = TRUE)
	dir.create(destination, recursive = TRUE, showWarnings = FALSE)
	stores <- new.env(parent = emptyenv())
	load <- function(directory) {
		if (!exists(directory, envir = stores, inherits = FALSE)) assign(directory, le8_table_load(directory), envir = stores)
		get(directory, envir = stores)
	}
	index <- load(file.path(root, 'shiny'))
	if (!length(index$files)) stop('Run final --index-only before creating a viewer package')
	read_index <- function(name) {
		entry <- index$files[[paste0(name, '.csv')]]
		if (is.null(entry)) stop('Missing viewer index: ', name)
		le8_table_materialize(entry, paste0(name, '.csv'))$data
	}
	tables <- read_index('tables') ; figures <- read_index('figures')
	resolve <- function(relative) {
		if (is.na(relative) || grepl('(^|/)[.][.]?(/|$)|^[A-Za-z]:|^/|\\\\', relative)) stop('Unsafe viewer path: ', relative)
		if (grepl('(^|/)([_.]|private|cache|input|neural|checkpoints)', relative, ignore.case = TRUE)) stop('Private viewer path: ', relative)
		path <- file.path(root, relative)
		directory <- normalizePath(dirname(path), winslash = '/', mustWork = TRUE)
		if (!startsWith(paste0(directory, '/'), paste0(root, '/'))) stop('Viewer path escapes the analysis root')
		file.path(directory, basename(path))
	}
	selected <- list(shiny = index$files)
	for (i in seq_len(nrow(tables))) {
		path <- resolve(tables$path[i])
		entry <- load(dirname(path))$files[[basename(path)]]
		if (is.null(entry) || !identical(entry$kind, 'table') || !identical(entry$sha256, tables$sha256[i])) stop('Viewer table is missing or changed: ', tables$path[i])
		directory <- dirname(tables$path[i])
		selected[[directory]][[basename(path)]] <- entry
	}
	for (i in seq_len(nrow(figures))) {
		path <- resolve(figures$path[i])
		if (!grepl('[.](png|jpg|jpeg)$', path, ignore.case = TRUE)) stop('Unexpected viewer image: ', path)
		if (!identical(digest::digest(file = path, algo = 'sha256'), figures$sha256[i])) stop('Viewer image changed: ', path)
		target <- file.path(destination, figures$path[i]) ; dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
		if (!file.copy(path, target, overwrite = TRUE)) stop('Cannot copy viewer image: ', path)
	}
	for (directory in names(selected)) {
		entries <- selected[[directory]]
		for (name in names(entries)) {
			entry <- le8_table_materialize(entries[[name]], name)
			if (!identical(entry$kind, 'table') || !is.data.frame(entry$data) || le8_table_private(entry$data)) stop('Non-aggregate source cannot be shared: ', name)
			entries[[name]] <- entry
		}
		target <- file.path(destination, directory) ; dir.create(target, recursive = TRUE, showWarnings = FALSE)
		le8_table_workbook(target, list(files = entries, private_files = list()))
		message('Shared aggregate tables: ', directory)
	}
	invisible(destination)
}

if (sys.nframe() == 0L && length(commandArgs(TRUE)) && commandArgs(TRUE)[1] %in% c('--tables-pack', '--tables-restore', '--tables-share')) {
	args <- commandArgs(TRUE)
	if (args[1] == '--tables-share') {
		if (length(args) != 3L) stop('Share requires an analysis root and destination')
		le8_tables_share(args[2], args[3])
		quit(save = 'no', status = 0)
	}
	if (length(args) != 2L) stop('Table storage commands require an analysis root')
	if (args[1] == '--tables-restore') {
		if (!startsWith(normalizePath(args[2], mustWork = TRUE), '/tmp/')) stop('Restore destination must be under /tmp')
		le8_tables_restore(args[2])
	} else le8_tables_pack(args[2])
	quit(save = 'no', status = 0)
}


# 🚩 Analysis settings and reusable results
# Shared R configuration, input preparation, caching and plotting for C1--C5/Final.
le8_y_date <- function(outcome = Y) { value <- Sys.getenv('LE8_Y_DATE', ''); if(nzchar(value)) value else paste0('fod_icd10_', outcome) }
le8_follow_end <- function() {
	value <- Sys.getenv("DATE_FOLLOW_END","2023-04-01")
	date <- tryCatch(as.Date(value),error=function(e) as.Date(NA))
	if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$",value) || is.na(date) || format(date,"%Y-%m-%d")!=value) stop("DATE_FOLLOW_END must be YYYY-MM-DD")
	date
}
le8_analysis_options <- function(outcome = Y) list(
	end_date=as.character(le8_follow_end()),
	group_column=Sys.getenv("LE8_GROUP_COLUMN",Sys.getenv("PGS_GROUP_COLUMN","")),
	group_hash=if(nzchar(Sys.getenv("LE8_GROUP_FILE",""))) unname(tools::md5sum(Sys.getenv("LE8_GROUP_FILE"))) else "participant IDs",
	outer_roster_hash=if(nzchar(Sys.getenv("LE8_OUTER_ROSTER",""))) unname(tools::md5sum(Sys.getenv("LE8_OUTER_ROSTER"))) else "not paired across modules",
	y_date = le8_y_date(outcome), vars_adj = unique(Filter(nzchar, trimws(strsplit(Sys.getenv("LE8_VARS_ADJ", ""), "[,[:space:]]+")[[1]]))),
	white_only = Sys.getenv('LE8_WHITE_ONLY', 'TRUE'), baseline_contract = get0('LE8_BASELINE_VERSION', ifnotfound = '2026-09-21.baseline-med-source-v3'),
	baseline_med_columns = Sys.getenv('LE8_BASELINE_MED_COLUMNS', ''),
	baseline_map = if (nzchar(Sys.getenv('LE8_BASELINE_MAP', ''))) tools::md5sum(Sys.getenv('LE8_BASELINE_MAP')) else ''
)

# Validate completed results without raw participant input or model fitting.
le8_code_fingerprint <- function() {
	code <- list.files(Sys.getenv("LE8_FDIR",file.path(Sys.getenv("DIRSCRIPT"),"f")),pattern="[.](R|py|sh)$",full.names=TRUE)
	f <- tempfile(); on.exit(unlink(f), add=TRUE)
	base::saveRDS(get0(".le8_loaded_code",ifnotfound=tools::md5sum(code)), f, version=2, compress=FALSE)
	unname(tools::md5sum(f))
}

# One policy for all locus evidence consumers; region evidence is never causal proof.
le8_evidence_policy <- function() {
	p <- list(version="locus-policy-v2",posterior=as.numeric(Sys.getenv("C3_H4","0.70")),
		mr_fdr=as.numeric(Sys.getenv("C3_MR_FDR","0.05")),prior_complete_required=TRUE,
		level=Sys.getenv("LE8_EVIDENCE_LEVEL","region_or_signal"))
	if (!is.finite(p$posterior) || p$posterior<=0 || p$posterior>=1 || !is.finite(p$mr_fdr) || p$mr_fdr<=0 || p$mr_fdr>=1 || !p$level %in% c("region_or_signal","signal_only")) stop("Invalid evidence policy")
	p$hash <- digest::digest(p,algo="sha256"); p
}

le8_module_policy <- function(module) {
	prefix <- switch(module,c1_correlate="^C1_",c2_cause="^(C2_|DANDELION_|MRLINK2_)",c3_coloc="^(C3_|COLOC_|LE8_REFERENCE_|LE8_GWAS_MANIFEST)",
		c4_connect="^C4_",c4_panel_validation="^C4_",final_prediction="^FINAL_","^$")
	env <- Sys.getenv(); settings <- env[grepl(prefix,names(env))]
	paths <- unname(settings[file.exists(settings)])
	list(settings=settings,files=if(length(paths)) tools::md5sum(paths) else character(),
		evidence=if(module %in% c("c2_cause","c3_coloc","c4_connect","c4_panel_validation","final_prediction")) le8_evidence_policy() else NULL)
}

le8_completed_results <- function(root, traits, layers, modules) {
	stores <- new.env(parent = emptyenv())
	table_entry <- function(path) {
		store <- dirname(path)
		if (!length(le8_table_archive_files(store))) return(NULL)
		if (!exists(store, envir = stores, inherits = FALSE)) assign(store, le8_table_load(store), envir = stores)
		entry <- get(store, envir = stores)$files[[basename(path)]]
		if (is.null(entry)) NULL else le8_table_materialize(entry, basename(path))
	}
	read_export <- function(path) {
		if (file.exists(path)) return(read.csv(path, stringsAsFactors = FALSE))
		table_entry(path)$data
	}
	required <- list(
		c1_correlate = c("association", "prevalent"),
		c2_cause = c("MR", "MR_reverse", "MR_best"), c3_coloc = c("summary", "variants"),
		c4_connect = c("scan", "membership", "mediation"), c4_panel_validation = "tables"
	)
	exported <- list(
		c1_correlate = c("c1.cohort.csv", "c1.directionality_triage.csv"),
		c2_cause = c("c2.MR_all.csv", "c2.evidence_grades.csv"),
		c3_coloc = c("c3.coloc_summary.csv", "c3.credible_set_audit.csv"),
		c4_connect = c("c4.LE8_feature_associations.csv", "c4.proxy_membership_YS_YSP_NS.csv", "c4.sex_interaction.csv"),
		c4_panel_validation = c("c4.focus.metrics.csv", "c4.focus.panel_members.csv")
	)
	rows <- list()
	for (y in traits) for (b in layers) for (m in intersect(modules, names(required))) {
		d <- file.path(root, y, b, if (m == 'c4_panel_validation') 'c4_connect' else m)
		f <- file.path(d, if (m == 'c4_panel_validation') 'c4.validation.rds' else paste0(substr(m, 1, 2), '.res.rds'))
		status <- "missing" ; detail <- "No completed result" ; generated <- version <- ""
		if (file.exists(f) && file.info(f)$size > 0) {
			x <- tryCatch(readRDS(f), error = function(e) stop("Unreadable completed cache: ", f, ": ", conditionMessage(e)))
			meta <- x$meta
			if (!identical(meta$trait, y) || !identical(meta$layer, if (b == "prot") "protein" else "metabolite") || !(identical(meta$module, m) || (m == "c4_panel_validation" && identical(meta$module, "c4_focus"))))
				stop("Cached trait/layer/module mismatch: ", f)
			settings_match <- identical(meta$module_policy,le8_module_policy(m)) && identical(meta$analysis_options, le8_analysis_options(y)) && identical(as.integer(meta$seed), as.integer(Sys.getenv("SEED", "2026")))
			generated <- meta$generated
			version <- if (is.null(meta$code_version)) "saved fitted result" else meta$code_version
			present <- function(paths) all(vapply(paths, function(path) {
				if (file.exists(path)) return(isTRUE(file.info(path)$size > 0))
				grepl('[.]csv$', path) && !is.null(table_entry(path))
			}, logical(1)))
			files <- file.path(d, exported[[m]])
			manifest <- file.path(d, "figure_manifest.csv")
			figures_ok <- FALSE
			if (present(manifest)) {
				mf <- read_export(manifest)
				if ('group' %in% names(mf) && m %in% c('c4_connect', 'c4_panel_validation'))
					mf <- mf[(mf$group %in% c('equal_budget','deployed_concept_fidelity')) == (m == 'c4_panel_validation'), , drop = FALSE]
				figures_ok <- "file" %in% names(mf) && nrow(mf) > 0 && present(file.path(d, mf$file))
			}
			reusable <- all(required[[m]] %in% names(x)) && present(files[grepl("[.]csv$", files)])
			complete <- reusable && present(files) && figures_ok
			if (m == "c4_connect") {
				shared <- d
				shared_ready <- present(file.path(shared, c(
					"c4.interactions.res.rds", "c4.nonlin.res.rds", "c4.penalty.res.rds",
					"c4.LE8_pairwise_interactions.csv", "c4.nonlin_tests.csv", "c4.penalty_summary.csv"
				)))
				reusable <- reusable && shared_ready
				complete <- complete && shared_ready
				shared_manifest <- file.path(shared, "figure_manifest.csv")
				if (present(shared_manifest)) {
					views <- read_export(shared_manifest)
					if ('group' %in% names(views)) views <- views[!views$group %in% c('equal_budget','deployed_concept_fidelity'), , drop = FALSE]
					complete <- complete && "file" %in% names(views) && nrow(views) > 0 && present(file.path(shared, views$file))
				} else complete <- FALSE
			}
			if (m == "c4_panel_validation") {
				reusable <- reusable && all(c("metrics", "panel_members", "design") %in% names(x$tables))
				complete <- complete && reusable
				budgets <- sort(unique(as.numeric(strsplit(Sys.getenv("C4_FOCUS_BUDGETS", "5,10,50"), ",")[[1]])))
				actual <- sort(unique(x$tables$metrics$budget[x$tables$metrics$budget > 0]))
				if (!identical(as.numeric(actual), budgets)) settings_match <- FALSE
			}
			if (m == "c2_cause") for (n in c("RUN_MRlink2", "RUN_Dandelion"))
				if (!identical(meta[[n]], Sys.getenv(n, "Top"))) settings_match <- FALSE
			if (!settings_match || !identical(meta$code_signature,le8_code_fingerprint())) { reusable <- FALSE; complete <- FALSE }
			status <- if (complete) "completed" else if (reusable) "cached" else "incomplete"
			detail <- if (complete) "Saved estimates and original scope retained; no refitting or relabelling" else if (reusable)
				"Numerical results complete; regenerate presentation without refitting" else "Missing result fields or numerical exports; stage resume required"
		}
		rows[[length(rows) + 1L]] <- data.frame(trait = y, layer = b, module = m, status, generated, version, detail, stringsAsFactors = FALSE)
	}
	if (length(rows)) do.call(rbind, rows) else data.frame()
}

if (sys.nframe() == 0L && length(commandArgs(TRUE)) && commandArgs(TRUE)[1] == "--check-completed") {
	args <- commandArgs(TRUE)
	if (length(args) != 6L) stop("--check-completed requires root, traits, layers, modules and audit path")
	split <- function(z) strsplit(z, ",", fixed = TRUE)[[1]]
	z <- le8_completed_results(args[2], split(args[3]), split(args[4]), split(args[5]))
	if (nrow(z)) {
		write.csv(z, args[6], row.names = FALSE)
		# Output completeness alone cannot prove unchanged inputs, groups or LD.
		# Each module revalidates its full stage signature before reusing its fits.
		# Explicit presentation-only entry points retain their frozen-results path.
	}
	quit(save = "no", status = 0)
}

# Shared prediction validation.
# Shared prediction validation for C4 panel validation and Final joint models.
# Explicit assay budgets, training-only feature selection and held-out evaluation.
# No GDF15 or NT-proBNP penalty is hard coded. Ablations are separately refitted.
le8_prepare_prediction_matrix <- function(train, test, vars) {
	train <- as.data.frame(train) ; test <- as.data.frame(test)
	vars <- intersect(vars, intersect(names(train), names(test))) ; a <- b <- list() ; audit <- list()
	for (v in vars) {
		x <- train[[v]] ; y <- test[[v]]
		if (is.numeric(x) || is.integer(x) || is.logical(x)) {
			x <- as.numeric(x) ; y <- as.numeric(y) ; m <- median(x[is.finite(x)], na.rm = TRUE)
			if (!is.finite(m)) next
			x[!is.finite(x)] <- m ; y[!is.finite(y)] <- m ; s <- sd(x)
			if (!is.finite(s) || s <= 0) next
			a[[v]] <- (x - m) / s ; b[[v]] <- (y - m) / s
			audit[[length(audit) + 1L]] <- tibble(variable = v, type = "numeric", center = m, scale = s, reference = "training median", unseen_test_levels = 0L)
		} else {
			x <- as.character(x) ; y <- as.character(y) ; x[x == ""] <- NA ; y[y == ""] <- NA
			tab <- table(x) ; if (!length(tab)) next
			mode <- names(tab)[which.max(tab)] ; lev <- sort(names(tab)) ; x[is.na(x)] <- mode
			unseen <- !is.na(y) & !y %in% lev ; y[is.na(y) | unseen] <- mode
			if (length(lev) < 2L) next
			# Reference and levels are learned exclusively from the training data.
			for (l in setdiff(lev, mode)) {
				n <- paste0(v, "__", l) ; a[[n]] <- as.numeric(x == l) ; b[[n]] <- as.numeric(y == l)
			}
			audit[[length(audit) + 1L]] <- tibble(variable = v, type = "categorical", center = NA_real_, scale = NA_real_, reference = mode, unseen_test_levels = sum(unseen))
		}
	}
	if (!length(a)) return(list(train = matrix(nrow = nrow(train), ncol = 0), test = matrix(nrow = nrow(test), ncol = 0), audit = bind_rows(audit)))
	xx <- as.matrix(as.data.frame(a, check.names = FALSE)) ; yy <- as.matrix(as.data.frame(b, check.names = FALSE))
	# Remove exactly redundant columns based on training predictors, not outcomes.
	qr0 <- qr(cbind(Intercept = 1, xx)) ; keep <- sort(setdiff(qr0$pivot[seq_len(qr0$rank)], 1L) - 1L)
	list(train = xx[, keep, drop = FALSE], test = yy[, keep, drop = FALSE], audit = bind_rows(audit))
}
le8_breslow_hazard <- function(time, event, lp) {
	# lp must use the same training centering as subsequent predictions.
	if (any(!is.finite(lp)) || max(abs(lp)) > 700) stop("Unstable linear predictor")
	z <- data.frame(time = time, event = event, risk = exp(lp))
	a <- stats::aggregate(z[c("event", "risk")], list(time = z$time), sum)
	a <- a[order(a$time), , drop = FALSE]
	denom <- rev(cumsum(rev(a$risk)))
	data.frame(time = a$time, hazard = cumsum(a$event / denom))
}

# 🚩 Shared participant groups and proxy statistics
le8_validate_ids <- function(d, label = "input") {
	if (!"eid" %in% names(d)) stop(label, " lacks eid")
	d$eid <- as.character(d$eid)
	if (anyNA(d$eid) || any(!nzchar(trimws(d$eid))) || anyDuplicated(d$eid)) stop(label, " contains missing/duplicate eid")
	d
}
le8_participant_groups <- function(d) {
	d <- le8_validate_ids(d)
	column <- Sys.getenv("LE8_GROUP_COLUMN", Sys.getenv("PGS_GROUP_COLUMN", ""))
	file <- Sys.getenv("LE8_GROUP_FILE", "")
	if (nzchar(file)) {
		m <- le8_validate_ids(as.data.frame(data.table::fread(file, colClasses="character")), "group mapping")
		if (!"group" %in% names(m)) stop("LE8_GROUP_FILE requires eid,group (kinship connected components)")
		if (anyNA(m$group) || any(!nzchar(trimws(m$group)))) stop("Incomplete group file")
		g <- as.character(m$group[match(d$eid, m$eid)])
		if (nzchar(column) && (!column %in% names(d) || anyNA(d[[column]]) || any(as.character(d[[column]]) != g, na.rm=TRUE))) stop("Conflicting family definitions in file and phenotype column")
	} else if (nzchar(column)) {
		if (!column %in% names(d)) stop("Required group column missing: ", column)
		g <- as.character(d[[column]])
	} else if (".group" %in% names(d)) g <- as.character(d$.group) else g <- d$eid
	if (anyNA(g) || any(!nzchar(trimws(g)))) stop("Incomplete participant group mapping")
	g
}
le8_outer_roles <- function(d, file = Sys.getenv("LE8_OUTER_ROSTER", "")) {
	if (!nzchar(file)) return(NULL)
	d <- le8_validate_ids(d)
	m <- le8_validate_ids(as.data.frame(data.table::fread(file, colClasses="character")), "outer roster")
	if (!"role" %in% names(m) || anyNA(m$role) || !all(m$role %in% c("training","test"))) stop("Outer roster requires eid,role with training/test")
	role <- m$role[match(d$eid,m$eid)]
	if (anyNA(role) || !setequal(role,c("training","test"))) stop("Outer roster must cover participants and both roles")
	g <- le8_participant_groups(d)
	if (any(vapply(split(role,g),function(x) length(unique(x))>1L,logical(1)))) stop("Families cross shared outer roles")
	ifelse(role=="test","validation","training")
}
le8_group_folds <- function(group, k = 5L, seed = 2026L) {
	group <- as.character(group); u <- sort(unique(group))
	if (anyNA(group) || any(!nzchar(group)) || length(u) < k) stop("Insufficient valid independent groups")
	set.seed(seed); f <- sample(rep(seq_len(k), length.out = length(u))); f[match(group, u)]
}
le8_group_bootstrap <- function(group) {
	blocks <- split(seq_along(group), as.character(group))
	unlist(blocks[sample(seq_along(blocks), length(blocks), replace = TRUE)], use.names = FALSE)
}
le8_proxy_map <- function(d, features, components, covars, adjustment = "basic_adjusted", min_n = 100L) {
	if (length(setdiff(covars, names(d)))) stop("Missing required proxy covariates: ", paste(setdiff(covars, names(d)), collapse = ","))
	components <- intersect(components, names(d)); features <- intersect(features, names(d))
	if (!adjustment %in% c("basic_adjusted", "conditional_specificity")) stop("Unknown proxy adjustment")
	rows <- lapply(components, function(cmp) {
		cv <- unique(c(covars, if (adjustment == "conditional_specificity") setdiff(components, cmp)))
		dd <- d[complete.cases(d[, unique(c(cmp, cv)), drop = FALSE]), , drop = FALSE]
		if (nrow(dd) < min_n || !is.finite(sd(dd[[cmp]])) || sd(dd[[cmp]]) <= 0) return(tibble())
		cv <- cv[vapply(dd[cv], function(x) length(unique(x)) > 1L, logical(1))]
		q <- qr(model.matrix(reformulate(if (length(cv)) cv else "1"), dd))
		y <- qr.resid(q, as.numeric(scale(dd[[cmp]]))); sy <- sqrt(sum(y^2)); df <- nrow(dd) - q$rank - 1L
		if (sy <= 0 || df < 10) return(tibble())
		bind_rows(lapply(split(features, ceiling(seq_along(features)/64)), function(bb) {
			x <- as.matrix(dd[, bb, drop = FALSE]); storage.mode(x) <- "double"
			for (j in seq_len(ncol(x))) {
				v <- x[,j]; med <- median(v[is.finite(v)], na.rm = TRUE)
				v[!is.finite(v)] <- if (is.finite(med)) med else 0; ss <- sd(v)
				x[,j] <- if (is.finite(ss) && ss > 0) (v-mean(v))/ss else 0
			}
			xr <- qr.resid(q,x); den <- sqrt(colSums(xr^2))*sy
			r <- as.numeric(crossprod(xr,y))/den; r <- pmax(-.999999,pmin(.999999,r))
			z <- r*sqrt(df/(1-r*r)); p <- 2*pt(abs(z),df,lower.tail=FALSE)
			tibble(feature=bb, component=cmp, r, z, p, N=nrow(dd), adjustment,
				se_r=sqrt((1-r*r)/df), status=ifelse(is.finite(p),"ok","constant/unavailable"))
		}))
	})
	z <- bind_rows(rows)
	if (!nrow(z)) return(z)
	z |> group_by(component) |> mutate(FDR=p.adjust(p,"BH",n=length(features))) |> ungroup() |>
		mutate(FDR_all=p.adjust(p,"BH",n=length(features)*length(components)))
}

le8_fit_budget_model <- function(train, test, clinical, features, tvar, evar,
			solver = "cox", seed = 2026) {
	if (!solver %in% c("cox", "ridge")) stop("solver must be cox or ridge")
	vars <- unique(c(clinical, features)) ; xx <- le8_prepare_prediction_matrix(train, test, vars)
	if (!ncol(xx$train)) return(list(status = "no usable predictors"))
	orig <- colnames(xx$train) ; safe <- paste0("x", seq_len(ncol(xx$train)))
	tr <- as.data.frame(xx$train) ; te <- as.data.frame(xx$test) ; names(tr) <- names(te) <- safe
	tr$.time <- train[[tvar]] ; tr$.event <- train[[evar]]
	if (any(!is.finite(tr$.time) | tr$.time <= 0 | !tr$.event %in% c(0, 1)))
		return(list(status = "invalid incident training outcomes"))
	condition <- kappa(scale(xx$train, center = TRUE, scale = FALSE), exact = TRUE)
	warn <- character() ; capture <- function(w) {
		warn <<- c(warn, conditionMessage(w)) ; invokeRestart("muffleWarning")
	}
	penalty <- as.numeric(vapply(orig, function(v) any(v == features | startsWith(v, paste0(features, "__"))), logical(1)))
	fit_method <- if (solver == "ridge" && any(penalty > 0) && ncol(xx$train) > 1L) "ridge Cox, inner training CV" else "unpenalized Cox"
	lambda <- NA_real_ ; lp_center <- 0
	if (fit_method == "unpenalized Cox") {
		if (!is.finite(condition) || condition > 1e6)
			return(list(status = "ill-conditioned training design; use prespecified ridge solver", condition_number = condition))
		f <- reformulate(safe, response = "survival::Surv(.time,.event)")
		fit <- tryCatch(withCallingHandlers(survival::coxph(f, tr,
			ties = "breslow", x = TRUE, y = TRUE,
			model = TRUE, singular.ok = FALSE
		), warning = capture), error = function(e) e)
		if (inherits(fit, "condition")) return(list(status = conditionMessage(fit)))
		beta <- as.numeric(coef(fit)) ; lp <- drop(as.matrix(te) %*% beta)
		bh <- survival::basehaz(fit, centered = FALSE)
	} else {
		if (!requireNamespace("glmnet", quietly = TRUE)) stop("ridge requires glmnet")
		if (min(table(factor(tr$.event, levels = 0 : 1))) < 5) return(list(status = "insufficient events for inner CV"))
		set.seed(seed) ; foldid <- integer(nrow(tr))
		for (e in 0 : 1) {
			ix <- which(tr$.event == e) ; foldid[ix] <- sample(rep(1 : 5, length.out = length(ix)))
		}
		if (".group" %in% names(train)) foldid <- le8_group_folds(train$.group, 5, seed)
		fit <- tryCatch(withCallingHandlers(glmnet::cv.glmnet(xx$train,
			survival::Surv(tr$.time, tr$.event),
			family = "cox", alpha = 0, foldid = foldid, cox.ties = "breslow",
			type.measure = "deviance", standardize = TRUE, penalty.factor = penalty
		), warning = capture), error = function(e) e)
		if (inherits(fit, "condition")) return(list(status = conditionMessage(fit)))
		lambda_rule <- Sys.getenv("LE8_RIDGE_LAMBDA", "lambda.min")
		if (!lambda_rule %in% c("lambda.min", "lambda.1se")) stop("Invalid LE8_RIDGE_LAMBDA")
		lambda <- fit[[lambda_rule]] ; beta <- as.numeric(coef(fit, s = lambda_rule))
		lp_train <- drop(xx$train %*% beta) ; lp_center <- mean(lp_train)
		lp <- drop(xx$test %*% beta) - lp_center
		bh <- tryCatch(le8_breslow_hazard(tr$.time, tr$.event, lp_train - lp_center), error = function(e) e)
		if (inherits(bh, "condition")) return(list(status = conditionMessage(bh)))
	}
	if (any(!is.finite(beta)) || any(!is.finite(lp)) || any(!is.finite(bh$hazard)))
		return(list(status = "non-finite fit/prediction"))
	fatal_warning <- any(grepl("converg|infinite|overflow|numerical", warn, ignore.case = TRUE))
	if (fatal_warning) return(list(
		status = paste("unstable fit:", paste(unique(warn), collapse = "; ")),
		condition_number = condition
	))
	count <- sum(vapply(features, function(v) any(orig == v | startsWith(orig, paste0(v, "__"))), logical(1)))
	# Downstream scoring uses coefficients, preprocessing and baseline hazard.
	# A coxph fit's formula environment captures this entire call (including
	# full train/test omics); serializing it can turn one small model into GBs.
	list(
		status = "ok", lp = lp, contributions=sweep(xx$test,2,beta,`*`), baseline_hazard = bh, N_selected = count,
		coefficient = tibble(variable = orig, beta = beta), preprocess = xx$audit,
		N_effective = sum(vapply(features, function(v)
			any(abs(beta[orig == v | startsWith(orig, paste0(v, "__"))]) > 1e-8), logical(1))),
		effective_threshold = 1e-8,
		effective_interpretation = "Numerical coefficient-use diagnostic; not variable importance or assay utility",
		fit_method = fit_method, tie_method = "breslow", condition_number = condition, lambda = lambda, lp_center = lp_center,
		N_train = nrow(train), events_train = sum(tr$.event), design_rank = qr(xx$train)$rank,
		penalty = tibble(variable = orig, penalty_factor = penalty),
		cv_curve = if (inherits(fit, "cv.glmnet")) tibble(lambda = fit$lambda, loss = fit$cvm, se = fit$cvsd) else tibble(),
		molecular_lp_sd_train = sd(drop(xx$train[, penalty > 0, drop = FALSE] %*% beta[penalty > 0])),
		molecular_lp_sd_test = sd(drop(xx$test[, penalty > 0, drop = FALSE] %*% beta[penalty > 0])),
		warnings = paste(unique(warn), collapse = "; ")
	)
}
le8_risk_at <- function(obj, horizon) {
	bh <- obj$baseline_hazard ; idx <- which(bh$time <= horizon)
	H <- if (length(idx)) bh$hazard[max(idx)] else 0
	pmin(1 - 1e-8, pmax(1e-8, - expm1( - H * exp(pmin(30, pmax( - 30, obj$lp))))))
}
le8_ipcw <- function(time, event, horizon) {
	sf <- survival::survfit(survival::Surv(time, 1 - event) ~ 1)
	at <- function(t, left = FALSE) {
		i <- findInterval(t, sf$time)
		if (left) {
			hit <- i > 0L & sf$time[pmax(i, 1L)] == t ; i[hit] <- i[hit] - 1L
		}
		z <- rep(1, length(t)) ; z[i > 0] <- sf$surv[i[i > 0]] ; z
	}
	y <- as.numeric(event == 1 & time <= horizon) ; w <- numeric(length(time))
	case <- y == 1 ; ctrl <- time > horizon
	gcase <- at(time[case], TRUE) ; gt <- at(horizon)
	valid <- length(gt) == 1L && is.finite(gt) && gt > .01 && sum(case) >= 20 && sum(ctrl) >= 20
	if (valid) {
		w[case] <- 1 / pmax(gcase, .01) ; w[ctrl] <- 1 / gt
	}
	list(
		y = y, w = w, status = if (valid) "ok" else "insufficient supported follow-up/cases/controls",
		N_case = sum(case), N_control = sum(ctrl), G_horizon = gt,
		assumption = "marginal independent censoring; death is censored, so predicted risk is not a competing-risk CIF"
	)
}
le8_weighted_auc <- function(p, y, w) {
	ii <- is.finite(p) & is.finite(y) & is.finite(w) & w > 0
	p <- p[ii] ; y <- y[ii] ; w <- w[ii] ; o <- order(p) ; p <- p[o] ; y <- y[o] ; w <- w[o]
	if (!length(p) || sum(w * y) <= 0 || sum(w * (1 - y)) <= 0) return(NA_real_)
	g <- cumsum(c(TRUE, diff(p) != 0)) ; a <- as.numeric(rowsum(w * y, g, reorder = FALSE))
	b <- as.numeric(rowsum(w * (1 - y), g, reorder = FALSE)) ; cum <- cumsum(b) - b
	sum(a * (cum + .5 * b)) / (sum(a) * sum(b))
}
le8_evaluate_risk <- function(time, event, p, horizon, model, budget, paradigm, ablation = "none", B = 0, seed = 2026, groups = seq_along(time)) {
	iw <- le8_ipcw(time, event, horizon)
	base <- tibble(model, budget, paradigm, ablation, horizon,
		N = length(time), events_by_horizon = iw$N_case,
		controls_at_horizon = iw$N_control, status = iw$status, AUC = NA_real_, AUC_lo = NA_real_, AUC_hi = NA_real_,
		Brier = NA_real_, Brier_lo = NA_real_, Brier_hi = NA_real_, calibration_intercept = NA_real_, calibration_slope = NA_real_,
		censoring_survival = iw$G_horizon, estimand = iw$assumption
	)
	if (iw$status != "ok") return(list(metrics = base, calibration = tibble(), decision = tibble(), bootstrap = tibble()))
	y <- iw$y ; w <- iw$w ; n <- length(y) ; base$AUC <- le8_weighted_auc(p, y, w) ; base$Brier <- sum(w * (y - p) ^ 2) / n
	lp <- qlogis(pmin(1 - 1e-8, pmax(1e-8, p))) ; ok <- w > 0
	fit <- tryCatch(glm(y[ok] ~ lp[ok], family = quasibinomial(), weights = w[ok]), error = function(e) NULL)
	ic <- tryCatch(glm(y[ok] ~ 1 + offset(lp[ok]), family = quasibinomial(), weights = w[ok]), error = function(e) NULL)
	if (!is.null(fit)) base$calibration_slope <- unname(coef(fit)[2])
	if (!is.null(ic)) base$calibration_intercept <- unname(coef(ic)[1])
	group <- pmin(10L, ceiling(rank(p, ties.method = "first") / n * 10L))
	cal <- tibble(group, p, y, w) |>
		group_by(group) |>
		summarise(
			N = n(), predicted = mean(p),
			observed_ipcw = sum(w * y) / sum(w), .groups = "drop"
		) |>
		mutate(model, horizon, budget, paradigm, ablation)
	dec <- map_dfr(seq(.01, .20, by = .01), function(th) {
		tr <- p >= th
		tibble(model, horizon,
			threshold = th, net_benefit = (sum(w * y * tr) - sum(w * (1 - y) * tr) * th / (1 - th)) / n,
			treat_all = (sum(w * y) - sum(w * (1 - y)) * th / (1 - th)) / n, treat_none = 0
		)
	})
	boot <- tibble()
	if (B > 0L) {
		set.seed(seed)
		boot <- map_dfr(seq_len(B), function(i) {
			ix <- le8_group_bootstrap(groups)
			z <- le8_ipcw(time[ix], event[ix], horizon)
			tibble(
				replicate = i, AUC = if (z$status == "ok") le8_weighted_auc(p[ix], z$y, z$w) else NA_real_,
				Brier = if (z$status == "ok") mean(z$w * (z$y - p[ix]) ^ 2) else NA_real_
			)
		}) |> mutate(model, horizon)
		for (metric in c("AUC", "Brier")) {
			v <- boot[[metric]] ; v <- v[is.finite(v)]
			if (length(v) >= max(20, ceiling(B * .8))) {
				base[[paste0(metric, "_lo")]] <- quantile(v, .025, names = FALSE) ; base[[paste0(metric, "_hi")]] <- quantile(v, .975, names = FALSE)
			}
		}
	}
	list(metrics = base, calibration = cal, decision = dec, bootstrap = boot)
}
le8_final_additions <- function(dat, biom_vars, ranked, training_screen, distalset, geneticset,
			clinical, tvar, evar, split, layer, outdir) {
	values <- suppressWarnings(as.numeric(le8_csv_env("FINAL_ASSAY_BUDGETS", paste(LE8_ASSAY_BUDGETS, collapse = ","))))
	if (!length(values) || any(!is.finite(values)) || any(values != floor(values)) || any(values < 1 | values > 500)) stop("FINAL_ASSAY_BUDGETS must contain integers 1..500")
	budgets <- sort(unique(as.integer(values)))
	horizons <- as.numeric(le8_csv_env("FINAL_REVIEW_HORIZONS", "2,5,10"))
	if (any(!is.finite(horizons) | horizons <= 0)) stop("Invalid FINAL_REVIEW_HORIZONS")
	B <- as.integer(le8_num_env("FINAL_REVIEW_BOOT", 100)) ; if (B < 0L) stop("FINAL_REVIEW_BOOT must be nonnegative")
	cv <- le8_csv_env("FINAL_REVIEW_CLINICAL_VARS", paste(clinical, collapse = ","))
	missing <- setdiff(cv, names(dat)) ; if (length(missing)) stop("Specified clinical predictors are missing: ", paste(missing, collapse = ","))
	train <- dat[split == "training", , drop = FALSE] ; test <- dat[split == "validation", , drop = FALSE]
	if (nrow(train) < 500 || nrow(test) < 100 || sum(train[[evar]]) < 40) return(list(status = tibble(status = "insufficient training/validation data")))
	# Near-term candidate ranking uses only training participants; early censoring
	# is retained by a time-truncated Cox model, not converted into non-cases.
	near <- train ; near$.near_time <- pmin(near[[tvar]], 2) ; near$.near_event <- as.integer(near[[evar]] == 1 & near[[tvar]] <= 2)
	near_scan <- cox_scan(near, head(ranked, as.integer(le8_num_env("FINAL_NEAR_SCREEN_MAX", 100))), cv, Y,
		time_var = ".near_time", event_var = ".near_event"
	)
	near_rank <- near_scan |>
		filter(is.finite(p.value)) |>
		arrange(p.value) |>
		pull(term)
	paradigms <- list(`Overall-selected` = ranked, `Distal-selected` = ranked[ranked %in% distalset], `Near-term-selected` = near_rank)
	if (truthy(Sys.getenv("FINAL_GENETIC_EVIDENCE_INDEPENDENT", unset = "FALSE")))
		paradigms[["Genetic-region-selected"]] <- ranked[ranked %in% geneticset]
	designs <- list(list(name = "Clinical", features = character(), budget = 0L, paradigm = "Clinical", ablation = "none"))
	for (pa in names(paradigms)) for (k in budgets) if (length(paradigms[[pa]]) >= k)
		designs[[length(designs) + 1L]] <- list(name = paste(pa, k, sep = "_"), features = head(paradigms[[pa]], k), budget = k, paradigm = pa, ablation = "none")
	# Fixed-panel ablation: do not fill the removed slot with a substitute assay.
	# This answers incremental contribution of the removed marker(s).
	full <- head(ranked, max(budgets))
	if (layer == "protein") for (nm in c("GDF15", "NTPROBNP", "GDF15+NTPROBNP")) {
		drop <- strsplit(nm, "+", fixed = TRUE)[[1]]
		designs[[length(designs) + 1L]] <- list(
			name = paste0("Overall_minus_", nm), features = setdiff(full, drop),
			budget = length(full), paradigm = "Overall-selected", ablation = nm
		)
	}
	rd <- le8_job_dir(outdir, "final_prediction") ; dir.create(file.path(rd, "review_models"), showWarnings = FALSE)
	metrics <- cal <- dca <- boot <- members <- prep <- coefs <- list() ; pred_cache <- list()
	for (i in seq_along(designs)) {
		ds <- designs[[i]] ; obj <- le8_fit_budget_model(train, test, cv, ds$features, tvar, evar)
		members[[i]] <- tibble(
			model = ds$name, feature = if (length(ds$features)) ds$features else NA_character_,
			requested_budget = ds$budget, actual_panel_size = length(ds$features), ablation = ds$ablation, status = obj$status
		)
		if (obj$status != "ok") next
		prep[[i]] <- obj$preprocess |> mutate(model = ds$name)
		coefs[[i]] <- obj$coefficient |> mutate(model = ds$name)
		saveRDS(obj, file.path(rd, "review_models", paste0(gsub("[^A-Za-z0-9_.-]", "_", ds$name), ".rds")))
		ev <- lapply(horizons, function(h) {
			p <- le8_risk_at(obj, h)
			# Identical seed at a given horizon gives paired validation resamples.
			z <- le8_evaluate_risk(test[[tvar]], test[[evar]], p, h, ds$name, ds$budget, ds$paradigm, ds$ablation,
				B = B, seed = SEED + as.integer(h * 100)
			)
			z$metrics$actual_n_assays <- length(ds$features) ; z$metrics$fitted_n_assays <- obj$N_selected
			z
		})
		metrics[[i]] <- bind_rows(lapply(ev, `[[`, "metrics")) ; cal[[i]] <- bind_rows(lapply(ev, `[[`, "calibration"))
		dca[[i]] <- bind_rows(lapply(ev, `[[`, "decision")) ; boot[[i]] <- bind_rows(lapply(ev, `[[`, "bootstrap"))
	}
	met <- bind_rows(metrics) ; ca <- bind_rows(cal) ; dc <- bind_rows(dca) ; bo <- bind_rows(boot)
	delta <- tibble()
	if (nrow(bo)) {
		ref <- bo |>
			filter(model == "Clinical") |>
			select(replicate, horizon, AUC_ref = AUC, Brier_ref = Brier)
		delta <- bo |>
			filter(model != "Clinical") |>
			inner_join(ref, by = c("replicate", "horizon")) |>
			mutate(delta_AUC = AUC - AUC_ref, delta_Brier = Brier - Brier_ref) |>
			group_by(model, horizon) |>
			summarise(
				valid = sum(is.finite(delta_AUC)), delta_AUC_lo = quantile(delta_AUC, .025, na.rm = TRUE),
				delta_AUC_hi = quantile(delta_AUC, .975, na.rm = TRUE), delta_Brier_lo = quantile(delta_Brier, .025, na.rm = TRUE),
				delta_Brier_hi = quantile(delta_Brier, .975, na.rm = TRUE), .groups = "drop"
			)
	}
	tables <- list(
		budget_metrics = met, panel_members = bind_rows(members), calibration = ca, decision_curves = dc,
		paired_delta_CI = delta, model_coefficients = bind_rows(coefs), training_preprocessing = bind_rows(prep), near_training_scan = near_scan,
		design = tibble(
			item = c("Clinical comparator", "Genetic evidence used in prediction", "Uncertainty", "Prediction risk"),
			value = c(
				paste(cv, collapse = ";"), Sys.getenv("FINAL_GENETIC_EVIDENCE_INDEPENDENT", unset = "FALSE"),
				"Validation bootstrap of frozen fits; not uncertainty of the entire training-selection procedure",
				"Cause-specific/net risk with deaths censored; not a competing-risk cumulative incidence"
			)
		)
	)
	for (nm in names(tables)) write_raw_csv(tables[[nm]], paste0("review_", nm, ".csv"), rd)
	# The row split is saved only in the local raw directory, not the workbook.
	data.table::fwrite(tibble(eid = dat$eid, split = split), file.path(rd, "review_split.csv.gz"), compress = "gzip")
	le8_plot_final_review(met, ca, dc, horizons, budgets, outdir)
	tables
}


le8_plot_final_review <- function(met, ca, dc, horizons, budgets, outdir) {
	if (nrow(met)) {
		h <- max(horizons) ; mm <- met |> filter(horizon == h, status == "ok")
		a <- mm |>
			filter(ablation == "none") |>
			ggplot(aes(actual_n_assays, AUC, color = paradigm)) +
			geom_line() +
			geom_point() +
			geom_errorbar(aes(ymin = AUC_lo, ymax = AUC_hi), width = .2) +
			labs(title = paste0("a. Assay budget and ", h, "-year discrimination"), x = "Actual measured assays", y = "IPCW AUC", color = NULL) +
			theme_5c(9)
		b <- ca |>
			filter(horizon == h, model %in% c("Clinical", paste0("Overall-selected_", max(budgets)))) |>
			ggplot(aes(predicted, observed_ipcw, color = model)) +
			geom_abline(slope = 1, intercept = 0, linetype = 2) +
			geom_line() +
			geom_point() +
			labs(title = "b. Frozen-model calibration", x = "Predicted risk", y = "Observed IPCW risk", color = NULL) +
			theme_5c(9)
		c <- mm |>
			filter(ablation != "none" | model == paste0("Overall-selected_", max(budgets))) |>
			ggplot(aes(AUC, model)) +
			geom_errorbarh(aes(xmin = AUC_lo, xmax = AUC_hi), height = .15) +
			geom_point() +
			labs(title = "c. Fixed-panel ablations, refitted in training", x = "IPCW AUC", y = NULL) +
			theme_5c(9)
		d <- dc |>
			filter(horizon == h, model %in% c("Clinical", paste0("Overall-selected_", max(budgets)))) |>
			ggplot(aes(threshold, net_benefit, color = model)) +
			geom_line() +
			geom_hline(yintercept = 0, linetype = 2) +
			labs(title = "d. Exploratory decision curves", x = "Risk threshold", y = "Net benefit", color = NULL) +
			theme_5c(9)
		save_plot((a | b) / (c | d) + plot_annotation(caption = "Same held-out participants. No validation outcomes are used to select features, tune budgets or recalibrate predictions. Genetic-region selection requires independent evidence declaration."),
			"Fig13.budget_ablation_calibration.png", 18, 12,
			outdir = outdir
		)
	}
}

# Relearn connection-based feature membership inside the outer training split.
# No validation LE8 measurements or outcomes enter this procedure.
le8_training_connection_set <- function(train, features, components, covars, rawdir) {
	d <- as.data.frame(train); components <- intersect(components,names(d))
	if (!length(components)) return(character())
	d$.group <- le8_participant_groups(d)
	half <- if (".le8_proxy_half" %in% names(d)) d$.le8_proxy_half else le8_group_folds(d$.group,2,SEED+411)
	if (anyNA(half) || !all(half %in% 1:2)) stop("Invalid proxy halves")
	if (any(vapply(split(half,d$.group),function(x) length(unique(x))>1,logical(1)))) stop("Family crosses proxy halves")
	scans <- bind_rows(lapply(c("basic_adjusted","conditional_specificity"),function(adj)
		bind_rows(lapply(1:2,function(h) le8_proxy_map(d[half==h,,drop=FALSE],features,components,covars,adj) |> mutate(half=h)))))
	if (!nrow(scans)) return(character())
	write_raw_csv(scans,"connection_proxy_maps_training_only.csv",rawdir)
	a <- scans |> filter(half==1,adjustment=="basic_adjusted") |> group_by(feature) |>
		mutate(specificity=abs(r)/sum(abs(r),na.rm=TRUE)) |> ungroup() |>
		select(feature,component,r1=r,FDR1=FDR,specificity)
	b <- scans |> filter(half==2,adjustment=="basic_adjusted") |> select(feature,component,r2=r,FDR2=FDR)
	joined <- left_join(a,b,by=c("feature","component")) |> mutate(
		selected=is.finite(FDR1)&is.finite(FDR2)&FDR1<.05&FDR2<.05&sign(r1)==sign(r2),
		proxy_strength=pmin(abs(r1),abs(r2)), scope="development only; disease-blind replicated basic-adjusted proxy",
		training_hash=le8_hash_object(sort(d$eid))) |> arrange(desc(proxy_strength),feature,component)
	# Greedy redundancy control is learned from development X only. Keep the full
	# qualified pool, but postpone near-duplicate assays when forming the prefix.
	out <- unique(joined$feature[joined$selected %in% TRUE])
	cutoff <- as.numeric(Sys.getenv("C4_PROXY_REDUNDANCY_R", ".95"))
	if (!is.finite(cutoff) || cutoff<=0 || cutoff>1) stop("Invalid C4_PROXY_REDUNDANCY_R")
	candidates <- head(out,500L); chosen <- deferred <- character()
	if (length(candidates)>1L && cutoff<1) {
		ids <- order(vapply(d$eid,function(id) digest::digest(paste(SEED,id),algo="xxhash64"),character(1)))
		ids <- head(ids,5000L)
		x <- as.matrix(d[ids,candidates,drop=FALSE]); storage.mode(x) <- "double"
		for(j in seq_len(ncol(x))) { v<-x[,j]; m<-median(v[is.finite(v)],na.rm=TRUE);v[!is.finite(v)]<-if(is.finite(m)) m else 0;x[,j]<-v }
		r <- suppressWarnings(cor(x));r[!is.finite(r)]<-0
		for(f in candidates) {
			if (length(chosen) && any(abs(r[f,chosen])>cutoff)) deferred<-c(deferred,f) else chosen<-c(chosen,f)
		}
		out <- c(chosen,setdiff(out,candidates),deferred)
	}
	joined$assay_rank <- match(joined$feature,out)
	joined$redundancy_rule <- paste0("Greedy |r|<=",cutoff," in top 500 proxy candidates; deterministic <=5000 development donors; qualified correlated assays postponed")
	write_raw_csv(joined,"connection_membership_training_only.csv",rawdir)
	attr(out,"membership") <- joined
	out
}

# Shared functions for the LE8-supervised 5C omics pipeline.
# The functions in this file deliberately keep protein and metabolite branches separate.

.required_pkgs <- c(
	"data.table", "dplyr", "tidyr", "purrr", "stringr", "tibble",
	"ggplot2", "ggrepel", "patchwork", "survival", "scales", "openxlsx", "forcats"
)
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(
	pacman::p_load(char = .required_pkgs)
)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
bt <- function(x) paste0("`", x, "`")

# Reproducible outcome-blind family folds shared by Final and C4 nonlinear analyses. Prediction-model
# functions themselves live in 0f/prediction.R and Final extensions in final.R.
make_folds <- function(dat, event_var, k = 5, seed = 2026) {
	group <- if (".group" %in% names(dat)) dat$.group else if ("eid" %in% names(dat)) le8_participant_groups(dat) else seq_len(nrow(dat))
	le8_group_folds(group,k,seed)
}

cap <- function(x, limit) pmax(pmin(x, limit), - limit)
safe_log <- function(x) ifelse(is.finite(x) & x > 0, log(x), NA_real_)
truthy <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
std_num <- function(x) {
	x <- suppressWarnings(as.numeric(x)) ; s <- stats::sd(x, na.rm = TRUE)
	if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
	as.numeric(scale(x))
}
inormal2 <- function(x) {
	x <- suppressWarnings(as.numeric(x)) ; n <- sum(is.finite(x))
	if (n < 5 || !is.finite(sd(x, na.rm = TRUE)) || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
	qnorm((rank(x, na.last = "keep", ties.method = "average") - 0.5) / n)
}


# 🚩 Paths and global settings
if (!exists("dir0")) dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
PHE_F_R <- Sys.getenv("PHE_F", unset = file.path(dir0, "scripts/0f/phenotype.sh"))
find_first_existing <- function(paths, default = paths[[1]]) {
	paths <- unique(paths[!is.na(paths) & nzchar(paths)])
	hit <- paths[file.exists(paths) | dir.exists(paths)]
	if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = FALSE) else default
}
Y <- if (exists("Y")) Y else Sys.getenv("Y", unset = "cvd_cad")
BIOM <- if (exists("BIOM")) BIOM else Sys.getenv("BIOM", unset = "prot,met")
LE8_JOB <- if (exists("LE8_JOB")) LE8_JOB else Sys.getenv("LE8_JOB", unset = "")
prot_DO <- truthy(Sys.getenv("PROT_DO", unset = "TRUE"))
met_DO <- truthy(Sys.getenv("MET_DO", unset = "TRUE"))
LE8_REPLACE <- truthy(Sys.getenv("LE8_REPLACE", unset = "FALSE"))
# Internal routing from le8.sh: completed results bypass analysis initialization.
LE8_REUSE_RESULTS <- !LE8_REPLACE && truthy(Sys.getenv("LE8_REUSE_RESULTS", unset = "FALSE"))
N_CORES <- max(1L, suppressWarnings(as.integer(Sys.getenv("N_CORES", unset = "2"))))
if (!is.finite(N_CORES)) N_CORES <- 1L
SEED <- as.integer(Sys.getenv("SEED", unset = "2026"))
set.seed(SEED)

indir <- find_first_existing(c(
	Sys.getenv("UKB_PHE", unset = ""), file.path(dir0, "data/ukb/phe"),
	"/mnt/ddata/ukb/phe", "D:/data/ukb/phe"
), file.path(dir0, "data/ukb/phe"))
common_dir <- find_first_existing(c(
	Sys.getenv("UKB_COMMON", unset = ""), file.path(indir, "common"),
	"/mnt/ddata/ukb/phe/common", "D:/data/ukb/phe/common"
), file.path(indir, "common"))
analysis_root <- Sys.getenv("LE8_ANALYSIS_ROOT", unset = file.path(dir0, "analysis/le8"))
out.base <- file.path(analysis_root, Y)
out.prot <- file.path(out.base, "prot")
out.met <- file.path(out.base, "met")
invisible(lapply(c(out.base, out.prot, out.met), dir.create, recursive = TRUE, showWarnings = FALSE))

dir.Y <- Sys.getenv("LE8_GWAS_DIR", unset = file.path("/mnt/d/data", "gwas/main"))
dir.X <- Sys.getenv("LE8_PQTL_IV_DIR", unset = file.path("/mnt/d/data", "gwas/prot"))
dir.met.gwas <- Sys.getenv("LE8_MQTL_IV_DIR", unset = file.path("/mnt/d/data", "gwas/met"))
gwas_trait_dir <- function(project_dir, trait, category = "common") file.path(project_dir, category, trait)
gwas_clean_dir <- function(project_dir, trait, category = "common") file.path(gwas_trait_dir(project_dir, trait, category), "gwas")
gwas_clean_file <- function(project_dir, trait, suffix = ".gz") file.path(gwas_clean_dir(project_dir, trait), paste0(trait, suffix))
gwas_magma_dir <- function(project_dir, trait, category = "common") file.path(gwas_trait_dir(project_dir, trait, category), "magma")
gwas_trait_dir_from_clean_file <- function(file) {
	gwas_dir <- dirname(file) ; trait_dir <- dirname(gwas_dir)
	if (basename(gwas_dir) != "gwas") return(NA_character_)
	trait_dir
}
gwas_magma_dir_from_clean_file <- function(file) {
	trait_dir <- gwas_trait_dir_from_clean_file(file)
	if (is.na(trait_dir)) NA_character_ else file.path(trait_dir, "magma")
}
prot_bed_file <- find_first_existing(
	c(
		Sys.getenv("LE8_PROT_BED", unset = ""),
		file.path(dir.X, "ppp_3k.b38.bed"),
		"/mnt/f/gwas/prot/ppp_3k.38.bed"
	),
	file.path(dir.X, "ppp_3k.b38.bed")
)
met_list_file <- find_first_existing(c(
	Sys.getenv("LE8_MET_LIST", unset = ""), file.path(common_dir, "met.lst"),
	"D:/data/ukb/phe/common/met.lst"
), file.path(common_dir, "met.lst"))
ukb_bgen_dir <- Sys.getenv("UKB_BGEN_DIR", unset = "/mnt/d/data/ukb/gen/imp")

# Source existing project helpers when present. The new code does not require association.R,
# but uses the user's existing t2e() definition from phenotype.R when available.
helper_names <- c("phenotype.R", "association.R", "plotting.R", "prediction.R")
helper_dirs <- unique(c(Sys.getenv("LE8_SHARED_HELPERS", unset = ""), file.path(dirname(Sys.getenv("DIRSCRIPT", unset = file.path(dir0, "scripts/le8"))), "0f"), file.path(dir0, "scripts/0f")))
for (f0 in unique(unlist(lapply(helper_dirs, function(d0) file.path(d0, helper_names))))) {
	if (file.exists(f0)) try(source(f0), silent = TRUE)
}
date_follow_end <- le8_follow_end()
if (!exists("vars.basic")) stop("vars.basic was not loaded from scripts/0f/phenotype.R.", call. = FALSE)
if (!exists("vars.le8")) stop("vars.le8 was not loaded from scripts/0f/phenotype.R.", call. = FALSE)
if (!exists("names.le8")) stop("names.le8 was not loaded from scripts/0f/phenotype.R.", call. = FALSE)
if (!exists("vars.adj2")) vars.adj2 <- unique(c(vars.basic, vars.le8))
# Analysis-wide covariate choice.  Keep vars.basic as the default; replace the
# next line with `covs_use <- vars.adj2` when downstream modules should consume
# the comprehensively adjusted C1 scan.
if (!exists("covs_use")) covs_use <- vars.adj2
covs_use_name <- if (identical(unique(covs_use), unique(vars.adj2))) "adj2" else "basic"
LE8_LABS <- c(
	diet = "Diet", pa = "Physical activity", smoke = "Smoking", bmi = "BMI",
	nonhdl = "Non-HDL-C", hba1c = "HbA1c", bp = "Blood pressure", sleep = "Sleep"
)
cols_le8 <- c(
	diet = "#E68613", pa = "#5B9BD5", smoke = "#CC79A7", bmi = "#E76F51",
	nonhdl = "#00A6B2", hba1c = "#009E73", bp = "#B79F00", sleep = "#4063D8"
)
cols_evidence <- c(
	Observational = "#4C78A8", MR = "#F58518", Colocalization = "#54A24B",
	Connection = "#B279A2", Prediction = "#E45756"
)


# 🚩 IO and output management
required_file <- function(path, label = "file") {
	if (is.na(path) || !nzchar(path) || !file.exists(path) || isTRUE(file.size(path) == 0))
		stop("Missing required ", label, ": ", path, call. = FALSE)
	normalizePath(path, winslash = "/", mustWork = FALSE)
}
read_all <- function(select_vars = NULL) {
	x <- readRDS(required_file(file.path(indir, "Rdata/all.rds"), "UKB phenotype all.rds"))
	x <- le8_select_phenotypes(x)
	x <- le8_rebuild_baseline(x, out.base)
	if (!is.null(select_vars)) x <- x[, intersect(unique(c(select_vars, le8_custom_covars, paste0("fod_icd10_",Y), Sys.getenv("LE8_GROUP_COLUMN",Sys.getenv("PGS_GROUP_COLUMN","")))), names(x)), drop = FALSE]
	x
}
read_prot <- function(required = TRUE) {
	f <- file.path(indir, "Rdata/prot.rds")
	if (!file.exists(f) || file.size(f) == 0) {
		if (required) stop("Missing proteomics RDS: ", f, call. = FALSE)
		return(NULL)
	}
	readRDS(f)
}
read_met <- function(required = TRUE) {
	f <- file.path(indir, "Rdata/met.rds")
	if (!file.exists(f) || file.size(f) == 0) {
		if (required) stop("Missing metabolomics RDS: ", f, call. = FALSE)
		return(NULL)
	}
	readRDS(f)
}
setwd2 <- function(path) {
	dir.create(path, recursive = TRUE, showWarnings = FALSE) ; setwd(path) ; invisible(path)
}
# All Final products share one stable tree, independent of native module outputs.
le8_final_dir <- function(layer = NULL, kind = NULL, trait = Y, root = analysis_root) {
	path <- file.path(root, "final", trait)
	if (!is.null(layer)) path <- file.path(path, layer)
	if (!is.null(kind)) path <- file.path(path, kind)
	path
}
le8_final_layer <- function(outdir) {
	parts <- strsplit(normalizePath(outdir, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]
	layers <- parts[parts %in% c("prot", "met")]
	if (!length(layers)) stop("Final layer output needs prot or met: ", outdir)
	tail(layers, 1)
}
le8_job_dir <- function(outdir = getwd(), job = LE8_JOB) {
	if (is.null(job) || !nzchar(job)) job <- "raw"
	if (job == "final_prediction") return(le8_final_dir(le8_final_layer(outdir), "prediction"))
	if (job == 'c4_panel_validation') job <- 'c4_connect'
	bn <- basename(normalizePath(outdir, winslash = "/", mustWork = FALSE))
	if (bn %in% c("prot", "met")) file.path(outdir, job) else outdir
}
# Module artifacts stay beside their numerical results; root is publication-only.
le8_artifact_path <- function(file, outdir = getwd(), module = LE8_JOB) {
	if (module == "final_prediction") outdir <- le8_job_dir(outdir, module)
	if (module == 'c4_panel_validation') module <- 'c4_connect'
	target <- if (grepl("^(/|[A-Za-z]:)", file)) file else file.path(outdir, file)
	parent <- dirname(target)
	if (basename(parent) %in% c("prot", "met") && nzchar(module) && module != "final")
		target <- file.path(parent, module, basename(target))
	dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
	target
}
safe_sheet <- function(x) {
	x <- gsub("[\\[\\]:*?/\\\\]", "_", x)
	x <- substr(x, 1, 31) ; make.unique(x, sep = "_")
}
write_raw_csv <- function(x, file, rawdir = le8_job_dir()) {
	dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
	if (is.null(x) || !is.data.frame(x) || ncol(x) == 0L) x <- data.frame(note = "No rows were produced for this optional output.")
	data.table::fwrite(as.data.frame(x), file.path(rawdir, file), sep = ",", na = "NA")
	invisible(file.path(rawdir, file))
}
write_raw_tsv <- function(x, file, rawdir = le8_job_dir()) {
	dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
	if (is.null(x) || !is.data.frame(x) || ncol(x) == 0L) x <- data.frame(note = "No rows were produced for this optional output.")
	data.table::fwrite(as.data.frame(x), file.path(rawdir, file), sep = "\t", na = "")
	invisible(file.path(rawdir, file))
}
write_xlsx2 <- function(x, file) {
	file <- le8_artifact_path(file)
	if (is.data.frame(x)) x <- list(data = x)
	x <- x[!vapply(x, is.null, logical(1))]
	if (!length(x)) x <- list(note = data.frame(note = "No table was produced."))
	names(x) <- safe_sheet(names(x))
	x <- lapply(x, function(z) {
		if (is.matrix(z)) z <- as.data.frame(z)
		if (!is.data.frame(z)) z <- data.frame(value = I(list(z)))
		# openxlsx::write.xlsx(asTable = TRUE) cannot serialize a data frame with
		# zero columns and fails internally with `df[[i]]: subscript out of bounds`.
		# Keep the worksheet, but make the absence of rows explicit.
		if (ncol(z) == 0L) z <- data.frame(note = "No table was produced.")
		z
	})
	openxlsx::write.xlsx(x, file, overwrite = TRUE, asTable = TRUE, freezePane = TRUE, autoFilter = TRUE)
	invisible(file)
}
module_cache <- function(outdir, module = LE8_JOB) file.path(le8_job_dir(outdir, module), paste0(module, ".res.rds"))
cache_valid <- function(path) {
	if (LE8_REPLACE || length(path) != 1L || is.na(path)) return(FALSE)
	info <- file.info(path)
	isTRUE(!info$isdir && info$size > 0)
}
# Reuse stage results only when their internal settings match the requested analysis.
.le8_loaded_code <- tools::md5sum(list.files(Sys.getenv("LE8_FDIR",file.path(Sys.getenv("DIRSCRIPT"),"f")),pattern="[.](R|py|sh)$",full.names=TRUE))
le8_stage_fingerprint <- function() {
	code <- list.files(Sys.getenv("LE8_FDIR",file.path(Sys.getenv("DIRSCRIPT"),"f")),pattern="[.](R|py|sh)$",full.names=TRUE)
	env <- Sys.getenv(); env <- env[grepl("^(C[1-5]_|PGS_|RUN_|LE8_(GWAS|PQTL|MQTL|LD|GROUP|GRCH|REFGEN|ENDPOINT|BASELINE|OUTER|Y_DATE|VARS_ADJ|VARIANT|REFERENCE|EVIDENCE)|FINAL_CONNECTION|DATE_FOLLOW_END)",names(env))]
	files <- unique(c(file.path(get0("indir",ifnotfound=""),"Rdata",c("all.rds","prot.rds","met.rds","prot.pgs.rds","met.pgs.rds")),
		unname(env[file.exists(env)]),get0(".le8_stage_source_files",ifnotfound=character()),file.path(get0("indir",ifnotfound=""),"rap/vip.tab.gz")))
	for (mf in env[grepl("MANIFEST$",names(env)) & file.exists(env)]) {
		m <- tryCatch(data.table::fread(mf,showProgress=FALSE),error=function(e) NULL)
		if (!is.null(m)) for (nm in intersect(c("file","path","weights_file","ld_file","normalization_proof","reference_fasta"),names(m))) {
			v <- as.character(m[[nm]]); relative <- !grepl("^/|^[A-Za-z]:",v)
			v[relative] <- file.path(dirname(mf),v[relative]); files <- unique(c(files,v))
		}
	}
	info <- file.info(files)
	# Transaction workspaces move between runs; identify upstream results by
	# their stable path within the analysis root, while inspecting actual files.
	root <- paste0(sub("/+$", "", get0("analysis_root", ifnotfound="")), "/")
	logical_files <- files
	inside <- !is.na(files) & nzchar(root) & root != "/" & startsWith(files, root)
	logical_files[inside] <- paste0("<analysis-root>/", substring(files[inside], nchar(root) + 1L))
	le8_hash_object(list(code=.le8_loaded_code,env=env,files=logical_files,size=info$size,mtime=as.numeric(info$mtime)))
}
read_stage_cache <- function(path, version = NULL) {
	if (!cache_valid(path)) return(NULL)
	z <- tryCatch(readRDS(path), error = function(e) NULL)
	if (is.null(z) || !is.list(z) || is.null(z$data)) return(NULL)
	if (!identical(z$source_signature,le8_stage_fingerprint())) return(NULL)
	if (!is.null(version) && !identical(z$version, version)) return(NULL)
	if (!is.null(z$analysis_options) && !identical(z$analysis_options, le8_analysis_options())) return(NULL)
	z$data
}
write_stage_cache <- function(data, path, version = NULL) {
	dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
	# Publish only after serialization finishes; interruptions leave the old cache intact.
	tmp <- tempfile(pattern = ".stage-", tmpdir = dirname(path))
	on.exit(unlink(tmp), add = TRUE)
	value <- list(generated = format(Sys.time(), "%F %T %z"), data = data, code_signature=le8_code_fingerprint(), source_signature=le8_stage_fingerprint(), analysis_options = le8_analysis_options())
	if (!is.null(version)) value$version <- version
	saveRDS(value, tmp, compress = "xz")
	if (!file.rename(tmp, path)) stop("Cannot publish stage cache: ", path, call. = FALSE)
	invisible(data)
}
cache_message <- function(label, path) {
	message(
		label, " already exists, skip ", label, " step: ", path,
		". Delete this existing file (or its step folder) to re-run."
	)
}
# Emit only major stage boundaries; routine package output stays in the file log.
le8_stage_start <- function(label, detail = "") {
	message("[LE8] START ", label, if (nzchar(detail)) paste0(" ", detail) else "")
	proc.time()[["elapsed"]]
}
le8_stage_done <- function(label, started, detail = "") {
	message(
		"[LE8] DONE ", label, " elapsed=", sprintf("%.1f", (proc.time()[["elapsed"]] - started) / 60),
		" min", if (nzchar(detail)) paste0(" ", detail) else ""
	)
	invisible(NULL)
}
le8_stage <- function(label, expr, detail = "") {
	started <- le8_stage_start(label, detail)
	value <- tryCatch(withVisible(force(expr)), error = function(e) {
		message("[LE8] FAIL ", label, ": ", conditionMessage(e))
		stop(e)
	})
	status <- ""
	if (is.list(value$value) && is.character(value$value$status) && length(value$value$status) == 1L)
		status <- paste0("status=", value$value$status)
	if (is.numeric(value$value) && length(value$value) == 1L && is.finite(value$value))
		status <- paste0("exit=", value$value)
	failed <- (is.numeric(value$value) && length(value$value) == 1L && is.finite(value$value) && value$value != 0) ||
		(nzchar(status) && grepl("status=failed", status, fixed = TRUE))
	if (failed) message("[LE8] FAIL ", label, " ", status) else le8_stage_done(label, started, status)
	if (value$visible) value$value else invisible(value$value)
}
parallel_map <- function(x, fun) {
	# PGS invokes six one-feature models; the enclosing feature loop owns GC.
	# Avoid repeatedly walking the entire resident score matrix for each model.
	if (length(x) <= 1L) return(lapply(x, fun))
	# Release unreachable model frames before forking so children inherit a
	# smaller heap. Collect between tasks after fun's local frame has returned;
	# otherwise each long-lived worker can retain several fits' worth of garbage.
	invisible(gc())
	run_one <- function(item) {
		value <- fun(item)
		invisible(gc())
		value
	}
	if (.Platform$OS.type != "windows" && N_CORES > 1L && length(x) > 1L) {
		message("Parallel scan: ", min(N_CORES, length(x)), " workers for ", length(x), " tasks")
		# Wrap successful results so an intentional NULL can be distinguished from
		# a child killed by OOM. Never silently save an incomplete association scan.
		result <- parallel::mclapply(x, function(item) list(value = run_one(item)),
			mc.cores = min(N_CORES, length(x)), mc.preschedule = TRUE,
			mc.allow.recursive = FALSE
		)
		failed <- vapply(result, function(z) is.null(z) || inherits(z, "try-error"), logical(1))
		if (any(failed)) stop("Parallel scan worker failed or was killed; refusing partial results. Reduce --cores and inspect the task log.", call. = FALSE)
		lapply(result, `[[`, "value")
	} else lapply(x, run_one)
}
# Third-party analysis packages often print one routine progress line per
# feature.  With thousands of omic traits that output overwhelms the useful
# stage-level messages.  Capture ordinary output and messages at the package
# boundary while deliberately leaving warnings and errors visible.
quiet_package_call <- function(expr) {
	value <- NULL
	invisible(utils::capture.output(
		value <- suppressMessages(force(expr)),
		type = "output"
	))
	value
}
module_meta <- function(layer, module = LE8_JOB, extra = list()) {
	c(list(
		module = module, layer = layer, trait = Y, generated = format(Sys.time(), "%F %T %z"),
		seed = SEED, R_runtime = R.version.string, biom = BIOM, configuration_sources=jsonlite::fromJSON(Sys.getenv("LE8_SHARED_CONFIG_SOURCES","{}")), module_policy=le8_module_policy(module), code_signature=le8_code_fingerprint(), source_signature=le8_stage_fingerprint(), analysis_options = le8_analysis_options()
	), extra)
}
finalize_outputs <- function(module, outdir = getwd()) {
	rawdir <- le8_job_dir(outdir, module) ; dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
	if (exists("le8_flush_figures", mode = "function")) le8_flush_figures(rawdir)
	png <- list.files(rawdir, pattern = "\\.png$", full.names = FALSE)
	xlsx <- list.files(rawdir, pattern = "\\.xlsx$", full.names = FALSE)
	index <- data.frame(
		module = module, trait = Y, directory = normalizePath(rawdir, winslash = "/", mustWork = FALSE),
		png = paste(png, collapse = "; "), workbook = paste(xlsx, collapse = "; "),
		generated = format(Sys.time(), "%F %T %z")
	)
	write_raw_csv(index, if (module == "final_prediction") "output_index.csv" else paste0(module, ".output_index.csv"), rawdir)
	message("Completed ", module, ": ", outdir)
	invisible(index)
}


# 🚩 Plotting
theme_5c <- function(base_size = 12) {
	theme_classic(base_size = base_size) +
		theme(
			plot.title = element_text(face = "bold", hjust = 0),
			plot.subtitle = element_text(color = "grey35"),
			axis.text = element_text(color = "black"),
			legend.title = element_text(face = "bold"),
			panel.grid.major.y = element_line(color = "grey91", linewidth = 0.25),
			strip.background = element_blank(), strip.text = element_text(face = "bold"),
			plot.margin = margin(8, 12, 8, 12)
		)
}
save_plot <- function(p, file, w = 8, h = 6, dpi = 320, outdir = getwd()) {
	target <- le8_artifact_path(file, outdir)
	if (exists("le8_queue_figure", mode = "function")) return(le8_queue_figure(p, target, w, h, dpi))
	ggplot2::ggsave(target, p, width = w, height = h, dpi = dpi, bg = "white", limitsize = FALSE)
	invisible(target)
}
blank_plot <- function(title, subtitle = "Required data were not available") {
	p <- ggplot() +
		theme_void(base_size = 13) +
		annotate("text", x = 0, y = .08, label = title, fontface = "bold", size = 5) +
		annotate("text", x = 0, y =  - .06, label = subtitle, color = "grey35", size = 3.7)
	attr(p, "le8_unavailable") <- paste(title, subtitle, sep = ": ") ; p
}
forest_theme <- function(base_size = 10) theme_5c(base_size) + theme(panel.grid.major.y = element_blank())


# 🚩 Phenotype and association models

# Add exact attained-age entry/exit columns for a delayed-entry Cox model.  A
# participant contributes risk time only after the baseline blood draw.  This
# avoids the immortal-time error that would arise from pretending the adult
# omic measurement had been observed continuously since birth.
add_attained_age_time <- function(dat, outcome = Y, entry_name = ".attained_entry",
			exit_name = ".attained_exit") {
	tvar <- paste0(outcome, ".t2e") ; bivar <- paste0(outcome, ".bi2e")
	if (!all(c(tvar, bivar) %in% names(dat))) stop("Attained-age time needs ", tvar, " and ", bivar, call. = FALSE)
	dat[[exit_name]] <- suppressWarnings(as.numeric(dat[[bivar]]))
	dat[[entry_name]] <- dat[[exit_name]] - suppressWarnings(as.numeric(dat[[tvar]]))
	invalid <- !is.finite(dat[[entry_name]]) | !is.finite(dat[[exit_name]]) |
		dat[[entry_name]] < 0 | dat[[exit_name]] <= dat[[entry_name]]
	dat[[entry_name]][invalid] <- NA_real_ ; dat[[exit_name]][invalid] <- NA_real_
	dat
}
filter_analysis_cohort <- function(dat) {
	if ("ethnic.c" %in% names(dat) && truthy(Sys.getenv("LE8_WHITE_ONLY", unset = "TRUE"))) {
		keep <- as.character(dat$ethnic.c) %in% c("White", "1")
		dat <- dat[keep, , drop = FALSE]
	}
	dat
}
cox_scan <- function(dat, xs, covars, outcome = Y, scale_x = TRUE, min_n = 500, min_event = 20,
			time_var = paste0(outcome, ".t2e"), event_var = paste0(outcome, ".Yt2e")) {
	tvar <- time_var ; evar <- event_var
	xs <- intersect(xs, names(dat)) ; covars <- intersect(covars, names(dat))
	bind_rows(parallel_map(xs, function(x) {
		need <- unique(c(tvar, evar, x, covars)) ; d <- dat[, need, drop = FALSE]
		d <- d[stats::complete.cases(d), , drop = FALSE]
		ne <- sum(d[[evar]] == 1, na.rm = TRUE)
		empty <- tibble(
			term = x, estimate = NA_real_, beta = NA_real_, std.error = NA_real_, conf.low = NA_real_,
			conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_, N_total = nrow(d), N_event = ne
		)
		if (nrow(d) < min_n || ne < min_event || length(unique(d[[evar]])) < 2) return(empty)
		d[[x]] <- suppressWarnings(as.numeric(d[[x]])) ; sx <- sd(d[[x]], na.rm = TRUE)
		if (!is.finite(sx) || sx == 0) return(empty)
		if (scale_x) d[[x]] <- as.numeric(scale(d[[x]]))
		fit <- tryCatch(coxph(as.formula(paste0(
			"Surv(", bt(tvar), ",", bt(evar), ") ~ ",
			paste(bt(c(x, covars)), collapse = " + ")
		)), d, ties = "efron"), error = function(e) NULL)
		if (is.null(fit)) return(empty)
		sm <- coef(summary(fit)) ; if (!x %in% rownames(sm)) return(empty)
		b <- sm[x, "coef"] ; se <- sm[x, "se(coef)"]
		tibble(
			term = x, estimate = exp(b), beta = b, std.error = se, conf.low = exp(b - 1.96 * se),
			conf.high = exp(b + 1.96 * se), statistic = b / se,
			p.value = 2 * pnorm(abs(b / se), lower.tail = FALSE), N_total = nrow(d), N_event = ne
		)
	})) |>
		mutate(FDR = p.adjust(p.value, "BH")) |>
		arrange(p.value)
}

cox_scan_delayed_entry <- function(dat, xs, covars, outcome = Y, scale_x = TRUE,
			min_n = 500, min_event = 20,
			entry_var = ".attained_entry", exit_var = ".attained_exit",
			event_var = paste0(outcome, ".Yt2e")) {
	xs <- intersect(xs, names(dat)) ; covars <- intersect(covars, names(dat))
	bind_rows(parallel_map(xs, function(x) {
		need <- unique(c(entry_var, exit_var, event_var, x, covars)) ; d <- dat[, need, drop = FALSE]
		d <- d[stats::complete.cases(d), , drop = FALSE]
		d <- d[d[[exit_var]] > d[[entry_var]], , drop = FALSE]
		ne <- sum(d[[event_var]] == 1, na.rm = TRUE)
		empty <- tibble(
			term = x, estimate = NA_real_, beta = NA_real_, std.error = NA_real_, conf.low = NA_real_,
			conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_, N_total = nrow(d), N_event = ne,
			time_scale = "attained age with delayed entry"
		)
		if (nrow(d) < min_n || ne < min_event || length(unique(d[[event_var]])) < 2) return(empty)
		d[[x]] <- suppressWarnings(as.numeric(d[[x]])) ; sx <- sd(d[[x]], na.rm = TRUE)
		if (!is.finite(sx) || sx == 0) return(empty)
		if (scale_x) d[[x]] <- as.numeric(scale(d[[x]]))
		f <- as.formula(paste0(
			"Surv(", bt(entry_var), ",", bt(exit_var), ",", bt(event_var), ") ~ ",
			paste(bt(c(x, covars)), collapse = " + ")
		))
		fit <- tryCatch(coxph(f, d, ties = "efron"), error = function(e) NULL)
		if (is.null(fit)) return(empty) ; sm <- coef(summary(fit)) ; if (!x %in% rownames(sm)) return(empty)
		b <- sm[x, "coef"] ; se <- sm[x, "se(coef)"]
		tibble(
			term = x, estimate = exp(b), beta = b, std.error = se, conf.low = exp(b - 1.96 * se), conf.high = exp(b + 1.96 * se),
			statistic = b / se, p.value = 2 * pnorm(abs(b / se), lower.tail = FALSE), N_total = nrow(d), N_event = ne,
			time_scale = "attained age with delayed entry"
		)
	})) |>
		mutate(FDR = p.adjust(p.value, "BH")) |>
		arrange(p.value)
}
lm_scan <- function(dat, xs, y, covars, scale_x = TRUE, min_n = 500) {
	xs <- intersect(xs, names(dat)) ; covars <- intersect(covars, names(dat))
	bind_rows(parallel_map(xs, function(x) {
		d <- dat[, unique(c(y, x, covars)), drop = FALSE] ; d <- d[complete.cases(d), , drop = FALSE]
		empty <- tibble(term = x, beta = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, N_total = nrow(d))
		if (nrow(d) < min_n) return(empty)
		d[[x]] <- suppressWarnings(as.numeric(d[[x]])) ; d[[y]] <- suppressWarnings(as.numeric(d[[y]]))
		if (!is.finite(sd(d[[x]], na.rm = TRUE)) || sd(d[[x]], na.rm = TRUE) == 0) return(empty)
		if (scale_x) d[[x]] <- std_num(d[[x]])
		fit <- tryCatch(lm(as.formula(paste0(bt(y), " ~ ", paste(bt(c(x, covars)), collapse = " + "))), d), error = function(e) NULL)
		if (is.null(fit)) return(empty)
		sm <- coef(summary(fit)) ; if (!x %in% rownames(sm)) return(empty)
		tibble(
			term = x, beta = sm[x, "Estimate"], std.error = sm[x, "Std. Error"], statistic = sm[x, "t value"],
			p.value = sm[x, "Pr(>|t|)"], N_total = nrow(d)
		)
	})) |>
		mutate(FDR = p.adjust(p.value, "BH")) |>
		arrange(p.value)
}
partial_r2 <- function(full, reduced) {
	r2f <- summary(full)$r.squared ; r2r <- summary(reduced)$r.squared
	pmax(0, pmin(1, r2f - r2r))
}


# 🚩 Protein and metabolite annotations
read_prot_bed <- function(proteins = NULL) {
	f <- required_file(prot_bed_file, "protein BED annotation")
	x <- data.table::fread(f, header = FALSE, fill = TRUE, showProgress = FALSE)
	if (ncol(x) < 5) stop("protein BED must have at least five columns (chr,start,end,protein,group): ", f, call. = FALSE)
	ans <- tibble(
		chr = as.character(x[[1]]), start = as.numeric(x[[2]]), end = as.numeric(x[[3]]),
		protein = as.character(x[[4]]), group = as.character(x[[5]])
	) |>
		mutate(
			chr = str_remove(chr, "^chr"), pos = floor((start + end) / 2),
			group = coalesce(na_if(group, ""), "Other")
		) |>
		distinct(protein, .keep_all = TRUE)
	if (!is.null(proteins)) tibble(protein = proteins) |> left_join(ans, by = "protein") else ans
}
met_super_group <- function(group, subgroup = NA_character_) {
	txt <- str_to_lower(paste(coalesce(group, ""), coalesce(subgroup, "")))
	case_when(
		str_detect(txt, "amino") ~ "Amino acids",
		str_detect(txt, "apolipoprotein|cholesterol|cholesteryl") ~ "Cholesterol and apolipoproteins",
		str_detect(txt, "fatty acid|omega|pufa|mufa|sfa") ~ "Fatty acids",
		str_detect(txt, "triglyceride|lipoprotein|hdl|ldl|vldl|lipid|phospholipid|sphingomyelin|ceramide") ~ "Lipoprotein lipids",
		str_detect(txt, "glycolysis|glucose|ketone|citrate|lactate") ~ "Energy metabolism",
		TRUE ~ "Other metabolites"
	)
}
read_met_annotation <- function(vars.met = NULL) {
	f <- required_file(met_list_file, "met.lst")
	# Let fread detect whether met.lst has a header.
	x <- data.table::fread(f, fill = TRUE, check.names = FALSE, showProgress = FALSE)
	nm0 <- names(x) ; nml <- tolower(nm0)
	pick <- function(pattern, fallback) {
		i <- grep(pattern, nml)[1] ; if (is.na(i)) i <- fallback ; nm0[[i]]
	}
	var_col <- pick("^var_name$|^variable$|^trait$|^field$|^name$", 1)
	if (!is.null(vars.met) && length(vars.met)) {
		vm <- unique(c(vars.met, str_remove(vars.met, "^met_"), paste0("met_", str_remove(vars.met, "^met_"))))
		overlap <- vapply(x, function(z) sum(unique(as.character(z)) %in% vm, na.rm = TRUE), numeric(1))
		if (length(overlap) && max(overlap, na.rm = TRUE) > 0) var_col <- nm0[[which.max(overlap)]]
	}
	label_col <- pick("^full_name$|biomarker|description|label|long_name", min(2, ncol(x)))
	# Prefer the explicit header in current met.lst files.  Only use the
	# use the penultimate column when no group header exists.
	group_i <- grep("^group$|^class$|^category$", nml)[1]
	group_col <- if (!is.na(group_i)) nm0[[group_i]] else nm0[[max(1, ncol(x) - 1)]]
	subgroup_col <- if (any(grepl("subgroup|sub_group|class", nml))) nm0[[grep("subgroup|sub_group|class", nml)[1]]] else nm0[[ncol(x)]]
	ans <- x |>
		transmute(
			trait0 = as.character(.data[[var_col]]), label = as.character(.data[[label_col]]),
			group = as.character(.data[[group_col]]), subgroup = as.character(.data[[subgroup_col]])
		) |>
		mutate(
			trait_prefixed = ifelse(str_detect(trait0, "^met_"), trait0, paste0("met_", trait0)),
			trait_stripped = str_remove(trait0, "^met_"),
			trait = case_when(
				!is.null(vars.met) & trait0 %in% vars.met ~ trait0,
				!is.null(vars.met) & trait_prefixed %in% vars.met ~ trait_prefixed,
				!is.null(vars.met) & trait_stripped %in% vars.met ~ trait_stripped,
				TRUE ~ trait_prefixed
			),
			label = coalesce(na_if(label, ""), trait), group = coalesce(na_if(group, ""), "Other"),
			subgroup = coalesce(na_if(subgroup, ""), group), super_group = met_super_group(group, subgroup)
		) |>
		select( - trait_prefixed, - trait_stripped) |>
		distinct(trait, .keep_all = TRUE)
	if (!is.null(vars.met)) tibble(trait = vars.met) |>
		left_join(ans, by = "trait") |>
		mutate(
			label = coalesce(label, trait), group = coalesce(group, "Other"), subgroup = coalesce(subgroup, group),
			super_group = coalesce(super_group, "Other metabolites")
		) else ans
}
le8_gene_coordinates <- function(layer, features) {
	if (layer == "protein") read_prot_bed(features) |> transmute(feature = protein, label = protein, group, chr, start, end, pos)
	else read_met_annotation(features) |> transmute(feature = trait, label, group, subgroup, super_group)
}

layer_annotation_audit <- function(layer, features) {
	features <- unique(as.character(features))
	if (layer == "protein") {
		ref <- read_prot_bed() |> transmute(feature = protein, annotation_group = group, annotation_key = protein)
		return(tibble(feature = features) |> left_join(ref, by = "feature") |>
			mutate(annotation_matched = !is.na(annotation_key), annotation_source = prot_bed_file))
	}
	# Pass the assayed feature names so read_met_annotation() can identify the
	# variable column even when met.lst is headerless or its variable-name
	# column is not first.  Calling it without this overlap information caused
	# a false 0/300 annotation-match audit in the current CAD output.
	ref <- read_met_annotation(features) |>
		transmute(
			annotation_key = trait, .key = str_remove(trait, "^met_"),
			annotation_label = label, annotation_group = group, annotation_subgroup = subgroup, annotation_super_group = super_group
		) |>
		distinct(.key, .keep_all = TRUE)
	tibble(feature = features, .key = str_remove(feature, "^met_")) |>
		left_join(ref, by = ".key") |>
		mutate(annotation_matched = !is.na(annotation_key), annotation_source = met_list_file) |>
		select( - .key)
}
qtl_name_candidates <- function(feature, layer) {
	z <- unique(c(feature, str_remove(feature, "^met_"), str_replace_all(str_remove(feature, "^met_"), "[^A-Za-z0-9_.-]", "_")))
	z[nzchar(z)]
}
find_qtl_directory <- function(feature, base_dir, layer) {
	cand <- qtl_name_candidates(feature, layer)
	hit <- cand[dir.exists(vapply(cand, function(x) gwas_clean_dir(base_dir, x), character(1)))][1]
	if (length(hit) && !is.na(hit)) list(name = hit, dir = gwas_clean_dir(base_dir, hit)) else list(name = cand[[1]], dir = gwas_clean_dir(base_dir, cand[[1]]))
}
find_qtl_files <- function(feature, base_dir, layer = c("protein", "metabolite")) {
	layer <- match.arg(layer) ; loc <- find_qtl_directory(feature, base_dir, layer) ; nm <- loc$name ; d <- loc$dir
	first <- function(x) {
		y <- x[file.exists(x) & file.size(x) > 0] ; if (length(y)) normalizePath(y[[1]], winslash = "/", mustWork = FALSE) else NA_character_
	}
	list(
		feature = feature, qtl_name = nm, dir = d,
		joint = first(c(file.path(d, paste0(nm, ".jma.cojo")), file.path(d, paste0(feature, ".jma.cojo")))),
		full = first(c(file.path(d, paste0(nm, ".gz")), file.path(d, paste0(feature, ".gz")), file.path(d, paste0(nm, ".4gcta")))),
		cis = first(c(file.path(d, paste0(nm, ".cis.gz")), file.path(d, paste0(feature, ".cis.gz")))),
		trans = first(c(file.path(d, paste0(nm, ".trans.gz")), file.path(d, paste0(nm, ".clump.assoc")), file.path(d, paste0(feature, ".clump.assoc"))))
	)
}


# 🚩 Standardised summary statistics, harmonisation, MR and R2
.norm_name <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
.norm_chr_value <- function(x) {
	z <- toupper(sub("^chr", "", as.character(x), ignore.case = TRUE))
	z <- sub("\\.0$", "", z)
	z[z %in% "X"] <- "23" ; z[z %in% "Y"] <- "24" ; z[z %in% c("M", "MT")] <- "25"
	sub("^0+", "", z)
}
.pick_col <- function(nms, candidates, prefer = NULL) {
	nn <- .norm_name(nms) ; cc <- .norm_name(candidates)
	if (!is.null(prefer)) {
		pp <- .norm_name(prefer) ; i <- match(pp, nn) ; i <- i[!is.na(i)] ; if (length(i)) return(nms[[i[[1]]]])
	}
	i <- match(cc, nn) ; i <- i[!is.na(i)] ; if (length(i)) nms[[i[[1]]]] else NA_character_
}
read_table_auto <- function(file, nrows =  - 1) {
	if (is.na(file) || !file.exists(file) || file.size(file) == 0) return(data.table())
	data.table::fread(file, nrows = nrows, fill = TRUE, showProgress = FALSE, check.names = FALSE)
}
match_GRCH_table <- function(query_file, reference_file, position_only = FALSE) {
	query_file <- required_file(query_file, "match_GRCH query")
	reference_file <- required_file(reference_file, "match_GRCH reference")
	phe_f <- required_file(PHE_F_R, "phenotype.sh")
	bash <- Sys.which("bash")
	if (!nzchar(bash)) stop("bash is required to run match_GRCH.", call. = FALSE)

	runner <- tempfile("le8_match_grch_", fileext = ".sh")
	output <- tempfile("le8_match_grch_", fileext = ".tsv")
	audit <- tempfile("le8_match_grch_audit_", fileext = ".tsv")
	on.exit(unlink(c(runner, output, audit), force = TRUE), add = TRUE)
	flag <- if (isTRUE(position_only)) " --position-only" else ""
	writeLines(
		c(
			"#!/usr/bin/env bash", "set -euo pipefail", 'source "$1"',
			paste0('match_GRCH --reference "$3" --output "$4" --audit "$5"', flag, ' "$2"')
		),
		runner,
		useBytes = TRUE
	)
	log <- suppressWarnings(system2(bash, shQuote(c(runner, phe_f, query_file, reference_file, output, audit)),
		stdout = TRUE, stderr = TRUE
	))
	status <- attr(log, "status") ; if (is.null(status)) status <- 0L
	if (status != 0L || !file.exists(output)) {
		detail <- paste(tail(log, 20), collapse = "\n")
		stop("match_GRCH failed for ", query_file, " against ", reference_file,
			if (nzchar(detail)) paste0(":\n", detail) else "",
			call. = FALSE
		)
	}
	d <- read_table_auto(output)
	attr(d, "match_GRCH_audit") <- if (file.exists(audit)) read_table_auto(audit) else data.table()
	d
}
read_sumstat_matched <- function(query_file, reference_file,
			N_default = suppressWarnings(as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "NA"))),
			joint = grepl("jma\\.cojo$", query_file %||% ""), position_only = FALSE) {
	standardize_sumstat(match_GRCH_table(query_file, reference_file, position_only = position_only),
		N_default,
		joint = joint, source_file = query_file
	)
}
le8_variant_key <- function(z) {
	pair <- paste(pmin(z$EA,z$NEA),pmax(z$EA,z$NEA),sep=":")
	snv <- nchar(z$EA)==1 & nchar(z$NEA)==1 & !is.na(z$EA) & !is.na(z$NEA)
	a <- chartr("ATCG","TAGC",z$EA); b <- chartr("ATCG","TAGC",z$NEA)
	pair[snv] <- pmin(pair[snv],paste(pmin(a[snv],b[snv]),pmax(a[snv],b[snv]),sep=":"))
	build <- if("BUILD" %in% names(z)) as.character(z$BUILD) else rep("unknown",nrow(z))
	build[is.na(build)|!nzchar(build)] <- "unknown"
	key <- paste(build,.norm_chr_value(z$CHR),format(z$POS,scientific=FALSE,trim=TRUE),pair,sep=":")
	missing <- is.na(z$CHR)|!is.finite(z$POS)|is.na(z$EA)|is.na(z$NEA)
	key[missing] <- paste(build[missing],z$SNP[missing],z$EA[missing],z$NEA[missing],sep=":")
	key
}
.le8_identity_hashes <- new.env(parent=emptyenv())
le8_file_sha256 <- function(file) {
	if (!file.exists(file)) stop("Missing identity/proof input: ",file)
	info <- file.info(file); key <- paste(normalizePath(file),info$size,as.numeric(info$mtime),sep="|")
	if (!exists(key,envir=.le8_identity_hashes,inherits=FALSE)) assign(key,digest::digest(file=file,algo="sha256"),envir=.le8_identity_hashes)
	get(key,envir=.le8_identity_hashes,inherits=FALSE)
}
le8_variant_identity <- function(z, metadata=le8_gwas_metadata(z$source_file[1])) {
	z$BUILD <- rep(as.character(metadata$build),nrow(z))
	z$variant_id <- rep(NA_character_,nrow(z)); z$normalization_status <- rep("orientation_key_only",nrow(z))
	z$normalization_proof_hash <- rep(NA_character_,nrow(z))
	if (!nrow(z)) return(z)
	sequence <- !is.na(z$REF)&!is.na(z$ALT)&grepl("^[ACGT]+$",z$REF)&grepl("^[ACGT]+$",z$ALT)&z$REF!=z$ALT
	indel <- nchar(z$EA)!=1L | nchar(z$NEA)!=1L
	z$normalization_status[indel %in% TRUE] <- "indel_normalization_unverified"
	proof <- metadata$normalization_proof
	if (!is.null(proof) && !is.na(proof) && nzchar(proof)) {
		if (is.na(metadata$build) || !nzchar(metadata$build)) stop("Normalization proof requires genome build")
		a <- jsonlite::fromJSON(proof)
		if (!all(c("file_sha256","reference_sha256","tool","tool_version","build","reference_checked","left_aligned","multiallelic_split") %in% names(a)) ||
			!isTRUE(a$reference_checked) || !isTRUE(a$left_aligned) || !isTRUE(a$multiallelic_split) ||
			!identical(as.character(a$build),as.character(metadata$build)) || !nzchar(a$tool) || !nzchar(a$tool_version) ||
			!grepl("^[a-fA-F0-9]{64}$",a$reference_sha256) || !identical(tolower(a$file_sha256),le8_file_sha256(z$source_file[1]))) stop("Invalid or stale GWAS normalization proof")
		alleles_match <- (z$EA==z$ALT & z$NEA==z$REF) | (z$EA==z$REF & z$NEA==z$ALT)
		if (any(sequence & !alleles_match,na.rm=TRUE)) stop("Effect alleles conflict with proven REF/ALT identity")
		ok <- (sequence & alleles_match & is.finite(z$POS) & !is.na(z$CHR)) %in% TRUE
		z$variant_id[ok] <- paste(z$BUILD[ok],z$CHR[ok],format(z$POS[ok],scientific=FALSE,trim=TRUE),z$REF[ok],z$ALT[ok],sep=":")
		z$normalization_status[ok] <- "reference_checked_left_aligned_split"
		z$normalization_proof_hash[ok] <- le8_file_sha256(proof)
	}
	z
}
le8_normalize_variants <- function(z, metadata) {
	if (!nrow(z)) return(z)
	fasta <- metadata$reference_fasta
	if (is.null(fasta) || is.na(fasta) || !nzchar(fasta)) fasta <- Sys.getenv(paste0("LE8_REFERENCE_FASTA_",metadata$build),"")
	if (!nzchar(fasta) || all(z$normalization_status=="reference_checked_left_aligned_split")) return(z)
	if (!file.exists(fasta) || !file.exists(paste0(fasta,".fai"))) stop("Variant normalization needs an indexed reference FASTA")
	if (!nzchar(Sys.which("bcftools"))) stop("Variant normalization needs bcftools")
	# Only explicit REF/ALT biallelic associations can be normalized. EA/NEA sorting is not REF inference.
	valid <- !is.na(z$REF)&!is.na(z$ALT)&grepl("^[ACGT]+$",z$REF)&grepl("^[ACGT]+$",z$ALT)&
		((z$EA==z$ALT&z$NEA==z$REF)|(z$EA==z$REF&z$NEA==z$ALT))&is.finite(z$POS)
	ix <- which(valid %in% TRUE & z$normalization_status!="reference_checked_left_aligned_split")
	if (!length(ix)) return(z)
	fai <- data.table::fread(paste0(fasta,".fai"),header=FALSE,colClasses=c(V1="character"))
	contig <- as.character(fai$V1[match(.norm_chr_value(z$CHR[ix]),.norm_chr_value(fai$V1))])
	if (anyNA(contig)) stop("Reference lacks requested chromosome")
	tmp <- tempfile("le8-normalize-",tmpdir="/tmp"); dir.create(tmp); on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
	input <- file.path(tmp,"input.vcf"); output <- file.path(tmp,"normalized.vcf")
	header <- c("##fileformat=VCFv4.2",paste0("##contig=<ID=",fai$V1,",length=",fai$V2,">"),"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	writeLines(header,input)
	data.table::fwrite(data.frame(contig,z$POS[ix],paste0("row",ix),z$REF[ix],z$ALT[ix],".",".","."),input,sep="\t",append=TRUE,col.names=FALSE,quote=FALSE)
	rc <- system2(Sys.which("bcftools"),c("norm","--check-ref","e","-f",shQuote(fasta),"-m","-any","-Ov","-o",shQuote(output),shQuote(input)),stdout=file.path(tmp,"stdout"),stderr=file.path(tmp,"stderr"))
	if (rc!=0L) stop("Reference validation/normalization failed: ",paste(readLines(file.path(tmp,"stderr"),warn=FALSE),collapse="; "))
	r <- data.table::fread(output,skip="#CHROM",check.names=FALSE)
	ids <- as.integer(sub("^row","",r$ID))
	if (!setequal(ids,ix) || anyDuplicated(ids)) stop("Normalization lost or ambiguously split scalar associations")
	effect_alt <- z$EA[ids]==z$ALT[ids]
	z$CHR[ids] <- .norm_chr_value(r[["#CHROM"]]); z$POS[ids] <- r$POS
	z$REF[ids] <- r$REF; z$ALT[ids] <- r$ALT
	z$EA[ids] <- ifelse(effect_alt,r$ALT,r$REF); z$NEA[ids] <- ifelse(effect_alt,r$REF,r$ALT)
	z$variant_id[ids] <- paste(z$BUILD[ids],z$CHR[ids],format(z$POS[ids],scientific=FALSE,trim=TRUE),z$REF[ids],z$ALT[ids],sep=":")
	z$normalization_status[ids] <- "reference_checked_left_aligned_split"
	z$normalization_proof_hash[ids] <- le8_hash_object(list(reference=le8_file_sha256(fasta),tool=system2(Sys.which("bcftools"),"--version-only",stdout=TRUE),build=metadata$build))
	z$variant_key <- le8_variant_key(z)
	z
}
.le8_sumstat_qc <- new.env(parent=emptyenv()); .le8_sumstat_qc$rows <- list()
standardize_sumstat <- function(d, N_default = suppressWarnings(as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "NA"))),
			joint = FALSE, source_file = NA_character_, conflict_policy=c("error","exclude")) {
	conflict_policy <- match.arg(conflict_policy)
	d <- as.data.frame(d, check.names = FALSE) ; nms <- names(d)
	col <- list(
		SNP = .pick_col(nms, c("SNP", "RSID", "RS_NUMBER", "VARIANT_ID", "ID", "MARKERNAME", "POS_NAME")),
		CHR = .pick_col(nms, c("CHR", "CHROM", "CHROMOSOME")),
		POS = .pick_col(nms, c("POS", "BP", "POSITION", "BASE_PAIR_LOCATION")),
		EA = .pick_col(nms, c("EA", "A1", "ALT", "ALLELE1", "EFFECT_ALLELE", "EFFECTALLELE")),
		NEA = .pick_col(nms, c("NEA", "A2", "REF", "ALLELE0", "OTHER_ALLELE", "REFERENCE_ALLELE", "NONEFFECTALLELE")),
		REF_A = .pick_col(nms, c("REFA")),
		REF = .pick_col(nms, c("REF", "REFERENCE_ALLELE")),
		ALT = .pick_col(nms, c("ALT", "ALTERNATE_ALLELE")),
		EAF = .pick_col(nms, c("EAF", "FREQ", "FREQ_GENO", "A1FREQ", "EFFECT_ALLELE_FREQUENCY")),
		BETA = if (joint) .pick_col(nms, c("BJ", "BETAJ", "BETA", "B", "EFFECT"), prefer = c("BJ", "BETAJ")) else .pick_col(nms, c("BETA", "B", "EFFECT", "LOG_ODDS", "LOGOR", "BJ")),
		SE = if (joint) .pick_col(nms, c("BJ_SE", "SEJ", "SE", "STDERR", "SEBETA"), prefer = c("BJ_SE", "SEJ")) else .pick_col(nms, c("SE", "STDERR", "SEBETA", "BJ_SE")),
		P = if (joint) .pick_col(nms, c("PJ", "P_J", "P", "PVAL", "PVALUE"), prefer = c("PJ", "P_J")) else .pick_col(nms, c("P", "PVAL", "P_VALUE", "PVALUE", "PJ")),
		N = .pick_col(nms, c("N", "SAMPLESIZE", "N_TOTAL", "NEFF", "N_IIDS"))
	)
	getv <- function(nm, default = NA) if (!is.na(col[[nm]]) && col[[nm]] %in% names(d)) d[[col[[nm]]]] else rep(default, nrow(d))
	ea <- toupper(as.character(getv("EA"))) ; nea <- toupper(as.character(getv("NEA"))) ; refa <- toupper(as.character(getv("REF_A")))
	if (joint) {
		# GCTA .jma.cojo supplies refA but not both alleles. Keep refA as a provisional effect allele;
		# read_qtl_instruments() replaces it using the corresponding full QTL file.
		ea[is.na(ea) | ea == ""] <- refa[is.na(ea) | ea == ""]
	}
	ans <- tibble(
		SNP = as.character(getv("SNP")), CHR = .norm_chr_value(getv("CHR")),
		POS = suppressWarnings(as.numeric(getv("POS"))), EA = ea, NEA = nea,
		REF=toupper(as.character(getv("REF"))), ALT=toupper(as.character(getv("ALT"))),
		EAF = suppressWarnings(as.numeric(getv("EAF"))), BETA = suppressWarnings(as.numeric(getv("BETA"))),
		SE = suppressWarnings(as.numeric(getv("SE"))), P = suppressWarnings(as.numeric(getv("P"))),
		N = suppressWarnings(as.numeric(getv("N"))), source_file = source_file, joint = joint
	)
	meta_N <- le8_gwas_metadata(source_file)$N
	if (is.finite(meta_N) && meta_N > 0) ans$N[!is.finite(ans$N)] <- meta_N
	# Unknown N stays NA; no fabricated default sample size.
	ans$P[!is.finite(ans$P) | ans$P <= 0 | ans$P > 1] <- 2 * pnorm(abs(ans$BETA[!is.finite(ans$P) | ans$P <= 0 | ans$P > 1] /
		ans$SE[!is.finite(ans$P) | ans$P <= 0 | ans$P > 1]), lower.tail = FALSE)
	ans <- ans |> filter(!is.na(SNP), SNP!="",is.finite(BETA),is.finite(SE),SE>0) |> distinct()
	ans <- le8_variant_identity(ans,le8_gwas_metadata(source_file))
	ans$variant_key <- le8_variant_key(ans)
	dup <- unique(ans$variant_key[duplicated(ans$variant_key)])
	remove <- integer(); conflicting <- character(); redundant <- 0L
	for (key in dup) {
		ix <- which(ans$variant_key==key); z <- ans[ix,,drop=FALSE]; a <- z[1,,drop=FALSE]
		comp <- function(v) ifelse(nchar(v)==1,chartr("ACGT","TGCA",v),NA_character_)
		plus <- (z$EA==a$EA & z$NEA==a$NEA) | (comp(z$EA)==a$EA & comp(z$NEA)==a$NEA)
		minus <- (z$EA==a$NEA & z$NEA==a$EA) | (comp(z$EA)==a$NEA & comp(z$NEA)==a$EA)
		direction <- ifelse(plus %in% TRUE,1,ifelse(minus %in% TRUE,-1,NA_real_))
		eq <- function(v,ref) all(!is.finite(v) & !is.finite(ref) | is.finite(v) & is.finite(ref) & abs(v-ref)<=1e-8*pmax(1,abs(ref)))
		consistent <- eq(z$BETA*direction,a$BETA) && eq(z$SE,a$SE) && eq(z$N,a$N) && eq(ifelse(direction==1,z$EAF,1-z$EAF),a$EAF)
		if (consistent) { remove<-c(remove,ix[-1]); redundant<-redundant+length(ix)-1L } else { conflicting<-c(conflicting,key);remove<-c(remove,ix) }
	}
	if (length(conflicting) && conflict_policy=="error") stop("Conflicting duplicate coordinate/allele statistics: ",source_file)
	if (length(remove)) ans <- ans[-remove,,drop=FALSE]
	if (length(conflicting)) warning("Excluded ",length(conflicting)," ambiguous coordinate/allele identities from ",source_file,call.=FALSE)
	qc <- tibble(source_file,input_rows=nrow(d),retained_rows=nrow(ans),redundant_rows_collapsed=redundant,
		ambiguous_variant_identities=length(conflicting),policy=conflict_policy)
	attr(ans,"variant_qc") <- qc
	if (exists(".le8_sumstat_qc",inherits=TRUE)) .le8_sumstat_qc$rows[[length(.le8_sumstat_qc$rows)+1L]] <- qc
	ans
}
read_sumstat <- function(file, N_default = suppressWarnings(as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "NA"))), joint = grepl("jma\\.cojo$", file %||% "")) {
	if (is.na(file) || !file.exists(file) || file.size(file) == 0) return(tibble())
	standardize_sumstat(read_table_auto(file), N_default, joint = joint, source_file = file, conflict_policy="exclude")
}
read_sumstat_header <- function(path) {
	con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else base::file(path, "rt")
	on.exit(close(con), add = TRUE)
	readLines(con, n = 1, warn = FALSE)
}
read_sumstat_region <- function(file, chr, start, end, N_default = suppressWarnings(as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "NA")))) {
	if (is.na(file) || !file.exists(file) || file.size(file) == 0) return(tibble())
	hdr <- tryCatch(read_sumstat_header(file), error = function(e) "")
	d <- NULL
	nms <- if (nzchar(hdr)) strsplit(hdr, "\t", fixed = TRUE)[[1]] else character()

	# Coordinate-indexed files are queried in O(region) time. The ordinary-gzip
	# streaming scan remains as a compatibility fallback for legacy inputs.
	tabix <- Sys.which("tabix")
	index_files <- c(paste0(file, ".tbi"), paste0(file, ".csi"))
	file_mtime <- suppressWarnings(file.info(file)$mtime[[1]])
	index_info <- suppressWarnings(file.info(index_files))
	fresh_index <- rownames(index_info)[!is.na(index_info$size) & index_info$size > 0 &
		!is.na(index_info$mtime) & index_info$mtime >= file_mtime]
	if (length(nms) && nzchar(tabix) && length(fresh_index)) {
		contigs <- tryCatch(system2(tabix, c("-l", shQuote(file)), stdout = TRUE, stderr = FALSE), error = function(e) character())
		source_chr <- contigs[match(.norm_chr_value(chr), .norm_chr_value(contigs))]

		if (length(source_chr) && !is.na(source_chr[[1]])) {
			region <- sprintf("%s:%d-%d", source_chr[[1]], floor(start), ceiling(end))
			rows <- tryCatch(system2(tabix, c(shQuote(file), shQuote(region)), stdout = TRUE, stderr = FALSE), error = function(e) NULL)
			if (!is.null(rows)) {
				if (length(rows)) {
					d <- tryCatch(data.table::fread(
						text = paste(rows, collapse = "\n"), header = FALSE,
						col.names = nms, fill = TRUE, showProgress = FALSE,
						check.names = FALSE
					), error = function(e) NULL)
				} else {
					d <- data.table::as.data.table(matrix(nrow = 0, ncol = length(nms)))
					data.table::setnames(d, nms)
				}
			}
		}
	}

	if (is.null(d) && length(nms) && Sys.info()[["sysname"]] != "Windows" && nzchar(Sys.which("awk"))) {
		cchr <- .pick_col(nms, c("CHR", "CHROM", "CHROMOSOME")) ; cpos <- .pick_col(nms, c("POS", "BP", "POSITION", "BASE_PAIR_LOCATION"))
		ichr <- match(cchr, nms) ; ipos <- match(cpos, nms)
		if (is.finite(ichr) && is.finite(ipos)) {
			dec <- if (grepl("\\.gz$", file)) paste("gzip -cd", shQuote(file)) else paste("cat", shQuote(file))
			awk <- sprintf(
				"awk -F '\\t' 'NR==1 || (($%d==\"%s\" || $%d==\"chr%s\") && $%d>=%d && $%d<=%d)'",
				ichr, chr, ichr, chr, ipos, floor(start), ipos, ceiling(end)
			)
			d <- tryCatch(data.table::fread(cmd = paste(dec, "|", awk), fill = TRUE, showProgress = FALSE, check.names = FALSE), error = function(e) NULL)
		}
	}
	if (is.null(d)) d <- read_table_auto(file)
	z <- standardize_sumstat(d, N_default, joint = grepl("jma\\.cojo$", file), source_file = file, conflict_policy="exclude")
	z |> filter(.norm_chr_value(CHR) == .norm_chr_value(chr), POS >= start, POS <= end)
}

read_sumstat_snps <- function(file, snps, N_default = suppressWarnings(as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "NA")))) {
	snps <- unique(as.character(snps)) ; snps <- snps[!is.na(snps) & nzchar(snps)]
	if (!length(snps) || is.na(file) || !file.exists(file)) return(tibble())
	hdr <- tryCatch(read_sumstat_header(file), error = function(e) "")
	d <- NULL
	if (nzchar(hdr) && Sys.info()[["sysname"]] != "Windows" && nzchar(Sys.which("awk"))) {
		nms <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
		csnp <- .pick_col(nms, c("SNP", "RSID", "RS_NUMBER", "VARIANT_ID", "ID", "MARKERNAME", "POS_NAME")) ; isnp <- match(csnp, nms)
		if (is.finite(isnp)) {
			tf <- tempfile("le8_snps_") ; writeLines(snps, tf)
			dec <- if (grepl("\\.gz$", file)) paste("gzip -cd", shQuote(file)) else paste("cat", shQuote(file))
			awk <- sprintf("awk -F '\\t' 'NR==FNR {a[$1]=1; next} FNR==1 || ($%d in a)' %s -", isnp, shQuote(tf))
			d <- tryCatch(data.table::fread(cmd = paste(dec, "|", awk), fill = TRUE, showProgress = FALSE, check.names = FALSE), error = function(e) NULL)
			unlink(tf)
		}
	}
	if (is.null(d)) d <- read_table_auto(file)
	standardize_sumstat(d, N_default, joint = grepl("jma\\.cojo$", file), source_file = file, conflict_policy="exclude") |> filter(SNP %in% snps)
}

get_y_gwas_file <- function(trait = Y, required = TRUE) {
	p <- Sys.getenv("Y_GWAS", unset = Sys.getenv("CAD_GWAS", unset = ""))
	if (!nzchar(p)) p <- gwas_clean_file(dir.Y, trait)
	if (file.exists(p) && file.size(p) > 0) return(normalizePath(p, winslash = "/", mustWork = FALSE))
	if (required) stop("Missing outcome GWAS: ", p, call. = FALSE)
	NA_character_
}
recover_qtl_alleles <- function(iv, full_file) {
	if (!nrow(iv) || is.na(full_file) || !file.exists(full_file)) return(iv)
	full <- read_sumstat_snps(full_file, iv$SNP)
	if (!nrow(full)) return(iv)
	a <- full |> transmute(SNP,
		EA_full = toupper(EA), NEA_full = toupper(NEA),
		EAF_full = EAF, CHR_full = CHR, POS_full = POS,
		BUILD_full=BUILD,REF_full=REF,ALT_full=ALT,variant_id_full=variant_id,
		normalization_status_full=normalization_status,normalization_proof_hash_full=normalization_proof_hash
	)
	out <- iv |>
		mutate(EA_joint = toupper(EA), NEA_joint = toupper(NEA)) |>
		left_join(a, by = "SNP") |>
		mutate(
			same_to_full = EA_joint == EA_full & (is.na(NEA_joint) | NEA_joint == "" | NEA_joint == NEA_full),
			reverse_to_full = EA_joint == NEA_full & (is.na(NEA_joint) | NEA_joint == "" | NEA_joint == EA_full),
			allele_ok = same_to_full | reverse_to_full
		) |>
		filter(allele_ok) |>
		mutate(
			BETA = ifelse(reverse_to_full, - BETA, BETA),
			EAF = ifelse(reverse_to_full & is.finite(EAF), 1 - EAF, EAF),
			EA = EA_full, NEA = NEA_full, EAF = coalesce(EAF, EAF_full),
			# The joint file may be GRCh37 while the full QTL has already been
			# lifted to GRCh38. Once the stable ID and alleles agree, the full
			# QTL is authoritative for downstream cis/local coordinates.
			CHR = CHR_full, POS = POS_full, BUILD=BUILD_full, REF=REF_full, ALT=ALT_full,
			variant_id=variant_id_full,normalization_status=normalization_status_full,normalization_proof_hash=normalization_proof_hash_full
		) |>
		select( - ends_with("_full"), - EA_joint, - NEA_joint, - same_to_full, - reverse_to_full, - allele_ok)
	ambiguous <- out$SNP[duplicated(out$SNP)|duplicated(out$SNP,fromLast=TRUE)]
	if(length(ambiguous)) warning("Excluded ",length(unique(ambiguous))," ambiguous multiallelic COJO identifiers",call.=FALSE)
	out <- out[!out$SNP %in% ambiguous,,drop=FALSE]
	out$variant_key <- le8_variant_key(out);out
}
is_palindromic <- function(a1, a2) paste0(a1, a2) %in% c("AT", "TA", "CG", "GC")
allele_complement <- function(a) chartr("ATCG", "TAGC", toupper(a))
harmonize_sumstats <- function(x, y, drop_palindromic = TRUE) {
	x$variant_key <- le8_variant_key(x); y$variant_key <- le8_variant_key(y)
	if (anyDuplicated(x$variant_key) || anyDuplicated(y$variant_key)) stop("Conflicting duplicate variant identities before harmonization")
	x$rsid_x <- x$SNP; y$rsid_y <- y$SNP
	x$SNP <- x$variant_key; y$SNP <- y$variant_key
	valid <- function(z) {
		normalized <- if("normalization_status" %in% names(z)) z$normalization_status %in% "reference_checked_left_aligned_split" else rep(FALSE,nrow(z))
		ok <- !is.na(z$EA)&!is.na(z$NEA)&grepl("^[ACGT]+$",z$EA)&grepl("^[ACGT]+$",z$NEA)&z$EA!=z$NEA &
			((nchar(z$EA)==1L & nchar(z$NEA)==1L) | normalized)
		z[ok %in% TRUE,,drop=FALSE]
	}
	x <- valid(x); y <- valid(y)
	d <- inner_join(x |> select(SNP, CHR_x = CHR, POS_x = POS, EA_x = EA, NEA_x = NEA, EAF_x = EAF, BETA_x = BETA, SE_x = SE, P_x = P, N_x = N, everything()),
		y |> select(SNP, CHR_y = CHR, POS_y = POS, EA_y = EA, NEA_y = NEA, EAF_y = EAF, BETA_y = BETA, SE_y = SE, P_y = P, N_y = N, rsid_y, any_of(c("BUILD","REF","ALT","variant_id","normalization_status","normalization_proof_hash"))),
		by = "SNP", suffix = c("", ".dup")
	)
	if (!nrow(d)) return(d)
	d <- d |>
		mutate(
			EA_x = toupper(EA_x), NEA_x = toupper(NEA_x), EA_y = toupper(EA_y), NEA_y = toupper(NEA_y),
			EA_y_comp = allele_complement(EA_y), NEA_y_comp = allele_complement(NEA_y),
			same = EA_x == EA_y & NEA_x == NEA_y,
			flip = EA_x == NEA_y & NEA_x == EA_y,
			strand_same = nchar(EA_x)==1 & nchar(NEA_x)==1 & EA_x == EA_y_comp & NEA_x == NEA_y_comp,
			strand_flip = nchar(EA_x)==1 & nchar(NEA_x)==1 & EA_x == NEA_y_comp & NEA_x == EA_y_comp,
			pal = is_palindromic(EA_x, NEA_x),
			pal_same_diff = abs(EAF_x - EAF_y),
			pal_flip_diff = abs(EAF_x - (1 - EAF_y)),
			pal_freq_ok = pal & is.finite(EAF_x) & is.finite(EAF_y) &
				pmin(pal_same_diff, pal_flip_diff) <= .10 & abs(pal_same_diff - pal_flip_diff) >= .05,
			pal_reverse = pal_freq_ok & pal_flip_diff < pal_same_diff,
			nonpal_keep = !pal & (same | flip | strand_same | strand_flip),
			keep = if (drop_palindromic) nonpal_keep | pal_freq_ok else (same | flip | strand_same | strand_flip),
			reverse = ifelse(pal, pal_reverse, flip | strand_flip),
			BETA_y = ifelse(reverse, - BETA_y, BETA_y), EAF_y = ifelse(reverse, 1 - EAF_y, EAF_y),
			harmonization = case_when(
				pal & pal_freq_ok & !pal_reverse ~ "palindrome_frequency_same",
				pal & pal_freq_ok & pal_reverse ~ "palindrome_frequency_swapped",
				!pal & same ~ "same", !pal & flip ~ "swapped",
				!pal & strand_same ~ "strand", !pal & strand_flip ~ "strand_swapped",
				TRUE ~ "mismatch"
			)
		) |>
		filter(keep) |>
		distinct(SNP, .keep_all = TRUE) |>
		select(
 - EA_y_comp, - NEA_y_comp, - same, - flip, - strand_same, - strand_flip, - pal,
 - pal_same_diff, - pal_flip_diff, - pal_freq_ok, - pal_reverse, - nonpal_keep, - keep, - reverse
		)
	d
}
calc_iv_metrics <- function(iv) {
	if (!nrow(iv)) return(tibble(
		n_iv = 0L,
		r2_median = NA_real_, r2_q25 = NA_real_, r2_q75 = NA_real_, r2_p90 = NA_real_, r2_max = NA_real_,
		mean_F = NA_real_, median_F = NA_real_, min_F = NA_real_, max_F = NA_real_
	))
	# Partial R2 derived from each SNP's t statistic is bounded and does not depend on phenotype scaling.
	t2 <- (iv$BETA / iv$SE) ^ 2 ; df <- pmax(iv$N - 2, 1)
	r2_i <- pmin(pmax(t2 / (t2 + df), 0), .999999)
	F_i <- t2
	r2_ok <- r2_i[is.finite(r2_i)] ; f_ok <- F_i[is.finite(F_i)]
	# Report per-IV r2 only; aggregate approximations inflate correlated instruments.
	tibble(
		n_iv = nrow(iv),
		r2_median = if (length(r2_ok)) median(r2_ok) else NA_real_,
		r2_q25 = if (length(r2_ok)) as.numeric(quantile(r2_ok, .25, names = FALSE)) else NA_real_,
		r2_q75 = if (length(r2_ok)) as.numeric(quantile(r2_ok, .75, names = FALSE)) else NA_real_,
		r2_p90 = if (length(r2_ok)) as.numeric(quantile(r2_ok, .90, names = FALSE)) else NA_real_,
		r2_max = if (length(r2_ok)) max(r2_ok) else NA_real_,
		mean_F = if (length(f_ok)) mean(f_ok) else NA_real_,
		median_F = if (length(f_ok)) median(f_ok) else NA_real_,
		min_F = if (length(f_ok)) min(f_ok) else NA_real_,
		max_F = if (length(f_ok)) max(f_ok) else NA_real_
	)
}
le8_fit_mr <- function(iv, ygwas, exposure, analysis) {
	empty <- tibble(
		exposure = exposure, analysis = analysis, method = NA_character_, n_IV = 0L,
		b = NA_real_, se = NA_real_, pval = NA_real_, Q = NA_real_, Q_p = NA_real_,
		egger_intercept = NA_real_, egger_intercept_p = NA_real_,
		egger_slope = NA_real_, egger_slope_se = NA_real_, egger_slope_p = NA_real_,
		steiger_r2_exposure = NA_real_, steiger_r2_outcome = NA_real_,
		steiger_support_fraction = NA_real_, steiger_n = 0L,
		steiger_support = NA, tsmr_verified = FALSE,
		tsmr_ivw_b = NA_real_, tsmr_ivw_p = NA_real_,
		tsmr_weighted_median_b = NA_real_, tsmr_weighted_median_p = NA_real_
	)
	if (!nrow(iv)) return(bind_cols(empty, calc_iv_metrics(iv)[, - 1, drop = FALSE]))
	d <- harmonize_sumstats(iv, ygwas)
	if (!nrow(d)) return(bind_cols(empty, calc_iv_metrics(iv)[, - 1, drop = FALSE]))
	d <- d |> filter(
		is.finite(BETA_x), BETA_x != 0, is.finite(SE_x), SE_x > 0,
		is.finite(BETA_y), is.finite(SE_y), SE_y > 0
	) |>
		# Orient every instrument to an exposure-increasing allele. Ratios are
		# unchanged, while the MR-Egger intercept now has the intended meaning.
		mutate(
			exposure_reverse = BETA_x < 0, BETA_y = ifelse(exposure_reverse, - BETA_y, BETA_y), BETA_x = abs(BETA_x),
			ratio = BETA_y / BETA_x, ratio_se = abs(SE_y / BETA_x), w = 1 / ratio_se ^ 2
		)
	if (!nrow(d)) return(bind_cols(empty, calc_iv_metrics(iv)[, - 1, drop = FALSE]))
	if (nrow(d) == 1) {
		b <- d$ratio[[1]] ; se <- d$ratio_se[[1]] ; method <- "Wald ratio" ; Q <- Qp <- NA_real_
	} else {
		b <- sum(d$w * d$ratio) / sum(d$w) ; se_fixed <- sqrt(1 / sum(d$w))
		Q <- sum(d$w * (d$ratio - b) ^ 2) ; Qp <- pchisq(Q, df = nrow(d) - 1, lower.tail = FALSE)
		# Multiplicative random-effects IVW is the primary multi-variant estimate.
		# It equals fixed-effect IVW when Q/(K-1)<=1 and inflates the uncertainty
		# when the instruments are heterogeneous.
		phi <- max(1, Q / (nrow(d) - 1)) ; se <- se_fixed * sqrt(phi)
		method <- "IVW multiplicative random effects"
	}
	ei <- eip <- es <- ese <- esp <- NA_real_
	if (nrow(d) >= 3) {
		eg <- tryCatch(lm(BETA_y ~ BETA_x, weights = 1 / SE_y ^ 2, data = d), error = function(e) NULL)
		if (!is.null(eg)) {
			sm <- coef(summary(eg)) ; ei <- sm["(Intercept)", "Estimate"] ; eip <- sm["(Intercept)", "Pr(>|t|)"]
			es <- sm["BETA_x", "Estimate"] ; ese <- sm["BETA_x", "Std. Error"] ; esp <- sm["BETA_x", "Pr(>|t|)"]
		}
	}
	metrics <- calc_iv_metrics(iv[le8_variant_key(iv) %in% d$variant_key,,drop=FALSE])
	tx <- (d$BETA_x / d$SE_x) ^ 2 ; ty <- (d$BETA_y / d$SE_y) ^ 2
	rxi <- tx / (tx + pmax(d$N_x - 2, 1)) ; ryi <- ty / (ty + pmax(d$N_y - 2, 1))
	sok <- is.finite(rxi) & is.finite(ryi) ; sf <- if (any(sok)) mean(rxi[sok] > ryi[sok]) else NA_real_
	# Sums are exposed for audit only.  Directional support is based on the
	# fraction of harmonized IVs for which r2(exposure)>r2(outcome), avoiding the
	# product-to-one artefact in very large distal instrument sets.
	r2x <- if (any(sok)) sum(pmin(pmax(rxi[sok], 0), .999999)) else NA_real_
	r2y <- if (any(sok)) sum(pmin(pmax(ryi[sok], 0), .999999)) else NA_real_
	# Optional package-level verification.  Failure never changes the primary
	# estimate; it is exposed as an audit field because TwoSampleMR APIs differ
	# across input layouts.
	tsmr <- list(ok = FALSE, ivw_b = NA_real_, ivw_p = NA_real_, wm_b = NA_real_, wm_p = NA_real_)
	if (requireNamespace("TwoSampleMR", quietly = TRUE) && nrow(d) >= 2) {
		td <- d |> transmute(SNP,
			beta.exposure = BETA_x, se.exposure = SE_x,
			beta.outcome = BETA_y, se.outcome = SE_y, effect_allele.exposure = EA_x,
			other_allele.exposure = NEA_x, effect_allele.outcome = EA_y,
			other_allele.outcome = NEA_y, pval.exposure = P_x, pval.outcome = P_y,
			exposure = .env$exposure, outcome = Y, id.exposure = .env$exposure, id.outcome = Y, mr_keep = TRUE
		)
		zz <- tryCatch(quiet_package_call(
			TwoSampleMR::mr(td, method_list = c("mr_ivw_mre", "mr_weighted_median"))
		), error = function(e) NULL)
		if (!is.null(zz) && nrow(zz)) {
			ivw <- zz |>
				filter(str_detect(method, "Inverse variance weighted")) |>
				slice(1)
			wm <- zz |>
				filter(str_detect(method, "Weighted median")) |>
				slice(1)
			tsmr <- list(
				ok = TRUE, ivw_b = if (nrow(ivw)) ivw$b[[1]] else NA_real_, ivw_p = if (nrow(ivw)) ivw$pval[[1]] else NA_real_,
				wm_b = if (nrow(wm)) wm$b[[1]] else NA_real_, wm_p = if (nrow(wm)) wm$pval[[1]] else NA_real_
			)
		}
	}
	bind_cols(
		tibble(
			exposure = exposure, analysis = analysis, method = method, n_IV = nrow(d), b = b, se = se,
			pval = 2 * pnorm(abs(b / se), lower.tail = FALSE), Q = Q, Q_p = Qp,
			egger_intercept = ei, egger_intercept_p = eip,
			egger_slope = es, egger_slope_se = ese, egger_slope_p = esp,
			steiger_r2_exposure = r2x, steiger_r2_outcome = r2y,
			steiger_support_fraction = sf, steiger_n = sum(sok),
			steiger_support = if (is.finite(sf)) sf > .5 else NA,
			tsmr_verified = tsmr$ok, tsmr_ivw_b = tsmr$ivw_b, tsmr_ivw_p = tsmr$ivw_p,
			tsmr_weighted_median_b = tsmr$wm_b, tsmr_weighted_median_p = tsmr$wm_p
		),
		metrics |> select( - n_iv)
	)
}


# 🚩 Reusable plots
stable_neglog10_p <- function(p, statistic = NULL) {
	ans <-  - log10(pmax(suppressWarnings(as.numeric(p)), .Machine$double.xmin))
	if (!is.null(statistic)) {
		s <- suppressWarnings(as.numeric(statistic))
		lp <- log(2) + pnorm(abs(s), lower.tail = FALSE, log.p = TRUE)
		z <-  - lp / log(10)
		use <- is.finite(z) & (!is.finite(ans) | !is.finite(p) | p <= .Machine$double.xmin * 10)
		ans[use] <- z[use]
	}
	ans
}

# Compress only the extreme display tail while preserving original-scale labels.
compress_extreme_tail <- function(x, threshold = 50, tail_scale = 8) {
	x <- suppressWarnings(as.numeric(x))
	tail <- is.finite(x) & x > threshold
	x[tail] <- threshold + tail_scale * log10(1 + (x[tail] - threshold) / tail_scale)
	x
}
compressed_tail_scale <- function(x, threshold = 50, tail_scale = 8) {
	mx <- suppressWarnings(max(x, na.rm = TRUE))
	candidates <- c(0, 2, 5, 10, 20, 50, 100, 200, 300, 500, 1000, 3000)
	originals <- candidates[candidates <= max(mx, threshold)]
	if (!length(originals) || tail(originals, 1) < mx)
		originals <- c(originals, signif(mx, 2))
	list(
		breaks = compress_extreme_tail(originals, threshold, tail_scale),
		labels = format(originals, trim = TRUE, scientific = FALSE)
	)
}

plot_pwas_manhattan <- function(res, proteins, title = "Proteome-wide association") {
	bed <- read_prot_bed(proteins)
	d <- res |>
		mutate(beta = coalesce(beta, safe_log(estimate))) |>
		transmute(protein = term, beta, p.value) |>
		left_join(bed, by = "protein") |>
		filter(!is.na(chr), is.finite(pos), is.finite(p.value)) |>
		mutate(chr_num = suppressWarnings(as.numeric(chr))) |>
		arrange(chr_num, pos) |>
		mutate(index = row_number(), sig = p.value < .05 / n(), y =  - log10(pmax(p.value, 1e-300)))
	labs <- d |>
		arrange(p.value) |>
		slice_head(n = 25)
	ggplot(d, aes(index, y)) +
		geom_hline(yintercept =  - log10(.05 / max(1, nrow(d))), linetype = 2, color = "#B2182B") +
		geom_point(aes(color = sig), alpha = .85, size = 1.5) +
		ggrepel::geom_text_repel(data = labs, aes(label = protein), size = 2.6, seed = 1, max.overlaps = Inf) +
		scale_color_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "grey78"), guide = "none") +
		labs(title = title, x = "Protein genomic order", y = expression( - log[10](P))) +
		theme_5c(12)
}
plot_volcano <- function(res, xvar = "beta", label_col = "term", title = "Volcano", pth = NULL,
			top_n = 25, x_quantile = NULL) {
	if (is.null(pth)) pth <- .05 / max(1, nrow(res))
	stat <- if ("statistic" %in% names(res)) res$statistic else NULL
	d <- res |> mutate(
		x_original = .data[[xvar]], y_original = stable_neglog10_p(p.value, stat),
		direction = case_when(p.value < pth & x_original > 0 ~ "Positive", p.value < pth & x_original < 0 ~ "Inverse", TRUE ~ "NS"),
		label = ifelse(p.value < pth & min_rank(p.value) <= top_n, .data[[label_col]], NA_character_)
	)
	clipped <- FALSE
	if (!is.null(x_quantile)) {
		x_cap <- suppressWarnings(as.numeric(quantile(abs(d$x_original), x_quantile, na.rm = TRUE, names = FALSE)))
		if (is.finite(x_cap) && x_cap > 0) {
			clipped <- any(abs(d$x_original) > x_cap, na.rm = TRUE)
			d$x <- pmax( - x_cap, pmin(x_cap, d$x_original))
		} else d$x <- d$x_original
	} else d$x <- d$x_original
	use_compression <- any(d$y_original > 150, na.rm = TRUE)
	threshold <- if (use_compression) 50 else Inf
	d$y <- if (use_compression) compress_extreme_tail(d$y_original, threshold) else d$y_original
	sc <- if (use_compression) compressed_tail_scale(d$y_original, threshold) else NULL
	p <- ggplot(d, aes(x, y)) +
		geom_hline(yintercept = if (use_compression) compress_extreme_tail( - log10(pth), threshold) else - log10(pth), linetype = 2, color = "grey45") +
		geom_vline(xintercept = 0, color = "grey65") +
		geom_point(aes(color = direction), alpha = .82, size = 1.8) +
		ggrepel::geom_text_repel(aes(label = label), size = 2.7, seed = 2, max.overlaps = Inf, na.rm = TRUE) +
		scale_color_manual(values = c(Positive = "#D7301F", Inverse = "#2C7FB8", NS = "grey82")) +
		labs(
			title = title, subtitle = if (clipped) paste0("Effect axis winsorized at the ", 100 * x_quantile, "th percentile; labels are Bonferroni-significant") else NULL,
			x = "Effect estimate",
			y = if (use_compression) expression( - log[10](P) ~ "(compressed above 50)") else expression( - log[10](P)),
			color = NULL
		) +
		theme_5c(12) +
		theme(legend.position = "bottom")
	if (use_compression) p <- p + scale_y_continuous(breaks = sc$breaks, labels = sc$labels)
	p
}
# Category-wise metabolite radial plotting is defined in 0.common.R below.

# Find PRS columns without reading a second file. Preference is CAD, then BMI, then other score_sum variables.
find_prs_vars <- function(dat, outcome = Y, max_n = 4) {
	z <- grep("score_sum$|\\.prs$|_prs$", names(dat), value = TRUE, ignore.case = TRUE)
	if (!length(z)) return(character())
	key <- c(outcome, str_remove(outcome, "^cvd_"), "cad", "bmi")
	ord <- order(!vapply(z, function(v) any(str_detect(tolower(v), fixed(tolower(key)))), logical(1)), z)
	head(z[ord], max_n)
}


# Final shared plotting and analysis helpers.
.le8_source_dir <- Sys.getenv("LE8_FDIR", unset = "")
if (!nzchar(.le8_source_dir)) {
	.le8_source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
	.le8_source_dir <- if (!is.null(.le8_source_file)) dirname(.le8_source_file) else file.path(dir0, "scripts/le8/f")
}
Sys.setenv(LE8_FDIR = .le8_source_dir)

# c1 met circle
# Metabolite associations by biochemical class, without genomic coordinates.
# Inspired by the category-wise radial bars in Li et al., Nature Medicine 2026,
# Fig2 (s41591-025-04105-8); all plotted observations come from the current run.
le8_met_circle_data <- function(d) {
	d <- d |> filter(is.finite(p.value), p.value >= 0, p.value <= 1)
	# Keep the derived columns even when every association is unavailable.
	if (!nrow(d)) return(d |> mutate(
		neglog10_FDR = numeric(), FDR = numeric(),
		circle_group = character(), significant = logical()
	))
	logp <- log(pmax(d$p.value, .Machine$double.xmin))
	if ('statistic' %in% names(d)) {
		recover <- d$p.value == 0 & is.finite(d$statistic)
		logp[recover] <- log(2) + pnorm(abs(d$statistic[recover]), lower.tail = FALSE, log.p = TRUE)
	}
	# BH family includes all tested metabolites, before selecting significant bars.
	ix <- order(logp) ; m <- length(logp)
	logq <- numeric(m) ; logq[ix] <- pmin(0, rev(cummin(rev(logp[ix] + log(m / seq_len(m))))))
	d$neglog10_FDR <-  - logq / log(10) ; d$FDR <- exp(logq)
	if (!'super_group' %in% names(d)) d$super_group <- met_super_group(d$group)
	d |> mutate(
		circle_group = coalesce(super_group, 'Other metabolites'),
		significant = is.finite(beta) & neglog10_FDR >  - log10(.05)
	)
}

plot_met_circle <- function(d, title, label_n = Inf, fdr_cap = 50, beta_limit = .3, show_legend = TRUE) {
	all <- le8_met_circle_data(d)
	if (!nrow(all)) return(blank_plot(title, 'No valid metabolite P values available'))
	d <- all |> filter(significant)
	if (!nrow(d)) return(blank_plot(title, 'No metabolites with BH FDR < 0.05'))
	stopifnot(is.finite(fdr_cap), fdr_cap > 0, is.finite(beta_limit), beta_limit > 0)
	group_order <- c(
		'Lipoprotein lipids', 'Cholesterol and apolipoproteins', 'Fatty acids',
		'Amino acids', 'Energy metabolism', 'Other metabolites'
	)
	group_names <- c(
		'Lipoprotein lipids' = 'Lipoprotein lipids',
		'Cholesterol and apolipoproteins' = 'Cholesterol / Apo', 'Fatty acids' = 'Fatty acids',
		'Amino acids' = 'Amino acids', 'Energy metabolism' = 'Energy', 'Other metabolites' = 'Other'
	)
	d <- d |> arrange(match(circle_group, group_order), term)
	# Blank angular slots separate categories and leave room for a radial scale.
	parts <- split(d, factor(d$circle_group, levels = unique(d$circle_group)))
	offset <- 6L
	for (i in seq_along(parts)) {
		parts[[i]]$id <- offset + seq_len(nrow(parts[[i]])) ; offset <- max(parts[[i]]$id) + 3L
	}
	z <- bind_rows(parts) ; slots <- offset + 3L
	z <- z |> mutate(
		angle0 = 90 - 360 * (id - .5) / slots, hjust = ifelse(angle0 <  - 90, 1, 0),
		angle = ifelse(angle0 <  - 90, angle0 + 180, angle0), height = pmin(neglog10_FDR, fdr_cap),
		lab = ifelse(rank( - neglog10_FDR, ties.method = 'first') <= label_n, sub('^met_', '', term), NA_character_)
	)
	bands <- z |>
		group_by(circle_group) |>
		summarise(x1 = min(id) - .45, x2 = max(id) + .45, .groups = 'drop') |>
		mutate(
			x = (x1 + x2) / 2, angle = (180 - 360 * (x - .5) / slots + 90) %% 180 - 90,
			label = unname(group_names[circle_group])
		)
	ticks <- data.frame(y = c(0, fdr_cap / 2, fdr_cap), label = c('0', format(fdr_cap / 2), format(fdr_cap)))
	p <- ggplot(z, aes(id, height)) +
		geom_col(aes(fill = beta), width = .88, na.rm = TRUE) +
		geom_segment(
			data = bands, aes(x = x1, xend = x2, y =  - fdr_cap * .02, yend =  - fdr_cap * .02),
			inherit.aes = FALSE, linewidth = .5, color = 'grey35'
		) +
		geom_text(
			data = bands, aes(x = x, y =  - fdr_cap * .15, label = label, angle = angle),
			inherit.aes = FALSE, size = 2.3, color = 'grey20'
		) +
		geom_text(aes(y = height + fdr_cap * .025, label = lab, angle = angle, hjust = hjust),
			size = 1.85,
			color = 'grey35', na.rm = TRUE
		) +
		annotate('segment', x = 3, xend = 3, y = 0, yend = fdr_cap, linewidth = .3, color = 'grey60') +
		geom_segment(data = ticks, aes(x = 2, xend = 4, y = y, yend = y), inherit.aes = FALSE, linewidth = .3, color = 'grey55') +
		geom_text(data = ticks, aes(x = 1, y = y, label = label), inherit.aes = FALSE, hjust = 1, size = 2.2, color = 'grey40') +
		annotate('text', x = 3, y = fdr_cap * 1.10, label = '-log10(FDR)', size = 2.8, hjust = .5) +
		scale_x_continuous(limits = c(.5, slots + .5), expand = c(0, 0), breaks = NULL) +
		scale_y_continuous(limits = c( - fdr_cap * .5, fdr_cap * 1.35), expand = c(0, 0), breaks = NULL) +
		scale_fill_gradient2(
			low = '#0878BF', mid = 'white', high = '#FF5A52', midpoint = 0,
			limits = c( - beta_limit, beta_limit), oob = scales::squish, name = 'Log effect per SD',
			breaks = c( - beta_limit, 0, beta_limit)
		) +
		coord_polar(clip = 'off') +
		labs(
			title = title, x = NULL, y = NULL,
			subtitle = paste0(
				nrow(d), ' / ', nrow(all), ' metabolites: FDR < 0.05; ', sum(z$neglog10_FDR > fdr_cap),
				' bars capped at ', fdr_cap
			)
		) +
		theme_void(11) +
		theme(
			plot.title = element_text(face = 'bold', hjust = .5), legend.position = if (show_legend) 'bottom' else 'none',
			axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
			plot.margin = margin(12, 15, 12, 15)
		)
	attr(p, 'le8_circle_data') <- all
	p
}


# 0.common.R
# Explicit, trait-independent settings. Disk phenotype data remain unchanged.
le8_custom_covars <- unique(Filter(nzchar, trimws(strsplit(Sys.getenv('LE8_VARS_ADJ', ''), '[,[:space:]]+')[[1]])))
le8_custom_adjustment <- function() length(le8_custom_covars) > 0L
le8_check_options <- function(obj) {
	current <- le8_analysis_options() ; old <- obj$meta$analysis_options
	default <- list(y_date = paste0('fod_icd10_', Y), vars_adj = character(), white_only = 'TRUE')
	if (is.null(old)) old <- default
	if (!identical(old, current)) stop('Existing results use different analysis options. Use --replace TRUE or a new output directory; changing options must not silently relabel cached estimates.')
}
le8_select_phenotypes <- function(x) {
	source <- le8_y_date() ; required <- unique(c(source, le8_custom_covars))
	if (length(setdiff(required, names(x)))) stop('Unknown --Y-date/--vars.adj fields: ', paste(setdiff(required, names(x)), collapse = ', '))
	# Alias only inside this analysis, so all existing consumers use the selected date.
	x[[paste0('fod_icd10_', Y)]] <- x[[source]]
	x
}
if (le8_custom_adjustment()) {
	covs_use <- le8_custom_covars ; covs_use_name <- 'adj2'
}

le8_final_request <- function() Sys.getenv('LE8_FINAL_REQUEST', '')
le8_check_final_request <- function(obj) {
	request <- le8_final_request()
	if (nzchar(request) && !identical(obj$training_request, request))
		stop('Saved Final training settings differ or are unavailable. Use --replace TRUE or a new output directory to fit the requested analysis.')
	invisible(NULL)
}

# Cache paths stay stable across disposable execution workspaces.
le8_cache_root <- function() {
	project <- Sys.getenv('LE8_PUBLISHED_ROOT', Sys.getenv('LE8_ANALYSIS_ROOT', '/mnt/d/analysis/le8'))
	file.path('/tmp/le8-cache', substr(digest::digest(project, algo = 'sha256', serialize = FALSE), 1, 12))
}
le8_cache_dir <- function(stage, layer, ...) {
	layer <- if (layer %in% c('protein', 'prot')) 'prot' else 'met'
	file.path(le8_cache_root(), Y, layer, stage, ...)
}

le8_guard_stage_options <- function(layer, module) {
	outdir <- if (layer == 'protein') out.prot else out.met
	d <- le8_job_dir(outdir, module) ; dir.create(d, recursive = TRUE, showWarnings = FALSE)
	name <- if (module == 'final_prediction') 'res.rds' else if (module == 'c4_panel_validation') 'c4.validation.rds' else paste0(substr(module, 1, 2), '.res.rds')
	path <- file.path(d, name)
	if (file.exists(path) && !LE8_REPLACE) le8_check_options(readRDS(path))
	invisible(d)
}


# 0.common.R
# Reconstruct time-safe LE8 variables in memory, before ANY phenotype selection.
# Unsuffixed source assays follow ukb/f/phe.R's baseline convention. Custom
# datasets must supply LE8_BASELINE_MAP (canonical,column,unit; see runbook).
LE8_BASELINE_VERSION <- "2026-09-21.baseline-med-source-v3"
.le8_med_cache <- new.env(parent = emptyenv())
le8_attach_baseline_medication <- function(x) {
	# Raw i0 responses retain -7 (none of the above), which phe.R may have
	# converted to NA. Never infer a negative response from that processed NA.
	root <- get0("indir", ifnotfound = Sys.getenv("UKB_PHE", "/mnt/d/data/ukb/phe"))
	path <- Sys.getenv("LE8_BASELINE_MED_FILE", file.path(root, "rap/vip.tab.gz"))
	if (!file.exists(path)) return(x)
	explicit <- nzchar(Sys.getenv("LE8_BASELINE_MED_FILE", ""))
	key <- paste(path, file.info(path)$size, as.numeric(file.info(path)$mtime), sep = "|")
	if (!identical(.le8_med_cache$key, key)) {
		hdr <- if (grepl("[.]rds$", path, ignore.case = TRUE)) names(readRDS(path)) else names(data.table::fread(path, nrows = 0, showProgress = FALSE))
		cols <- grep("^(p)?(6153|6177)(_i0(_a[0-9]+)?|[.]0[.][0-9]+|-0[.][0-9]+)$", hdr, value = TRUE)
		ids <- intersect(c("eid", "f.eid", "IID"), hdr)
		if (!length(cols) || !length(ids)) {
			if (explicit) stop("LE8_BASELINE_MED_FILE needs eid and explicit baseline p6153/p6177 fields")
			return(x)
		}
		med <- if (grepl("[.]rds$", path, ignore.case = TRUE)) as.data.frame(readRDS(path))[, c(ids[1], cols), drop = FALSE] else
			as.data.frame(data.table::fread(path, select = c(ids[1], cols), showProgress = FALSE))
		names(med)[1] <- "eid" ; med$eid <- as.character(med$eid)
		if (anyNA(med$eid) || anyDuplicated(med$eid)) stop("Raw baseline medication input has missing/duplicate eid")
		names(med)[ - 1] <- paste0(".le8_raw_med_i0_", seq_along(cols))
		.le8_med_cache$data <- med ; .le8_med_cache$key <- key ; .le8_med_cache$source <- paste(path, paste(cols, collapse = ","), sep = " : ")
	}
	med <- .le8_med_cache$data ; i <- match(as.character(x$eid), med$eid)
	for (v in setdiff(names(med), "eid")) x[[v]] <- med[[v]][i]
	attr(x, "baseline_med_source") <- .le8_med_cache$source ; x
}
le8_rebuild_baseline <- function(dat, audit_dir = NULL) {
	x <- le8_attach_baseline_medication(as.data.frame(dat)) ; n <- nrow(x)
	getnum <- function(v) if (v %in% names(x)) suppressWarnings(as.numeric(as.character(x[[v]]))) else rep(NA_real_, n)
	dt <- function(z) if (inherits(z, "Date")) as.Date(z) else if (is.numeric(z)) as.Date(z, origin = "1970-01-01") else as.Date(as.character(z))
	if (!"date_attend" %in% names(x)) stop("Baseline date_attend is required")
	baseline <- dt(x$date_attend)
	old_names <- intersect(c("drug.lipid", "drug.dm", "drug.htn", "dm.yes", "htn.yes", "nonhdl.pts", "hba1c.pts", "bp.pts"), names(x))
	old <- x[, old_names, drop = FALSE]
	mapfile <- Sys.getenv("LE8_BASELINE_MAP", "")
	if (nzchar(mapfile)) {
		m <- read.csv(mapfile, stringsAsFactors = FALSE)
		if (!all(c("canonical", "column", "unit") %in% names(m)) || anyDuplicated(m$canonical)) stop("Invalid baseline map")
		if (any(!m$column %in% names(x))) stop("Baseline mapping refers to missing columns")
		expected <- c(bmi = "kg/m2", bb_TC = "mmol/L", bb_HDL = "mmol/L", bb_HBA1C = "mmol/mol", sbp = "mmHg", dbp = "mmHg")
		if (any(!m$canonical %in% names(expected)) || any(m$unit != expected[m$canonical])) stop("Convert baseline units before mapping")
		for (j in seq_len(nrow(m))) x[[m$canonical[j]]] <- x[[m$column[j]]]
	}
	# Only explicitly tagged first-visit responses, or explicitly supplied columns.
	medcols <- trimws(strsplit(Sys.getenv("LE8_BASELINE_MED_COLUMNS", ""), ",", fixed = TRUE)[[1]])
	medcols <- medcols[nzchar(medcols)]
	if (!length(medcols)) medcols <- grep("^([.]le8_raw_med_i0_|drug[.]big3.*[_.]i0([_.]|$)|(p)?(6153|6177)_i0(_a[0-9]+)?$)", names(x), value = TRUE)
	if (any(!medcols %in% names(x))) stop("Unknown baseline medication column")
	meds <- matrix(NA_integer_, n, 3, dimnames = list(NULL, c("drug.lipid", "drug.htn", "drug.dm")))
	if (length(medcols)) {
		patterns <- do.call(paste, c(lapply(x[medcols], as.character), sep = " "))
		keys <- unique(patterns)
		lookup <- t(vapply(keys, function(s) {
			tokens <- suppressWarnings(as.integer(regmatches(s, gregexpr("-?[0-9]+", s))[[1]]))
			known <- any(tokens %in% c( - 7L, 1 : 5)) && !any(tokens %in% c( - 1L, - 3L)) &&
				!(any(tokens ==  - 7L, na.rm = TRUE) && any(tokens %in% 1 : 5))
			if (known) as.integer(1 : 3 %in% tokens) else rep(NA_integer_, 3)
		}, integer(3)))
		meds[, ] <- lookup[match(patterns, keys), , drop = FALSE]
	}
	for (v in colnames(meds)) x[[v]] <- meds[, v]
	baseline_diagnosis <- function(suffix) {
		cc <- grep(paste0("^fod_(srd|icd10|ref)_", suffix, "$"), names(x), value = TRUE)
		hit <- rep(FALSE, n)
		for (v in cc) {
			d <- dt(x[[v]]) ; hit <- hit | (!is.na(d) & !is.na(baseline) & d <= baseline)
		}
		# No recorded diagnosis is not proof of absence; this indicator describes records.
		hit
	}
	x$nonhdl <- (getnum("bb_TC") - getnum("bb_HDL")) * 38.67
	x$nonhdl[!is.finite(x$nonhdl) | x$nonhdl < 0] <- NA_real_
	x$hba1c_ngsp <- getnum("bb_HBA1C") * .0915 + 2.15
	x$hba1c_ngsp[getnum("bb_HBA1C") <= 0] <- NA_real_
	for (v in c("sbp", "dbp")) {
		cc <- grep(paste0("^", v, "_.*_i0([_.]|$)"), names(x), value = TRUE)
		if (length(cc)) {
			z <- as.matrix(x[, cc, drop = FALSE]) ; storage.mode(z) <- "double" ; z[z <= 0] <- NA_real_
			x[[v]] <- rowMeans(z, na.rm = TRUE)
		}
		if (!v %in% names(x)) x[[v]] <- rep(NA_real_, n)
		x[[v]][!is.finite(x[[v]]) | x[[v]] <= 0] <- NA_real_
	}
	x$dm.yes <- ifelse(baseline_diagnosis("t2dm") | x$drug.dm == 1 | x$hba1c_ngsp >= 6.5, 1, 0)
	# Unknown medication + no positive evidence remains unknown.
	x$htn.yes <- ifelse(baseline_diagnosis("cvd_htn") | x$drug.htn == 1, 1, 0)
	cutscore <- function(v, br, sc) as.numeric(as.character(cut(v, br, labels = sc, right = FALSE)))
	x$bmi.pts <- cutscore(getnum("bmi"), c(0, 25, 30, 35, 40, Inf), c(100, 70, 30, 15, 0))
	z <- cutscore(x$nonhdl, c(0, 130, 160, 190, 220, Inf), c(100, 60, 40, 20, 0))
	x$nonhdl.pts <- ifelse(z == 0, 0, pmax(z - 20 * x$drug.lipid, 0))
	h <- x$hba1c_ngsp
	x$hba1c.pts <- ifelse(x$dm.yes == 0 & h < 5.7, 100, ifelse(x$dm.yes == 0 & h < 6.5, 60,
		ifelse(h >= 6.5 | x$dm.yes == 1, cutscore(h, c(0, 7, 8, 9, 10, Inf), c(40, 30, 20, 10, 0)), NA_real_)
	))
	s <- x$sbp ; d <- x$dbp
	z <- ifelse(s >= 160 | d >= 100, 0, ifelse(s >= 140 | d >= 90, 25, ifelse(s >= 130 | d >= 80, 50, ifelse(s >= 120, 75, 100))))
	x$bp.pts <- ifelse(z == 0, 0, pmax(z - 20 * x$drug.htn, 0))
	x$.le8_baseline_version <- LE8_BASELINE_VERSION
	audit <- do.call(rbind, lapply(old_names, function(v) {
		a <- as.character(old[[v]]) ; b <- as.character(x[[v]])
		data.frame(
			variable = v, N = n, changed = sum((is.na(a) != is.na(b)) | (!is.na(a) & !is.na(b) & a != b)),
			old_missing = sum(is.na(a)), new_missing = sum(is.na(b)), version = LE8_BASELINE_VERSION
		)
	}))
	if (is.null(audit)) audit <- data.frame(variable = character(), N = integer())
	attr(x, "baseline_audit") <- audit
	if (!is.null(audit_dir)) {
		dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
		dir.create(file.path(analysis_root, "logs", Y), recursive = TRUE, showWarnings = FALSE)
		write.csv(audit, file.path(analysis_root, "logs", Y, "baseline_rebuild_audit.csv"), row.names = FALSE)
		writeLines(c(
			LE8_BASELINE_VERSION, paste("baseline medication columns:", paste(medcols, collapse = ",")),
			paste("raw medication source:", if (is.null(attr(x, "baseline_med_source"))) "not available" else attr(x, "baseline_med_source")),
			paste("mapping:", if (nzchar(mapfile)) mapfile else "ukb/f/phe.R baseline raw-assay convention"),
			"Diagnosis requires date <= date_attend; unresolved medication status stays NA.",
			"Raw assays and other LE8 components still require source visit provenance review."
		), file.path(analysis_root, "logs", Y, "baseline_provenance.txt"))
	}
	x
}


# 0.common.R
# Publication-style companions. Always plot this run's results, never paper data.
le8_mock_volcano <- function(a, layer = 'protein') {
	z <- a |>
		filter(is.finite(beta), is.finite(p.value)) |>
		mutate(HR = exp(beta), significance = case_when(p.value < .05 / nrow(a) & beta > 0 ~ 'Risk', p.value < .05 / nrow(a) & beta < 0 ~ 'Protective', TRUE ~ 'Not significant'))
	lab <- z |>
		arrange(p.value) |>
		slice_head(n = 15)
	ggplot(z, aes(HR, - log10(pmax(p.value, 1e-300)), color = significance)) +
		geom_point(size = 1.5, alpha = .7) +
		geom_hline(yintercept =  - log10(.05 / nrow(a)), linetype = 2) +
		geom_vline(xintercept = 1, color = 'grey75') +
		ggrepel::geom_text_repel(data = lab, aes(label = term), seed = SEED, max.overlaps = Inf, size = 3) +
		scale_color_manual(values = c(Risk = '#CE4965', Protective = '#7463AB', `Not significant` = 'grey80')) +
		labs(
			title = paste('A. Plasma', if (layer == 'protein') 'proteins' else 'metabolites', 'and incident', Y),
			x = paste('Hazard ratio per 1-SD', if (layer == 'protein') 'protein' else 'metabolite'), y = expression( - log[10](P)), color = NULL
		) +
		theme_5c(10)
}
le8_mock_trajectories <- function(dat, a, covars, tvar, evar, layer = 'protein') {
	features <- a |>
		filter(is.finite(p.value), p.value < .05 / nrow(a)) |>
		arrange(p.value) |>
		pull(term) |>
		intersect(names(dat))
	empty <- list(trajectories = tibble(), clusters = tibble(), status = tibble(status = 'No significant biomarkers or insufficient matched observations', layer = layer))
	if (!length(features)) return(empty)
	if (!requireNamespace('MatchIt', quietly = TRUE)) {
		empty$status$status <- 'Install MatchIt for matched-control trajectories' ; return(empty)
	}
	z <- as.data.frame(dat[, unique(c('eid', covars, tvar, evar, features)), drop = FALSE])
	z <- z[complete.cases(z[, c(covars, tvar, evar), drop = FALSE]) & z[[tvar]] > 0, , drop = FALSE]
	z$.case <- z[[evar]] ; if (sum(z$.case == 1) < 20 || sum(z$.case == 0) < 20) return(empty)
	set.seed(SEED)
	match_covars <- covars[vapply(z[covars], function(v) length(unique(v)) > 1L, logical(1))]
	form <- if (length(match_covars)) reformulate(match_covars, response = '.case') else as.formula('.case ~ 1')
	m <- tryCatch(MatchIt::matchit(form, z, method = 'nearest', ratio = max(1L, min(10L, floor(sum(z$.case == 0) / sum(z$.case == 1)))), replace = FALSE), error = identity)
	if (inherits(m, 'error')) {
		empty$status$status <- paste('Matching unavailable:', conditionMessage(m)) ; return(empty)
	}
	md <- MatchIt::match.data(m) ; ca <- md[md$.case == 1, , drop = FALSE] ; co <- md[md$.case == 0, , drop = FALSE]
	grid <- seq( - 15, 0, by = .1)
	tr <- bind_rows(lapply(features, function(x) {
		mu <- mean(co[[x]], na.rm = TRUE) ; s <- sd(co[[x]], na.rm = TRUE)
		if (!is.finite(s) || s <= 0) return(tibble())
		d <- data.frame(years =  - ca[[tvar]], z = (ca[[x]] - mu) / s) ; d <- d[complete.cases(d) & d$years >=  - 15, , drop = FALSE]
		if (nrow(d) < 20 || length(unique(d$years)) < 5) return(tibble())
		# Only fitted means are used; omit unused standard-error calculations.
		fit <- tryCatch(loess(z ~ years, d, span = .75, control = loess.control(statistics = 'none')), error = function(e) NULL) ; if (is.null(fit)) return(tibble())
		tibble(feature = x, years = grid, z = as.numeric(predict(fit, data.frame(years = grid))))
	}))
	if (!nrow(tr)) return(empty)
	eligible <- tr |>
		group_by(feature) |>
		summarise(cross = any(abs(z) > .25, na.rm = TRUE), .groups = 'drop') |>
		filter(cross) |>
		pull(feature)
	mat <- tr |>
		filter(feature %in% eligible) |>
		pivot_wider(names_from = years, values_from = z)
	cl <- tibble(feature = character(), cluster = integer())
	if (nrow(mat) >= 3) {
		v <- as.matrix(mat[, - 1]) ; ok <- colSums(is.finite(v)) == nrow(v)
		if (sum(ok) > 5) {
			hc <- hclust(dist(v[, ok, drop = FALSE]), method = 'ward.D2') ; cl <- tibble(feature = mat$feature, cluster = unname(cutree(hc, k = 3)))
		}
	}
	list(trajectories = tr, clusters = cl, status = tibble(
		status = 'ok', layer = layer, matching_covariates = paste(match_covars, collapse = ','), cases = nrow(ca), controls = nrow(co), biomarkers = length(features),
		proteins = if (layer == 'protein') length(features) else NA_integer_, metabolites = if (layer == 'metabolite') length(features) else NA_integer_,
		interpretation = 'One baseline sample per person; diagnosis-timed LOESS, not longitudinal within-person change', reference = 'Mean/SD of matched controls; absolute Z > 0.25 for clustering'
	))
}
le8_plot_mock2 <- function(z, outdir, layer = 'protein') {
	tr <- z$trajectories ; cl <- z$clusters
	omic <- if (layer == 'protein') 'protein' else 'metabolite'
	write_xlsx2(z, le8_artifact_path('c1.Fig7.profiles.xlsx', outdir))
	if (!nrow(tr)) {
		save_plot(blank_plot(paste('Diagnosis-timed', omic, 'trajectories'), z$status$status[1]), 'c1.Fig7.diagnosis_timed_profiles.png', 14, 9, outdir = outdir) ; return()
	}
	ord <- if (nrow(cl)) c(cl$feature[order(cl$cluster)], setdiff(unique(tr$feature), cl$feature)) else unique(tr$feature)
	h <- tr |> mutate(feature = factor(feature, levels = rev(ord)), display = ifelse(abs(z) > .25, z, 0))
	p <- ggplot(h, aes(years, feature, fill = display)) +
		geom_tile() +
		scale_fill_gradient2(low = '#3878B9', mid = 'white', high = '#CF4050', midpoint = 0, na.value = 'grey90') +
		guides(fill = guide_colourbar(barwidth = grid::unit(6, 'cm'), barheight = grid::unit(.3, 'cm'), title.position = 'top')) +
		labs(title = paste0('A. ', tools::toTitleCase(omic), ' levels before diagnosis'), x = 'Years before diagnosis', y = NULL, fill = 'Z score') +
		theme_5c(8) +
		theme(panel.grid = element_blank())
	if (length(ord) > 150) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(y = paste(length(ord), paste0(omic, 's; row names in CSV')))
	lines <- lapply(1 : 3, function(k) {
		d <- inner_join(tr, cl, by = 'feature') |> filter(cluster == k)
		if (!nrow(d)) return(blank_plot(paste('Cluster', k), 'Fewer than three estimable clusters'))
		avg <- d |>
			group_by(years) |>
			summarise(z = mean(z, na.rm = TRUE), .groups = 'drop')
		ggplot(d, aes(years, z, group = feature)) +
			geom_line(color = c('#4C8CB5', '#D87463', '#7C67A4')[k], alpha = .3, na.rm = TRUE) +
			geom_line(data = avg, aes(group = 1), linewidth = 1, color = 'black', na.rm = TRUE) +
			geom_hline(yintercept = 0, linetype = 2, color = 'grey65') +
			labs(title = paste0(LETTERS[k + 1], '. Cluster ', k, ' (n=', n_distinct(d$feature), ')'), x = 'Years before diagnosis', y = 'Z score') +
			theme_5c(9)
	})
	caption <- paste0(
		z$status$interpretation[1], '. Heatmap: ', length(ord), ' incident-associated ', omic, 's; ', nrow(cl),
		' enter Ward clustering (k = 3; |Z| > 0.25 at any supported time). Thin lines: individual ', omic, 's; black: cluster mean. White heatmap cells: |Z| ≤ 0.25; gray: unsupported LOESS range.'
	)
	save_plot(wrap_plots(c(list(p), lines), design = 'AB\nAC\nAD', widths = c(1.3, 1), heights = c(1, 1, 1)) + plot_annotation(caption = caption),
		'c1.Fig7.diagnosis_timed_profiles.png', 22, 12,
		outdir = outdir
	)
}
le8_mock_bar <- function(z, title) {
	if (!nrow(z)) return(blank_plot(title, 'No FDR-significant terms / annotation unavailable'))
	z <- z |>
		filter(is.finite(adjusted_p), adjusted_p < .05) |>
		arrange(adjusted_p) |>
		slice_head(n = 12)
	if (!nrow(z)) return(blank_plot(title, 'No FDR-significant terms'))
	ggplot(z, aes( - log10(pmax(adjusted_p, 1e-300)), reorder(term_name, - adjusted_p), fill = source)) +
		geom_col(width = .7) +
		scale_fill_manual(values = c('GO:BP' = '#BC5366', 'GO:CC' = '#286B8B', 'GO:MF' = '#6D629A', KEGG = '#388B78', MGI = '#C58D37', TF = '#687784')) +
		scale_y_discrete(labels = function(x) stringr::str_wrap(x, 45)) +
		labs(title = title, x = expression( - log[10]('FDR')), y = NULL, fill = NULL) +
		theme_5c(8) +
		theme(axis.text.y = element_text(size = 8))
}
le8_mock_network <- function(edges, title, values = NULL) {
	if (!nrow(edges)) return(blank_plot(title, 'No supported edges / annotation unavailable'))
	edges <- edges |> distinct(from, to, .keep_all = TRUE) ; g <- igraph::graph_from_data_frame(as.data.frame(edges[, c('from', 'to')]), directed = FALSE)
	set.seed(SEED) ; xy <- igraph::layout_with_fr(g) ; n <- tibble(node = igraph::V(g)$name, x = xy[, 1], y = xy[, 2], value = as.numeric(igraph::degree(g)))
	if (!is.null(values)) n$value <- values[n$node]
	label_nodes <- names(sort(igraph::degree(g), decreasing = TRUE))[seq_len(min(18, nrow(n)))]
	n$label <- ifelse(n$node %in% label_nodes, n$node, '')
	if (!'score' %in% names(edges)) edges$score <- 1
	e <- edges |>
		left_join(n |> select(from = node, x, y), by = 'from') |>
		left_join(n |> select(to = node, xend = x, yend = y), by = 'to')
	ggplot() +
		geom_segment(data = e, aes(x, y, xend = xend, yend = yend, linewidth = score), color = 'grey75') +
		geom_point(data = n, aes(x, y, color = value), size = 3) +
		ggrepel::geom_text_repel(data = n, aes(x, y, label = label), max.overlaps = Inf, seed = SEED, size = 2.3) +
		scale_linewidth_continuous(range = c(.2, 1.1), guide = 'none') +
		scale_color_gradient(low = '#AEC9DB', high = '#9E2847', na.value = 'grey60') +
		labs(title = title, subtitle = 'Labels: up to 18 highest-degree nodes; all edges retained', color = if (is.null(values)) 'Degree' else '-log10(P)') +
		theme_void() +
		theme(plot.title = element_text(face = 'bold'))
}
# Public annotation dictionaries and atlas polygons live outside the script tree.
le8_go_dir <- function() Sys.getenv("LE8_GO_DIR", unset = if (Sys.info()[["sysname"]] == "Windows") "F:/annot/go" else "/mnt/f/annot/go")
le8_ggseg_dir <- function() Sys.getenv("LE8_GGSEG_DIR", unset = if (Sys.info()[["sysname"]] == "Windows") "F:/annot/ggseg" else "/mnt/f/annot/ggseg")

le8_mock_enrich <- function(a, enrich, outdir, enrich_prev = tibble()) {
	rawdir <- le8_job_dir(outdir, 'c1_correlate') ; genes <- a$term[a$p.value < .05 / nrow(a) & is.finite(a$p.value)]
	sf <- file.path(rawdir, 'c1.mock_gene_universe.csv') ; write_raw_csv(a |> transmute(gene = term, selected = term %in% genes), 'c1.mock_gene_universe.csv', rawdir)
	cached <- all(file.exists(file.path(rawdir, c('c1.mock_function_terms.csv', 'c1.mock_tf_edges.csv', 'c1.mock_ppi_edges.csv'))))
	# Download dictionaries are reusable inputs; only derived analyses are published.
	annotation_cache <- le8_cache_dir('annotations',basename(outdir))
	py <- Sys.which('python3') ; status <- if (LE8_REUSE_RESULTS && cached) 0L else if (nzchar(py)) system2(py, shQuote(c(file.path(Sys.getenv('LE8_FDIR'), 'c1.abm.py'), 'annotations', sf, rawdir, annotation_cache)), timeout = 600) else 1L
	rd <- function(n) {
		f <- file.path(rawdir, n) ; if (status == 0L && file.exists(f)) as_tibble(data.table::fread(f)) else tibble()
	}
	terms <- rd('c1.mock_function_terms.csv') ; tf <- rd('c1.mock_tf_edges.csv') ; ppi <- rd('c1.mock_ppi_edges.csv')
	if (!'source' %in% names(terms)) terms <- tibble(source = character(), term_name = character(), adjusted_p = numeric())
	go <- enrich |> filter(source %in% c('GO:BP', 'GO:CC', 'GO:MF'))
	label_file <- file.path(le8_go_dir(), 'go_terms.tsv')
	if (nrow(go) && file.exists(label_file)) {
		labels <- data.table::fread(label_file) ; ix <- match(go$term_name, labels$term)
		hit <- !is.na(ix) ; go$term_name[hit] <- labels$TERM[ix[hit]]
	}
	pathway <- bind_rows(go, terms |> filter(source == 'KEGG')) |>
		group_by(source) |>
		slice_min(adjusted_p, n = 3, with_ties = FALSE) |>
		ungroup()
	values <- setNames( - log10(pmax(a$p.value, 1e-300)), a$term)
	# Six readable panels integrate prevalent and incident enrichment with annotations.
	gp <- enrich_prev
	if (nrow(gp) && file.exists(label_file)) {
		labels <- data.table::fread(label_file) ; ix <- match(gp$term_name, labels$term)
		hit <- !is.na(ix) ; gp$term_name[hit] <- labels$TERM[ix[hit]]
	}
	top_sources <- function(d) if (nrow(d)) d |>
		group_by(source) |>
		slice_min(adjusted_p, n = 3, with_ties = FALSE) |>
		ungroup() else d
	panels <- list(
		le8_mock_bar(top_sources(go), 'A. Incident GO enrichment'),
		le8_mock_bar(top_sources(gp), 'B. Prevalent enrichment'),
		le8_mock_bar(terms |> filter(source == 'KEGG'), 'C. Incident pathways'),
		le8_mock_bar(terms |> filter(source == 'MGI'), 'D. Mammalian phenotype'),
		le8_mock_network(tf, 'E. Transcription-factor targets', values),
		le8_mock_network(ppi, 'F. Physical protein interactions')
	)
	save_plot(wrap_plots(panels, ncol = 2) + plot_annotation(caption = 'Bonferroni-selected associations; library-wise FDR and assayed background. Networks are annotation support, not causal evidence.'),
		'c1.Fig8.enrich_sig.png', 20, 18,
		outdir = outdir
	)
	write_xlsx2(list(incident = enrich, prevalent = enrich_prev, annotations = terms, TF_edges = tf, PPI_edges = ppi, annotation_status = data.frame(exit_code = status)), le8_artifact_path('c1.Fig8.enrichment.xlsx', outdir))
	list(terms = terms, tf_edges = tf, ppi_edges = ppi, annotation_exit = status)
}
le8_mock_c1 <- function(dat, a, enrich, layer, covars, tvar, evar, outdir, enrich_prev = tibble()) {
	save_plot(le8_mock_volcano(a, layer), 'c1.Fig17.measured_volcano.png', 10, 7, outdir = outdir)
	z <- le8_mock_trajectories(dat, a, covars, tvar, evar, layer)
	for (n in names(z)) write_raw_csv(z[[n]], paste0('c1.mock_', n, '.csv'), le8_job_dir(outdir, 'c1_correlate'))
	le8_plot_mock2(z, outdir, layer)
	e <- if (layer == 'protein') le8_mock_enrich(a, enrich, outdir, enrich_prev) else list()
	list(temporal = z, functional = e)
}

le8_mock_final <- function(obj, outdir) {
	imp <- as_tibble(obj$scores$Yu$importance %||% tibble())
	if (!nrow(imp) && !is.null(obj$scores$Yu$fit)) imp <- tryCatch(as_tibble(lightgbm::lgb.importance(obj$scores$Yu$fit)), error = function(e) tibble())
	pa <- if (nrow(imp)) imp |>
		arrange(desc(Gain)) |>
		slice_head(n = 15) |>
		ggplot(aes(reorder(Feature, Gain), Gain)) +
		geom_col(fill = '#BA5272') +
		coord_flip() +
		labs(title = 'A. Training LightGBM feature importance', x = NULL, y = 'Split gain') +
		theme_5c(8) else blank_plot('A. Feature importance', 'LightGBM importance unavailable')
	seq <- obj$sequential$log %||% tibble()
	curve <- if (nrow(seq)) ggplot(seq, aes(step, auc)) +
		geom_line(color = '#356C96') +
		geom_point() +
		labs(title = 'B. Training-only forward selection', x = 'Number of predictors', y = 'Inner holdout AUC') +
		theme_5c(9) else blank_plot('B. Forward selection', 'No estimable training selection curve')
	pred <- obj$prediction ; label <- unique(pred$biom_set[grepl('LightGBM', pred$biom_set)])
	if (!length(label)) label <- unique(pred$biom_set)[1]
	d <- pred |> filter(biom_set == label[1]) ; if ('split' %in% names(d)) d <- d |> filter(split == 'validation')
	full_h <- as.numeric(Sys.getenv('FINAL_MOCK_HORIZON', '15')) ; windows <- list(c(0, full_h), c(0, 10), c(10, full_h))
	curves <- list() ; metrics <- list() ; plots <- lapply(seq_along(windows), function(i) {
		lo <- windows[[i]][1] ; hi <- windows[[i]][2] ; z <- d |>
			filter(time > lo) |>
			mutate(time = time - lo)
		rr <- z |>
			group_by(model) |>
			group_modify( ~ {
				q <- ipcw_roc_curve(.x$time, .x$event, .x$score, hi - lo)
				if (nrow(q)) q$AUC <- weighted_time_auc(.x$time, .x$event, .x$score, hi - lo)
				q
			}) |>
			ungroup()
		title <- paste0(LETTERS[i + 2], '. ', if (lo == 0) paste0('Within ', hi, ' years') else paste0('Years ', lo, '–', hi, '; event-free at ', lo))
		if (!nrow(rr)) return(blank_plot(title, 'Insufficient cases/controls with valid predictions'))
		rr$window <- paste(lo, hi, sep = '-') ; curves[[i]] <<- rr
		mm <- rr |>
			group_by(model) |>
			summarise(AUC = first(AUC), .groups = 'drop') |>
			mutate(window = paste(lo, hi, sep = '-')) ; metrics[[i]] <<- mm
		text <- paste(sprintf('%s: %.3f', mm$model, mm$AUC), collapse = '\n')
		ggplot(rr, aes(fpr, tpr, color = model)) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = 'grey70') +
			geom_step(linewidth = .8) +
			annotate('text', x = .98, y = .05, label = text, hjust = 1, vjust = 0, size = 2.7) +
			coord_equal() +
			scale_color_manual(values = c(`Model 0` = 'grey45', Biomarkers = '#B64C69', Combined = '#367F9B')) +
			labs(title = title, x = 'False positive rate', y = 'True positive rate', color = NULL) +
			theme_5c(9)
	})
	save_plot(
		(pa | curve) / wrap_plots(plots, nrow = 1) +
			plot_annotation(caption = 'A: LightGBM importance; B: separate forward-selection path. Held-out predictions; IPCW ROC accounting for censoring. Fixed baseline models, no validation-based feature selection. The >10-year panel is a landmark evaluation.'),
		'Fig14.heldout_ROC.png', 20, 11.5,
		outdir = outdir
	)
	rawdir <- le8_job_dir(outdir, 'final_prediction')
	write_raw_csv(imp, 'mock_importance.csv', rawdir) ; write_raw_csv(bind_rows(metrics), 'mock_auc.csv', rawdir)
	write_raw_csv(bind_rows(curves), 'mock_roc.csv', rawdir)
	write_xlsx2(list(importance = imp, selection = seq, auc = bind_rows(metrics), ROC = bind_rows(curves)), le8_artifact_path('Fig14.heldout_ROC.xlsx', outdir))
	list(importance = imp, auc = bind_rows(metrics))
}

le8_mock_brain <- function(counts, atlas_name, title) {
	if (!requireNamespace('sf', quietly = TRUE)) return(blank_plot(title, 'Install sf to render atlas polygons'))
	env <- new.env() ; load(file.path(le8_ggseg_dir(), paste0(atlas_name, '.rda')), envir = env)
	a <- get(atlas_name, env) ; d <- sf::st_as_sf(as.data.frame(a$data))
	d$n_sig <- counts$n_sig[match(d$label, counts$atlas_label)]
	if (atlas_name == 'dk') {
		# Reposition the four real atlas views in a compact 2x2 grid (translation only).
		group <- interaction(d$hemi, d$side, drop = TRUE)
		for (g in levels(group)) {
			ii <- which(group == g) ; bb <- sf::st_bbox(d[ii, ]) ; dx <- if (d$hemi[ii[1]] == 'right') 400 else 0 ; dy <- if (d$side[ii[1]] == 'lateral') 240 else 0
			sf::st_geometry(d)[ii] <- sf::st_geometry(d)[ii] + c(dx - bb['xmin'], dy - bb['ymin'])
		}
	} else d <- d[d$side == 'coronal', ]
	# These cached MULTIPOLYGON atlases are planar drawings. Render their native
	# vertices without coord_sf's geographic transformations or PROJ database.
	# Keep separate polygons and rings (including holes) within each atlas region.
	vertices <- as.data.frame(sf::st_coordinates(d))
	vertices$n_sig <- d$n_sig[vertices$L3]
	vertices$polygon <- interaction(vertices$L3, vertices$L2, drop = TRUE)
	vertices$ring <- interaction(vertices$L3, vertices$L2, vertices$L1, drop = TRUE)
	ggplot(vertices, aes(X, Y, group = polygon, subgroup = ring, fill = n_sig)) +
		geom_polygon(rule = 'evenodd', color = 'grey45', linewidth = .12) +
		coord_equal() +
		scale_fill_gradient(low = '#FFF3DC', high = '#BD3651', na.value = 'grey88', limits = c(0, max(1, counts$n_sig, na.rm = TRUE))) +
		labs(title = title, fill = 'Significant\nbiomarkers') +
		theme_void() +
		theme(plot.title = element_text(face = 'bold', size = 10), legend.position = 'right')
}
le8_mock_c4 <- function(a, outdir) {
	if (!nrow(a)) return(invisible(NULL))
	global <- a |>
		filter(measure %in% c('img_total_cortical_gm', 'img_p26517', 'img_p24486')) |>
		mutate(structure = recode(measure, img_total_cortical_gm = 'CGV', img_p26517 = 'SGV', img_p24486 = 'WMH'), mark = case_when(bonferroni < .05 ~ '**', FDR_all < .05 ~ '*', TRUE ~ ''))
	top_global <- global |>
		group_by(feature) |>
		summarise(pmin = min(p, na.rm = TRUE), .groups = 'drop') |>
		arrange(pmin) |>
		slice_head(n = 35) |>
		pull(feature)
	global <- global |>
		filter(feature %in% top_global) |>
		mutate(feature = factor(feature, levels = rev(top_global)))
	pa <- ggplot(global, aes(structure, feature, fill = beta)) +
		geom_tile() +
		geom_text(aes(label = mark), size = 2) +
		scale_fill_gradient2(low = '#3678B0', mid = 'white', high = '#CC4656', midpoint = 0, na.value = 'grey88') +
		labs(title = 'A. Global structure: top 35 features', x = NULL, y = NULL, fill = 'Standardized\nbeta') +
		theme_5c(7) +
		theme(panel.grid = element_blank())
	if (n_distinct(global$feature) > 150) pa <- pa + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + labs(y = paste(n_distinct(global$feature), 'features; row names in CSV'))
	counts <- a |>
		group_by(family, label, measure) |>
		summarise(n_sig = sum(FDR_all < .05, na.rm = TRUE), tested = sum(is.finite(p)), .groups = 'drop')
	counts <- counts |> mutate(atlas_label = case_when(
		family %in% c('Cortical area', 'Cortical volume') ~ sub('^aparc-Desikan_([lr]h)_(area|volume)_', '\\1_', label),
		family == 'Subcortical volume' ~ sub('^aseg_rh_volume_', 'Right-', sub('^aseg_lh_volume_', 'Left-', label)), TRUE ~ NA_character_
	))
	pb <- le8_mock_brain(counts |> filter(family == 'Cortical area'), 'dk', 'B. Cortical surface area')
	pc <- le8_mock_brain(counts |> filter(family == 'Cortical volume'), 'dk', 'C. Cortical gray matter volume')
	pd <- le8_mock_brain(counts |> filter(family == 'Subcortical volume'), 'aseg', 'D. Subcortical gray matter volume')
	wm <- counts |>
		filter(family == 'White matter FA/MD') |>
		mutate(metric = ifelse(grepl('_FA_', label), 'FA', 'MD'), tract = sub('^.*_(FA|MD)_', '', label))
	pe <- ggplot(wm, aes(metric, tract, fill = n_sig)) +
		geom_tile(color = 'white') +
		geom_text(aes(label = n_sig), size = 2.4) +
		scale_fill_gradient(low = '#FFF3DC', high = '#BD3651') +
		labs(title = 'E. White matter tracts', x = NULL, y = NULL, fill = 'Significant\nbiomarkers') +
		theme_5c(8) +
		theme(panel.grid = element_blank())
	save_plot(
		(pa | wrap_plots(pb, pc, pd, pe, ncol = 2)) + plot_layout(widths = c(.7, 2)) +
			plot_annotation(caption = 'A: * global BH < 0.05; ** Bonferroni < 0.05. B–E: number of biomarkers passing global BH. Gray atlas regions: not mapped/tested. DK/aseg atlas; FA/MD shown by measured tract.'),
		'c4.Fig14.imaging_atlas.png', 21, 12,
		outdir = outdir
	)
	write_raw_csv(counts, 'c4.mock_brain_region_counts.csv', le8_job_dir(outdir, 'c4_connect'))
	write_xlsx2(list(global_top35 = global, region_counts = counts), le8_artifact_path('c4.Fig14.imaging_atlas.xlsx', outdir))
	invisible(counts)
}


# 0.common.R
# Final LE8 5C methods. UKB data are read from external inputs.
# Shared integrity, outcome and reproducibility functions.
LE8_CODE_VERSION <- "2026-10-05.integrated-5c-v2"
.le8_analysis_state <- new.env(parent = emptyenv())
.le8_analysis_state$hashes <- new.env(parent = emptyenv())

le8_num_env <- function(name, default) {
	z <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
	if (length(z) != 1L || !is.finite(z)) stop("Invalid numeric setting: ", name, call. = FALSE)
	z
}
le8_csv_env <- function(name, default = "") {
	z <- trimws(strsplit(Sys.getenv(name, unset = default), ",", fixed = TRUE)[[1]])
	unique(z[nzchar(z)])
}
le8_hash_object <- function(x) {
	f <- tempfile() ; on.exit(unlink(f), add = TRUE)
	base::saveRDS(x, f, version = 2, compress = FALSE)
	unname(tools::md5sum(f))
}
le8_file_hash <- function(path) {
	if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) return("missing")
	path <- normalizePath(path, winslash = "/", mustWork = TRUE)
	fi <- file.info(path)
	if (isTRUE(fi$isdir)) return("directory")
	key <- paste(path, fi$size, as.numeric(fi$mtime), as.numeric(fi$ctime), sep = "|")
	h <- .le8_analysis_state$hashes
	if (exists(key, envir = h, inherits = FALSE)) return(get(key, envir = h))
	# Persist full-content digests; file ctime is part of the invalidation key.
	hd <- file.path(tempdir(), "le8-input-hashes") ; dir.create(hd, recursive = TRUE, showWarnings = FALSE)
	cf <- file.path(hd, paste0(le8_hash_object(path), ".rds"))
	old <- if (file.exists(cf)) tryCatch(readRDS(cf), error = function(e) NULL) else NULL
	if (!is.null(old) && identical(old$key, key)) value <- old$value else {
		value <- unname(tools::md5sum(path))
		if (is.na(value)) stop("Cannot fingerprint input: ", path, call. = FALSE)
		base::saveRDS(list(key = key, value = value), cf)
	}
	assign(key, value, envir = h) ; value
}
le8_table <- function(x) if (is.null(x)) tibble::tibble() else tibble::as_tibble(x)
le8_append_workbook <- function(file, tables, prefix = "review_") {
	if (!length(tables)) return(invisible(NULL))
	file <- le8_artifact_path(file)
	wb <- if (file.exists(file)) openxlsx::loadWorkbook(file) else openxlsx::createWorkbook()
	for (nm in names(tables)) {
		d <- as.data.frame(tables[[nm]])
		if (!ncol(d)) d <- data.frame(status = "No estimable result / unavailable")
		sh <- substr(paste0(prefix, nm), 1, 31)
		if (sh %in% names(wb)) openxlsx::removeWorksheet(wb, sh)
		openxlsx::addWorksheet(wb, sh)
		if (nrow(d)) openxlsx::writeDataTable(wb, sh, d) else openxlsx::writeData(wb, sh, d)
		openxlsx::freezePane(wb, sh, firstRow = TRUE)
		openxlsx::setColWidths(wb, sh, cols = seq_len(ncol(d)), widths = 18)
	}
	openxlsx::saveWorkbook(wb, file, overwrite = TRUE) ; invisible(file)
}

# This function is also inserted verbatim into 0f/phenotype.R by the installer.
t2e <- function(dat, domain, Y_date, birth_date, date_attend, date_lost, date_death,
			date_end = date_follow_end, prefix = NA, time_unit = "year") {
	n <- nrow(dat)
	getdate <- function(x) {
		if (inherits(x, "Date")) {
			if (length(x) == 1L) return(rep(x, n))
			if (length(x) != n) stop("Date vector length mismatch")
			return(x)
		}
		if (is.null(x) || length(x) != 1L || is.na(x) || !x %in% names(dat)) return(rep(as.Date(NA), n))
		z <- dat[[x]]
		if (is.numeric(z) && !inherits(z, "Date")) as.Date(z, origin = "1970-01-01") else as.Date(z)
	}
	ba <- getdate(date_attend) ; bd <- getdate(birth_date) ; yd <- getdate(Y_date)
	lost <- getdate(date_lost) ; dead <- getdate(date_death)
	end <- as.Date(date_end) ; if (length(end) == 1L) end <- rep(end, n)
	if (length(end) != n || anyNA(end)) stop("Administrative end date must be valid", call. = FALSE)
	censor <- pmin(lost, dead, end, na.rm = TRUE)
	valid <- !is.na(ba) & !is.na(censor) & censor >= ba
	# Same-day diagnosis is baseline disease (never an incident non-case).
	prevalent <- !is.na(yd) & !is.na(ba) & yd <= ba
	incident <- valid & !is.na(yd) & yd > ba & yd <= censor
	beyond <- valid & !is.na(yd) & yd > censor
	dirty <- rep(FALSE, n)
	if (!is.null(domain) && length(domain) == 1L && !is.na(domain)) {
		dv <- if (domain %in% names(dat)) domain else paste0("icd10Ct_", domain)
		if (dv %in% names(dat)) {
			zz <- dat[[dv]]
			evidence <- if (inherits(zz, "Date")) !is.na(zz) else suppressWarnings(as.numeric(zz)) > 0
			dirty <- is.na(yd) & !is.na(evidence) & evidence
		}
	}
	yt <- ifelse(valid & !prevalent & !dirty, as.numeric(incident), NA_real_)
	yr <- ifelse(valid & !dirty & !incident, as.numeric(prevalent), NA_real_)
	stopdate <- censor ; stopdate[incident] <- yd[incident]
	forward <- ifelse(!is.na(yt), as.numeric(stopdate - ba), NA_real_)
	# Backward clock remains descriptive; non-cases' r2e equals baseline age.
	reverse <- ifelse(!is.na(yr), ifelse(prevalent, as.numeric(ba - yd), as.numeric(ba - bd)), NA_real_)
	signed <- ifelse(prevalent, as.numeric(yd - ba), ifelse(incident, forward, NA_real_))
	exitdate <- stopdate ; exitdate[prevalent] <- yd[prevalent]
	ageexit <- ifelse((prevalent | !is.na(yt)) & !is.na(bd), as.numeric(exitdate - bd), NA_real_)
	divisor <- if (time_unit == "year") 365.25 else if (time_unit %in% c("day", "days")) 1 else
		stop("time_unit must be year or day", call. = FALSE)
	vals <- list(
		Yt2e = yt, Yr2e = yr, t2e = forward / divisor, r2e = reverse / divisor,
		b2e = signed / divisor, bi2e = ageexit / divisor
	)
	nm <- names(vals) ; if (length(prefix) == 1L && !is.na(prefix)) nm <- paste0(prefix, ".", nm)
	dat[nm] <- vals
	attr(dat, "le8_outcome_audit") <- data.frame(
		metric = c(
			"N", "prevalent_including_same_day",
			"same_day", "incident_within_followup", "diagnosis_after_censor", "invalid_followup", "unknown_date_with_disease_evidence"
		),
		value = c(
			n, sum(prevalent), sum(!is.na(yd) & !is.na(ba) & yd == ba), sum(incident),
			sum(beyond), sum(!valid), sum(dirty)
		)
	)
	attr(dat, "le8_t2e_version") <- "2026-09-05.censor-inclusive-baseline-v1"
	dat
}

make_outcome <- function(dat, outcome = Y) {
	req <- c("birth_date", "date_attend", paste0("fod_icd10_", outcome))
	if (length(setdiff(req, names(dat)))) stop("Missing outcome fields: ", paste(setdiff(req, names(dat)), collapse = ", "))
	# Always rebuild: old phenotype RDS may contain pre-fix t2e variables.
	for (v in c("date_lost", "date_death")) if (!v %in% names(dat)) dat[[v]] <- as.Date(NA)
	ans <- t2e(
		dat, NA, paste0("fod_icd10_", outcome), "birth_date", "date_attend",
		"date_lost", "date_death", date_follow_end, outcome, "year"
	)
	if (!is.null(.le8_analysis_state$rawdir)) {
		au <- attr(ans, "le8_outcome_audit")
		au$outcome <- outcome
		data.table::fwrite(au, file.path(
			.le8_analysis_state$rawdir,
			paste0("outcome_audit.", nrow(ans), ".csv")
		))
	}
	ans
}
make_prevalent_status <- function(dat, outcome = Y) {
	y <- as.Date(dat[[paste0("fod_icd10_", outcome)]]) ; b <- as.Date(dat$date_attend)
	ifelse(is.na(b), NA_real_, as.numeric(!is.na(y) & y <= b))
}

le8_begin_analysis <- function(layer, module) {
	.le8_analysis_state$rawdir <- le8_guard_stage_options(layer, module)
	# Legacy drug.lipid may include post-baseline visits: do not silently use it.
	if (!truthy(Sys.getenv("LE8_TREATMENT_BASELINE_CONFIRMED", unset = "FALSE"))) {
		for (nm in c("C1_TREATMENT_VARS")) if (exists(nm, .GlobalEnv, inherits = FALSE)) {
			z <- get(nm, .GlobalEnv)
			if ("drug.lipid" %in% z) {
				warning("drug.lipid excluded: baseline-only provenance unconfirmed; use a baseline-only field or LE8_TREATMENT_BASELINE_CONFIRMED=TRUE after verification", call. = FALSE)
				assign(nm, setdiff(z, "drug.lipid"), .GlobalEnv)
			}
		}
		Sys.setenv(C2_TREATMENT_VARS = paste(setdiff(le8_csv_env("C2_TREATMENT_VARS", Sys.getenv("C1_TREATMENT_VARS")), "drug.lipid"), collapse = ","))
	}
	invisible(NULL)
}

# An assay can differ from its encoding-gene symbol. User mappings override
# this small nomenclature map. Gene coordinates still come from annotation.
le8_assay_genes <- function(features) {
	out <- setNames(as.character(features), as.character(features))
	alias <- c(NTPROBNP = "NPPB", NTproBNP = "NPPB", `NT-proBNP` = "NPPB")
	hit <- intersect(names(out), names(alias)) ; out[hit] <- alias[hit]
	f <- Sys.getenv("LE8_ASSAY_GENE_MAP", unset = "")
	if (nzchar(f)) {
		if (!file.exists(f)) stop("LE8_ASSAY_GENE_MAP is missing: ", f)
		m <- data.table::fread(f)
		if (!all(c("assay", "gene") %in% names(m))) stop("Assay map needs assay and gene columns")
		if (anyDuplicated(m$assay)) stop("Assay map has duplicate assay identifiers")
		hit <- intersect(names(out), m$assay) ; out[hit] <- m$gene[match(hit, m$assay)]
	}
	out
}
le8_coloc_pass <- function(z, policy=le8_evidence_policy()) {
	if (!nrow(z)) return(logical())
	if ("policy_hash" %in% names(z) && any(!is.na(z$policy_hash) & nzchar(z$policy_hash) & z$policy_hash!=policy$hash)) stop("Stored C3 evidence uses a different policy; use its original settings or refresh affected C3 evidence before consolidation")
	if (!all(c("status","PP.H4_robust_min","prior_complete") %in% names(z))) return(rep(FALSE,nrow(z)))
	pass <- z$status %in% "ok" & z$prior_complete %in% TRUE & is.finite(z$PP.H4_robust_min) & z$PP.H4_robust_min>=policy$posterior
	if ("prior_status" %in% names(z)) pass <- pass & !z$prior_status %in% c("prior_incomplete","coverage_incomplete")
	if (policy$level=="signal_only") pass <- pass & if("MR_signal_status" %in% names(z)) z$MR_signal_status %in% "all_instrument_signals_supported" else FALSE
	pass
}
le8_mr_class <- function(layer, scope="primary") {
	if (!layer %in% c("protein","prot","metabolite","met")) stop("Unknown MR layer")
	if (!scope %in% c("primary","secondary")) stop("Unknown MR scope")
	if(layer %in% c("protein","prot")) { if(scope=="primary") "cis" else "trans" } else if(scope=="primary") "local" else "distal"
}
le8_select_mr <- function(mr, layer, scope="primary") {
	if (!nrow(mr)) return(mr)
	if (!all(c("exposure","analysis") %in% names(mr))) stop("MR results lack exposure/analysis")
	z <- dplyr::distinct(mr[mr$analysis %in% le8_mr_class(layer,scope),,drop=FALSE])
	if (anyDuplicated(z$exposure)) stop("Multiple MR estimates in one prespecified analysis class; resolve method identity, do not select by P")
	z
}
le8_bidirectional_evidence <- function(forward, reverse, layer, policy=le8_evidence_policy()) {
	rv <- if(nrow(reverse) && all(c("feature","b","pval","FDR_reverse","n_IV") %in% names(reverse))) reverse |>
		dplyr::transmute(feature,reverse_beta=b,reverse_p=pval,reverse_FDR=FDR_reverse,n_reverse_IV=n_IV) else
		tibble::tibble(feature=character(),reverse_beta=numeric(),reverse_p=numeric(),reverse_FDR=numeric(),n_reverse_IV=integer())
	if (anyDuplicated(rv$feature)) stop("Duplicate reverse MR estimates")
	features <- union(forward$exposure,rv$feature)
	dplyr::bind_rows(lapply(c("primary","secondary"),function(scope) {
		fw <- le8_select_mr(forward,layer,scope)
		if(nrow(fw)) fw <- fw |> dplyr::transmute(feature=exposure,forward_beta=b,forward_p=pval,forward_FDR=FDR_all,forward_present=TRUE) else
			fw <- tibble::tibble(feature=character(),forward_beta=numeric(),forward_p=numeric(),forward_FDR=numeric(),forward_present=logical())
		dplyr::left_join(tibble::tibble(feature=features),fw,by="feature") |> dplyr::left_join(rv,by="feature") |>
			dplyr::mutate(scope=.env$scope,analysis=le8_mr_class(layer,.env$scope),policy_fdr=policy$mr_fdr,
				forward_status=dplyr::case_when(is.na(forward_present)~"not_tested",!is.finite(forward_p)~"failed_or_not_estimable",TRUE~"tested"),
				reverse_status=dplyr::case_when(is.na(n_reverse_IV)~"not_tested",!is.finite(reverse_p)~"failed_or_not_estimable",TRUE~"tested"),
				forward_support=is.finite(forward_FDR)&forward_FDR<policy$mr_fdr,
				reverse_support=is.finite(reverse_FDR)&reverse_FDR<policy$mr_fdr,
				class=dplyr::case_when(forward_support & reverse_support~"Bidirectional genetic support",forward_support~"Omic -> disease only",
					reverse_support~"Disease liability -> omic only",forward_status!="tested"|reverse_status!="tested"~"Incomplete/untested evidence",TRUE~"No FDR support"))
	}))
}
le8_same_locus_evidence <- function(mr, coloc, layer, policy=le8_evidence_policy()) {
	empty <- tibble::tibble(
		feature = character(), locus = character(), locus_class = character(),
		MR_FDR = numeric(), coloc_robust_min = numeric(), instrument_overlap = logical(),
		eligible = logical(), claim = character()
	)
	if (!is.data.frame(mr) || !nrow(mr) || !is.data.frame(coloc) || !nrow(coloc)) return(empty)
	need_m <- c("exposure", "analysis", "FDR_all", "instrument_chr", "instrument_pos_min", "instrument_pos_max")
	need_c <- c("feature", "locus", "chr", "start", "end", "status", "PP.H4_robust_min", "locus_class")
	if (!all(need_m %in% names(mr)) || !all(need_c %in% names(coloc))) return(empty)
	primary <- le8_mr_class(layer)
	m <- le8_select_mr(mr,layer)
	c <- coloc[coloc$locus_class == primary & coloc$status == "ok", , drop = FALSE]
	if (!nrow(m) || !nrow(c)) return(empty)
	z <- dplyr::inner_join(c, m, by = c("feature" = "exposure"))
	# The interval must overlap actual retained MR instruments, not merely the
	# same feature. Single-signal ABF remains region-level, not signal-resolved.
	z$instrument_overlap <- vapply(seq_len(nrow(z)), function(i) {
		same <- as.character(z$chr[i]) %in% strsplit(as.character(z$instrument_chr[i]), ";", fixed = TRUE)[[1]]
		if (!same) return(FALSE)
		if ("instrument_positions" %in% names(z)) {
			keys <- strsplit(as.character(z$instrument_positions[i]), ";", fixed = TRUE)[[1]]
			pp <- strsplit(keys, ":", fixed = TRUE)
			return(any(vapply(pp, function(q) length(q) == 2L && q[1] == as.character(z$chr[i]) &&
				is.finite(suppressWarnings(as.numeric(q[2]))) && as.numeric(q[2]) >= z$start[i] && as.numeric(q[2]) <= z$end[i], logical(1))))
		}
		# Old min/max ranges do not prove a particular retained SNP overlaps.
		FALSE
	}, logical(1))
	z$region_IV_coverage <- vapply(seq_len(nrow(z)),function(i) {
		if (!"instrument_positions" %in% names(z)) return(NA_real_)
		pp <- strsplit(strsplit(as.character(z$instrument_positions[i]),";",fixed=TRUE)[[1]],":",fixed=TRUE)
		if (!length(pp)) return(NA_real_)
		mean(vapply(pp,function(q) length(q)==2 && q[1]==as.character(z$chr[i]) && is.finite(suppressWarnings(as.numeric(q[2]))) && as.numeric(q[2])>=z$start[i] && as.numeric(q[2])<=z$end[i],logical(1)))
	},numeric(1))
	z$IV_coverage <- vapply(seq_len(nrow(z)),function(i) {
		if (!all(c("instrument_snps","aligned_MR_IVs") %in% names(z))) return(NA_real_)
		ids <- strsplit(as.character(z$instrument_snps[i]),";",fixed=TRUE)[[1]]
		kept <- strsplit(as.character(z$aligned_MR_IVs[i]),";",fixed=TRUE)[[1]]
		if (!length(ids) || anyNA(ids) || any(!nzchar(ids))) return(NA_real_)
		mean(ids %in% kept)
	},numeric(1))
	z$MR_signal_status <- if ("MR_signal_status" %in% names(z)) z$MR_signal_status else NA_character_
	z$MR_signal_coverage <- if ("MR_signal_coverage" %in% names(z)) z$MR_signal_coverage else NA_real_
	z$resolved_signal_pass <- is.na(z$MR_signal_status) | z$MR_signal_status=="all_instrument_signals_supported"
	z$input_scale_verified <- if ("beta_scale_status" %in% names(z)) z$beta_scale_status %in% "verified" else FALSE
	z$posterior_pass <- le8_coloc_pass(z,policy)
	z |> dplyr::transmute(feature, locus, locus_class, region_IV_coverage, IV_coverage,input_scale_verified,
		signal_support=ifelse(!is.na(MR_signal_status),MR_signal_status,ifelse(IV_coverage<1,"partial_signal_support","region_only; independent signals unresolved")),
		MR_signal_coverage, prior_complete=if("prior_complete" %in% names(z)) prior_complete else FALSE,
		policy_h4=policy$posterior,policy_fdr=policy$mr_fdr,policy_hash=policy$hash,policy_level=policy$level,
		MR_FDR = FDR_all,
		coloc_robust_min = PP.H4_robust_min, instrument_overlap,
		eligible = resolved_signal_pass & input_scale_verified & instrument_overlap & is.finite(IV_coverage) & IV_coverage==1 & is.finite(FDR_all) & FDR_all < policy$mr_fdr & posterior_pass,
		claim = "Same cis/local region MR + ABF support; not a resolved causal signal or causal proof"
	)
}
le8_same_locus_candidates <- function(mr, coloc, layer, policy=le8_evidence_policy()) {
	z <- le8_same_locus_evidence(mr, coloc, layer, policy)
	unique(z$feature[z$eligible %in% TRUE])
}
# Optional additions fail locally and leave the already-produced core outputs.
le8_optional <- function(label, expr) {
	tryCatch(force(expr), error = function(e) {
		warning(label, ": ", conditionMessage(e), call. = FALSE)
		rd <- .le8_analysis_state$rawdir
		if (!is.null(rd)) data.table::fwrite(data.frame(
			module = label, status = "failed",
			message = conditionMessage(e)
		), file.path(rd, paste0(gsub("[^A-Za-z0-9_]", "_", label), ".status.csv")))
		list(status = tibble::tibble(status = "failed", message = conditionMessage(e)))
	})
}
le8_gwas_metadata <- function(file) {
	out <- list(
		N = NA_real_, N_case = NA_real_, N_control = NA_real_, sdY = NA_real_,
		type = NA_character_, beta_scale = NA_character_, build = NA_character_, ancestry = NA_character_,
		source = "unprovided", discovery_overlap = "unknown",
		normalization_proof=NA_character_, reference_fasta=NA_character_
	)
	if (length(file) != 1L || is.na(file) || !nzchar(file)) return(out)
	f <- Sys.getenv("LE8_GWAS_MANIFEST", unset = "")
	stem <- sub("[.](cis[.])?gz$", "", basename(file))
	bf <- file.path(dirname(dirname(file)), "qc", paste0(stem, ".grch"))
	if (file.exists(bf)) {
		b <- trimws(readLines(bf, warn = FALSE)[1])
		if (b %in% c("37","38")) { out$build <- b; out$source <- bf }
	}
	if (!nzchar(f) || !file.exists(f)) return(out)
	m <- data.table::fread(f, showProgress = FALSE)
	if (!"file" %in% names(m)) stop("LE8_GWAS_MANIFEST needs a file column")
	normal <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
	idx <- which(vapply(as.character(m$file), normal, character(1)) == normal(file))
	if (length(idx) > 1) stop("Duplicate exact file in GWAS manifest: ", file)
	if (!length(idx)) return(out)
	for (n in intersect(names(out), names(m))) out[[n]] <- m[[n]][idx]
	for (n in c("N", "N_case", "N_control", "sdY")) out[[n]] <- suppressWarnings(as.numeric(out[[n]]))
	if (!is.finite(out$N) && is.finite(out$N_case) && is.finite(out$N_control)) out$N <- out$N_case + out$N_control
	for (nm in c("normalization_proof","reference_fasta")) if (!is.na(out[[nm]]) && nzchar(out[[nm]]) && !grepl("^/|^[A-Za-z]:",out[[nm]])) out[[nm]] <- file.path(dirname(f),out[[nm]])
	out$source <- f ; out
}
get_case_fraction <- function(outcome = Y) {
	explicit <- suppressWarnings(as.numeric(Sys.getenv("COLOC_CASE_FRAC", unset = "NA")))
	if (is.finite(explicit) && explicit > 0 && explicit < 1) return(explicit)
	f <- get_y_gwas_file(outcome, FALSE)
	if (is.na(f)) return(NA_real_)
	m <- le8_gwas_metadata(f)
	if (is.finite(m$N_case) && is.finite(m$N_control) && m$N_case > 0 && m$N_control > 0)
		m$N_case / (m$N_case + m$N_control) else NA_real_
}

layer_annotation <- function(layer, features) {
	z <- le8_gene_coordinates(layer, features)
	if (layer != "protein") return(z)
	genes <- le8_assay_genes(features) ; z$encoding_gene <- unname(genes[z$feature])
	aliases <- names(genes)[genes != names(genes)]
	if (length(aliases)) {
		a <- le8_gene_coordinates(layer, unique(unname(genes[aliases])))
		for (i in which(z$feature %in% aliases & (!is.finite(z$start) | !is.finite(z$end) | is.na(z$chr)))) {
			j <- match(z$encoding_gene[i], a$feature)
			if (!is.na(j)) for (nm in intersect(c("chr", "start", "end", "pos"), names(z))) z[[nm]][i] <- a[[nm]][j]
		}
	}
	z
}
le8_finish_analysis <- function(layer, module, env) {
	if (exists(".le8_sumstat_qc",inherits=TRUE) && length(.le8_sumstat_qc$rows)) {
		qc <- dplyr::distinct(dplyr::bind_rows(.le8_sumstat_qc$rows))
		write_raw_csv(qc,paste0(substr(module,1,2),".variant_identity_qc.csv"),le8_job_dir(if(layer=="protein") out.prot else out.met,module))
	}
	# An early cache return has already-complete additions. A failed core run
	# must not be disguised by an expensive on-exit analysis.
	if (!exists("out", envir = env, inherits = FALSE)) return(invisible(NULL))
	obj <- get("out", env) ; outdir <- if (layer == "protein") out.prot else out.met
	get0local <- function(n, default = NULL) if (exists(n, envir = env, inherits = FALSE)) get(n, env) else default
	tables <- le8_optional(paste0(module, "_review_additions"), switch(module,
		c1_correlate = le8_c1_additions(get0local("dat"), layer, get0local("covs_adj2"), get0local("tvar"), get0local("evar")),
		c2_cause = {
			dan <- obj$DANDELION %||% list() ; de <- obj$individual_decomposition %||% list()
			le8_dandelion_plot_bundle(dan, outdir)
			list(
				dan_targets_all = dan$targets_all %||% tibble(), dan_input_audit = dan$input_audit %||% tibble(),
				decomp_folds = de$folds %||% tibble(), decomp_summary = de$summary %||% tibble()
			)
		},
		c3_coloc = le8_c3_additions(get0local("res", obj$summary), get0local("mr", tibble()), layer, outdir, obj$susie %||% list()),
		c4_connect = le8_c4_additions(obj, outdir),
		final_prediction = le8_final_additions(
			get0local("dat"), get0local("biom_vars"), get0local("ranked"),
			get0local("training_screen"), get0local("distalset", character()), get0local("geneticset", character()),
			get0local("clinical"), get0local("tvar"), get0local("evar"), get0local("split"), layer, outdir
		),
		list()
	))
	workbook <- file.path(le8_job_dir(outdir, module), if (module == "final_prediction") "prediction_panels.xlsx" else paste0(sub("_.*$", "", module), ".out.xlsx"))
	le8_optional(paste0(module, "_review_workbook"), le8_append_workbook(workbook, tables))
	# Attach aggregate additions to the same result RDS so they are not orphaned
	# from C1-C5. Do not overwrite unknown files or serialized individual data.
	cache <- get0local("selected_cache", get0local("final_cache", get0local("cache")))
	if (is.character(cache) && length(cache) == 1L && file.exists(cache) && grepl("(res[.]rds)$", cache)) {
		obj$review <- tables ; obj$meta$common_code_version <- LE8_CODE_VERSION
		base::saveRDS(obj, cache, compress = "xz")
	}
	le8_optional(paste0(module, "_output_index"), finalize_outputs(module, outdir))
	invisible(NULL)
}


# Shared marginal-effect MR and interval-specific Cox methods
le8_interval_cox <- function(dat, x, tvar, evar, covars, lo, hi, scale_x = TRUE, min_event = 20L) {
	if(length(setdiff(covars,names(dat)))) stop("Missing required interval-Cox covariates: ",paste(setdiff(covars,names(dat)),collapse=","))
	cols <- unique(c(x, tvar, evar, covars,intersect(".group",names(dat))))
	d <- as.data.frame(dat[, cols, drop = FALSE])
	d <- d[complete.cases(d), , drop = FALSE]
	d <- d[is.finite(d[[tvar]]) & d[[tvar]] > lo & d[[evar]] %in% c(0, 1), , drop = FALSE]
	d$.stop <- pmin(d[[tvar]], hi) - lo
	d$.event <- as.integer(d[[evar]] == 1 & d[[tvar]] <= hi)
	ans <- tibble(
		term = x, beta = NA_real_, std.error = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_,
		N_total = nrow(d), N_event = sum(d$.event), N_case = sum(d$.event), N_control = sum(d$.event == 0), person_years = sum(d$.stop),
		window_lo = lo, window_hi = hi, time = (lo + hi) / 2, effect_measure = "log HR; interval-specific Cox", status = "insufficient events or variation"
	)
	if (!nrow(d) || sum(d$.event) < min_event || !is.finite(sd(d[[x]])) || sd(d[[x]]) <= 0)
		return(ans)
	# Fixed baseline reference SD, not a different SD in every late risk set.
	if (scale_x) {
		ref <- as.numeric(dat[[x]])
		s <- sd(ref, na.rm = TRUE)
		m <- mean(ref, na.rm = TRUE)
		d[[x]] <- (d[[x]] - m) / s
	}
	ff <- as.formula(paste0("survival::Surv(.stop,.event) ~ ", paste(bt(c(x, covars)), collapse = " + ")))
	fit <- tryCatch(if(".group" %in% names(d) && anyDuplicated(d$.group)) survival::coxph(ff,d,ties="efron",cluster=d$.group,robust=TRUE) else survival::coxph(ff, d, ties = "efron"), error = function(e) NULL)
	if (is.null(fit))
		return(ans)
	sm <- coef(summary(fit))
	if (!x %in% rownames(sm))
		return(ans)
	b <- sm[x, "coef"]
	se <- sqrt(vcov(fit)[x,x])
	ans$beta <- b
	ans$std.error <- se
	ans$conf.low <- b - 1.96 * se
	ans$conf.high <- b + 1.96 * se
	ans$p.value <- sm[x, "Pr(>|z|)"]
	ans$status <- "ok"
	ans
}

le8_ld_read <- function(path) {
	if (length(path) != 1L || is.na(path) || !file.exists(path))
		return(NULL)
	z <- tryCatch(if (grepl("\\.rds$", path))
		readRDS(path) else read.table(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
	if (is.list(z) && !is.data.frame(z) && !is.matrix(z))
		z <- z$R
	if (is.data.frame(z)) {
		if (ncol(z) == nrow(z) + 1L) {
			rn <- as.character(z[[1]])
			z <- as.matrix(z[, - 1, drop = FALSE])
			rownames(z) <- rn
		} else z <- as.matrix(z)
	}
	if (!is.matrix(z) || nrow(z) != ncol(z) || is.null(rownames(z)) || is.null(colnames(z)))
		return(NULL)
	storage.mode(z) <- "double"
	if (anyDuplicated(rownames(z)) || anyDuplicated(colnames(z)) || !setequal(rownames(z), colnames(z)))
		return(NULL)
	z <- z[, rownames(z), drop = FALSE]
	if (any(!is.finite(z)) || max(abs(z - t(z))) > 1e-05 || max(abs(diag(z) - 1)) > 0.01 || max(abs(z)) > 1.001)
		return(NULL)
	z
}

le8_ld_greedy <- function(iv, R, r2 = 0.001) {
	stopifnot(is.matrix(R), r2 >= 0, r2 < 1)
	iv <- iv[order(iv$P, iv$SNP), , drop = FALSE]
	candidates <- intersect(iv$SNP, rownames(R))
	selected <- character()
	for (s in candidates) if (!length(selected) || all(R[s, selected] ^ 2 <= r2))
		selected <- c(selected, s)
	iv[match(selected, iv$SNP), , drop = FALSE]
}

le8_ld_for_iv <- function(iv, feature) {
	explicit <- Sys.getenv("C2_MR_LD_DIR", unset = "")
	sf <- unique(iv$source_file)
	sf <- sf[!is.na(sf) & nzchar(sf)]
	dirs <- unique(c(explicit, dirname(sf)))
	paths <- unique(unlist(lapply(dirs, function(d) file.path(d, c(paste0(feature, ".ld.rds"), paste0(
		feature,
		".ldr.cojo"
	), paste0(feature, ".jma.ldr"), paste0(feature, ".jma.cojo.ldr"))))))
	for (p in paths) {
		R <- le8_ld_read(p)
		if (!is.null(R) && sum(iv$SNP %in% rownames(R)) >= 2L)
			return(list(R = R, source = p))
	}
	list(R = NULL, source = "unavailable")
}

le8_prepare_mr_iv <- function(iv, feature) {
	iv <- iv[is.finite(iv$BETA) & is.finite(iv$SE) & iv$SE > 0 & is.finite(iv$P), , drop = FALSE]
	if (!nrow(iv))
		return(list(iv = iv, status = "no eligible variants", source = "none", n_input = 0L))
	# Even COJO-conditional discoveries must satisfy marginal relevance here.
	iv <- iv[iv$P <= le8_num_env("C2_MR_P", 5e-08) & (iv$BETA / iv$SE) ^ 2 >= le8_num_env("C2_MR_MIN_F", 10), , drop = FALSE]
	if (!nrow(iv))
		return(list(iv = iv, status = "no marginally strong variants", source = "none", n_input = 0L))
	n0 <- nrow(iv)
	ld <- le8_ld_for_iv(iv, feature)
	if (!is.null(ld$R)) {
		known <- iv[iv$SNP %in% rownames(ld$R), , drop = FALSE]
		kept <- le8_ld_greedy(known, ld$R, le8_num_env("C2_MR_LD_R2", 0.001))
		list(iv = kept, status = "LD-pruned marginal IVs", source = ld$source, n_input = n0, n_without_LD = n0 -
			nrow(known), max_retained_r2 = if (nrow(kept) > 1) max(ld$R[kept$SNP, kept$SNP][upper.tri(ld$R[
			kept$SNP,
			kept$SNP
		])] ^ 2) else 0)
	} else {
		# Never assume COJO conditional signals are mutually uncorrelated.
		kept <- iv[order(iv$P, iv$SNP)[1L], , drop = FALSE]
		list(
			iv = kept, status = if (n0 == 1) "single marginal IV" else "single-IV fallback; LD unavailable", source = "none",
			n_input = n0, n_without_LD = n0, max_retained_r2 = NA_real_
		)
	}
}

read_qtl_instruments <- function(feature, base_dir, layer = c("protein", "metabolite"), annotation = NULL) {
	layer <- match.arg(layer)
	fs <- find_qtl_files(feature, base_dir, layer)
	empty <- list(instruments = tibble(), score_instruments = tibble(), files = fs)
	if (is.na(fs$joint) || is.na(fs$full))
		return(empty)
	j <- read_sumstat(fs$joint, joint = TRUE)
	if (!nrow(j))
		return(empty)
	q <- read_sumstat_snps(fs$full, j$SNP)
	if (nrow(q) < nrow(j)) {
		jm <- tryCatch(read_sumstat_matched(fs$joint, fs$full, joint = TRUE), error = function(e) tibble())
		if (nrow(jm)) {
			j <- jm
			q <- read_sumstat_snps(fs$full, j$SNP)
		}
	}
	# Alleles, BETA, SE, P and N all come from the same marginal QTL record.
	score <- recover_qtl_alleles(j,fs$full)
	iv <- q[q$variant_key %in% score$variant_key & !is.na(q$EA) & nzchar(q$EA) & !is.na(q$NEA) & nzchar(q$NEA), , drop = FALSE]
	if (!nrow(iv))
		return(empty)
	classify <- function(z) {
		if (layer == "protein") {
			a <- if (!is.null(annotation))
				annotation[annotation$feature == feature, , drop = FALSE] else NULL
			if (!is.null(a) && nrow(a) && all(is.finite(c(a$start[1], a$end[1]))) && !is.na(a$chr[1])) {
				pad <- le8_num_env("C2_CIS_WINDOW_BP", 1e+06)
				z$analysis <- ifelse(as.character(z$CHR) == as.character(a$chr[1]) & z$POS >= a$start[1] - pad &
					z$POS <= a$end[1] + pad, "cis", "trans")
				z$cis_annotation_status <- "gene coordinate defined"
			} else {
				z$analysis <- "unknown"
				z$cis_annotation_status <- "missing annotation; never inferred from lead SNP"
			}
		} else {
			lead <- iv[which.min(iv$P), , drop = FALSE]
			z$analysis <- ifelse(z$CHR == lead$CHR & abs(z$POS - lead$POS) <= le8_num_env(
				"C2_LOCAL_WINDOW_BP",
				1e+06
			), "local", "distal")
			z$cis_annotation_status <- "metabolite local locus; not gene cis"
		}
		z
	}
	iv <- classify(iv)
	iv$effect_type <- "marginal"
	score <- classify(score)
	score$effect_type <- "COJO_joint_PGS_only"
	list(instruments = iv, score_instruments = score, files = fs)
}

run_mr <- function(iv, ygwas, exposure, analysis) {
	# Also repairs reverse-MR callers that passed a raw COJO disease set.
	if (nrow(iv) && "joint" %in% names(iv) && any(iv$joint %in% TRUE)) {
		sources <- unique(iv$source_file)
		sources <- sources[!is.na(sources)]
		full <- if (length(sources))
			sub("\\.jma\\.cojo$", ".gz", sources[1]) else NA_character_
		marginal <- if (!is.na(full) && file.exists(full))
			read_sumstat_snps(full, iv$SNP) else tibble()
		if (!nrow(marginal)) {
			warning("MR ", exposure, ": joint effects cannot be repaired; no estimate", call. = FALSE)
			iv <- iv[0, , drop = FALSE]
		} else iv <- marginal
	}
	# Harmonize before LD pruning/fallback, so a missing top outcome SNP does not hide another valid
	# overlapping instrument.
	if (nrow(iv)) {
		h <- harmonize_sumstats(iv, ygwas)
		iv <- iv[le8_variant_key(iv) %in% h$variant_key, , drop = FALSE]
	}
	pp <- le8_prepare_mr_iv(iv, exposure)
	if (analysis == "unknown")
		pp$iv <- pp$iv[0, , drop = FALSE]
	ans <- le8_fit_mr(pp$iv, ygwas, exposure, analysis)
	ans$effect_type <- "marginal/marginal"
	ans$LD_status <- pp$status
	ans$LD_source <- pp$source
	ans$n_IV_before_LD <- pp$n_input
	ans$n_without_LD <- pp$n_without_LD %||% NA_integer_
	ans$max_retained_r2 <- pp$max_retained_r2 %||% NA_real_
	ans$steiger_method <- "per-IV t-statistic heuristic; not formal liability-scale Steiger"
	ans$steiger_formal_tested <- FALSE
	ans$steiger_r2_sum_exposure_audit <- ans$steiger_r2_exposure
	ans$steiger_r2_sum_outcome_audit <- ans$steiger_r2_outcome
	ans$steiger_r2_exposure <- NA_real_
	ans$steiger_r2_outcome <- NA_real_
	ans$steiger_support <- NA # Do not let a heuristic masquerade as a formal test.
	ans$tsmr_verified <- is.finite(ans$tsmr_ivw_b) & is.finite(ans$b) & abs(ans$tsmr_ivw_b - ans$b) < 1e-08 & is.finite(ans$tsmr_ivw_p) &
		abs(ans$tsmr_ivw_p - ans$pval) < 1e-05
	if (nrow(pp$iv)) {
		ans$instrument_snps <- paste(pp$iv$SNP, collapse = ";")
		ans$instrument_chr <- paste(unique(pp$iv$CHR), collapse = ";")
		ans$instrument_positions <- paste(paste(pp$iv$CHR, pp$iv$POS, sep = ":"), collapse = ";")
		ans$instrument_pos_min <- min(pp$iv$POS)
		ans$instrument_pos_max <- max(pp$iv$POS)
	} else {
		ans$instrument_snps <- ""
		ans$instrument_chr <- ""
		ans$instrument_positions <- ""
		ans$instrument_pos_min <- NA_real_
		ans$instrument_pos_max <- NA_real_
	}
	ans
}


# 0.common.R
# Rebuild presentation outputs without entering analysis/cache invalidation paths.
# Core result RDS remains unchanged. C1 companion plots need selected-feature
# data and matching; C4 independently validates/updates its imaging-only cache.
le8_figure_result <- function(path, fields = character()) {
	required_file(path, "result for output regeneration")
	z <- tryCatch(readRDS(path), error = function(e)
		stop("Cannot read result for output regeneration: ", path, ": ", conditionMessage(e), call. = FALSE))
	if (!is.list(z) || length(setdiff(fields, names(z))))
		stop("Incomplete result for output regeneration: ", path,
			"; required fields: ", paste(fields, collapse = ", "),
			call. = FALSE
		)
	z
}

le8_figure_csv <- function(rawdir, file) {
	as_tibble(data.table::fread(required_file(
		file.path(rawdir, file),
		"table for output regeneration"
	), data.table = FALSE, check.names = FALSE))
}

# Older C1 results did not embed this table. Recover it from the same plotting
# cohort and cached analysis metadata when the optional CSV has been removed.
le8_figure_c1_cohort <- function(obj, dat, layer, rawdir) {
	if (!is.null(obj$cohort)) return(obj$cohort)
	path <- file.path(rawdir, "c1.cohort.csv")
	if (file.exists(path)) return(le8_figure_csv(rawdir, "c1.cohort.csv"))
	message("C1: rebuilding missing cohort summary from trajectory inputs and cached metadata")
	meta <- obj$meta
	outcome <- meta$trait
	events <- sum(dat[[paste0(outcome, ".Yt2e")]] == 1, na.rm = TRUE)
	if ((!is.null(meta$N) && nrow(dat) != meta$N) ||
		(!is.null(meta$events) && events != meta$events))
		stop("C1 cohort inputs no longer match the cached sample/event counts; rerun C1 with --replace TRUE", call. = FALSE)
	dat <- add_attained_age_time(dat, outcome)
	audit <- obj$input_feature_annotation_audit
	removed <- meta$le8_covariates_removed %||% character()
	cohort <- tibble(
		layer = layer, N_omics = nrow(dat), incident_events = events,
		prevalent_cases = sum(make_prevalent_status(dat, outcome) == 1, na.rm = TRUE),
		features = length(audit$feature %||% unique(obj$association_adj2$term)),
		annotation_matched = if (is.null(audit$annotation_matched)) NA_integer_ else sum(audit$annotation_matched, na.rm = TRUE),
		annotation_unmatched = if (is.null(audit$annotation_matched)) NA_integer_ else sum(!audit$annotation_matched, na.rm = TRUE),
		attained_age_N = sum(is.finite(dat$.attained_entry) & is.finite(dat$.attained_exit)),
		bi2e_available = sum(is.finite(dat[[paste0(outcome, ".bi2e")]])),
		PGS_matched = meta$pgs_matched %||% NA_integer_,
		PGS_file_signature = meta$pgs_signature %||% NA_character_,
		LE8_components = length(intersect(vars.le8, names(dat))),
		covs_use = meta$covs_use %||% NA_character_,
		LE4_components = length(meta$le4_covariates),
		le8_covariates_removed = paste(removed, collapse = ";"),
		overlap_covariates_removed = paste(meta$overlap_covariates_removed %||% removed, collapse = ";"),
		treatment_covariates = paste(meta$treatment_covariates, collapse = ";")
	)
	write_raw_csv(cohort, "c1.cohort.csv", rawdir)
	cohort
}

# Explicit workbook mappings retain the original sheet names and avoid exporting
# fitted model objects or individual-level data that were never in the workbook.
le8_figure_tables <- function(obj, mapping) {
	lapply(mapping, function(path) {
		z <- obj
		for (key in strsplit(path, "$", fixed = TRUE)[[1]]) z <- z[[key]]
		z %||% tibble()
	})
}

le8_figure_workbook <- function(tables, file, review = list()) {
	# Write once: reopening a 200+ MB workbook just to append review sheets can
	# cost more time and memory than all of the figures combined.
	if (length(review)) {
		names(review) <- substr(paste0("review_", names(review)), 1, 31)
		tables <- c(tables[!names(tables) %in% names(review)], review)
	}
	message("Writing workbook: ", file)
	write_xlsx2(tables, file)
}

le8_restore_outputs <- function(layer, module) {
	outdir <- if (layer == "protein") out.prot else out.met
	rawdir <- le8_job_dir(outdir, module)
	oldwd <- getwd() ; on.exit(setwd(oldwd), add = TRUE) ; setwd2(if (module == "final_prediction") rawdir else outdir)
	cache <- file.path(rawdir, if (module == "final_prediction") "res.rds" else paste0(sub("_.*$", "", module), ".res.rds"))
	fields <- switch(module,
		c1_correlate = c("association", "clusters", "cluster_selection", "enrichment"),
		c2_cause = c("MR", "MR_reverse", "MR_best", "observational"),
		c3_coloc = c("summary", "regional", "variants"),
		c4_connect = c("scan", "membership", "modules", "mediation"),
		final_prediction = c("scores", "prediction", "summary", "score_vs_topN")
	)
	message("Rebuilding ", module, "/", layer, " PNG/XLSX; reading cached result: ", cache)
	obj <- le8_figure_result(cache, c("meta", fields))
	le8_check_options(obj)
	if (!is.null(obj$meta$trait) && !identical(obj$meta$trait, Y))
		stop("Result trait does not match --trait: ", cache, call. = FALSE)
	if (!is.null(obj$meta$layer) && !identical(obj$meta$layer, layer))
		stop("Result layer does not match --biom: ", cache, call. = FALSE)
	renderer <- get(paste0("le8_render_", sub("_.*$", "", module)), mode = "function")
	renderer(obj, layer, outdir, rawdir)
	finalize_outputs(module, outdir)
	invisible(obj)
}

le8_render_c1 <- function(obj, layer, outdir, rawdir) {
	needed <- c("association_adj2", "prevalent_adj2", "birthline_adj2", "pgs_incident", "pgs_prevalent", "pgs_attained_age")
	if (length(setdiff(needed, names(obj)))) stop("C1 result lacks fields needed for figures: ", paste(setdiff(needed, names(obj)), collapse = ", "))
	aliases <- c(
		assoc = "association", assoc_prevalent = "prevalent", assoc_basic = "association_basic",
		assoc_adj2 = "association_adj2", prevalent_basic = "prevalent_basic", prevalent_adj2 = "prevalent_adj2",
		birthline_basic = "birthline_basic", birthline_adj2 = "birthline_adj2",
		assoc_adj2_full_le8 = "association_adj2_full_le8_sensitivity", prevalent_adj2_full_le8 = "prevalent_adj2_full_le8_sensitivity",
		birthline_adj2_full_le8 = "birthline_adj2_full_le8_sensitivity", reverse_adj2 = "reverse_prevalent",
		duration_adj2 = "prevalent_duration", landmark_adj2 = "landmark_incident", risk_window_adj2 = "diagnosis_window_riskset",
		attenuation_sameN = "attenuation", enrich = "enrichment", enrich_prev = "enrichment_prevalent",
		pgs_incident = "pgs_incident", pgs_prevalent = "pgs_prevalent", pgs_attained_age = "pgs_attained_age",
		pgs_incident_same = "pgs_incident_same_omic", pgs_prevalent_same = "pgs_prevalent_same_omic",
		pgs_attained_age_same = "pgs_attained_age_same_omic", pgs_concordance = "pgs_actual_concordance",
		input_feature_audit = "input_feature_annotation_audit"
	)
	list2env(le8_figure_tables(obj, aliases), envir = environment())
	features_all <- input_feature_audit$feature %||% unique(assoc_adj2$term)
	covs_adj2 <- obj$meta$covariates
	if (is.null(covs_adj2)) stop("C1 cache does not record plotting covariates; rerun C1 explicitly first")
	fig_features <- unique(c(
		head(obj$top_features, YY_TOP), obj$clusters$feature,
		assoc$term[order(assoc$p.value)][seq_len(min(nrow(assoc), max(TOP_N, CLUSTER_TOP)))],
		assoc_adj2$term[order(assoc_adj2$p.value)][seq_len(min(nrow(assoc_adj2), GRADIENT_TOP))],
		prevalent_adj2$term[order(prevalent_adj2$p.value)][seq_len(min(nrow(prevalent_adj2), max(GRADIENT_TOP, CLUSTER_TOP)))]
	))
	message("C1: reading selected-feature trajectory inputs; association/PGS scans and enrichment are reused")
	biom <- if (layer == "protein") read_prot() else read_met()
	fig_features <- unique(c(fig_features, assoc_adj2$term[is.finite(assoc_adj2$p.value) & assoc_adj2$p.value < .05 / nrow(assoc_adj2)]))
	missing <- setdiff(fig_features, names(biom))
	if (length(missing)) stop("C1 plotting inputs lack cached features: ", paste(missing, collapse = ", "))
	biom <- biom[, unique(c("eid", fig_features)), drop = FALSE] ; invisible(gc())
	need <- unique(c("eid", "ethnic.c", vars.basic, vars.le8, covs_adj2, "birth_date", "date_attend", "date_lost", "date_death", paste0("fod_icd10_", Y)))
	dat <- read_all(need) |>
		filter_analysis_cohort() |>
		inner_join(biom, by = "eid") |>
		make_outcome(Y)
	rm(biom) ; invisible(gc())
	if (length(setdiff(covs_adj2, names(dat)))) stop("C1 plotting inputs lack recorded covariates")
	cohort <- le8_figure_c1_cohort(obj, dat, layer, rawdir)
	covs_basic <- intersect(vars.basic, names(dat))
	tvar <- paste0(Y, ".t2e") ; evar <- paste0(Y, ".Yt2e") ; bvar <- paste0(Y, ".b2e")
	pgs <- list(status = obj$pgs_status)
	vldl_deep <- if (layer == "metabolite") build_vldl_tg_deep_dive(
		list(incident = pgs_incident, prevalent = pgs_prevalent, attained_age = pgs_attained_age),
		list(incident = pgs_incident_same, prevalent = pgs_prevalent_same, attained_age = pgs_attained_age_same),
		list(incident = assoc_adj2, prevalent = prevalent_adj2, attained_age = birthline_adj2),
		list(incident = assoc_basic, prevalent = prevalent_basic, attained_age = birthline_basic),
		list(incident = assoc_adj2_full_le8, prevalent = prevalent_adj2_full_le8, attained_age = birthline_adj2_full_le8),
		risk_window_adj2, obj$L_VLDL_TG_pct_pgs_conditional, obj$L_VLDL_TG_pct_measured_conditional
	) else
		list(data = tibble(), trajectory = tibble(), pgs_conditional = tibble(), measured_conditional = tibble())
	top <- assoc |>
		filter(is.finite(p.value)) |>
		slice_min(p.value, n = TOP_N, with_ties = FALSE) |>
		pull(term)
	top6 <- head(top, YY_TOP) ; top_inc <- assoc_adj2 |>
		filter(is.finite(p.value)) |>
		slice_min(p.value, n = GRADIENT_TOP, with_ties = FALSE) |>
		pull(term)
	top_prev <- prevalent_adj2 |>
		filter(is.finite(p.value)) |>
		slice_min(p.value, n = GRADIENT_TOP, with_ties = FALSE) |>
		pull(term)
	cluster_features <- obj$clusters |>
		filter(analysis == "incident_behavioral_LE4") |>
		pull(feature)

	# Fig1 + Fig2: the third row is an attained-age sensitivity analysis with
	# delayed entry at baseline. The left column is a fixed-at-conception omic
	# PGS in the full genetic cohort; the right column is the measured adult omic
	# level with behavioral-LE4 adjustment. Basic observed-level analyses remain in the
	# workbook and raw files but are no longer displayed here.
	plot_c1_fig12(
		layer, features_all, pgs_incident, pgs_prevalent, pgs_attained_age,
		assoc_adj2, prevalent_adj2, birthline_adj2, outdir
	)

	# Fig3: same participants in basic and behavioral-LE4 residualized YY panels.
	yy_need <- unique(c(bvar, top6, covs_adj2)) ; yy_dat <- dat[complete.cases(dat[, intersect(yy_need, names(dat)), drop = FALSE]) & is.finite(dat[[bvar]]), intersect(yy_need, names(dat)), drop = FALSE]
	omic_name <- ifelse(layer == "protein", "protein", "metabolite")
	p3a <- make_yy_panel(yy_dat, top6, bvar, covs_basic, paste0("a. Yin–Yang ", omic_name, " patterns — basic adjustment"),
		ifelse(layer == "protein", "Mean protein z-score", "Mean metabolite z-score"),
		smooth_lines = FALSE
	)
	adj2_label <- if (le8_custom_adjustment()) paste0("adjusted for ", paste(le8_custom_covars, collapse = ", ")) else "basic + behavioral LE4 adjustment"
	p3b <- make_yy_panel(yy_dat, top6, bvar, covs_adj2, paste0("b. Yin–Yang trajectories — ", adj2_label),
		ifelse(layer == "protein", "Mean protein z-score", "Mean metabolite z-score"),
		smooth_lines = TRUE
	)
	save_plot(p3a / p3b + plot_layout(heights = c(1, 1)), "c1.Fig3.yy_top.png", 16, 13.5, outdir = outdir)

	# Fig4 shows adjusted diagnosis-timed cross-sectional omics gradients.
	p4i <- make_gradient_panel(dat, top_inc, bvar, "incident", GRADIENT_STEP, MIN_BIN_N,
		paste0("a. Incident (Yin): top ", length(top_inc), if (le8_custom_adjustment()) " adjusted biomarkers" else " behavioral-LE4 biomarkers"),
		covars = covs_adj2
	)
	p4p <- make_gradient_panel(dat, top_prev, bvar, "prevalent", GRADIENT_STEP, MIN_BIN_N,
		paste0("b. Baseline-prevalent (Yang): top ", length(top_prev), if (le8_custom_adjustment()) " adjusted biomarkers" else " behavioral-LE4 biomarkers"),
		covars = covs_adj2
	)
	save_plot((p4i / p4p + plot_layout(guides = "collect")) & theme(legend.position = "right"),
		"c1.Fig4.gradient.png", 17, 9.6,
		outdir = outdir
	)
	gradient_rank <- bind_rows(
		assoc_adj2 |> filter(term %in% top_inc) |> transmute(analysis = "incident_behavioral_LE4", feature = term, p_value = p.value) |> arrange(p_value) |> mutate(rank = row_number()),
		prevalent_adj2 |> filter(term %in% top_prev) |> transmute(analysis = "prevalent_behavioral_LE4", feature = term, p_value = p.value) |> arrange(p_value) |> mutate(rank = row_number())
	)

	# Fig5: separately selected and clustered incident and prevalent adj2 sets.
	cluster_prev <- obj$clusters |>
		filter(analysis == "prevalent_behavioral_LE4") |>
		pull(feature)
	cli <- make_cluster_figure(dat, cluster_features, bvar, CLUSTER_STEP, YY_MAX_YEAR, "incident", "a",
		saved_membership = filter(obj$clusters, analysis == "incident_behavioral_LE4"),
		saved_metrics = filter(obj$cluster_selection, analysis == "incident_behavioral_LE4")
	)
	clp <- make_cluster_figure(dat, cluster_prev, bvar, CLUSTER_STEP, YY_MAX_YEAR, "prevalent", "b",
		saved_membership = filter(obj$clusters, analysis == "prevalent_behavioral_LE4"),
		saved_metrics = filter(obj$cluster_selection, analysis == "prevalent_behavioral_LE4")
	)
	# Keep the five top-level patchwork rows explicit.  Nesting the three-row
	# cluster_main object here is flattened by patchwork, so a three-entry outer
	# heights vector leaves room for only three of the resulting five panels.
	cluster_figure <- cli$cluster_plot / plot_spacer() / clp$cluster_plot / plot_spacer() /
		(cli$diagnostic_plot | clp$diagnostic_plot) +
		plot_layout(heights = c(1, .10, 1, .05, .38))
	save_plot(cluster_figure,
		"c1.Fig5.gradient_cluster.png", 18, 10.5,
		outdir = outdir
	)
	cl_members <- bind_rows(cli$cluster |> mutate(analysis = "incident_behavioral_LE4"), clp$cluster |> mutate(analysis = "prevalent_behavioral_LE4"))
	cl_metrics <- bind_rows(cli$metrics |> mutate(analysis = "incident_behavioral_LE4"), clp$metrics |> mutate(analysis = "prevalent_behavioral_LE4"))

	# Fig6: adjusted Q5-vs-Q1 HRs on the same adj2 complete-case sample; no subtitle.
	qplots <- plot_quantile_top(dat, top6, tvar, evar, covs_adj2, paste0("Adjusted for ", paste(covs_adj2, collapse = ", ")))
	save_plot(wrap_plots(qplots, ncol = 3), "c1.Fig6.quantile_top.png", 16, 10, outdir = outdir)

	# Fig8: functional coherence of incident and baseline-prevalent adj2 scans.
	n_sig <- sum(is.finite(assoc_adj2$p.value) & assoc_adj2$p.value * nrow(assoc_adj2) < .05)
	n_sig_prev <- sum(is.finite(prevalent_adj2$p.value) & prevalent_adj2$p.value * nrow(prevalent_adj2) < .05)
	pe_i <- plot_functional_enrichment(enrich, n_sig, assoc_adj2, layer, if (le8_custom_adjustment()) "a. Incident (Yin), selected covariates" else "a. Incident (Yin), behavioral LE4")
	pe_p <- plot_functional_enrichment(enrich_prev, n_sig_prev, prevalent_adj2, layer, if (le8_custom_adjustment()) "b. Baseline-prevalent (Yang), selected covariates" else "b. Baseline-prevalent (Yang), behavioral LE4")
	if (layer != "protein") save_plot((pe_i / pe_p + plot_layout(guides = "collect")) & theme(legend.position = "right"),
		"c1.Fig8.enrich_sig.png", 18, 10.5,
		outdir = outdir
	)

	# Fig9 separates distal antecedent prediction from disease-state evidence.
	# GDF15 and PCSK9 are anchors by default, but the table is generated for all
	# assayed proteins/metabolites.
	directionality <- build_directionality_table(assoc_adj2, prevalent_adj2, duration_adj2, landmark_adj2, birthline_adj2, reverse_adj2)
	save_plot(plot_directionality_triage(directionality, C1_DIRECTION_ANCHORS),
		"c1.Fig9.directionality_triage.png", 16, 9,
		outdir = outdir
	)
	save_plot(plot_landmark_and_attained_age(landmark_adj2, assoc_adj2, birthline_adj2, C1_DIRECTION_ANCHORS),
		"c1.Fig10.landmark_birthline_sensitivity.png", 16, 9,
		outdir = outdir
	)
	save_plot(plot_risk_window_scan(risk_window_adj2, C1_DIRECTION_ANCHORS),
		"c1.Fig11.diagnosis_window_riskset.png", 18, 11,
		outdir = outdir
	)
	save_plot(plot_directionality_supplement(directionality, C1_DIRECTION_ANCHORS),
		"c1.Fig12.directionality_detail.png", 14, 7.5,
		outdir = outdir
	)
	save_plot(plot_reverse_time_exploratory(reverse_adj2, prevalent_adj2, duration_adj2, C1_DIRECTION_ANCHORS),
		"c1.Fig13.reverse_time_exploratory.png", 15.5, 8,
		outdir = outdir
	)
	save_plot(plot_pgs_actual_concordance(pgs_concordance, C1_DIRECTION_ANCHORS),
		"c1.Fig14.pgs_actual_concordance.png", 18, 10.5,
		outdir = outdir
	)
	if (layer == "metabolite") save_plot(vldl_deep$figure,
		"c1.Fig15.L_VLDL_TG_pct_deep_dive.png", 19, 26,
		outdir = outdir
	)

	le8_mock_c1(dat, assoc_adj2, enrich, layer, covs_adj2, tvar, evar, outdir, enrich_prev)
	gradient_rank <- obj$gradient_top10_provenance ; cl_members <- obj$clusters ; cl_metrics <- obj$cluster_selection
	le8_figure_workbook(list(
		cohort = cohort, input_feature_audit = input_feature_audit, association = assoc, prevalent = assoc_prevalent, incident_basic = assoc_basic, incident_adj2 = assoc_adj2,
		birthline_basic = birthline_basic, birthline_adj2 = birthline_adj2,
		prevalent_basic = prevalent_basic, prevalent_adj2 = prevalent_adj2,
		association_adj2_full_le8_sensitivity = assoc_adj2_full_le8,
		prevalent_adj2_full_le8_sensitivity = prevalent_adj2_full_le8,
		birthline_adj2_full_le8_sensitivity = birthline_adj2_full_le8,
		attenuation_sameN = attenuation_sameN,
		reverse_prevalent = reverse_adj2, prevalent_duration = duration_adj2, incident_landmark = landmark_adj2,
		diagnosis_window_riskset = risk_window_adj2, directionality = directionality,
		pgs_status = pgs$status, pgs_incident = pgs_incident, pgs_prevalent = pgs_prevalent,
		pgs_attained_age = pgs_attained_age, pgs_incident_same_omic = pgs_incident_same,
		pgs_prevalent_same_omic = pgs_prevalent_same, pgs_attained_same_omic = pgs_attained_age_same,
		pgs_actual_concordance = pgs_concordance,
		L_VLDL_TG_pct_deep_dive = vldl_deep$data,
		L_VLDL_TG_pct_riskset_trajectory = vldl_deep$trajectory,
		L_VLDL_TG_pct_pgs_conditional = vldl_deep$pgs_conditional,
		L_VLDL_TG_pct_measured_conditional = vldl_deep$measured_conditional,
		gradient_top10_provenance = gradient_rank,
		cluster_membership = cl_members, cluster_selection = cl_metrics,
		enrichment_incident = enrich, enrichment_prevalent = enrich_prev
	), "c1.out.xlsx", obj$review %||% list())
	review <- obj$review %||% list()
	if (length(review)) le8_plot_c1_temporal(
		review$paired_associations %||% tibble(),
		review$time_heterogeneity %||% tibble(), review$paired_status %||% tibble(status = "Unavailable in saved result"),
		features_all, layer
	)
}

le8_render_c2 <- function(obj, layer, outdir, rawdir) {
	aliases <- c(
		mr = "MR", reverse_mr = "MR_reverse", reverse_audit = "MR_reverse_audit", assoc = "observational",
		method_scope = "method_scope", evidence_grades = "evidence_grades", top_candidates = "top_candidates",
		genetic_score_manifest = "genetic_score_manifest", best = "MR_best", availability = "instrument_availability",
		top_r2 = "R2_QTL", arch = "architecture", directionality_integration = "directionality_integration",
		jobs_audit = "MRLink2_job_audit", jobs = "MRLink2_jobs"
	)
	list2env(le8_figure_tables(obj, aliases), envir = environment())
	individual_decomposition <- obj$individual_decomposition %||% list()
	dandelion <- obj$DANDELION %||% list() ; mrlink2 <- obj$MRLink2 %||% read_mrlink2_results(rawdir)
	fig1 <- plot_c2_fig1(mr, assoc, layer) ; save_c2_plot(fig1$plot, "c2.Fig1.prots.top.png", 17, 13.5, outdir = outdir)
	save_c2_plot(plot_qtl_variance(mr, layer), "c2.Fig2.pQTL_R2.png", 15.5, 10.5, outdir = outdir)
	fig3 <- plot_c2_fig2(mr, assoc, layer) ; save_plot(fig3$plot, "c2.Fig3.effect_concordance.png", 16, 12.5, outdir = outdir)
	write_raw_csv(fig3$wide, "c2.cis_trans_comparison.csv", le8_job_dir(outdir, "c2_cause"))
	le8_emit_restored_c2(mr, assoc, layer, outdir)
	save_plot(plot_c2_fig4(mr, layer), "c2.Fig4.sensitivity_architecture.png", 15.5, 11.5, outdir = outdir)
	plot_mrlink2_results(mrlink2, outdir)
	if (layer == "protein") {
		le8_dandelion_plot_bundle(dandelion, outdir)
		plot_dandelion_mr_integration(dandelion, mr, assoc, outdir)
	} else {
		suffix <- c("dandelion", "dandelion_evidence", "dandelion_mr_integration", "dandelion_network")
		for (i in seq_along(suffix)) save_plot(
			blank_plot(
				paste0("C2 Figure ", i + 5),
				"DANDELION is a gene/protein regulatory-network analysis and is not defined for metabolites"
			),
			paste0("c2.Fig", i + 5, ".", suffix[i], ".png"), 10, 6,
			outdir = outdir
		)
	}
	save_plot(plot_c2_directionality(directionality_integration), "c2.Fig10.directionality_causal.png", 16, 12, outdir = outdir)
	save_plot(plot_bidirectional_mr(mr, reverse_mr, layer), "c2.Fig11.bidirectional_mr.png", 16, 8.5, outdir = outdir)
	save_plot(plot_individual_genetic_decomposition(individual_decomposition$summary %||% tibble()), "c2.Fig12.genetic_decomposition.png", 17, 13, outdir = outdir)
	save_plot(plot_component_leadtime(individual_decomposition$trajectory %||% tibble()), "c2.Fig13.genetic_component_leadtime.png", 17, 11, outdir = outdir)
	save_plot(plot_c2_evidence_grades(evidence_grades, layer), "c2.Fig14.evidence_grades.png", 17, 9, outdir = outdir)
	le8_figure_workbook(list(
		MR_all = mr, MR_reverse = reverse_mr, MR_reverse_audit = reverse_audit,
		method_scope = method_scope, evidence_grades = evidence_grades, top_candidates = top_candidates,
		genetic_score_manifest = genetic_score_manifest,
		genetic_decomp_status = individual_decomposition$status %||% tibble(),
		heritability_status = individual_decomposition$heritability_status %||% tibble(),
		genetic_decomp_summary = individual_decomposition$summary %||% tibble(),
		genetic_leadtime = individual_decomposition$trajectory %||% tibble(),
		MR_best = best, instrument_availability = availability, evidence_overlap = fig1$bars, effect_forest = fig1$forest, cis_local_vs_trans_distal = fig3$wide, QTL_R2 = top_r2, instrument_architecture = arch,
		DANDELION_input_audit = dandelion$input_audit %||% tibble(), DANDELION_QTL_coverage = dandelion$qtl_audit %||% tibble(), DANDELION_exposure_QC = dandelion$exposure_qc %||% tibble(),
		DANDELION_lead_snps = dandelion$lead_snps %||% tibble(), DANDELION_snp_gene_map = dandelion$snp_gene_map %||% tibble(), DANDELION_pairs = dandelion$pairs %||% tibble(), DANDELION_gene_pairs = dandelion$gene_pairs %||% tibble(), DANDELION_targets = dandelion$targets %||% tibble(), DANDELION_MR_integration = dandelion$integration %||% tibble(),
		directionality_causal = directionality_integration,
		MRLink2_job_audit = jobs_audit, MRLink2_jobs = jobs,
		MRLink2_results = mrlink2$results %||% tibble(),
		MRLink2_status = mrlink2$status %||% tibble()
	), "c2.out.xlsx", obj$review %||% list())
}

le8_render_c3 <- function(obj, layer, outdir, rawdir) {
	# Cell-type analyses are owned exclusively by c5_cellulation.
	plot_coloc_results(obj$summary, obj$regional, obj$variants, layer, outdir)
	gpu <- obj$GPU_coloc %||% read_gpu_coloc_results(rawdir)
	plot_gpu_coloc_validation(gpu, obj$summary, outdir)
	aud <- obj$credible_set_audit %||% credible_set_audit(obj$summary, obj$variants)
	tri <- obj$pgs_triangulation %||% read_c3_pgs_integration(layer, outdir, obj$summary)
	plot_c3_pgs_integration(tri, outdir)
	lists <- obj$causal_lists %||% list()
	le8_figure_workbook(list(
		coloc_summary = obj$summary, credible_set_audit = aud$overall, credible_set_by_locus = aud$by_locus,
		GPU_results = gpu$results, GPU_status = gpu$status,
		GPU_manifest = obj$manifest, causal_sets = if (length(lists)) stack(lists) else tibble(), pgs_observed_coloc = tri
	), "c3.out.xlsx", obj$review %||% list())
}

le8_render_c4 <- function(obj, layer, outdir, rawdir) {
	if (nrow(obj$modules$membership %||% tibble()) && nrow(obj$modules$metrics %||% tibble())) {
		obj$modules$metrics$selected_k <- n_distinct(obj$modules$membership$module)
		write_raw_csv(obj$modules$metrics, "c4.supervised_module_selection.csv", rawdir)
	}
	sets <- list(primary = obj$primary, membership = obj$membership, YS_edges = obj$YS_edges)
	if (identical(obj$meta$status, "unavailable")) {
		suffix <- c(
			"proxy_heatmap", "group_pillar_flow", "connection_bridge", "connection_evidence", "mediation_forest",
			"mediation_diagnostics", "supervised_atlas", "network_globe", "selection_mediation", "state_network_remodeling", "state_network_edges"
		)
		if (layer == "metabolite") suffix[3 : 5] <- c("module_network", "relationship_globe", "mediation_wheel")
		reason <- paste(obj$status$detail %||% "Unavailable in saved result", collapse = "; ")
		for (i in seq_along(suffix)) save_plot(blank_plot(paste0("C4 Figure ", i), reason), paste0("c4.Fig", i, ".", suffix[i], ".png"), 10, 6, outdir = outdir)
	} else {
		plot_c4(obj$scan, sets, obj$mediation, obj$genetic_edges, layer, outdir)
		plot_supervised_atlas(obj$modules, sets, obj$mediation, layer, outdir)
		plot_module_globe(obj$modules, outdir)
		plot_c4_state_network(obj$state_network, layer, outdir)
	}
	tables <- le8_figure_tables(obj, c(
		primary_assignment = "primary", proxy_membership = "membership", YS_edges = "YS_edges",
		supervised_modules = "modules$membership", module_selection = "modules$metrics", all_LE8_associations = "scan",
		PRS_associations = "genetic_scan", genetic_omic_bridges = "genetic_edges", mediation = "mediation",
		state_network_status = "state_network$status", state_counts = "state_network$state_counts",
		state_network_edges = "state_network$edges", state_network_hubs = "state_network$hubs"
	))
	if (!is.null(obj$status)) tables$status <- obj$status
	# Reuse C4 core fits. The separate imaging cache checks its own input/option signature.
	imaging <- run_c4_imaging(layer, outdir = outdir)
	if (!is.null(imaging)) {
		associations <- imaging$associations %||% tibble()
		if (nrow(associations)) plot_c4_imaging(associations, unique(associations$feature), outdir)
		tables <- c(tables, le8_figure_tables(list(imaging = imaging), c(imaging_status = "imaging$status", imaging_associations = "imaging$associations", imaging_fields = "imaging$fields")))
	}
	le8_figure_workbook(tables, "c4.out.xlsx", obj$review %||% list())
}

le8_render_final <- function(obj, layer, outdir, rawdir) {
	le8_mock_final(obj, outdir)
	all_glmnet_label <- if (layer == "protein") "Pradeep-style / glmnet" else "All-metabolite / glmnet"
	lightgbm_label <- if (layer == "protein") "Yu-style / LightGBM" else "MWAS-ranked / LightGBM"
	# The final model cache can predate prevalent-row augmentation. The exported
	# person_scores table is the input actually used for the completed figures.
	score_rows <- le8_figure_csv(rawdir, "person_scores.csv")
	pred <- obj$prediction ; psum <- obj$summary
	pairs_native <- obj$preclinical_pairs %||% list() ; boot_by_method <- obj$yy_auc_bootstrap %||% list()
	orders <- list(
		c(all_glmnet_label, lightgbm_label, "User specified"),
		c("Distal antecedent", "Genetic-region evidence", "Hybrid triangulated", "Evidence-selected compact"),
		c("C4 NS", "C4 YS", "C4 YSplus")
	)
	titles <- c("Prediction paradigms", "Evidence-screened prediction", "C4 connection-guided prediction")
	files <- c("Fig1.simple_pred.png", "Fig2.screened_pred.png", "Fig3.connection_pred.png")
	for (i in seq_along(orders)) save_plot(final_make_5row_grid(orders[[i]], pred, score_rows, pairs_native, boot_by_method, titles[i]), files[i], 24, 13.5, outdir = outdir)
	save_plot(plot_leadtime_prediction(obj$leadtime, obj$leadtime_windows), "Fig4.leadtime_prediction.png", 20, 11.25, outdir = outdir)
	save_plot(plot_topn_mechanism(obj$score_vs_topN), "Fig5.score_vs_topN_mechanism.png", 19, 10.7, outdir = outdir)
	save_plot(plot_evidence_matrix(obj$evidence) | plot_causal_reactive_map(obj$evidence), "Fig6.evidence_matrix.png", 18, 10.5, outdir = outdir)
	save_plot(plot_performance_benchmark(psum), "Fig7.performance_benchmark.png", 13, 9, outdir = outdir)
	save_plot(plot_incremental_performance(psum), "Fig8.incremental_performance.png", 13, 9, outdir = outdir)
	save_plot(plot_complexity_performance(psum), "Fig9.parsimony_performance.png", 10.5, 8, outdir = outdir)
	save_plot(plot_score_correlation(obj$score_correlations), "Fig10.score_concordance.png", 11, 9.5, outdir = outdir)
	save_plot(plot_subgroup_auc(obj$subgroup_AUC), "Fig11.subgroup_discrimination.png", 12, 10, outdir = outdir)
	save_plot(plot_attained_age_score_sensitivity(obj$attained_age_sensitivity), "Fig12.attained_age_sensitivity.png", 14, 8.5, outdir = outdir)
	review <- obj$review %||% list()
	if (nrow(review$budget_metrics %||% tibble())) {
		le8_plot_final_review(
			review$budget_metrics, review$calibration, review$decision_curves,
			unique(review$budget_metrics$horizon), unique(review$budget_metrics$budget), outdir
		)
	}
	tables <- le8_figure_tables(obj, c(
		upstream_availability = "upstream_availability", prediction_summary = "summary",
		input_set_sizes = "set_sizes", candidate_sets = "candidate_sets", mechanism_weight_prior = "mechanism_weight_prior",
		distal_landmark_screen = "distal_screen", score_coverage_audit = "score_coverage_audit", sequential_forward = "sequential$log",
		leadtime_discrimination = "leadtime", leadtime_headline = "leadtime_headline", leadtime_windows = "leadtime_windows",
		score_vs_topN = "score_vs_topN$summary", topN_individuals = "score_vs_topN$individuals", topN_individual_range = "score_vs_topN$individual_range",
		attained_age_sensitivity = "attained_age_sensitivity", evidence_consolidation = "evidence", score_correlations = "score_correlations",
		subgroup_AUC = "subgroup_AUC", nested_cv_summary = "nested_cv$summary"
	))
	tables$training_screen <- le8_figure_csv(rawdir, "training_association_screen.csv")
	tables$score_method_summary <- score_rows |>
		group_by(method, kind, split) |>
		summarise(rows = n(), finite_scores = sum(is.finite(score_z)), mean_score = mean(score_z, na.rm = TRUE), sd_score = sd(score_z, na.rm = TRUE), .groups = "drop")
	tables$preclinical_pairs <- bind_rows(imap(pairs_native, ~ .x |> mutate(method = .y)))
	tables$yy_auc_bootstrap <- bind_rows(imap(boot_by_method, ~ .x |> mutate(method = .y)))
	le8_figure_workbook(tables, "prediction_panels.xlsx", obj$review %||% list())
}


# 0.common.R
# Shared user-editable assay budgets for C4 focus and newly fitted Final reviews.
LE8_ASSAY_BUDGETS <- c(5, 10, 50)


# 0.common.R
# Preserve scientific views while grouping panels and retaining one final
# PNG/workbook pair per page. Temporary source layouts stay under /tmp.
if (!exists(".le8_figure_queue", inherits = FALSE)) .le8_figure_queue <- new.env(parent = emptyenv())

le8_figure_rule <- function(file) {
	key <- sub("[.]png$", "", sub("^[^.]+[.]Fig[0-9]+[.]", "", basename(file)))
	module <- if (grepl("^Fig", basename(file))) "final" else sub("[.].*$", "", basename(file))
	if (module == "final") key <- sub("^Fig[0-9]+[.]", "", sub("[.]png$", "", basename(file)))
	omit <- character()
	group <- key
	groups <- switch(module,
		c1 = list(temporal_profiles = c("yy_top", "gradient"), temporal_evidence = c("directionality_triage", "directionality_detail")),
		c2 = list(instrument_diagnostics = c("pQTL_R2"), dandelion_sensitivity = c("dandelion", "dandelion_evidence")),
		c3 = list(posterior_evidence = c("evidence_triage", "posterior_diagnostics", "prior_sensitivity")),
		c4 = list(supervised_connections = c("group_pillar_flow", "supervised_atlas"), mediation = c("mediation_forest", "mediation_diagnostics"), state_remodeling = c("state_network_remodeling", "state_network_edges"), imaging_context = c("imaging_associations", "imaging_overview")),
		final = list(prediction_benchmark = c("performance_benchmark", "incremental_performance", "parsimony_performance"), prediction_sensitivity = c("subgroup_discrimination", "attained_age_sensitivity")),
		list()
	)
	for (nm in names(groups)) if (key %in% groups[[nm]]) group <- nm
	list(module = module, key = key, group = group, omit = key %in% omit)
}

le8_queue_figure <- function(p, target, w, h, dpi) {
	rd <- dirname(target) ; q <- .le8_figure_queue[[rd]]
	if (is.null(q)) q <- list()
	rule <- le8_figure_rule(target)
	q[[basename(target)]] <- list(plot = if (rule$omit) NULL else p, rule = rule, width = w, height = h, dpi = dpi)
	.le8_figure_queue[[rd]] <- q
	invisible(target)
}

le8_plot_leaves <- function(p) {
	if (is.null(p)) return(list())
	if (inherits(p, "patchwork")) return(unlist(lapply(seq_len(length(p)), function(i) le8_plot_leaves(p[[i]])), recursive = FALSE))
	if (inherits(p, "spacer")) return(list())
	list(p)
}

le8_clone_layer <- function(x) {
	force(x) ; ggplot2::ggproto(NULL, x)
}

# Shared themes precede layer-specific additions. Presentation numbers are
# independent of legacy source-file numbers and remain consecutive.
le8_figure_order <- function(module, groups) {
	preferred <- switch(module,
		c1 = c(
			'mh', 'circular', 'vc', 'temporal_profiles', 'diagnosis_timed_profiles',
			'gradient_cluster', 'quantile_top', 'enrich_sig', 'temporal_evidence',
			'landmark_birthline_sensitivity', 'diagnosis_window_riskset', 'paired_temporal_validation'
		),
		c2 = c(
			'instrument_diagnostics', 'effect_concordance', 'mr_incident_prevalent', 'mrlink2', 'bidirectional_mr',
			'genetic_decomposition', 'genetic_component_leadtime', 'evidence_grades',
			'dandelion_sensitivity', 'dandelion_mr_integration'
		),
		c3 = c('posterior_evidence', 'regional_top_loci', 'credible_sets', 'gpu_coloc_validation', 'pgs_coloc_triangulation'),
		c4 = c('supervised_connections', 'mediation', 'imaging_context', 'lifestyle_omics_risk', 'state_remodeling', 'sex_interaction', 'LE8_component_interactions', 'spline_patterns', 'pass_fail_penalty'),
		final = c(
			'prediction_benchmark', 'score_concordance', 'prediction_sensitivity',
			'budget_ablation_calibration', 'heldout_ROC', 'leadtime_prediction'
		),
		character()
	)
	order(match(groups, preferred, nomatch = length(preferred) + 1L), seq_along(groups))
}

# Merge freshly rendered themes with existing PNGs, then renumber via the same
# policy as a full module run. This permits aggregate-only presentation updates.
le8_refresh_figure_files <- function(rawdir, incoming = NULL) {
	mf <- file.path(rawdir, 'figure_manifest.csv')
	old <- if (file.exists(mf)) as.data.frame(data.table::fread(mf)) else data.frame()
	fresh <- if (!is.null(incoming)) as.data.frame(data.table::fread(file.path(incoming, 'figure_manifest.csv'))) else old[FALSE, ]
	normalize_group <- function(x) {
		if (nrow(x)) x$group <- sub('^(c[1-5][.])?Fig[0-9]+[.]', '', x$group)
		x
	}
	old <- normalize_group(old) ; fresh <- normalize_group(fresh)
	retained <- if (nrow(old)) old[!old$group %in% fresh$group & file.exists(file.path(rawdir, old$file)), , drop = FALSE] else old
	# PGS manifests use `source`; grouped figures add `sources` and `policy`.
	# Align by name and retain all metadata when refreshing mixed/older outputs.
	revised <- as.data.frame(data.table::rbindlist(list(retained, fresh), use.names = TRUE, fill = TRUE))
	if (!nrow(revised)) {
		data.table::fwrite(revised, mf)
		if (!is.null(incoming)) stopifnot(file.copy(file.path(incoming, 'figure_omission_audit.csv'),
			file.path(rawdir, 'figure_omission_audit.csv'),
			overwrite = TRUE
		))
		return(invisible(revised))
	}
	sources <- c(file.path(rawdir, retained$file), if (nrow(fresh)) file.path(incoming, fresh$file) else character())
	stopifnot(all(file.exists(sources)))
	# An earlier PGS refresh registered the same PNG under two figure numbers.
	# Deduplicate only byte-identical images in the same theme, retaining pages.
	signatures <- paste(revised$group, vapply(sources, function(path) digest::digest(file = path, algo = 'sha256'), character(1)))
	keep <- !duplicated(signatures) ; revised <- revised[keep, , drop = FALSE] ; sources <- sources[keep]
	prefix <- if (grepl('^Fig', revised$file[1])) 'final' else sub('[.].*$', '', revised$file[1])
	ix <- le8_figure_order(prefix, revised$group) ; revised <- revised[ix, , drop = FALSE] ; sources <- sources[ix]
	revised$file <- paste0(if (prefix == 'final') '' else paste0(prefix, '.'), 'Fig', seq_len(nrow(revised)), '.', revised$group, '.png')
	stage <- tempfile('.renumber-', tmpdir = '/tmp') ; dir.create(stage)
	on.exit(unlink(stage, recursive = TRUE), add = TRUE)
	stopifnot(all(file.copy(sources, file.path(stage, revised$file))))
	old_tables <- character()
	for (i in seq_along(sources)) {
		stem <- sub('[.]png$', '', basename(sources[i]))
		tables <- list.files(dirname(sources[i]), pattern = '[.]csv$', full.names = TRUE)
		tables <- tables[startsWith(basename(tables), paste0(stem, '.panel_'))]
		if (!length(tables)) next
		names <- paste0(sub('[.]png$', '', revised$file[i]), substring(basename(tables), nchar(stem) + 1L))
		stopifnot(all(file.copy(tables, file.path(stage, names), overwrite = TRUE)))
		if (dirname(sources[i]) == rawdir) old_tables <- c(old_tables, tables)
	}
	# Stage every source before overwriting names that may be reused by another theme.
	staged <- list.files(stage, full.names = TRUE)
	stopifnot(all(file.copy(staged, file.path(rawdir, basename(staged)), overwrite = TRUE)))
	unlink(setdiff(old_tables, file.path(rawdir, basename(staged))))
	stale <- if (nrow(old)) setdiff(old$file, revised$file) else character()
	if (length(stale)) unlink(file.path(rawdir, stale))
	data.table::fwrite(revised, mf)
	if (!is.null(incoming)) {
		af <- file.path(rawdir, 'figure_omission_audit.csv')
		a <- if (file.exists(af)) as.data.frame(data.table::fread(af)) else data.frame()
		b <- as.data.frame(data.table::fread(file.path(incoming, 'figure_omission_audit.csv')))
		if (nrow(a)) a <- a[!a$group %in% b$group, , drop = FALSE]
		data.table::fwrite(rbind(a, b), af)
	}
	invisible(revised)
}

le8_layer_figure_correspondence <- function(traitdir) {
	read_layer <- function(layer) {
		files <- list.files(file.path(traitdir, layer), pattern = '^figure_manifest[.]csv$', recursive = TRUE, full.names = TRUE)
		files <- files[!grepl('/_history/', files, fixed = TRUE)]
		ans <- lapply(files, function(f) {
			d <- as.data.frame(data.table::fread(f)) ; if (!nrow(d)) return(NULL)
			if (!'group' %in% names(d)) d$group <- sub('[.]png$', '', d$file)
			d$group[is.na(d$group)] <- sub('[.]png$', '', d$file[is.na(d$group)])
			d$module <- basename(dirname(f)) ; d$page <- ave(seq_len(nrow(d)), d$group, FUN = seq_along)
			# Genomic protein order and metabolite biochemical order share one theme.
			d$group[d$group %in% c('mh', 'circular')] <- 'association_overview'
			d$file <- file.path(d$module, d$file) ; d[, c('module', 'group', 'page', 'file', 'panels')]
		})
		ans <- do.call(rbind, ans)
		if (is.null(ans)) ans <- data.frame(module = character(), group = character(), page = integer(), file = character(), panels = integer())
		names(ans)[4 : 5] <- paste0(layer, c('_file', '_panels')) ; ans
	}
	d <- merge(read_layer('prot'), read_layer('met'), by = c('module', 'group', 'page'), all = TRUE)
	d$status <- ifelse(!is.na(d$prot_file) & !is.na(d$met_file), 'paired theme',
		ifelse(is.na(d$prot_file), 'met only / no usable prot panel', 'prot only / no usable met panel')
	)
	number <- function(x) suppressWarnings(as.integer(sub('.*[.]Fig([0-9]+).*', '\\1', x)))
	d$same_number <- ifelse(d$status == 'paired theme', number(d$prot_file) == number(d$met_file), NA)
	dest <- le8_final_dir(trait = basename(traitdir), root = dirname(traitdir))
	dir.create(dest, recursive = TRUE, showWarnings = FALSE)
	data.table::fwrite(d, file.path(dest, 'figure_layer_correspondence.csv'))
	invisible(d)
}

le8_readable_plot <- function(p) {
	# Limit labels only; every observation remains in points/lines and tables.
	for (i in seq_along(p$layers)) if (inherits(p$layers[[i]]$geom, "GeomTextRepel") || inherits(p$layers[[i]]$geom, "GeomLabelRepel")) {
		layer <- le8_clone_layer(p$layers[[i]])
		d <- if (is.data.frame(layer$data)) layer$data else p$data
		label <- layer$mapping$label %||% p$mapping$label
		if (is.data.frame(d) && !is.null(label)) {
			labs <- tryCatch(rlang::eval_tidy(label, data = d), error = function(e) rep(NA_character_, nrow(d)))
			keep <- head(which(!is.na(labs) & nzchar(as.character(labs))), 12L)
			layer$data <- d[keep, , drop = FALSE]
			layer$aes_params$colour <- "grey20" ; layer$geom_params$max.overlaps <- 12L
			p$layers[[i]] <- layer
		}
	}
	for (nm in c("title", "subtitle", "caption")) if (is.character(p$labels[[nm]]) && length(p$labels[[nm]]) == 1) {
		value <- p$labels[[nm]]
		if (nm == "title") value <- sub("^[A-Za-z][.] +", "", value)
		p <- p + do.call(ggplot2::labs, setNames(list(stringr::str_wrap(value, if (nm == "title") 60 else 85)), nm))
	}
	p + ggplot2::theme(
		plot.title = ggplot2::element_text(size = 11, face = "bold"),
		plot.subtitle = ggplot2::element_text(size = 9), plot.caption = ggplot2::element_text(size = 8, hjust = 0),
		axis.text = ggplot2::element_text(size = 9), axis.title = ggplot2::element_text(size = 10),
		legend.text = ggplot2::element_text(size = 8), legend.title = ggplot2::element_text(size = 9),
		legend.position = "bottom", legend.box = "vertical", legend.box.just = "center",
		plot.margin = ggplot2::margin(10, 15, 10, 12)
	)
}

le8_expand_facets <- function(p) {
	if (!is.null(attr(p, "le8_unavailable"))) return(list())
	if (inherits(p$facet, "FacetNull")) return(list(le8_readable_plot(p)))
	layout <- ggplot2::ggplot_build(p)$layout$layout
	vars <- setdiff(names(layout), c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y", "COORD"))
	if (nrow(layout) <= 1 || !length(vars)) return(list(le8_readable_plot(p)))
	# Each facet becomes an independent axis, so the six-panel limit is real.
	lapply(seq_len(nrow(layout)), function(j) {
		values <- layout[j, vars, drop = FALSE]
		subset_data <- function(d) {
			if (!is.data.frame(d)) return(d)
			common <- intersect(names(d), vars) ; ok <- rep(TRUE, nrow(d))
			for (v in common) ok <- ok & (if (is.na(values[[v]])) is.na(d[[v]]) else !is.na(d[[v]]) & as.character(d[[v]]) == as.character(values[[v]]))
			droplevels(d[ok, , drop = FALSE])
		}
		z <- p ; z$data <- subset_data(z$data)
		for (i in seq_along(z$layers)) {
			layer <- le8_clone_layer(z$layers[[i]]) ; layer$data <- subset_data(layer$data) ; z$layers[[i]] <- layer
		}
		z <- z + ggplot2::facet_null() + ggplot2::labs(title = paste(p$labels$title, paste(unlist(values), collapse = " | "), sep = " — "))
		le8_readable_plot(z)
	})
}

le8_figure_composition <- function(plots, group) {
	n <- length(plots) ; nc <- if (n == 1) 1 else 2 ; nr <- ceiling(n / nc)
	if (group == "diagnosis_timed_profiles" && n == 4L) {
		# Restore the original c1.MOCK.Fig2 geometry: one full-height heatmap
		# alongside three aligned cluster axes. Do not flatten this into a 2x2.
		return(list(
			plot = patchwork::wrap_plots(plots, design = "AB\nAC\nAD", widths = c(1.3, 1), heights = c(1, 1, 1)),
			width = 22, height = 12, policy = "full-height trajectory heatmap + three aligned cluster axes; all numeric rows retained"
		))
	}
	if (group == 'circular') return(list(
		plot = patchwork::wrap_plots(plots, ncol = nc, guides = 'collect') & ggplot2::theme(legend.position = 'bottom'),
		width = if (nc == 1) 12 else 24, height = nr * 12 + .8,
		policy = 'biochemical-category radial associations; shared effect scale; <=6 axes; all significant features labelled'
	))
	if (group == 'effect_concordance' && n == 2L) return(list(
		plot = patchwork::wrap_plots(plots, ncol = 2, widths = c(1.05, 1)), width = 22, height = 10,
		policy = 'paired genetic effects with 95% CIs; fixed-anchor forest; no fitted regression line'
	))
	if (group == 'mr_incident_prevalent' && n == 5L) return(list(
		plot = patchwork::wrap_plots(plots, design = 'AAABBB\nCCDDEE', heights = c(.55, 1.5)), width = 22, height = 14,
		policy = 'two support summaries above three row-aligned MR/incident/prevalent forests; five axes'
	))
	if (group == 'lifestyle_omics_risk' && n == 1L) return(list(
		plot = plots[[1]], width = 20, height = 14,
		policy = 'LE8-selected proxy association ribbons; original nodes and signed links retained; one axis'
	))
	if (group == 'prots.top' && n == 6L) return(list(
		plot = patchwork::wrap_plots(plots, ncol = 3, heights = c(.65, 1.6)), width = 22, height = 17,
		policy = 'restored three overlap panels above three aligned effect forests; six axes'
	))
	if (group == 'directionality_causal' && n == 3L) return(list(
		plot = patchwork::wrap_plots(plots, design = 'AA\nBC', heights = c(1.25, 1)), width = 20, height = 15,
		policy = 'restored full-width temporal/genetic evidence matrix above two diagnostic views'
	))
	if (group %in% c('mrlink2', 'genetic_decomposition') && n == 3L) return(list(
		plot = patchwork::wrap_plots(plots, design = 'AB\nAC', widths = c(1, 1.1)), width = 20, height = 12,
		policy = 'full-height estimate panel + two aligned diagnostic panels; three axes'
	))
	list(
		plot = patchwork::wrap_plots(plots, ncol = nc, widths = rep(1, nc), heights = rep(1, nr)),
		width = if (nc == 1) 10 else 19, height = nr * 6 + .6,
		policy = "equal panel cells; <=6 axes; max 12 labels in source order; all numeric rows retained"
	)
}

le8_flush_figures <- function(rawdir) {
	q <- .le8_figure_queue[[rawdir]] ; if (is.null(q)) return(invisible(NULL))
	nms <- names(q) ; ord <- order(as.integer(sub(".*[.]Fig([0-9]+).*", "\\1", nms)), nms, na.last = TRUE)
	q <- q[ord] ; groups <- list() ; audit <- list() ; notes <- list()
	for (src in names(q)) {
		item <- q[[src]] ; rule <- item$rule
		audit[[src]] <- data.frame(source = src, group = rule$group, status = if (rule$omit) "omitted: redundant, obsolete or misleading view" else "queued")
		if (rule$omit) next
		leaves <- le8_plot_leaves(item$plot)
		# Keep the approved first six C1 figures stable while retaining the old
		# cluster diagnostics on additional pages, rather than discarding them.
		extra <- if (rule$key == 'gradient_cluster' && length(leaves) > 2L) leaves[ - c(1, 2)] else list()
		if (length(extra)) leaves <- leaves[1 : 2]
		missing <- vapply(leaves, function(p) !is.null(attr(p, "le8_unavailable")), logical(1))
		if (any(missing)) notes[[rule$group]] <- c(notes[[rule$group]], vapply(leaves[missing], attr, character(1), "le8_unavailable"))
		pp <- unlist(lapply(leaves[!missing], le8_expand_facets), recursive = FALSE)
		if (!length(pp)) {
			audit[[src]]$status <- "omitted: no usable panel" ; next
		}
		caption <- if (inherits(item$plot, "patchwork")) item$plot$patches$annotation$caption else NULL
		if (is.character(caption)) notes[[rule$group]] <- c(notes[[rule$group]], caption)
		if (rule$group %in% c("temporal_profiles", "gradient_cluster", "diagnosis_window_riskset"))
			notes[[rule$group]] <- c(notes[[rule$group]], "Each participant contributes one baseline omic measurement. Diagnosis-time bins compare different people; these are not within-person longitudinal trajectories.")
		for (p in pp) groups[[rule$group]] <- c(groups[[rule$group]], list(list(plot = p, source = src)))
		if (length(extra)) {
			diagnostics <- unlist(lapply(extra, le8_expand_facets), recursive = FALSE)
			for (p in diagnostics) groups[['gradient_cluster_diagnostics']] <- c(groups[['gradient_cluster_diagnostics']], list(list(plot = p, source = src)))
			notes[['gradient_cluster_diagnostics']] <- 'Restored cluster diagnostics are descriptive heuristics; they do not establish causal direction.'
		}
		audit[[src]]$status <- "rendered in grouped figures"
	}
	staging <- tempfile(".figures-", tmpdir = "/tmp") ; dir.create(staging)
	on.exit(unlink(staging, recursive = TRUE), add = TRUE)
	manifest <- list() ; number <- 0L ; prefix <- q[[1]]$rule$module
	dpi <- as.integer(Sys.getenv("LE8_FIGURE_DPI", unset = "220"))
	for (group in names(groups)[le8_figure_order(prefix, names(groups))]) {
		pp <- groups[[group]]
		for (start in seq(1L, length(pp), by = 6L)) {
			page <- pp[start : min(length(pp), start + 5L)] ; number <- number + 1L
			file <- sprintf("%sFig%d.%s.png", if (prefix == "final") "" else paste0(prefix, "."), number, group)
			n <- length(page) ; composition <- le8_figure_composition(lapply(page, `[[`, "plot"), group)
			caption <- paste(unique(notes[[group]]), collapse = "; ")
			p <- composition$plot +
				patchwork::plot_annotation(
					caption = if (nzchar(caption)) stringr::str_wrap(caption, 150) else NULL,
					tag_levels = "A", theme = ggplot2::theme(plot.caption = ggplot2::element_text(size = 9, hjust = 0), plot.tag = ggplot2::element_text(face = "bold"))
				)
			ggplot2::ggsave(file.path(staging, file), p, width = composition$width, height = composition$height, dpi = dpi, bg = "white", limitsize = FALSE)
			le8_table_plot_exports(lapply(page, `[[`, 'plot'), staging, sub('[.]png$', '', file))
			manifest[[length(manifest) + 1L]] <- data.frame(
				file = file, group = group, panels = n, sources = paste(unique(vapply(page, `[[`, character(1), "source")), collapse = ";"),
				policy = composition$policy
			)
		}
	}
	# Preserve source geometry separately from the aligned/paginated figures.
	originals <- tempfile('le8-source-figures-', tmpdir = '/tmp') ; dir.create(originals)
	on.exit(unlink(originals, recursive = TRUE), add = TRUE)
	for (src in names(q)) {
		item <- q[[src]] ; target <- file.path(originals, src)
		if (is.null(item$plot)) next
		tryCatch(ggplot2::ggsave(target, item$plot,
			width = item$width, height = item$height,
			dpi = dpi, bg = 'white', limitsize = FALSE
		), error = function(e) {
			audit[[src]]$status <- paste(audit[[src]]$status, '; source-layout rendering failed:', conditionMessage(e))
		})
	}
	mf <- if (length(manifest)) do.call(rbind, manifest) else data.frame(file = character(), group = character(), panels = integer(), sources = character(), policy = character())
	data.table::fwrite(mf, file.path(staging, "figure_manifest.csv"))
	data.table::fwrite(do.call(rbind, audit), file.path(staging, "figure_omission_audit.csv"))
	mf <- le8_refresh_figure_files(rawdir, staging)
	.le8_figure_queue[[rawdir]] <- NULL
	invisible(mf)
}
