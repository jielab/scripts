# Whole-genome IBDmix evidence: retained native calls, grouped by the existing
# segment catalog. No GWAS overlap requirement or additional score threshold.
gu_ibdmix_evidence_loci <- function(Q, dataset, build) {
	Q("WITH supported AS (
		SELECT s.segment_code, COUNT(DISTINCT s.sample_id) AS n_carriers,
			COUNT(DISTINCT s.source) AS n_references, GROUP_CONCAT(DISTINCT s.source) AS reference_names,
			MAX(s.score) AS max_lod, COUNT(*) AS n_calls
		FROM segments_evidence s
		WHERE s.dataset_id=? AND s.genome_build=? AND s.method='ibdmix'
			AND s.source_class IN ('Neanderthal','Denisovan')
			AND s.reference_role='scientific_reference'
			AND s.sample_id IS NOT NULL AND s.sample_id!=''
			AND s.chr IN (SELECT chr FROM method_runs WHERE dataset_id=? AND genome_build=?
				AND method='ibdmix' AND status='complete' AND evidence_eligible=1)
		GROUP BY s.segment_code
	)
	SELECT c.segment_code,c.chr,c.start,c.end,c.length_bp,c.source_class,c.genome_build,
		s.n_carriers,s.n_references,s.reference_names,s.max_lod,s.n_calls
	FROM supported s JOIN segment_catalog c ON c.segment_code=s.segment_code
	WHERE c.dataset_id=? AND c.genome_build=?
	ORDER BY s.n_carriers DESC,s.max_lod DESC,c.chr,c.start,c.segment_code",
		list(dataset, build, dataset, build, dataset, build))
}

gu_ibdmix_evidence_cached <- function(Q, dataset, build, db_path,
		cache_dir = file.path("/tmp", "gu-shiny-evidence")) {
	stamp <- function() {
		info <- file.info(db_path)
		list(path = normalizePath(db_path), size = info$size,
			mtime = as.numeric(info$mtime), ctime = as.numeric(info$ctime))
	}
	before <- stamp()
	key <- digest::digest(list(version = 1L, database = before, dataset = dataset, build = build))
	path <- file.path(cache_dir, paste0(key, ".rds"))
	if(file.exists(path)) {
		d <- tryCatch(readRDS(path), error = function(e)NULL)
		if(is.data.frame(d)) return(d)
	}
	d <- gu_ibdmix_evidence_loci(Q, dataset, build)
	# A cache is disposable; never replace the scientific database or reuse a
	# result if final replaced it while the aggregation was running.
	if(identical(before, stamp())) {
		tryCatch({
			dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
			tmp <- tempfile(pattern = key, tmpdir = cache_dir)
			on.exit(unlink(tmp), add = TRUE)
			saveRDS(d, tmp)
			file.rename(tmp, path)
		}, error = function(e)warning("Cannot save IBDmix table cache: ", conditionMessage(e)))
	}
	d
}

gu_ibdmix_evidence_server <- function(input, output, session, Q, dataset, build,
		stamp, db_path, select_record, open_record) {
	all_loci <- reactive({
		stamp()
		withProgress(message = "读取 IBDmix 全基因组证据区间", {
			gu_ibdmix_evidence_cached(Q, dataset(), build(), db_path)
		})
	})
	loci <- reactive({
		d <- all_loci()
		lineage <- input$density_lineage %||% "all"
		chr <- input$ibdmix_evidence_chr %||% "ALL"
		minimum <- suppressWarnings(as.numeric(input$ibdmix_evidence_carriers))
		if(length(minimum) != 1L || !is.finite(minimum)) minimum <- 1
		d[(lineage == "all" | d$source_class == lineage) & (chr == "ALL" | d$chr == chr) &
			d$n_carriers >= max(1, minimum), , drop = FALSE]
	})
	output$ibdmix_evidence_note <- renderText({
		paste0("显示 ", format(nrow(loci()), big.mark = ","), " / ",
			format(nrow(all_loci()), big.mark = ","), " 条全基因组谱系区间（", dataset(), " · ", build(),
			"）；沿用原生片段过滤，仅纳入已完成且可用于证据评估的科学参考结果，不限 GWAS loci。")
	})
	output$ibdmix_evidence_loci <- renderDT({
		d <- loci()
		show <- data.frame(
			ID = d$segment_code, Chr = d$chr,
			region = paste0(format(d$start + 1, scientific = FALSE, trim = TRUE), "–",
				format(d$end, scientific = FALSE, trim = TRUE)),
			lineage = d$source_class, kb = round(d$length_bp / 1000, 2),
			carriers = d$n_carriers, references = d$n_references, names = d$reference_names,
			LOD = round(d$max_lod, 3), calls = d$n_calls, check.names = FALSE)
		names(show) <- c("ID", "Chr", paste0("代表区间 (", build(), ")"), "谱系", "长度 (kb)",
			"携带者人数", "参考数", "参考来源", "最高 LOD", "片段记录数")
		DT::datatable(show, rownames = FALSE, class = "compact stripe hover nowrap",
			selection = list(mode = "single", selected = NULL),
			options = list(pageLength = 10, lengthMenu = c(10, 25, 50), scrollX = TRUE,
				order = list(list(5, "desc"), list(8, "desc")),
				columnDefs = list(list(targets = 0, visible = FALSE)),
				language = list(emptyTable = "当前筛选条件下无符合证据资格的 IBDmix 区间")),
			# Stable IDs are required: server-side pages/sorting change row indices.
			callback = DT::JS("table.on('click','tbody tr',function(){var r=table.row(this).data();if(r)Shiny.setInputValue('ibdmix_evidence_select',r[0],{priority:'event'});});table.on('dblclick','tbody tr',function(){var r=table.row(this).data();if(r)Shiny.setInputValue('ibdmix_evidence_open',r[0],{priority:'event'});});"))
	}, server = TRUE)
	lookup <- function(code) {
		d <- loci(); i <- match(code, d$segment_code)
		req(length(i) == 1L, !is.na(i))
		r <- d[i, , drop = FALSE]
		r$super_population <- "ALL"
		r
	}
	observeEvent(input$ibdmix_evidence_select, select_record(lookup(input$ibdmix_evidence_select), "ibdmix"))
	observeEvent(input$ibdmix_evidence_open, open_record(lookup(input$ibdmix_evidence_open), "ibdmix"))
	invisible(loci)
}
