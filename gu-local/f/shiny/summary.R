# Current reference summaries use exact interval unions and certified sample rosters.
gu_reference_order <- function(references) {
	current <- c("Altai", "Chagyr", "Denisova", "Denisova25", "Vindija")
	actual <- unique(as.character(references))
	c(intersect(current, actual), sort(setdiff(actual, current)))
}

gu_ibdmix_summary <- function(samples, burden, reference) {
	d <- merge(samples, burden[burden$reference == reference, , drop = FALSE], by = "sample_id", all.x = TRUE, sort = FALSE)
	good <- is.finite(d$coordinate_coverage_pct) & d$tested_bp > 0 & d$n_chromosomes > 0
	measured <- d[good, , drop = FALSE]
	pops <- lapply(split(d, d$population), function(p) {
		x <- p[is.finite(p$coordinate_coverage_pct) & p$tested_bp > 0, , drop = FALSE]
		data.frame(population = p$population[[1]], super_population = p$super_population[[1]],
			n = nrow(x), n_total = nrow(p), pct = if(nrow(x))mean(x$coordinate_coverage_pct) else NA_real_,
			mb = if(nrow(x))mean(x$archaic_bp) / 1e6 else NA_real_)
	})
	list(reference = reference, samples = d, measured = measured,
		populations = if(length(pops))do.call(rbind, pops) else data.frame(),
		complete = nrow(d) > 0 && all(good) && all(d$n_chromosomes == 22),
		comparable_scope = nrow(measured) > 0 && length(unique(measured$chromosomes)) == 1L,
		overall_pct = if(nrow(measured))mean(measured$coordinate_coverage_pct) else NA_real_,
		overall_diploid_pct = if(nrow(measured))mean(measured$coverage_pct) else NA_real_,
		overall_mb = if(nrow(measured))mean(measured$archaic_bp) / 1e6 else NA_real_)
}

gu_ibdmix_map <- function(summary, land, locations) {
	old <- par(mar = c(1, 0, 2, 0), family = "sans", bg = "white");on.exit(par(old))
	plot.new();plot.window(xlim = c(-168, 181), ylim = c(-63, 93), asp = 1, xaxs = "i", yaxs = "i")
	for(d in split(land, land$polygon)) polygon(d$longitude, d$latitude, col = "#f4f5f5", border = "#bbc3c8", lwd = .65)
	title(main = paste(summary$reference, "coordinate coverage by population"), adj = 0, cex.main = 1.05, col.main = "#253746", line = .3)
	p <- merge(locations, summary$populations, by = c("population", "super_population"), sort = FALSE)
	if(!nrow(p)) {text(0, 0, "No mapped 1KG populations in this dataset", col = "#607080");return(invisible(NULL))}
	segments(p$longitude, p$latitude, p$label_x, p$label_y, col = "#aab3ba", lwd = .7)
	points(p$longitude, p$latitude, pch = 20, col = "#7f8a92", cex = .45)
	for(i in seq_len(nrow(p))) {
		r <- p[i, ];radius <- 4.5;theta <- seq(0, 2 * pi, length.out = 100)
		polygon(r$label_x + radius * cos(theta), r$label_y + radius * sin(theta), col = if(is.finite(r$pct))"#e0e3e5" else "white", border = "#89959e", lwd = .8)
		if(is.finite(r$pct) && r$pct > 0) {
			theta <- seq(pi / 2, pi / 2 - 2 * pi * min(r$pct, 100) / 100, length.out = 40)
			polygon(c(r$label_x, r$label_x + radius * cos(theta)), c(r$label_y, r$label_y + radius * sin(theta)), col = "#b3222b", border = NA)
		}
		label <- paste0(r$population, "\n", if(is.finite(r$pct))sprintf("%.2f%%", r$pct) else "N/A")
		text(r$label_x, r$label_y - 6, label, adj = c(.5, 1), cex = .78, col = "#253746", font = 2)
	}
	legend("bottomleft", inset = c(.015, 0), legend = c("Interval union / physical span", "Remaining span", "Not evaluated"),
		fill = c("#b3222b", "#e0e3e5", "white"), border = "#89959e", bty = "n", horiz = TRUE, cex = .75, text.col = "#52616d")
}

gu_ibdmix_summary_text <- function(s, manifest) {
	pct <- function(x)if(is.finite(x))sprintf("%.2f%%", x) else "N/A"
	mb <- function(x)if(is.finite(x))sprintf("%.1f", x) else "—"
	d <- s$measured
	scope <- if(s$complete)"22 / 22 autosomes · X excluded" else if(nrow(d))paste0("Partial coverage · ", paste(sort(unique(d$n_chromosomes)), collapse = ", "), " / 22 autosomes per person") else "No certified whole-autosome results"
	shiny::tagList(
		shiny::tags$div(class = "gu-summary-kicker", paste(s$reference, "· AUTOSOMAL COORDINATE COVERAGE U / L")),
		shiny::tags$div(class = "gu-summary-value", pct(s$overall_pct)),
		shiny::tags$p(class = "gu-summary-subvalue", paste0(mb(s$overall_mb), " Mb / person")),
		shiny::tags$p(paste0(format(nrow(d), big.mark = ","), " / ", format(nrow(s$samples), big.mark = ","), " individuals · ", sum(s$populations$n > 0), " populations")),
		shiny::tags$p(class = "gu-summary-scope", scope),
		shiny::tags$p(class = "text-muted", paste0("对应二倍体物理分母 U/(2L)：", pct(s$overall_diploid_pct), "。这两种统计都不等于已校准的祖源剂量比例或 callable-base 比例。")),
		shiny::tags$p(class = "text-muted", "只统计当前所选参考，先按个体合并重叠片段再等权平均。参考匹配不单独证明渗入来源，实际过滤和未校准项见运行记录。"),
		if(length(manifest$diploid_autosome_bp))shiny::tags$p(class = "text-muted", sprintf("完整常染色体二倍体物理长度：%.3f Gb（2 × 参考序列长度；未定相片段并集不能确定纯合双份剂量）。", manifest$diploid_autosome_bp / 1e9)),
		if(!s$complete)shiny::tags$p(class = "text-muted", "仅统计已确认完成的整条常染色体；未分析的染色体、人群或缺少目标样本记录的运行记为 N/A，不外推全基因组比例。"),
		if(!s$comparable_scope && nrow(d))shiny::tags$p(class = "text-muted", "个体的已分析染色体范围不同，人群比例不能直接比较。")
	)
}

gu_ibdmix_summary_table <- function(s, burden, x_male = FALSE) {
	expected_chr <- if(x_male)1L else 22L
	refs <- gu_reference_order(burden$reference[is.finite(burden$n_chromosomes) & burden$n_chromosomes > 0])
	if(!length(refs))return(shiny::tags$p(class = "gu-summary-caption", if(x_male)"暂无已完成且有目标样本记录的男性非 PAR X 参考结果。" else "暂无已完成且有目标样本记录的整条常染色体参考结果。"))
	afr_pops <- c("ESN", "GWD", "LWK", "MSL", "YRI")
	d <- merge(s$samples[, c("sample_id", "population", "super_population")], burden, by = "sample_id", sort = FALSE)
	d <- d[is.finite(d$archaic_bp) & d$tested_bp > 0 & d$n_chromosomes > 0, , drop = FALSE]
	present_super <- unique(d$super_population)
	present_pop <- intersect(afr_pops, unique(d$population))
	groups <- c("ALL", if(length(present_pop))c("AFR (5)", present_pop), intersect(c("AFR", "EUR", "EAS", "SAS", "AMR"), present_super), sort(setdiff(present_super, c("AFR", "EUR", "EAS", "SAS", "AMR"))))
	group_rows <- function(d, group) {
		if(group == "ALL")d else if(group == "AFR (5)")d[d$population %in% afr_pops, , drop = FALSE] else if(group %in% afr_pops)d[d$population %in% group, , drop = FALSE] else d[d$super_population %in% group, , drop = FALSE]
	}
	rows <- lapply(groups, function(group) {
		p <- group_rows(d, group)
		cells <- lapply(refs, function(ref) {
			x <- p[p$reference == ref, , drop = FALSE]
			if(!nrow(x))return(shiny::tags$td("N/A"))
			scope <- paste(sort(unique(x$n_chromosomes)), collapse = ",")
			partial <- any(x$n_chromosomes != expected_chr)
			detail <- paste0("n = ", nrow(x), "; chromosomes: ", paste(unique(x$chromosomes), collapse = "; "),
				"; ", if(x_male)"U/L_X = " else "U/(2L_A) = ", sprintf("%.3f%%", mean(x$coverage_pct)))
			shiny::tags$td(title = detail, sprintf(if(x_male)"%.3f Mb" else "%.1f Mb", mean(x$archaic_bp) / 1e6),
				shiny::tags$small(class = "gu-summary-cell-scope", paste0(sprintf("U/L %.3f%%", mean(x$coordinate_coverage_pct)), " · n = ", nrow(x), if(partial)paste0(" · partial ", scope, "/", expected_chr, " chr"))))
		})
		shiny::tags$tr(class = if(group %in% afr_pops)"gu-summary-subpopulation" else "gu-summary-group", shiny::tags$td(group), cells)
	})
	headers <- c(list(shiny::tags$th("Group")), lapply(refs, function(ref)shiny::tags$th(ref, shiny::tags$br(), "Mb / person · U/L (%)")))
	shiny::tagList(
		shiny::tags$h6(if(x_male)"chrX · male non-PAR · tested references" else "Population summary · autosomes · tested references"),
		shiny::tags$div(class = "gu-summary-table-scroll", shiny::tags$table(class = "table table-sm gu-summary-table", shiny::tags$thead(shiny::tags$tr(headers)), shiny::tags$tbody(rows))),
		shiny::tags$p(class = "gu-summary-caption", "只显示实际检测的参考。U 为个体区间并集长度，L 为已检测染色体的单倍体物理长度。各参考按个体合并重叠后再等权平均；已检测但零片段者计为 0，未检测者记 N/A。鼠标悬停可查看染色体范围和物理分母统计。"),
		shiny::tags$p(class = "gu-summary-caption", if(x_male)"男性 X 仅使用非 PAR 区域，女性不进入分母。伪二倍体 IBDmix 的 X LOD 尚未校准；与常染色体的过滤也可能不同，因此 X/常染色体比例仅作描述。" else "常染色体统计排除 X。U/L 和 U/(2L) 是不同物理分母下的覆盖统计；未分相片段不能确定真实拷贝剂量。不同参考或人群的检测范围不一致时，Mb/人不能直接比较。"),
		shiny::tags$p(class = "gu-summary-caption", "Altai、Chagyr、Vindija 属于 Neanderthal 参考；Denisova、Denisova25 属于 Denisovan 参考。实际参考版本和过滤记录保存在原生运行元数据中。AFR (5) 为已检测的 ESN、GWD、LWK、MSL、YRI 合计。")
	)
}


# Saved aggregate viewing requires neither an individual database nor raw sequences.
gu_summary_read_workbook <- function(path, table, optional = FALSE) {
	if (!file.exists(path)) return(data.frame())
	directory <- tempfile("gu-summary-", tmpdir = "/tmp")
	dir.create(directory)
	on.exit(unlink(directory, recursive = TRUE), add = TRUE)
	utils::unzip(path, files = "results/manifest.json", exdir = directory)
	manifest <- jsonlite::read_json(file.path(directory, "results/manifest.json"))
	entries <- Filter(function(x) identical(x$name, table), manifest$tables)
	if (length(entries) == 0L && optional) return(data.frame())
	if (length(entries) != 1L) stop("Missing or ambiguous saved summary table: ", table)
	entry <- manifest$files[[entries[[1L]]$source]]
	if (!grepl("^[a-f0-9]{64}$", entry$sha256) || !identical(entry$part, paste0("results/exports/", entry$sha256, ".bin"))) stop("Invalid saved summary source")
	utils::unzip(path, files = entry$part, exdir = directory)
	source <- file.path(directory, entry$part)
	if (!identical(digest::digest(file = source, algo = "sha256"), entry$sha256)) stop("Saved summary checksum mismatch")
	as.data.frame(data.table::fread(source, showProgress = FALSE, na.strings = c("", "NA")))
}

gu_summary_only_app <- function(path, phyml_path) {
	burden <- gu_summary_read_workbook(path, "population_burden")
	paired <- gu_summary_read_workbook(path, "same_males_X_autosomes", optional = TRUE)
	phyml <- gu_summary_read_workbook(phyml_path, "locus")
	ui <- bslib::page_fluid(
		shiny::tags$h2("GU · 保存的科学汇总"),
		shiny::tags$p("汇总查看模式：读取正式分析保存的 XLSX，不重新计算、不补零或外推未检测染色体。"),
		shiny::tags$p("个体查询、片段浏览、序列比对和交互树图不可用：需要完整数据库及原生结果。"),
		shiny::selectInput("scope", "范围", choices = unique(burden$scope)),
		shiny::selectInput("reference", "参考", choices = unique(burden$reference)),
		bslib::navset_card_tab(
			bslib::nav_panel("IBDmix 人群负荷", DT::DTOutput("burden")),
			bslib::nav_panel("同一男性 X / 常染色体", DT::DTOutput("paired")),
			bslib::nav_panel("全部 PhyML 位点状态", DT::DTOutput("phyml"))))
	server <- function(input, output, session) {
		output$burden <- DT::renderDT(DT::datatable(burden[burden$scope == input$scope & burden$reference == input$reference, , drop = FALSE], rownames = FALSE, options = list(scrollX = TRUE)))
		output$paired <- DT::renderDT(DT::datatable(if(nrow(paired))paired else data.frame(状态 = "没有保存可比的完整同一男性队列"), rownames = FALSE, options = list(scrollX = TRUE)))
		output$phyml <- DT::renderDT(DT::datatable(if(nrow(phyml))phyml else data.frame(状态 = "未提供 PhyML 汇总"), rownames = FALSE, options = list(scrollX = TRUE)))
	}
	shiny::shinyApp(ui, server)
}
