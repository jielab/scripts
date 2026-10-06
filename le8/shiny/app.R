#!/usr/bin/env Rscript


# 🚩 app
# Standalone Shiny explorer: data access, question views, UI and server.
arg <- grep("^--file=", commandArgs(), value = TRUE)
app_dir <- if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)) else normalizePath(".")
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Aggregate data access
# Read-only aggregate access through result workbooks; no individual/model loading or arbitrary paths.
suppressPackageStartupMessages({
	library(shiny)
	library(DT)
	library(data.table)
	library(ggplot2)
	library(digest)
})
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
	i <- match(name, args)
	if (!is.na(i)) {
		if (i == length(args)) stop(name, ' needs a value')
		return(args[i + 1L])
	}
	default
}
bundled <- file.path(app_dir, '../results')
default_root <- if (dir.exists(bundled)) bundled else '/mnt/d/analysis/le8'
LE8_ROOT <- normalizePath(option('--analysis-root', Sys.getenv('LE8_ANALYSIS_ROOT', default_root)), winslash = '/', mustWork = TRUE)
LE8_INDEX <- Sys.getenv("LE8_BROWSER_INDEX", file.path(LE8_ROOT, "shiny"))
LE8_MAX_BYTES <- as.numeric(Sys.getenv("LE8_SHINY_MAX_TABLE_MB", "128")) * 1024 ^ 2
if (!is.finite(LE8_MAX_BYTES) || LE8_MAX_BYTES <= 0) stop("Invalid LE8_SHINY_MAX_TABLE_MB")
private_columns <- c(
	"eid", "iid", "fid", "participant_id", "individual_id", "subject_id", "patient_id", "query_id",
	"target_id", "reference_id", "donor_id", "date_birth", "birth_date",
	"id", "id1", "id2", "sampleid", "ukbid", "ukbeid"
)

has_private_columns <- function(columns) {
	names <- tolower(columns)
	# Match the indexer's S7 aggregate schema: reference_id names a model.
	if (all(c('run_id', 'model_id', 'primary_id', 'reference_id', 'n', 'uno_c_horizon') %in% names))
		names <- setdiff(names, 'reference_id')
	any(gsub('[^a-z0-9]', '', names) %in% gsub('[^a-z0-9]', '', private_columns) |
		grepl('(^|[_. #])(eid|iid|fid)($|[_. ])', names))
}


# 🚩 Aggregate workbook access
# Read only the recorded aggregate exports; participant/model RDS files are never loaded.
for (expr in parse(file.path(app_dir, '../f/0.common.R'))) {
	if (is.call(expr) && identical(expr[[1]], as.name('<-')) && is.symbol(expr[[2]]) &&
		startsWith(as.character(expr[[2]]), 'le8_table')) eval(expr)
}
.table_stores <- new.env(parent = emptyenv())
read_stored_table <- function(path) {
	directory <- dirname(path)
	workbooks <- le8_table_archive_files(directory)
	if (!length(workbooks)) return(NULL)
	key <- normalizePath(directory, winslash = '/', mustWork = TRUE)
	stamp <- file.info(workbooks)$mtime
	cached <- .table_stores[[key]]
	if (is.null(cached) || !identical(cached$stamp, stamp)) {
		cached <- list(stamp = stamp, store = le8_table_load(directory))
		.table_stores[[key]] <- cached
	}
	entry <- cached$store$files[[basename(path)]]
	if (is.null(entry)) return(NULL)
	if (!identical(entry$kind, 'table')) stop('Requested source is not an aggregate table')
	if (entry$size > LE8_MAX_BYTES) stop('Source exceeds the configured table size limit')
	entry <- le8_table_materialize(entry, basename(path))
	if (has_private_columns(names(entry$data))) stop('Participant identifiers found')
	entry
}

read_index <- function(name) {
	path <- file.path(LE8_INDEX, paste0(name, ".csv"))
	if (file.exists(path)) {
		d <- as.data.frame(data.table::fread(path, na.strings = c("", "NA", "NaN"), showProgress = FALSE))
	} else {
		entry <- read_stored_table(path)
		if (is.null(entry)) return(data.frame())
		d <- entry$data
	}
	if (has_private_columns(names(d))) stop('Participant identifiers found in prepared index')
	bool <- grep("^measured_.*supported$|^measured_.*lt_005$|^same_people$|^same_N$|_is_MR$", names(d), value = TRUE)
	for (k in bool) if (is.character(d[[k]]))
		d[[k]] <- tolower(d[[k]]) %in% c("true", "t", "1")
	d
}
load_catalogue <- function() {
	if (!length(le8_table_archive_files(LE8_INDEX)) && !file.exists(file.path(LE8_INDEX, "tables.csv")))
		stop("No prepared evidence index. Run: ./le8.sh final --Y cvd_cad,ra --biom prot,met --index-only")
	ans <- setNames(lapply(c(
		"candidates", "effects", "loci", "proxy_comparisons", "discovery_counts", "tables",
		"figures", "status"
	), read_index), c(
		"candidates", "effects", "loci", "proxy_comparisons", "discovery_counts",
		"tables", "figures", "status"
	))
	for (name in question_names) ans[[paste0("question_", name)]] <- read_index(paste0("question_", name))
	ans
}
resolve_source <- function(relative, registry) {
	if (length(relative) != 1 || is.na(relative) || !relative %in% registry$path)
		stop("Source is not in the prepared allowlist")
	if (grepl("(^|/)([_.]|private|cache|input|neural|checkpoints)", relative, ignore.case = TRUE))
		stop("Private source path")
	path <- file.path(normalizePath(dirname(file.path(LE8_ROOT, relative)), winslash = "/", mustWork = TRUE), basename(relative))
	root <- paste0(LE8_ROOT, "/")
	if (!startsWith(path, root))
		stop("Source escapes the analysis root")
	entry <- if (!file.exists(path)) read_stored_table(path) else NULL
	if (!file.exists(path) && is.null(entry)) stop("Missing aggregate source")
	if (file.exists(path) && file.info(path)$size > LE8_MAX_BYTES)
		stop("Source exceeds the configured table/image size limit")
	expected <- registry$sha256[match(relative, registry$path)]
	actual <- if (is.null(entry)) digest::digest(file = path, algo = "sha256") else entry$sha256
	if (length(expected) && !is.na(expected) && nzchar(expected) && actual != expected)
		stop("Source changed after indexing. Run final --index-only, then Reload prepared index.")
	path
}
read_aggregate <- function(relative, registry) {
	path <- resolve_source(relative, registry)
	if (!grepl("[.](csv|tsv)([.]gz)?$", path, ignore.case = TRUE))
		stop("Only aggregate CSV/TSV tables are allowed")
	if (!file.exists(path)) return(read_stored_table(path)$data)
	header <- data.table::fread(path, nrows = 0, showProgress = FALSE)
	if (has_private_columns(names(header)))
		stop("Participant identifiers found; this table is not exposed")
	as.data.frame(data.table::fread(path, showProgress = FALSE, na.strings = c("", "NA", "NaN")))
}
feature_name <- function(value) {
	v <- gsub("^(prot|met)__", "", as.character(value))
	v <- sub("[.]pgs$", "", v, ignore.case = TRUE)
	v[tolower(v) %in% c("lactat", "lactate")] <- "Lactate"
	v
}
feature_column <- function(d) {
	key <- intersect(c(
		"feature", "assay", "exposure", "protein", "Protein", "metabolite", "Metabolite", "gene",
		"term"
	), names(d))
	if (length(key))
		key[1] else NULL
}
subset_context <- function(d, Y, layer) {
	if (!nrow(d))
		return(d)
	if ("Y" %in% names(d))
		d <- d[!is.na(d$Y) & d$Y %in% c(Y, "all"), , drop = FALSE]
	if ("layer" %in% names(d) && length(layer) == 1 && !is.na(layer) && layer != "both")
		d <- d[!is.na(d$layer) & d$layer %in% c(layer, "joint", "all"), , drop = FALSE]
	d
}
make_dt <- function(d, select = "single", page = 12) {
	if (!ncol(d))
		d <- data.frame(status = "No available aggregate rows for this view")
	tab <- DT::datatable(d,
		rownames = FALSE, filter = "top", selection = select, escape = TRUE, extensions = "Buttons",
		options = list(pageLength = page, scrollX = TRUE, autoWidth = FALSE, dom = "Bfrtip", buttons = c(
			"copy",
			"csv"
		), lengthMenu = c(10, 25, 50, 100), deferRender = TRUE)
	)
	num <- names(d)[vapply(d, is.numeric, logical(1))]
	# DT numeric sort uses the underlying full-precision value, not formatted text.
	if (length(num))
		tab <- DT::formatSignif(tab, columns = num, digits = 4)
	tab
}
empty_plot <- function(message) {
	ggplot() +
		annotate("text", x = 0, y = 0, label = message, size = 4) +
		theme_void()
}
plot_effects <- function(d, title) {
	if (!nrow(d) || !all(c("beta", "series") %in% names(d)))
		return(empty_plot("No estimate in this scope; missing evidence is not zero."))
	d <- d[is.finite(d$beta), , drop = FALSE]
	if (!nrow(d))
		return(empty_plot("No finite estimate."))
	if ("model" %in% names(d))
		d$label <- paste(d$series, d$model, sep = " | ") else d$label <- d$series
	d$label <- factor(d$label, levels = rev(unique(d$label)))
	ggplot(d, aes(beta, label)) +
		geom_vline(xintercept = 0, linetype = 2) +
		geom_errorbar(aes(xmin = lo, xmax = hi),
			orientation = "y", width = 0.16, na.rm = TRUE
		) +
		geom_point(size = 2.5) +
		labs(title = title, x = paste(unique(d$units),
			collapse = "; "
		), y = NULL, caption = "95% intervals are as reported/derived from the reported SE. G/R intervals may condition on fitted calibration.") +
		theme_bw(base_size = 12) +
		theme(plot.caption = element_text(hjust = 0, size = 9), axis.title.x = element_text(size = 10))
}


# Research questions
# Question-led Final dashboard. All controls filter frozen aggregate estimates.
question_names <- c("overview", "prediction", "contrasts", "proxy", "members", "pillars",
	"heterogeneity", "inflammation_definition", "design", "fit", "abm_metrics", "abm_coverage",
	"abm_paired", "abm_support", "abm_training", "abm_gate", "abm_risk_gain", "abm_audit", "abm_decision", "genetic", "temporal", "cohort", "mediation", "modules",
	"nonlinear", "nonlinear_curves", "same_locus", "susie_pairs", "susie_diagnostics", "mr_signal_evidence", "cell", "cell_status", "cell_contrasts", "cell_evidence", "concept_coefficients", "concept_status", "concept_fold_panels", "mr_scope", "dandelion", "dandelion_lolo", "dandelion_native", "state_projection", "decomposition", "age_models", "sources")

question_ui <- function(id) {
	ns <- NS(id)
	tabPanel("研究问题 · Final", value = "questions",
		div(class = "q-intro",
			h2("从研究假设到可核查的结果"),
			div(class = "q-path", "LE8 → omics → disease ← omics ← genetics"),
			p("箭头表示待检验的机制假设。下方分别检验预测、可解释性、个体分层与遗传证据。")),
		tabsetPanel(id = ns("section"),
			tabPanel("结果总览", value = "overview",
				uiOutput(ns("cards")), DTOutput(ns("overview")),
				h4("研究人群"), DTOutput(ns("cohort")),
				p("prot 与 met 的样本量可能不同；以下比较只在各自原生验证方案内进行。顶部 scope / adjustment 只过滤 biomarker 证据，不改变冻结的预测模型。")),
			tabPanel("1 · LE8 supervision", value = "supervision",
				p("NS：训练集多变量疾病模型选择；YS：只按重复验证的 LE8 代理强度选择；YSplus：总预算内按 80% YS 加 20% 疾病候选；YSbalanced：按代理强度均衡可用 pillars。YSconcept 用同一最终 assay panel 预测 LE8 概念，再进入可加的风险模型。实际组成见 biom 清单，不强制填满缺乏支持的pillar。"),
				fluidRow(column(2, selectInput(ns("budget"), "相同 biom 数量", choices = c(5, 10, 50), selected = 10)),
					column(3, selectInput(ns("family"), "Panel 方法", c("YS" = "YS", "YS + plus（YSP 扩展）" = "YSplus", "LE8 pillar 均衡" = "YSbalanced", "LE8 概念模型" = "YSconcept", "LE8 概念替代" = "YSconceptReplacement", "LE8 概念＋补充指标" = "YSconceptPlus"))),
					column(2, selectInput(ns("donors"), "Proxy discovery 人群", c("Yin + Yang" = "YinYang", "Yin" = "Yin"))),
					column(2, selectInput(ns("landmark"), "Landmark（年）", c(0, 2, 5))),
					column(3, selectInput(ns("stratum"), "验证人群", c("All", "Low baseline inflammation", "High baseline inflammation")))),
				uiOutput(ns("panel_note")),
				fluidRow(column(6, plotOutput(ns("prediction_plot"), height = 380)), column(6, plotOutput(ns("delta_plot"), height = 380))),
				DTOutput(ns("contrasts")),
				h4("解释性：相同预算能否更好地重建 LE8？"),
				plotOutput(ns("proxy_plot"), height = 400), DTOutput(ns("proxy")),
				p("ΔR² = 当前 panel − 同预算 NS。配对 bootstrap 条件于冻结 panel 和训练拟合；FDR 在每个疾病×omics内覆盖全部 panel×LE8 比较。提高重建保真度不等于证明干预响应。"),
				tabsetPanel(tabPanel("8 pillars 与 biom 清单", DTOutput(ns("pillars")), DTOutput(ns("members"))),
					tabPanel("YinYang − Yin", plotOutput(ns("yin_plot"), height = 360), DTOutput(ns("yin"))),
					tabPanel("炎症分层与异质性", DTOutput(ns("inflammation")), DTOutput(ns("heterogeneity")),
						p("低/高基线炎症由训练数据冻结阈值定义。低炎症组不等于已确立的 non-inflammatory CAD；看组间 ΔAUC 的直接差异及其 FDR。")),
					tabPanel("模型与研究设计", DTOutput(ns("design")), DTOutput(ns("fit"))),
					tabPanel("LE8 概念风险模型", DTOutput(ns("concept_status")), DTOutput(ns("concept_coefficients")), DTOutput(ns("concept_fold_panels"))))),
			tabPanel("2 · ABM 分层", value = "abm",
				fluidRow(column(3, selectInput(ns("backend"), "独立运行方案", c("Reference / Selective" = "reference", "TabICLv2" = "tabicl"))),
					column(4, selectInput(ns("abm_model"), "主模型", "elasticnet_weighted")),
					column(3, selectInput(ns("abm_metric"), "覆盖率曲线指标", c("IPCW AUC（越高越好）" = "AUC_IPCW", "IPCW Brier（越低越好）" = "Brier_IPCW", "相对 elastic net 的 AUC 差" = "delta_AUC_vs_elasticnet", "相对 elastic net 的 Brier 改善" = "Brier_gain_vs_elasticnet")))),
				uiOutput(ns("abm_note")),
				fluidRow(column(6, plotOutput(ns("abm_plot"), height = 380)), column(6, plotOutput(ns("coverage_plot"), height = 380))),
				DTOutput(ns("abm_metrics")),
				p("supported / rejected 由开发集冻结的规则分层，再对测试人群评价；它与外层验证互补。图中是固定时间窗的 AUC，不是 survival C-index。覆盖率曲线展示所有阈值，不自动选测试集上的最佳点。"),
				h4("同一人群上的配对比较"), DTOutput(ns("abm_paired")),
				h4("是否只是筛出了低风险人群？"), DTOutput(ns("abm_support")), DTOutput(ns("abm_risk_gain")),
				selectInput(ns("abm_selector"), "覆盖率曲线的冻结规则", c("预测增益"="gain", "仅支持度"="support_only", "绝对误差"="absolute_error", "临床低风险"="clinical_lowrisk", "随机对照"="random")),
				tabsetPanel(tabPanel("筛选训练与匹配对照", DTOutput(ns("abm_training"))),
					tabPanel("Gate 与独立审计", DTOutput(ns("abm_gate")), DTOutput(ns("abm_audit"))),
					tabPanel("决策曲线数据", DTOutput(ns("abm_decision")))) ,
				p("同时查看事件率、绝对误差及相对于同组 elastic net 的 Brier 改善。基线风险降低本身会降低 Brier。筛选后训练应与 target_tuned 和随机筛选对照比较，不能只比较入选者与全人群。")),
			tabPanel("3 · 遗传 / 实测 / 时间", value = "genetic",
				fluidRow(column(4, selectInput(ns("anchor"), "预设案例（点击同步 biomarker 明细）", choices = character())),
					column(5, selectInput(ns("genetic_screen"), "证据筛选", c("全部同人群比较" = "all", "实测 P≥0.05，PGS FDR<0.05" = "discordant", "G−R 差异 FDR<0.05" = "difference")))),
				p("PGS 来源：COJO 条件效应 bJ × 基因型 dosage 的加权和（gwas_format.sh pgs）。它捕获遗传倾向；校准 G 与剩余 R 的差异不等同于出生前/出生后或因果/非因果的纯分解。"),
				plotOutput(ns("genetic_plot"), height = 420), DTOutput(ns("genetic")),
				h4("Yin / Yang 与远期关联"), DTOutput(ns("temporal")),
				p("Incident、prevalent、病程及 landmark 是互补证据。原表的 birthline 是成年实测 biom 的 attained-age Cox 分析，并非出生时浓度；reverse 是探索性反向时间分析，并非下方的 reverse MR。患病后治疗、存活选择和亚临床病变均会影响解释；不能把 log OR 与 log HR 直接相减证明反向因果。"),
				detail_ui("questions_detail")),
			tabPanel("5C 机制与扩展", value = "mechanisms",
				p("选择顶部 biomarker 查看其 LE8 pillar、描述性中介与同位点遗传证据。Cell 注释、coloc 与 MR 保持不同证据层级。"),
				tabsetPanel(
					tabPanel("LE8 → omics → Y", DTOutput(ns("modules")), DTOutput(ns("mediation")),
						p("当前 mediation 是基线线性×Cox系数乘积的关联分解；尚未识别因果中介比例。")),
					tabPanel("MR + 同位点 coloc", DTOutput(ns("same_locus")), DTOutput(ns("mr_scope")), DTOutput(ns("dandelion"))),
					tabPanel("SuSiE 多信号", DTOutput(ns("susie_pairs")), DTOutput(ns("susie_diagnostics")), DTOutput(ns("mr_signal_evidence"))),
					tabPanel("DANDELION 诊断", DTOutput(ns("dandelion_native")), DTOutput(ns("dandelion_lolo"))),
					tabPanel("状态与年龄", DTOutput(ns("state_projection")), DTOutput(ns("age_models"))),
					tabPanel("C2 G/R 分解", DTOutput(ns("decomposition"))),
					tabPanel("非线性 / breakpoints", DTOutput(ns("nonlinear")),
						p("当前表检验 spline 非线性。nadir 是拟合曲线的最低点，不是断点；年龄阈值与出生队列变化须另行拟合、给出阈值区间并处理年龄–时期–队列不可识别性。")),
					tabPanel("Cellulation", DTOutput(ns("cell_status")), DTOutput(ns("cell")), DTOutput(ns("cell_contrasts")), DTOutput(ns("cell_evidence")),
						p("Cell-expression 富集使用完整 assay 背景；不等于蛋白释放来源或细胞衰老时钟。CIGMA 状态以完成的结果为准。")))),
			tabPanel("下载 / 来源", value = "sources",
				selectInput(ns("download_kind"), "当前疾病 / omics 的完整证据表", choices = question_names),
				downloadButton(ns("download"), "下载 CSV"), DTOutput(ns("sources")),
				p("下载保留全部比较及原始来源，不随图中预算或方法过滤。source_row 为原CSV数据行（不含表头）；source_sha256 用于核查来源版本。"),
				h4("方法与文献对应"),
				tags$ul(
					tags$li(tags$a(href = "https://pmc.ncbi.nlm.nih.gov/articles/11634769/", target = "_blank", rel = "noopener", "Schuermans / Natarajan：cardiac proteomics"), "：作为关联与预测参照，外部论文的AUC不能当作本队列的配对基线。"),
					tags$li(tags$a(href = "https://academic.oup.com/proteincell/article/17/3/231/8250438", target = "_blank", rel = "noopener", "Systematic CVD proteomics"), "：蛋白与临床信息的互补性需要同测试人群比较。"),
					tags$li(tags$a(href = "https://www.nature.com/articles/s41591-025-04105-8", target = "_blank", rel = "noopener", "Metabolites、genetics、lifestyle 与 T2D"), "：本页将分段证据与冻结panel验证同时展示。"),
					tags$li(tags$a(href = "https://pmc.ncbi.nlm.nih.gov/articles/PMC13496006/", target = "_blank", rel = "noopener", "DANDELION：trans-regulatory gene mapping"), "：本项目适配是否符合原方法输入要求，见C2输入审计。"),
					tags$li(tags$a(href = "https://www.nature.com/articles/s41586-026-10577-6", target = "_blank", rel = "noopener", "CIGMA"), "：需配对单细胞表达与基因型；现有cell-expression注释不替代CIGMA拟合。"),
					tags$li(tags$a(href = "https://pmc.ncbi.nlm.nih.gov/articles/PMC13279268/", target = "_blank", rel = "noopener", "Plasma proteomic signatures of cellular aging"), "：当前注释并未训练细胞年龄时钟。")))
		)
	)
}

question_server <- function(id, catalogue, context, select_feature) {
	moduleServer(id, function(input, output, session) {
		data <- function(name) {
			c <- context()
			subset_context(catalogue()[[paste0("question_", name)]], c$Y, c$layer)
		}
		filtered <- function(d, key, values) {
			if (!nrow(d) || !key %in% names(d)) return(d[FALSE, , drop = FALSE])
			d[!is.na(d[[key]]) & d[[key]] %in% values, , drop = FALSE]
		}
		front <- function(d, cols) d[, unique(c(intersect(cols, names(d)), names(d))), drop = FALSE]
		selected_model <- reactive(paste(input$family, input$donors, input$budget, sep = "_"))
		panel_filter <- function(d, scope = TRUE) {
			d <- filtered(d, "model", c("Clinical", paste0("NS_", input$budget), selected_model()))
			if (scope) {
				d <- filtered(d, "landmark", as.numeric(input$landmark))
				d <- filtered(d, "stratum", input$stratum)
			}
			d
		}
		observeEvent(list(context()$Y, context()$layer, catalogue()), {
			d <- data("prediction")
			b <- sort(unique(d$budget[d$budget > 0]))
			updateSelectInput(session, "budget", choices = b, selected = if (10 %in% b) 10 else head(b, 1))
			updateSelectInput(session, "stratum", choices = unique(d$stratum), selected = "All")
			g <- data("genetic")
			anchors <- c("GDF15", "MMP12", "PCSK9", "IL6", "ApoB", "L_VLDL_TG.pct")
			a <- intersect(anchors, g$feature)
			updateSelectInput(session, "anchor", choices = c("选择案例" = "", setNames(a, a)), selected = "")
		}, ignoreNULL = FALSE)
		observeEvent(input$anchor, {
			if (nzchar(input$anchor)) {
				d <- data("genetic")
				z <- which(d$feature == input$anchor)
				if (length(z)) select_feature(d, z[1])
			}
		})
		output$cards <- renderUI({
			d <- data("overview")
			labels <- c("LE8 supervision", "ABM 分层", "遗传与实测差异")
			fluidRow(lapply(labels, function(label) {
				z <- d[d$question %in% if (label == "LE8 supervision") c("LE8 预测增益", "LE8 可解释性") else label, , drop = FALSE]
				column(4, div(class = "q-card", h4(label),
					if (!nrow(z)) p("暂无可用结果") else lapply(seq_len(nrow(z)), function(i)
						tagList(tags$strong(paste(z$layer[i], "·", z$question[i], "·", z$status[i])), p(z$result[i])))))
			}))
		})
		output$overview <- renderDT(make_dt(front(data("overview"), c("question", "layer", "status", "result", "limitation")), select = "none", page = 12))
		output$cohort <- renderDT(make_dt(front(data("cohort"), c("Y", "layer", "N_omics", "incident_events", "prevalent_cases", "features", "PGS_matched")), select = "none"))
		prediction <- reactive(panel_filter(data("prediction")))
		contrasts <- reactive({
			d <- panel_filter(data("contrasts"))
			filtered(d, "reference", paste0("NS_", input$budget))
		})
		proxy <- reactive(panel_filter(data("proxy"), FALSE))
		output$panel_note <- renderUI({
			d <- prediction()
			if (!nrow(d)) return(p("此预算/分层没有可用模型。"))
			bad <- d[d$effective_assays == 0 & d$budget > 0 & !is.na(d$effective_assays), , drop = FALSE]
			tagList(p(paste("当前模型：", selected_model(), "；Clinical、NS与当前panel保持相同验证设计。")),
				if (nrow(bad)) div(class = "q-note", paste("分子系数趋近零：", paste(unique(bad$model), collapse = ", "), "。名义 biom 数量不能解释为有效分子贡献。")))
		})
		output$prediction_plot <- renderPlot({
			d <- prediction()
			if (!nrow(d)) return(empty_plot("No matched-budget prediction estimates"))
			d$label <- paste(d$layer, d$model, sep = " | ")
			ggplot(d, aes(AUC, reorder(label, AUC))) + geom_errorbar(aes(xmin = AUC_lo, xmax = AUC_hi), orientation = "y", width = .15, na.rm = TRUE) +
				geom_point(aes(color = layer), size = 3) + theme_bw(base_size = 12) + labs(x = "IPCW AUC (95% conditional interval)", y = NULL, title = "同预算预测性能")
		})
		contrast_plot <- function(d, title) {
			if (!nrow(d) || !"delta_AUC" %in% names(d)) return(empty_plot("No paired contrast available"))
			d <- d[!is.na(d$comparison_valid) & d$comparison_valid, , drop = FALSE]
			if (!nrow(d)) return(empty_plot("Comparison conditions could not be verified"))
			d$label <- paste(d$layer, d$model, sep = " | ")
			ggplot(d, aes(delta_AUC, label)) + geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
				geom_errorbar(aes(xmin = delta_lo, xmax = delta_hi), orientation = "y", width = .15) + geom_point(size = 3, color = "#176b73") +
				theme_bw(base_size = 12) + labs(x = "Paired ΔAUC (nominal 95% interval)", y = NULL, title = title)
		}
		output$delta_plot <- renderPlot(contrast_plot(contrasts(), "当前 panel − NS"))
		output$contrasts <- renderDT(make_dt(front(contrasts(), c("layer", "model", "reference", "delta_AUC", "delta_lo", "delta_hi", "interval_status", "N", "actual_assays", "effective_assays", "identical_panel")), select = "none"))
		output$proxy_plot <- renderPlot({
			d <- proxy()
			if (!nrow(d)) return(empty_plot("No matched-budget LE8 reconstruction"))
			d$label <- paste(d$layer, d$pillar, sep = " | ")
			p <- ggplot(d, aes(delta_R2_vs_NS, label)) + geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
				geom_point(aes(color = layer), size = 3) + theme_bw(base_size = 12) + labs(x = "Held-out ΔR² vs NS", y = NULL, title = "LE8 proxy 保真度")
			if (all(c("lo", "hi") %in% names(d))) p <- p + geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = .15, na.rm = TRUE)
			p
		})
		output$proxy <- renderDT(make_dt(front(proxy(), c("layer", "pillar", "model", "N", "delta_R2_vs_NS", "lo", "hi", "FDR_all_proxy_contrasts", "uncertainty")), select = "none"))
		output$pillars <- renderDT(make_dt(panel_filter(data("pillars"), FALSE), select = "none"))
		members <- reactive(panel_filter(data("members"), FALSE))
		output$members <- renderDT(make_dt(members()))
		observeEvent(input$members_rows_selected, select_feature(members(), input$members_rows_selected))
		yin <- reactive({
			d <- data("contrasts")
			d <- filtered(d, "model", paste(input$family, "YinYang", input$budget, sep = "_"))
			d <- filtered(d, "reference", paste(input$family, "Yin", input$budget, sep = "_"))
			d <- filtered(d, "stratum", input$stratum)
			filtered(d, "landmark", as.numeric(input$landmark))
		})
		output$yin_plot <- renderPlot(contrast_plot(yin(), "YinYang − Yin：相同 panel 方法与预算"))
		output$yin <- renderDT(make_dt(yin(), select = "none"))
		output$heterogeneity <- renderDT(make_dt(filtered(filtered(data("heterogeneity"), "model", selected_model()), "landmark", as.numeric(input$landmark)), select = "none"))
		output$inflammation <- renderDT(make_dt(data("inflammation_definition"), select = "none"))
		output$design <- renderDT(make_dt(data("design"), select = "none"))
		for (key in c("concept_status","concept_coefficients","concept_fold_panels","cell_contrasts","cell_evidence","abm_training","abm_gate","abm_risk_gain","abm_audit","abm_decision","susie_pairs","susie_diagnostics","mr_signal_evidence","dandelion_native","dandelion_lolo","state_projection","decomposition","age_models")) local({
			k <- key; output[[k]] <- renderDT(make_dt(data(k),select="none"))
		})
		output$fit <- renderDT(make_dt(panel_filter(data("fit"), FALSE), select = "none"))

		abm <- reactive(filtered(data("abm_metrics"), "backend", input$backend))
		observeEvent(list(input$backend, context()$Y, context()$layer, catalogue()), {
			m <- unique(abm()$model)
			default <- if (input$backend == "reference") "elasticnet_weighted" else "tabicl_finetuned"
			updateSelectInput(session, "abm_model", choices = m, selected = if (default %in% m) default else head(m, 1))
		}, ignoreNULL = FALSE)
		pm <- reactive(filtered(abm(), "model", c(input$abm_model, "elasticnet", "clinical")))
		output$abm_note <- renderUI({
			d <- pm()
			if (!nrow(d)) return(p("该方案暂无完成结果。"))
			p(paste("原生运行：", paste(unique(d$run), collapse = "；")))
		})
		output$abm_plot <- renderPlot({
			d <- pm()
			if (!nrow(d)) return(empty_plot("No ABM metrics"))
			d$subset <- factor(d$subset, c("all", "supported", "rejected"))
			lines <- d[ave(seq_len(nrow(d)), interaction(d$layer, d$model), FUN = length) > 1, , drop = FALSE]
			ggplot(d, aes(subset, AUC_IPCW, color = model, group = model)) + geom_point(size = 3) + geom_line(data = lines) +
				facet_wrap( ~ layer) + theme_bw(base_size = 12) + labs(x = "Frozen support strata", y = "IPCW AUC", title = "所有人 / supported / rejected") + theme(legend.position = "bottom")
		})
		output$coverage_plot <- renderPlot({
			d <- filtered(filtered(filtered(data("abm_coverage"), "backend", input$backend), "model", c(input$abm_model, "elasticnet", "clinical")), "selector", input$abm_selector)
			k <- input$abm_metric
			if (!nrow(d) || !k %in% names(d)) return(empty_plot("No frozen coverage curve for this backend"))
			d$value <- d[[k]]
			ggplot(d, aes(coverage, value, color = model, group = model)) + geom_line() + geom_point() + facet_wrap( ~ layer) +
				theme_bw(base_size = 12) + labs(x = "Fraction of native test cohort", y = k, title = "覆盖率与性能：全部冻结阈值") + theme(legend.position = "bottom")
		})
		output$abm_metrics <- renderDT(make_dt(front(pm(), c("layer", "model", "subset", "coverage", "n", "cases", "observed_IPCW", "AUC_IPCW", "Brier_IPCW", "delta_AUC_vs_elasticnet", "Brier_gain_vs_elasticnet")), select = "none"))
		output$abm_paired <- renderDT(make_dt(filtered(filtered(data("abm_paired"), "backend", input$backend), "model", input$abm_model), select = "none"))
		output$abm_support <- renderDT(make_dt(filtered(data("abm_support"), "backend", input$backend), select = "none"))

		genetic <- reactive({
			d <- data("genetic")
			if (!nrow(d)) return(d)
			d <- d[tolower(as.character(d$same_people)) %in% c("true", "1"), , drop = FALSE]
			c <- context()
			if (!is.null(c$scope) && c$scope != "all") d <- filtered(d, "scope", c$scope)
			if (!is.null(c$adjustment) && c$adjustment != "all") d <- filtered(d, "adjustment", c$adjustment)
			if (identical(input$genetic_screen, "discordant")) d <- d[which(d$measured_p >= .05 & d$pgs_FDR < .05), , drop = FALSE]
			if (identical(input$genetic_screen, "difference")) d <- d[which(d$GR_difference_FDR < .05), , drop = FALSE]
			front(d, c("feature", "layer", "measured_beta", "measured_p", "pgs_beta", "pgs_FDR", "pgs_measured_r", "GR_beta_difference", "GR_difference_FDR", "N", "events"))
		})
		output$genetic <- renderDT(make_dt(genetic()))
		observeEvent(input$genetic_rows_selected, select_feature(genetic(), input$genetic_rows_selected))
		output$genetic_plot <- renderPlot({
			d <- genetic()
			if (!nrow(d)) return(empty_plot("No matched genetic/measured rows in this scope"))
			d$pattern <- ifelse(d$measured_p >= .05 & d$pgs_FDR < .05, "Measured P≥0.05 / PGS FDR<0.05", "Other matched assays")
			z <- d[d$feature %in% c("GDF15", "MMP12", "ApoB", "L_VLDL_TG.pct", context()$feature), , drop = FALSE]
			ggplot(d, aes(measured_beta, pgs_beta, color = pattern)) + geom_hline(yintercept = 0, color = "grey80") + geom_vline(xintercept = 0, color = "grey80") +
				geom_point(alpha = .45) + geom_text(data = z, aes(label = feature), vjust =  - 1, check_overlap = TRUE, size = 3) + facet_wrap( ~ layer) +
				theme_bw(base_size = 12) + labs(x = "Measured association (per own SD)", y = "PGS association (per own SD)", title = "同人群实测与 PGS：不把显著性差异当作效应差异检验") + theme(legend.position = "bottom")
		})
		feature_data <- function(name) {
			d <- data(name)
			key <- feature_column(d)
			if (is.null(key)) return(d[FALSE, , drop = FALSE])
			d[!is.na(d[[key]]) & feature_name(d[[key]]) == context()$feature, , drop = FALSE]
		}
		for (name in c("temporal", "modules", "mediation", "same_locus")) local({
			nm <- name
			output[[nm]] <- renderDT(make_dt(feature_data(nm), select = "none"))
		})
		for (name in c("mr_scope", "dandelion", "nonlinear", "cell_status", "sources")) local({
			nm <- name
			output[[nm]] <- renderDT(make_dt(data(nm), select = "none"))
		})
		output$cell <- renderDT({
			d <- data("cell")
			if ("FDR_all_panel_cell_tests" %in% names(d)) d <- d[order(d$FDR_all_panel_cell_tests, d$p), , drop = FALSE]
			make_dt(d, select = "none")
		})
		output$download <- downloadHandler(filename = function() paste0("LE8_", context()$Y, "_", context()$layer, "_", input$download_kind, ".csv"),
			content = function(file) {
				req(input$download_kind %in% question_names)
				data.table::fwrite(data(input$download_kind), file)
			})
	})
}


# Application UI
# Six scientific entries: C1--C5 plus Final; no C5 prediction menu remains.
detail_ui <- function(id) {
	ns <- NS(id)
	tagList(h4(textOutput(ns("title"))), textOutput(ns("warning")), tabsetPanel(tabPanel("Measured / PGS", plotOutput(ns("paired"),
		height = "380px"
	), DTOutput(ns("effect_rows"))), tabPanel("G / remainder", plotOutput(ns("components"),
		height = "380px"
	), plotOutput(ns("contrast"), height = "250px")), tabPanel(
		"Time / MR", plotOutput(ns("temporal"),
			height = "400px"
		), plotOutput(ns("prevalent"), height = "330px"), plotOutput(ns("mr"), height = "400px"),
		plotOutput(ns("reverse"), height = "330px")
	), tabPanel(
		"Colocalization", plotOutput(ns("coloc"), height = "380px"),
		DTOutput(ns("locus_rows"))
	), tabPanel("Related source rows", uiOutput(ns("source_picker")), DTOutput(ns("source_rows")))))
}
module_ui <- function(id, label) {
	ns <- NS(id)
	tabPanel(label, value = id, tabsetPanel(tabPanel(
		"Evidence table", uiOutput(ns("controls")), DTOutput(ns("evidence")),
		uiOutput(ns("comparison_ui")), detail_ui(paste0(id, "_detail"))
	), tabPanel(
		"Figures", uiOutput(ns("figure_select")),
		textOutput(ns("figure_path")), imageOutput(ns("figure"), height = "auto"), uiOutput(ns("pdf_note")), downloadButton(
			ns("download_figure"),
			"Download selected figure"
		)
	), tabPanel("Source tables", uiOutput(ns("table_select")), checkboxInput(
		ns("selected_only"),
		"Filter source rows to selected biomarker", FALSE
	), DTOutput(ns("raw")), textOutput(ns("raw_note")), downloadButton(
		ns("download_table"),
		"Download aggregate table"
	))))
}
ui <- fluidPage(
	tags$head(tags$style(HTML("body{max-width:1800px;margin:auto;padding:18px;color:#193142;background:#fafcfd;} .well{border-radius:5px;} .shiny-output-error-validation{color:#805600;} .dataTables_wrapper{margin-bottom:18px;} h2{font-weight:650;} .tabbable{margin-top:12px;} img{max-width:100%;} .q-intro{background:#eaf4f5;border-left:5px solid #176b73;padding:14px 24px;margin:18px 0;} .q-intro h2{font-size:25px;margin-top:6px;} .q-path{font-size:23px;letter-spacing:.3px;margin:12px 0;} .q-card{background:white;border:1px solid #d3e1e7;border-radius:8px;padding:15px 20px;margin:18px 0;min-height:180px;} .q-card h4{color:#176b73;} .q-note{background:#fff3d6;border-left:4px solid #bb8316;padding:12px;margin:12px 0;} .nav-tabs>li.active>a{font-weight:650;color:#176b73;}"))),
	titlePanel("LE8 | Molecular evidence explorer"), fluidRow(
		column(2, selectInput("Y", "Outcome", choices = character())),
		column(2, selectInput("layer", "Molecular layer", choices = c(Proteins = "prot", Metabolites = "met", Both = "both"))),
		column(4, selectizeInput("feature", "Selected biomarker", choices = NULL, options = list(placeholder = "Search an assay, e.g. Lactate or L_VLDL_TG.pct"))),
		column(2, selectInput("scope", "Comparison scope", choices = "existing_full")), column(2, selectInput("adjustment",
			"Adjustment",
			choices = "basic"
		))
	), fluidRow(
		column(10, p("Click a table row to update the evidence plots. Column filters and sorting do not refit models. PGS association, MR and colocalization remain separate evidence types.")),
		column(2, actionButton("reload", "Reload prepared index"))
	), tabsetPanel(id = "module", selected = "questions", question_ui("questions"), module_ui(
		"c1_correlate",
		"C1 Correlation"
	), module_ui("c2_cause", "C2 Causation"), module_ui("c3_coloc", "C3 Colocalization"), module_ui(
		"c4_connect",
		"C4 Connection"
	), module_ui("c5_cellulation", "C5 Cellulation"), module_ui("final", "Final"), tabPanel(
		"Data status",
		DTOutput("status"), DTOutput("counts"), verbatimTextOutput("root_info")
	))
)


# Application server
server <- function(input, output, session) {
	catalogue <- reactiveVal(load_catalogue())
	observeEvent(input$reload, {
		catalogue(load_catalogue())
		showNotification("Prepared aggregate index reloaded. Run final --index-only to incorporate newly added source files.",
			type = "message"
		)
	})
	chosen <- reactiveVal(NULL)
	observe({
		ys <- sort(unique(c(catalogue()$tables$Y, catalogue()$candidates$Y)))
		ys <- setdiff(ys, "all")
		updateSelectInput(session, "Y", choices = ys, selected = if (!is.null(isolate(input$Y)) && isolate(input$Y) %in%
			ys)
			isolate(input$Y) else head(ys, 1))
	})
	observeEvent(list(input$Y, input$layer, catalogue()),
		{
			d <- subset_context(catalogue()$candidates, input$Y, input$layer)
			fs <- sort(unique(d$feature))
			prev <- input$feature
			anchor <- intersect(c(if (identical(input$layer, "met")) "ApoB" else "GDF15", "L_VLDL_TG.pct", "PCSK9"), fs)
			updateSelectizeInput(session, "feature", choices = fs, selected = if (!is.null(prev) && prev %in% fs)
				prev else if (length(anchor)) anchor[1] else head(fs, 1), server = TRUE)
			scopes <- unique(d$scope)
			adj <- unique(d$adjustment)
			updateSelectInput(session, "scope",
				choices = c(`All available scopes` = "all", setNames(scopes, scopes)),
				selected = if ("existing_full" %in% scopes)
					"existing_full" else "all"
			)
			updateSelectInput(session, "adjustment", choices = c(`All adjustments` = "all", setNames(adj, adj)), selected = if ("basic" %in%
				adj)
				"basic" else "all")
		},
		ignoreNULL = FALSE
	)
	context <- reactive(list(Y = input$Y, layer = input$layer, feature = input$feature, scope = input$scope, adjustment = input$adjustment))
	dataset_context <- reactive(list(Y = input$Y, layer = input$layer))
	comparison_context <- reactive(list(Y = input$Y, layer = input$layer, scope = input$scope, adjustment = input$adjustment))
	filter_scope <- function(d, ctx) {
		if (!nrow(d))
			return(d)
		if ("scope" %in% names(d) && !is.null(ctx$scope) && ctx$scope != "all")
			d <- d[!is.na(d$scope) & d$scope == ctx$scope, , drop = FALSE]
		if ("adjustment" %in% names(d) && !is.null(ctx$adjustment) && ctx$adjustment != "all")
			d <- d[!is.na(d$adjustment) & d$adjustment == ctx$adjustment, , drop = FALSE]
		d
	}
	select_feature <- function(d, i) {
		if (length(i) != 1 || i < 1 || i > nrow(d))
			return()
		key <- feature_column(d)
		if (is.null(key))
			return()
		f <- feature_name(d[[key]][i])
		if (is.na(f) || !nzchar(f))
			return()
		chosen(f)
		dc <- dataset_context()
		available <- subset_context(catalogue()$candidates, dc$Y, dc$layer)
		# With server-side Selectize, an assay need not be in the browser's loaded
		# option page. Supply the full catalogue when selecting a clicked row/anchor.
		if (!f %in% available$feature) return()
		updateSelectizeInput(session, "feature", choices = sort(unique(available$feature)), selected = f, server = TRUE)
		# Locus/MR rows often have no association scope. Never clear a valid matched-comparison selection.
		if ("scope" %in% names(d) && !is.na(d$scope[i]) && as.character(d$scope[i]) %in% available$scope)
			updateSelectInput(session, "scope", selected = as.character(d$scope[i]))
		if ("adjustment" %in% names(d) && !is.na(d$adjustment[i]) && as.character(d$adjustment[i]) %in% available$adjustment)
			updateSelectInput(session, "adjustment", selected = as.character(d$adjustment[i]))
	}
	detail_server <- function(id) {
		force(id)
		moduleServer(id, function(input, output, session) {
			ctx <- context
			all_effects <- reactive({
				c <- ctx()
				req(c$Y, c$layer, c$feature)
				d <- subset_context(catalogue()$effects, c$Y, c$layer)
				d[d$feature == c$feature & !is.na(d$feature), , drop = FALSE]
			})
			paired <- reactive({
				d <- all_effects()
				d <- d[d$kind == "measured_vs_PGS", , drop = FALSE]
				filter_scope(d, ctx())
			})
			output$title <- renderText(paste(ctx()$Y, ctx()$layer, "|", ctx()$feature))
			output$warning <- renderText("Evidence is exploratory. A null measured association does not prove no effect; a significant PGS association is not MR. Clinical deployment and intervention effects are not established here.")
			output$paired <- renderPlot({
				plot_effects(paired(), "Matched measured and PGS associations")
			})
			output$effect_rows <- renderDT(make_dt(paired(), select = "none"))
			output$components <- renderPlot({
				d <- all_effects()
				d <- d[d$kind == "G_R", , drop = FALSE]
				if ("model" %in% names(d) && any(d$model == "joint"))
					d <- d[d$model == "joint", , drop = FALSE]
				# G/R is only defined at the adjustment saved by that analysis; do not silently substitute a
				# different adjustment.
				if (!is.null(ctx()$scope) && ctx()$scope != "all")
					d <- d[d$scope == ctx()$scope, , drop = FALSE]
				plot_effects(d, "Calibrated genetic component and remaining measured component")
			})
			output$contrast <- renderPlot({
				d <- all_effects()
				d <- d[d$kind == "G_minus_R", , drop = FALSE]
				if (!is.null(ctx()$scope) && ctx()$scope != "all")
					d <- d[d$scope == ctx()$scope, , drop = FALSE]
				plot_effects(d, "Direct G minus R contrast")
			})
			output$temporal <- renderPlot({
				d <- all_effects()
				d <- d[d$kind %in% c("landmark", "event_window"), , drop = FALSE]
				if (nrow(d)) {
					d$series <- paste(d$kind, d$series, sep = " | ")
					if (nrow(d) > 60)
						d <- head(d, 60)
				}
				plot_effects(d, "Diagnosis-time windows / landmark associations (distinct analyses)")
			})
			output$prevalent <- renderPlot({
				d <- all_effects()
				d <- d[d$kind == "prevalent_window", , drop = FALSE]
				plot_effects(head(d, 40), "Baseline prevalent associations: log odds, not prospective hazards")
			})
			output$mr <- renderPlot({
				d <- all_effects()
				d <- d[d$kind == "MR", , drop = FALSE]
				plot_effects(head(d, 40), "MR estimates: original exposure scales")
			})
			output$reverse <- renderPlot({
				d <- all_effects()
				d <- d[d$kind == "MR_reverse", , drop = FALSE]
				plot_effects(head(d, 40), "Disease-liability to biomarker MR (not post-diagnosis treatment)")
			})
			locus <- reactive({
				c <- ctx()
				d <- subset_context(catalogue()$loci, c$Y, c$layer)
				d[d$feature == c$feature & !is.na(d$feature), , drop = FALSE]
			})
			output$locus_rows <- renderDT(make_dt(locus(), select = "none"))
			output$coloc <- renderPlot({
				d <- locus()
				if (!nrow(d))
					return(empty_plot("No available colocalization result for this assay."))
				d <- head(d, 30)
				g <- rbind(data.frame(locus = d$locus, value = d$PP_H4, prior = "Reported default"), data.frame(
					locus = d$locus,
					value = d$PP_H4_robust_min, prior = "Minimum across all required priors"
				))
				g$locus <- factor(g$locus, levels = rev(unique(g$locus)))
				ggplot(g, aes(value, locus, shape = prior)) +
					geom_point(
						position = position_dodge(width = 0.5),
						size = 3
					) +
					xlim(0, 1) +
					labs(
						title = "Regional shared-signal posterior", x = "PP(H4)", y = NULL,
						caption = "Rows may be overlapping loci for correlated assays. No independent causal-locus count is inferred."
					) +
					theme_bw(base_size = 12)
			})
			relevant <- reactive({
				c <- dataset_context()
				r <- subset_context(catalogue()$tables, c$Y, c$layer)
				r <- r[r$role != "browse_aggregate" | grepl("pgs|coloc|mediat|proxy|membership|MR|direction|temporal|genetic|cell",
					r$path,
					ignore.case = TRUE
				), , drop = FALSE]
				r
			})
			output$source_picker <- renderUI({
				r <- relevant()
				selectInput(session$ns("source"), "Evidence source", choices = setNames(r$source_id, r$path))
			})
			detail_rows <- reactive({
				r <- relevant()
				req(input$source)
				hit <- r[r$source_id == input$source, , drop = FALSE]
				req(nrow(hit) == 1)
				d <- read_aggregate(hit$path, catalogue()$tables)
				k <- feature_column(d)
				if (is.null(k))
					return(data.frame(status = "This is a model-level table with no assay key."))
				d <- d[feature_name(d[[k]]) == ctx()$feature & !is.na(d[[k]]), , drop = FALSE]
				d
			})
			output$source_rows <- renderDT({
				make_dt(detail_rows(), select = "none")
			})
		})
	}
	module_server <- function(id, module) {
		force(id)
		force(module)
		moduleServer(id, function(input, output, session) {
			registry <- reactive({
				dc <- dataset_context()
				r <- subset_context(catalogue()$tables, dc$Y, dc$layer)
				r[r$module == module, , drop = FALSE]
			})
			output$controls <- renderUI({
				if (module %in% c("c1_correlate", "c2_cause", "final"))
					selectInput(session$ns("screen"), "Discovery filter", choices = c(
						`All assays — no measured-association gate` = "all",
						`Measured P ≥ 0.05; matched PGS FDR < 0.05` = "pgs_nominal", `Measured FDR ≥ 0.05; matched PGS FDR < 0.05` = "pgs_fdr",
						`Measured P ≥ 0.05; cis/local MR FDR < 0.05` = "mr_null"
					)) else if (module == "c4_connect")
					p("YS−NS held-out reconstruction contrasts at the same assay budget. CI is unavailable unless paired prediction-level bootstrap was run.") else if (module == "c5_cellulation")
					p("External cell-expression annotation is not secretion tracing or a cellular-aging clock. Select an annotation table below for detailed results.") else p("Sort PP(H4) or its minimum across priors; click an assay to see its matched associations and genetic-component results.")
			})
			evidence <- reactive({
				c <- comparison_context()
				if (module == "c3_coloc")
					return(subset_context(catalogue()$loci, c$Y, c$layer))
				if (module == "c4_connect")
					return(subset_context(catalogue()$proxy_comparisons, c$Y, c$layer))
				if (module == "c5_cellulation") {
					r <- registry()
					a <- r[grepl("annotation|enrichment", r$path, ignore.case = TRUE), , drop = FALSE]
					if (!nrow(a)) {
						s <- r[grepl("c5[.]cellulation_status[.]csv$", r$path), , drop = FALSE]
						if (nrow(s))
							return(read_aggregate(s$path[1], catalogue()$tables))
						return(data.frame(status = "No available cell-expression result for this molecular layer."))
					}
					# Prefer the current per-trait/per-layer C5 export over historical joint-panel annotations.
					a <- a[order(!grepl("/c5_cellulation/", a$path), !grepl("[.]enrichment[.]csv$", a$path)), , drop = FALSE]
					return(read_aggregate(a$path[1], catalogue()$tables))
				}
				d <- filter_scope(subset_context(catalogue()$candidates, c$Y, c$layer), c)
				if (!nrow(d))
					return(d)
				col <- switch(input$screen %||% "all",
					pgs_nominal = "measured_P_ge_005_PGS_FDR_lt_005",
					pgs_fdr = "measured_FDR_ge_005_PGS_FDR_lt_005",
					mr_null = "measured_P_ge_005_cis_MR_FDR_lt_005",
					NULL
				)
				if (!is.null(col) && col %in% names(d))
					d <- d[!is.na(d[[col]]) & d[[col]], , drop = FALSE]
				front <- intersect(c(
					"feature", "Y", "layer", "measured_p", "measured_FDR", "pgs_p", "pgs_FDR",
					"pgs_measured_r", "best_reported_cis_MR_FDR_all", "GR_beta_difference", "GR_difference_FDR",
					"PP_H4_robust_min", "best_reported_coloc_region", "N", "events", "same_people", "evidence_pattern",
					"scope", "adjustment"
				), names(d))
				d[, c(front, setdiff(names(d), front)), drop = FALSE]
			})
			output$comparison_ui <- renderUI({
				if (module == "c4_connect")
					plotOutput(session$ns("comparison_plot"), height = "420px")
			})
			output$comparison_plot <- renderPlot({
				req(module == "c4_connect")
				d <- evidence()
				if (!nrow(d) || !all(c("model", "component", "delta_R2_vs_NS") %in% names(d)))
					return(empty_plot("No paired-budget reconstruction table available."))
				i <- input$evidence_rows_selected
				if (length(i) == 1 && i <= nrow(d))
					selected <- d$model[i] else selected <- if ("YS_YinYang_10" %in% d$model)
					"YS_YinYang_10" else d$model[1]
				z <- d[d$model == selected, , drop = FALSE]
				z$label <- paste(z$layer, z$component, sep = " | ")
				z$label <- factor(z$label, levels = rev(unique(z$label)))
				p <- ggplot(z, aes(delta_R2_vs_NS, label)) +
					geom_vline(xintercept = 0, linetype = 2) +
					geom_point(size = 3) +
					labs(
						title = paste(selected, "minus same-budget NS"), x = "Held-out R² difference: basic covariates + panel",
						y = NULL, caption = "All available LE8 targets, including negative differences. No CI means point estimates only; this is reconstruction, not intervention response."
					) +
					theme_bw(base_size = 12) +
					theme(plot.caption = element_text(hjust = 0, size = 9))
				if (all(c("lo", "hi") %in% names(z)))
					p <- p + geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.16, na.rm = TRUE)
				p
			})
			output$evidence <- renderDT(make_dt(evidence()))
			# DT returns indexes in the server data, not the visual order. This stays correct after sort/filter.
			observeEvent(input$evidence_rows_selected, select_feature(evidence(), input$evidence_rows_selected))
			output$table_select <- renderUI({
				r <- registry()
				selectInput(session$ns("table"), "Aggregate table", choices = setNames(r$source_id, paste(r$role,
					r$path,
					sep = " | "
				)))
			})
			raw_full <- reactive({
				r <- registry()
				req(input$table)
				hit <- r[r$source_id == input$table, , drop = FALSE]
				req(nrow(hit) == 1)
				read_aggregate(hit$path, catalogue()$tables)
			})
			raw <- reactive({
				d <- raw_full()
				k <- feature_column(d)
				if (isTRUE(input$selected_only) && !is.null(k))
					d <- d[feature_name(d[[k]]) == context()$feature & !is.na(d[[k]]), , drop = FALSE]
				d
			})
			output$raw <- renderDT(make_dt(raw()))
			observeEvent(input$raw_rows_selected, select_feature(raw(), input$raw_rows_selected))
			output$raw_note <- renderText("Numeric column sorting uses unrounded values. Downloads contain aggregate rows only. An empty filtered table is not a negative scientific result.")
			output$download_table <- downloadHandler(filename = function() paste0(module, "_aggregate.csv"), content = function(file) {
				data.table::fwrite(raw(), file)
			})
			figs <- reactive({
				dc <- dataset_context()
				r <- subset_context(catalogue()$figures, dc$Y, dc$layer)
				r[r$module == module, , drop = FALSE]
			})
			output$figure_select <- renderUI({
				r <- figs()
				selectInput(session$ns("figure_file"), "Saved figure", choices = setNames(r$path, r$name))
			})
			output$figure_path <- renderText(input$figure_file %||% "No saved figure in this view")
			output$figure <- renderImage(
				{
					req(input$figure_file)
					validate(need(grepl("[.](png|jpg|jpeg)$", input$figure_file, ignore.case = TRUE), "PDF selected: use Download selected figure."))
					f <- resolve_source(input$figure_file, catalogue()$figures)
					list(
						src = f, contentType = if (grepl("[.]png$", f)) "image/png" else "image/jpeg", width = "100%",
						alt = basename(f)
					)
				},
				deleteFile = FALSE
			)
			output$pdf_note <- renderUI({
				if (!is.null(input$figure_file) && grepl("[.]pdf$", input$figure_file, ignore.case = TRUE))
					p("PDF is available by download; the result directory is not published as a web resource.")
			})
			output$download_figure <- downloadHandler(filename = function() basename(input$figure_file), content = function(file) {
				req(input$figure_file)
				file.copy(resolve_source(input$figure_file, catalogue()$figures), file, overwrite = TRUE)
			})
		})
		detail_server(paste0(id, "_detail"))
	}
	for (m in c("c1_correlate", "c2_cause", "c3_coloc", "c4_connect", "c5_cellulation", "final")) module_server(
		m,
		m
	)
	question_server("questions", catalogue, context, select_feature)
	detail_server("questions_detail")
	output$status <- renderDT(make_dt(catalogue()$status, select = "none"))
	output$counts <- renderDT(make_dt(catalogue()$discovery_counts, select = "none"))
	output$root_info <- renderText(paste(
		"Read-only analysis root:", LE8_ROOT, "\nPrepared index:", LE8_INDEX,
		"\nNo participant-level data or analysis execution is exposed."
	))
}

if ("--review" %in% commandArgs(trailingOnly = TRUE)) {
	c <- load_catalogue()
	stopifnot(all(c("candidates", "effects", "loci", "tables", "figures") %in% names(c)))
	for (path in c$tables$path) resolve_source(path, c$tables)
	for (path in c$figures$path) resolve_source(path, c$figures)
	message("Index and table allowlist checked; no participant data loaded.")
	quit(status = 0)
}
port <- as.integer(option("--port", Sys.getenv("LE8_SHINY_PORT", "3839")))
if (is.na(port) || port < 1024 || port > 65535) stop("Invalid LE8_SHINY_PORT")
host <- option("--host", Sys.getenv("LE8_SHINY_HOST", "127.0.0.1"))
if (!host %in% c("127.0.0.1", "localhost", "::1") && Sys.getenv("LE8_SHINY_ALLOW_NETWORK", "FALSE") != "TRUE") stop("Network exposure requires LE8_SHINY_ALLOW_NETWORK=TRUE and your own authentication/reverse proxy.")
shiny::runApp(shiny::shinyApp(ui, server), host = host, port = port, launch.browser = FALSE)
