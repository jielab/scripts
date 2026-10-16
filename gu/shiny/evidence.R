# Window-based IBDmix support candidates. Raw calls remain in the method tab.
gu_ibdmix_evidence_cached <- function(Q, dataset, build, db_path,
        cache_dir = file.path("/tmp", "gu-shiny-evidence")) {
    helper <- normalizePath(file.path(app_dir, "..", "f", "ibdmix.py"))
    log <- tempfile("gu-evidence-stderr-")
    on.exit(unlink(log), add = TRUE)
    path <- system2("python3", c(shQuote(helper), "loci", "--database", shQuote(db_path),
        "--dataset", shQuote(dataset), "--build", shQuote(build),
        "--cache-dir", shQuote(cache_dir)), stdout = TRUE, stderr = log)
    status <- attr(path, "status")
    if ((!is.null(status) && status != 0L) || length(path) != 1L || !file.exists(path)) {
        detail <- paste(tail(readLines(log, warn = FALSE), 4), collapse = " ")
        stop("IBDmix 窗口汇总失败：", detail)
    }
    result <- jsonlite::fromJSON(path)
    if (!is.data.frame(result$loci)) {
        chars <- c("segment_code", "chr", "source_class", "genome_build", "reference_names")
        nums <- c("rank", "start", "end", "length_bp", "n_tested", "peak_start", "peak_end",
            "peak_fraction", "peak_raw_fraction", "peak_effective_copies", "mean_fraction",
            "background_fraction", "extension_threshold", "core_threshold", "core_windows",
            "n_windows", "excess_support_kb", "n_carriers", "n_references", "max_lod", "n_calls")
        result$loci <- as.data.frame(c(setNames(rep(list(character()), length(chars)), chars),
            setNames(rep(list(numeric()), length(nums)), nums)))
    }
    result$directory <- dirname(path)
    result
}

gu_ibdmix_evidence_filter <- function(d, lineage = "all", chr = "ALL", minimum_pct = 2) {
    minimum_pct <- suppressWarnings(as.numeric(minimum_pct))
    if(length(minimum_pct) != 1L || !is.finite(minimum_pct)) minimum_pct <- 2
    d <- d[(lineage == "all" | d$source_class == lineage) & (chr == "ALL" | d$chr == chr) &
        100 * d$peak_fraction >= max(2, minimum_pct), , drop = FALSE]
    # Mirror the scientific aggregation's ranking before taking the display shortlist.
    d[order(-d$excess_support_kb, -d$peak_fraction, as.integer(d$chr),
        d$start, d$source_class), , drop = FALSE]
}

gu_ibdmix_evidence_server <- function(input, output, session, Q, dataset, build,
        stamp, db_path, select_record, open_record) {
    result <- reactive({
        stamp()
        withProgress(message = "计算 IBDmix 窗口支持峰（首次需要汇总片段）", {
            gu_ibdmix_evidence_cached(Q, dataset(), build(), db_path)
        })
    })
    candidates <- reactive({
        gu_ibdmix_evidence_filter(result()$loci, input$density_lineage %||% "all",
            input$ibdmix_evidence_chr %||% "ALL", input$ibdmix_evidence_peak %||% 2)
    })
    loci <- reactive(head(candidates(), 50L))
    output$ibdmix_evidence_note <- renderText({
        meta <- result()$metadata
        warning <- if(length(meta$warnings)) paste0(" 数据范围提示：", paste(meta$warnings, collapse = "; ")) else ""
        paste0("Top 50 supported loci：显示当前筛选范围内支持最强的 ", nrow(loci()), " 个候选位点；",
            format(nrow(candidates()), big.mark = ","), " 个区间通过窗口规则（",
            dataset(), " · ", build(), "）。50 kb 窗口，3 窗口中位数平滑；背景门槛 = max(1%, 该谱系全基因组支持率)，",
            "核心 ≥ 2 × 背景且连续 ≥ 2 窗口。支持强度 = Σ(平滑支持率 − 背景门槛) × 窗口长度(kb)，",
            "从高到低排序，同分时按平滑峰值排序；完整候选集可下载。", warning)
    })
    output$ibdmix_evidence_loci <- renderDT({
        d <- loci()
        region <- if(nrow(d)) paste0(format(d$start + 1, scientific = FALSE, trim = TRUE), "–",
            format(d$end, scientific = FALSE, trim = TRUE)) else character()
        show <- data.frame(ID = d$segment_code, rank = d$rank,
            strength = round(d$excess_support_kb, 2), Chr = d$chr, region = region,
            lineage = d$source_class, kb = round(d$length_bp / 1000, 1),
            peak = round(100 * d$peak_fraction, 2), mean = round(100 * d$mean_fraction, 2),
            threshold = round(100 * d$core_threshold, 2), copies = round(d$peak_effective_copies, 1),
            tested = d$n_tested,
            core = d$core_windows, carriers = d$n_carriers,
            names = d$reference_names, LOD = round(d$max_lod, 2), check.names = FALSE)
        names(show) <- c("ID", "全局排名", "IBDmix 支持强度 (支持率×kb)", "Chr", paste0("候选区间 (", build(), ")"), "谱系", "长度 (kb)",
            "平滑峰值 (%)", "区间均值 (%)", "核心门槛 (%)", "峰窗等效拷贝数", "检测人数",
            "核心窗口数", "区间携带者", "参考来源", "最高 LOD")
        DT::datatable(show, rownames = FALSE, class = "compact stripe hover nowrap",
            selection = list(mode = "single", selected = NULL),
            options = list(pageLength = 10, lengthMenu = c(10, 25, 50), scrollX = TRUE,
                order = list(list(1, "asc")),
                columnDefs = list(list(targets = 0, visible = FALSE)),
                language = list(emptyTable = "当前范围没有通过窗口支持规则的候选位点")),
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
    windows <- reactive({
        as.data.frame(data.table::fread(file.path(result()$directory, "windows.tsv"), showProgress = FALSE))
    })
    output$ibdmix_evidence_profile <- renderPlotly({
        d <- loci(); req(nrow(d))
        i <- match(input$ibdmix_evidence_select %||% "", d$segment_code)
        r <- d[if(is.na(i)) 1L else i, , drop = FALSE]
        w <- windows()
        w <- w[w$chr == r$chr & w$source_class == r$source_class &
            w$end > r$start - 150000 & w$start < r$end + 150000, , drop = FALSE]
        p <- plotly::plot_ly(w, x = ~((start + end) / 2 / 1e6), y = ~(100 * raw_fraction),
            type = "scatter", mode = "lines", name = "原始窗口支持率", line = list(color = "#a5b1ba"))
        p <- plotly::add_trace(p, y = ~(100 * smooth_fraction), name = "3 窗口中位数", line = list(color = "#de7b20"))
        plotly::layout(p,
            title = list(text = paste0(r$source_class, " · chr", r$chr, " · 窗口支持曲线"), font = list(size = 14)),
            xaxis = list(title = "位置 (Mb)"), yaxis = list(title = "保守拷贝支持率 (%)", rangemode = "tozero"),
            margin = list(t = 35, b = 45), legend = list(orientation = "h", y = -0.3),
            shapes = list(
                list(type = "rect", x0 = r$start / 1e6, x1 = r$end / 1e6, y0 = 0, y1 = 1,
                    yref = "paper", fillcolor = "rgba(222,123,32,0.08)", line = list(width = 0), layer = "below"),
                list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 100 * r$extension_threshold,
                    y1 = 100 * r$extension_threshold, line = list(color = "#999999", dash = "dot")),
                list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = 100 * r$core_threshold,
                    y1 = 100 * r$core_threshold, line = list(color = "#a74b1b", dash = "dash"))))
    })
    output$ibdmix_evidence_download <- downloadHandler(
        filename = function()paste0("ibdmix_top50_supported_loci_", dataset(), "_", build(), ".tsv"),
        content = function(file)write.table(loci(), file, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
    output$ibdmix_candidates_download <- downloadHandler(
        filename = function()paste0("ibdmix_all_window_candidates_", dataset(), "_", build(), ".tsv"),
        content = function(file)file.copy(file.path(result()$directory, "candidates.tsv"), file, overwrite = TRUE))
    invisible(loci)
}
