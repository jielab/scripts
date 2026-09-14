server <- function(input, output, session) {
  state <- reactiveValues(build = opt$grch, chr = "All", range = c(1, sum(bplot_lengths(opt$grch))),
    ready = FALSE, loading = FALSE, error = NULL, progress = "准备缓存…",
    page = "overview", selection = NULL, overview = NULL, block_view = NULL)
  job <- NULL; target <- opt$grch; progress <- NULL; log <- NULL
  set_region <- function(chr, region = NULL) {
    lengths <- bplot_lengths(state$build)
    if (!chr %in% c("All", names(lengths))) return(invisible(NULL))
    limit <- if (chr == "All") sum(lengths) else unname(lengths[chr])
    if (is.null(region)) region <- c(1, limit)
    region <- round(pmax(1, pmin(limit, as.numeric(region))))
    if (length(region) != 2L || any(!is.finite(region)) || region[2] <= region[1]) return(invisible(NULL))
    state$chr <- chr; state$range <- region
  }
  finish <- function(build) {
    state$build <- build; state$loading <- FALSE; state$ready <- TRUE; state$error <- NULL
    state$selection <- NULL; state$overview <- NULL; state$block_view <- NULL
    state$page <- "overview"; bslib::nav_select("nav", "overview", session = session)
    set_region("All"); updateSelectInput(session, "build", selected = build)
  }
  prepare <- function(build) {
    if (!is.null(job) && job$is_alive()) job$kill_tree()
    target <<- build; state$error <- NULL
    if (all(vapply(seq_len(nrow(opt$tracks)), function(i) !is.null(bplot_cache_meta(opt, i, build)), logical(1)))) {
      finish(build); return(invisible(NULL))
    }
    state$loading <- TRUE; state$progress <- paste0("准备 GRCh", build, " 缓存…")
    session$sendCustomMessage("bplot_busy", TRUE)
    directory <- file.path(opt$output_dir, "cache", paste0("grch", build))
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    stamp <- paste0(session$token, "-", build, "-", format(Sys.time(), "%H%M%OS6"))
    progress <<- file.path(run_dir, paste0(stamp, ".json"))
    log <<- file.path(run_dir, paste0(stamp, ".log"))
    tryCatch({
      if (build != opt$source_build && !file.exists(bplot_chain(opt, build))) stop("Missing chain: ", bplot_chain(opt, build))
      job <<- processx::process$new("flock", c(file.path(directory, "prepare.lock"),
        file.path(R.home("bin"), "Rscript"), "--vanilla", file.path(app_dir, "app.R"),
        "--prepare-config", config, "--build", build, "--progress", progress),
        stdout = log, stderr = "2>&1", cleanup_tree = TRUE)
    }, error = function(e) {
      state$loading <- FALSE; state$error <- conditionMessage(e)
      updateSelectInput(session, "build", selected = state$build)
    })
  }
  observeEvent(input$build, {
    if (!input$build %in% c("37", "38")) return()
    if (state$ready && !state$loading && identical(input$build, state$build)) return()
    prepare(input$build)
  }, ignoreInit = FALSE)
  observe({
    if (!state$loading) return()
    invalidateLater(750, session)
    if (!is.null(progress) && file.exists(progress)) {
      msg <- tryCatch(jsonlite::read_json(progress), error = function(e) NULL)
      if (!is.null(msg$stage)) state$progress <- msg$stage
    }
    if (!is.null(job) && !job$is_alive()) {
      ok <- identical(job$get_exit_status(), 0L) &&
        all(vapply(seq_len(nrow(opt$tracks)), function(i) !is.null(bplot_cache_meta(opt, i, target)), logical(1)))
      if (ok) finish(target) else {
        state$loading <- FALSE; state$error <- paste0("缓存准备失败。日志：", log)
        updateSelectInput(session, "build", selected = state$build)
      }
    }
  })
  session$onSessionEnded(function() if (!is.null(job) && job$is_alive()) job$kill_tree())
  observeEvent(input$chr, { if (!state$loading) set_region(input$chr) })
  observeEvent(input$reset, { if (!state$loading) set_region(state$chr) })
  zoom <- function(factor) {
    if (state$loading || state$chr == "All") return()
    center <- mean(state$range); radius <- diff(state$range) * factor / 2
    set_region(state$chr, center + c(-radius, radius))
  }
  observeEvent(input$zoom_in, zoom(0.5))
  observeEvent(input$zoom_out, zoom(2))
  observeEvent(input$go, {
    match <- regmatches(gsub(",", "", trimws(input$region)), regexec("^(?:chr)?([0-9]+|X|x):([0-9]+)-([0-9]+)$", gsub(",", "", trimws(input$region))))[[1]]
    if (length(match) != 4L) { showNotification("输入区域，例如 chr1:10000000-12000000", type = "message"); return() }
    if (!state$loading) set_region(bplot_chr(match[2]), as.numeric(match[3:4]))
  })
  observeEvent(input$plot_event, {
    event <- input$plot_event
    if (state$loading || !identical(event$page, state$page) || !identical(event$build, state$build) || !identical(event$chr, state$chr)) return()
    if (identical(event$type, "block")) {
      if (!event$race %in% opt$tracks$race) return()
      selected <- bplot_triplet(opt$block_dir, event$race, state$build, bplot_chr(event$chrom), event$id)
      if (is.null(selected)) return()
      if (state$page == "overview") state$overview <- list(chr = state$chr, range = state$range)
      state$selection <- selected; state$page <- "blockview"
      set_region(selected$chr, c(selected$start, selected$end))
      state$block_view <- list(chr = state$chr, range = state$range)
      bslib::nav_select("nav", "blockview", session = session)
      return()
    }
    if (identical(event$type, "click")) {
      set_region(bplot_chr(event$chrom), as.numeric(event$pos) + c(-500000, 500000)); return()
    }
    if (isTRUE(event$reset)) { set_region(state$chr); return() }
    region <- as.numeric(c(event$start, event$end))
    if (length(region) != 2L || any(!is.finite(region))) return()
    if (state$chr == "All") {
      lengths <- bplot_lengths(state$build); offsets <- c(0, head(cumsum(lengths), -1))
      i <- findInterval(mean(region), offsets)
      i <- max(1L, min(length(lengths), i))
      set_region(names(lengths)[i], region - offsets[i])
    } else set_region(state$chr, region)
  })
  observeEvent(input$igv_event, {
    event <- input$igv_event
    if (state$loading || !identical(event$build, state$build)) return()
    if (identical(event$chr, "All")) set_region("All") else
      set_region(bplot_chr(event$chr), c(event$start, event$end))
  })
  observe({
    req(state$ready, !state$loading)
    build <- state$build; chr <- state$chr; region <- state$range
    session$sendCustomMessage("bplot_view", list(build = build, chr = chr, start = region[1], end = region[2]))
    genome <- if (build == "37") "hg19" else "hg38"
    ref <- refs$references[[genome]]
    if (is.null(ref)) {
      session$sendCustomMessage("bplot_igv", list(error = paste0("IGV 需要有效的 GRCH", build, ".fasta 和 .fai；GWAS 视图可继续使用。")))
    } else {
      ref$id <- genome; ref$name <- paste0("GRCh", build)
      session$sendCustomMessage("bplot_igv", list(build = build, reference = ref,
        locus = bplot_locus(chr, region),
        genes = if (file.exists(file.path(genes, paste0(build, ".bed")))) paste0("bplot_genes/", build, ".bed") else NULL))
    }
  })
  output$status <- renderText({
    if (state$loading) return(paste0(state$progress, " · 首次读取或转换需要一些时间"))
    if (!is.null(state$error)) return(state$error)
    paste0("GRCh", state$build, " · ", if (state$chr == "All") "All chromosomes" else paste0("chr", state$chr, ":",
      format(state$range[1], big.mark = ",", scientific = FALSE), "–", format(state$range[2], big.mark = ",", scientific = FALSE)),
      " · 每轨最多 ", format(opt$max_points, big.mark = ","), " 个显示点")
  })
  output$panel <- renderUI({
    if (state$loading || !state$ready) return(div(class = "loading-panel", if (!is.null(state$error)) state$error else "正在准备 GWAS 数据…"))
    plotly::plotlyOutput("plot", height = paste0(max(350L, nrow(opt$tracks) * 175L + 65L), "px"))
  })
  output$plot <- plotly::renderPlotly({
    req(state$ready, !state$loading)
    valid <- vapply(seq_len(nrow(opt$tracks)), function(i)
      !is.null(bplot_cache_meta(opt, i, state$build)), logical(1))
    validate(need(all(valid), "GWAS 或缓存已改变，请重新启动 bplot 以更新视图。"))
    views <- bplot_view(opt, state$build, state$chr, state$range, isTRUE(input$blocks))
    bplot_plot(views, state$build, state$chr, state$range, opt$p_threshold)
  })
  observeEvent(input$nav, {
    if (identical(input$nav, state$page)) return()
    current <- list(chr = state$chr, range = state$range)
    if (state$page == "overview") state$overview <- current else state$block_view <- current
    state$page <- input$nav
    next_view <- if (input$nav == "overview") state$overview else state$block_view
    if (!is.null(next_view)) set_region(next_view$chr, next_view$range)
  })
  observeEvent(input$back, bslib::nav_select("nav", "overview", session = session))
  observeEvent(input$block_home, {
    req(state$selection, !state$loading)
    set_region(state$selection$chr, c(state$selection$start, state$selection$end))
  })
  observeEvent(input$block_in, zoom(.5))
  observeEvent(input$block_out, zoom(2))
  output$block_title <- renderText({
    if (is.null(state$selection)) return("在 Overview 中双击一个 block")
    paste0("block: ", state$selection$race, " ", state$selection$id, " · GRCh", state$build,
      " · blocks ", paste(state$selection$ids, collapse = ", "))
  })
  output$block_status <- renderText({
    if (is.null(state$selection)) return("先选择有 BED 边界的轨道；没有对应 BED 的轨道不会推断 block。")
    paste0(bplot_locus(state$chr, state$range), " · 黄色阴影为双击选中的 block · 每轨最多 ", opt$max_points, " 个显示点")
  })
  output$block_panel <- renderUI({
    if (is.null(state$selection)) return(div(class = "loading-panel", "等待选择 block…"))
    if (state$loading || !state$ready) return(div(class = "loading-panel", "正在准备 GWAS…"))
    plotly::plotlyOutput("block_plot", height = paste0(max(350L, nrow(opt$tracks) * 175L + 65L), "px"))
  })
  output$block_plot <- plotly::renderPlotly({
    req(state$ready, !state$loading, state$selection)
    views <- bplot_view(opt, state$build, state$chr, state$range, TRUE)
    selected <- if (state$selection$chr == state$chr) state$selection$selected else NULL
    bplot_plot(views, state$build, state$chr, state$range, opt$p_threshold, "blockview", selected)
  })

  # LD I/O runs outside the Shiny event loop. A changed viewport cancels stale work.
  ld <- reactiveValues(result = NULL, loading = FALSE, error = NULL)
  ld_job <- NULL; ld_output <- NULL; ld_log <- NULL; ld_serial <- 0L
  ld_request <- debounce(reactive({
    if (!state$ready || state$loading || state$page != "blockview" || is.null(state$selection)) return(NULL)
    list(build = state$build, chr = state$chr, range = state$range)
  }), 450)
  observeEvent(ld_request(), {
    request <- ld_request()
    if (!is.null(ld_job) && ld_job$is_alive()) ld_job$kill_tree()
    ld$result <- NULL; ld$error <- NULL; ld$loading <- FALSE
    if (is.null(request)) return()
    ld_serial <<- ld_serial + 1L
    ld_output <<- file.path(run_dir, paste0(session$token, "-ld-", ld_serial, ".json"))
    ld_log <<- sub("[.]json$", ".log", ld_output)
    ld$loading <- TRUE
    tryCatch({
      ld_job <<- processx::process$new(opt$ld_python, c(file.path(app_dir, "ld.py"),
        "--root", opt$ld_dir, "--cache", file.path(opt$output_dir, "cache", "ld"),
        "--chr", request$chr, "--grch", request$build,
        "--start", format(request$range[1], scientific = FALSE, trim = TRUE),
        "--end", format(request$range[2], scientific = FALSE, trim = TRUE),
        "--max-snps", as.character(opt$ld_max_snps), "--output", ld_output,
        "--chain-dir", opt$chain_dir, "--liftover-bin", opt$liftover_bin),
        stdout = ld_log, stderr = "2>&1", cleanup_tree = TRUE)
    }, error = function(e) { ld$loading <- FALSE; ld$error <- conditionMessage(e) })
  }, ignoreNULL = FALSE)
  observe({
    if (!ld$loading) return()
    invalidateLater(400, session)
    if (!is.null(ld_job) && !ld_job$is_alive()) {
      ld$loading <- FALSE
      if (identical(ld_job$get_exit_status(), 0L) && file.exists(ld_output)) {
        ld$result <- tryCatch(jsonlite::read_json(ld_output, simplifyVector = FALSE),
          error = function(e) { ld$error <- conditionMessage(e); NULL })
      } else ld$error <- paste0("LD 读取失败；日志：", ld_log)
    }
  })
  session$onSessionEnded(function() if (!is.null(ld_job) && ld_job$is_alive()) ld_job$kill_tree())
  output$ld_status <- renderText({
    if (!is.null(ld$error)) return(ld$error)
    if (ld$loading) return("正在读取各人群 LD；首次建立 SNP 坐标索引需要一些时间…")
    if (is.null(ld$result)) return("选择 block 后展示该区域的 LD。")
    paste0("GRCh", ld$result$build, " · chr", ld$result$chr, " · ", length(ld$result$snps),
      " shown / ", ld$result$total, " reference SNPs · 每个人群共用同一坐标轴（最多 ", opt$ld_max_snps, " SNPs）")
  })
  output$ld_panel <- renderUI({
    if (is.null(ld$result)) return(div(class = "loading-panel", if (ld$loading) "读取 LD…" else "等待区域 LD"))
    tagList(lapply(ld$result$populations, function(item) div(class = "ld-track",
      h5(paste0(item$race, " · ", item$shown, " shown / ", item$available, " available")),
      if (nzchar(item$note)) div(class = "view-status", item$note),
      if (!is.null(item$r2)) plotly::plotlyOutput(paste0("ld_", item$race), height = "310px") else
        div(class = "ld-empty", "该人群在当前区域没有可显示的 LD"))))
  })
  for (race in c("EUR", "AFR", "EAS", "SAS", "AMR")) local({
    pop <- race
    output[[paste0("ld_", pop)]] <- plotly::renderPlotly({ req(ld$result); ld_plot(ld$result, pop) })
  })

}
