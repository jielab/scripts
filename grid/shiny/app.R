#!/usr/bin/env Rscript
script <- grep("^--file=",commandArgs(),value=TRUE)
app_dir <- Sys.getenv("GRID_SHINY_APP_DIR",if(length(script)) dirname(normalizePath(sub("^--file=","",script[1]))) else "/mnt/d/scripts/grid/shiny")
source(file.path(app_dir,"data.R"))
suppressPackageStartupMessages({library(shiny);library(bslib);library(DT);library(plotly);library(ggplot2)})

grid_app <- function(final_dir) {
  initial <- grid_read_snapshot(final_dir)
  labels <- c(pca="PCA",csx="CSx",disco="Disco",prsformer="PRSformer",abm="ABM")
  descriptions <- c(pca="参考 PCA 空间、群体分布与祖源分配质量。",
    csx="各祖源及组合评分的分布、覆盖情况与来源。评分大小不等同于预测准确性。",
    disco="基于遗传距离插值的评分分布，以及与 CSx 评分的相关性。",
    prsformer="训练过程与该模型自身测试集表现。跨方法比较请查看综合比较页。",
    abm="ABM 的训练审计、匹配支持和预测置信度。置信度结论属于内部测试集验证。")
  card <- function(...) div(class="grid-card",...)
  display_table <- function(d, n=10L) {
    d <- as.data.frame(d)
    widget <- datatable(d,rownames=FALSE,options=list(pageLength=n,scrollX=TRUE))
    numeric <- names(d)[vapply(d,is.numeric,logical(1))]
    continuous <- numeric[vapply(d[numeric],function(x) any(is.finite(x) & x!=round(x)),logical(1))]
    if(length(continuous)) widget <- formatSignif(widget,continuous,digits=4)
    widget
  }
  module_ui <- function(m) tabPanel(labels[[m]],value=m,
    div(class="page-heading",h2(labels[[m]]),p(descriptions[[m]])),
    fluidRow(column(3,card(selectInput(paste0(m,"_trait"),"性状",c("全部 / 通用"="all",initial$traits)),
      selectInput(paste0(m,"_table"),"结果表",character()),selectInput(paste0(m,"_figure"),"结果图",character()),
      downloadButton(paste0(m,"_download"),"下载当前表"))),
      column(9,card(uiOutput(paste0(m,"_status")),DTOutput(paste0(m,"_data"))),
        card(imageOutput(paste0(m,"_image"),height="auto")))) )
  ui <- navbarPage(title=div(class="brand",span("GRID"),tags$small("RESULTS ATLAS")),id="page",
    theme=bs_theme(version=5, bg="#f3f6f5",fg="#203b41",primary="#13796d",base_font="system-ui"),
    header=tagList(tags$head(tags$style(HTML('
      .navbar{background:#fff!important;border-bottom:1px solid #dce7e2;padding:14px 24px}
      .brand{display:flex;gap:12px;align-items:baseline}.brand span{font-weight:800;letter-spacing:2px}.brand small{font-size:10px;letter-spacing:2px;color:#68817b}
      .container-fluid{max-width:1600px}.page-heading{padding:28px 8px 10px}.page-heading h2{font-weight:700;margin-bottom:10px}.page-heading p{color:#5b7470}
      .grid-card{background:white;border:1px solid #deE7e3;border-radius:14px;padding:24px;margin-bottom:20px;box-shadow:0 3px 12px #193b3210}
      .kpi strong{display:block;font-size:32px;color:#13796d}.kpi span{color:#617a75}.snapshot-line{display:flex;justify-content:space-between;align-items:center;padding:18px 8px 0;color:#617a75;font-size:13px}
      .notice{border-left:4px solid #d69835;background:#fff8ea;padding:12px 16px;margin:12px 0}.finding{padding:10px 0;border-bottom:1px solid #edf0ed;line-height:1.7}
      .shiny-image-output img{max-width:100%;height:auto}.grid-card .dataTables_wrapper{overflow-x:auto}.btn{border-radius:8px}
    '))),div(class="snapshot-line",textOutput("snapshot_time",inline=TRUE),actionButton("refresh","刷新汇总"))),
    tabPanel("总览",value="overview",
      div(class="page-heading",h2("从群体结构到预测表现"),p("已有分析结果汇总，以及在共同测试队列上的方法比较。")),
      fluidRow(column(3,card(div(class="kpi",textOutput("module_count",container=strong),span("分析模块")))),
        column(3,card(div(class="kpi",textOutput("trait_count",container=strong),span("性状")))),
        column(3,card(div(class="kpi",textOutput("available_count",container=strong),span("已有结果的模块")))),
        column(3,card(downloadButton("summary_download","下载综合汇总"),p("状态、结论与性能表")))),
      fluidRow(column(7,card(h4("主要发现与待完成事项"),selectInput("finding_trait","性状",c("全部"="all",initial$traits),selected=initial$traits[1]),uiOutput("findings"))),
        column(5,card(h4("模块状态"),DTOutput("inventory"))))),
    tabPanel("综合比较",value="comparison",
      div(class="page-heading",h2("同一队列上的方法比较"),p("保留 final 的指标定义、样本量与区间；按性状和祖源分别解读。")),
      fluidRow(column(3,card(selectInput("comparison_trait","性状",initial$traits),
        selectInput("comparison_targets","祖源",character(),multiple=TRUE),
        selectInput("comparison_methods","方法",character(),multiple=TRUE),
        selectInput("comparison_detail","审计与补充结果",c("配对差值"="paired_comparison","距离分层"="distance","测试队列"="test_cohort","方法与来源"="metadata")),
        downloadButton("comparison_download","下载当前比较"))),
        column(9,card(uiOutput("comparison_notice"),plotlyOutput("comparison_plot",height="520px"),DTOutput("comparison_data")),
          card(h4("配对比较、距离分层与队列审计"),DTOutput("comparison_audit"))))),
    module_ui("pca"),module_ui("csx"),module_ui("disco"),
    if(any(initial$inventory$module=="prsformer" & initial$inventory$status=="available")) module_ui("prsformer"),module_ui("abm")
  )
  server <- function(input,output,session) {
    snapshot <- reactiveVal(initial)
    prs_visible <- reactiveVal(any(initial$inventory$module=="prsformer" & initial$inventory$status=="available"))
    observe({
      available <- any(snapshot()$inventory$module=="prsformer" & snapshot()$inventory$status=="available")
      if(available != isolate(prs_visible())) {
        if(available) insertTab("page",module_ui("prsformer"),target="abm",position="before")
        else removeTab("page","prsformer")
        prs_visible(available)
      }
    })
    observeEvent(input$refresh,{
      tryCatch({snapshot(grid_read_snapshot(final_dir));showNotification("已读取最新 final 汇总",type="message")},
        error=function(e) showNotification(conditionMessage(e),type="error",duration=NULL))
    })
    output$snapshot_time <- renderText(paste("数据快照",snapshot()$generated_at,"· 来源：final"))
    output$trait_count <- renderText(length(snapshot()$traits))
    output$module_count <- renderText(length(unique(snapshot()$inventory$module)))
    output$available_count <- renderText(length(unique(snapshot()$inventory$module[snapshot()$inventory$status=="available"])))
    output$summary_download <- downloadHandler(filename=function() "grid-summary.xlsx",content=function(file)
      file.copy(grid_snapshot_path(final_dir,snapshot()$workbook),file))
    output$inventory <- renderDT(datatable(snapshot()$inventory,rownames=FALSE,options=list(pageLength=15,dom="t",scrollX=TRUE)))
    output$findings <- renderUI({
      d <- snapshot()$findings
      if(!is.null(input$finding_trait) && input$finding_trait!="all") d <- d[d$trait %in% c("all",input$finding_trait),]
      tagList(lapply(seq_len(nrow(d)),function(i) div(class="finding",strong(paste0(d$trait[i]," · ")),d$text[i])))
    })
    observe({
      d <- grid_performance(snapshot(),input$comparison_trait)
      updateSelectInput(session,"comparison_targets",choices=unique(d$target),selected=unique(d$target))
      updateSelectInput(session,"comparison_methods",choices=unique(d$method),selected=unique(d$method))
    })
    compared <- reactive(grid_performance(snapshot(),input$comparison_trait,input$comparison_targets,input$comparison_methods))
    output$comparison_notice <- renderUI({
      d <- compared()
      if(!nrow(d)) return(div(class="notice","尚无所选条件的最终比较结果。"))
      if(any(d$report_status=="stale")) return(div(class="notice","上游结果已更新。这里显示的是此前的比较，请重新运行 final 后刷新。"))
      p(paste("指标：",paste(unique(d$metric),collapse=" / "),"。误差线表示保存的 95% 区间；区间缺失时仅显示点估计。"))
    })
    output$comparison_plot <- renderPlotly({
      d <- compared(); validate(need(nrow(d)>0,"暂无可比较的结果"))
      d$hover <- paste(d$method,d$target,paste0("N = ",d$N),paste0(d$metric," = ",signif(d$estimate,5)),sep="\n")
      p <- ggplot(d,aes(estimate,method,colour=method,text=hover)) + geom_point(size=3) +
        geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation="y",width=.15,na.rm=TRUE) +
        facet_wrap(~target,scales="free_x") + theme_minimal(base_size=12) +
        theme(legend.position="none",panel.grid.minor=element_blank()) + labs(y=NULL,x=paste(unique(d$metric),collapse=" / "))
      ggplotly(p,tooltip="text")
    })
    output$comparison_data <- renderDT(display_table(compared()))
    output$comparison_audit <- renderDT({
      z <- snapshot()$tables[[paste("final",input$comparison_trait,input$comparison_detail,sep="/")]]
      if(is.null(z)) z <- data.frame(message="暂无此项汇总")
      display_table(z,8L)
    })
    output$comparison_download <- downloadHandler(filename=function() paste0(input$comparison_trait,"-comparison.csv"),
      content=function(file) write.csv(compared(),file,row.names=FALSE,na=""))
    for (module in names(labels)) local({
      m <- module
      selected_assets <- reactive({
        z <- snapshot()$artifacts
        z[z$module==m & z$trait %in% c("all",input[[paste0(m,"_trait")]]),]
      })
      observe({
        trait <- input[[paste0(m,"_trait")]]
        keys <- names(snapshot()$tables)
        keys <- keys[startsWith(keys,paste0(m,"/"))]
        if(!is.null(trait) && trait!="all") keys <- keys[vapply(strsplit(keys,"/",fixed=TRUE),function(x) x[2] %in% c("all",trait),logical(1))]
        updateSelectInput(session,paste0(m,"_table"),choices=setNames(keys,sub(paste0("^",m,"/"),"",keys)))
        z <- selected_assets()
        if(is.null(trait) || trait=="all") z <- snapshot()$artifacts[snapshot()$artifacts$module==m,]
        updateSelectInput(session,paste0(m,"_figure"),choices=setNames(z$path,paste(z$trait,z$title,sep=" / ")))
      })
      output[[paste0(m,"_status")]] <- renderUI({
        z <- snapshot()$inventory[snapshot()$inventory$module==m,]
        if(all(z$status=="missing")) div(class="notice","该模块尚无可用结果。完成运行后执行 final，再刷新页面。")
      })
      current_table <- reactive({
        key <- input[[paste0(m,"_table")]]
        if(is.null(key)||!nzchar(key)) return(data.frame(message="暂无此项汇总"))
        snapshot()$tables[[key]]
      })
      output[[paste0(m,"_data")]] <- renderDT(display_table(current_table(),12L))
      output[[paste0(m,"_image")]] <- renderImage({
        p <- input[[paste0(m,"_figure")]]
        req(!is.null(p),nzchar(p),p %in% snapshot()$artifacts$path)
        list(src=grid_snapshot_path(final_dir,p),contentType="image/png",alt=paste(labels[[m]],"结果图"))
      },deleteFile=FALSE)
      output[[paste0(m,"_download")]] <- downloadHandler(filename=function() paste0(m,"-summary.csv"),
        content=function(file) write.csv(current_table(),file,row.names=FALSE,na=""))
    })
  }
  app <- shinyApp(ui,server)
  app$grid_ui <- ui
  app
}

if(sys.nframe()==0L) {
  app <- grid_app(Sys.getenv("GRID_FINAL_DIR","/mnt/d/analysis/grid/final"))
  if("--check" %in% commandArgs(TRUE)) {
    invisible(htmltools::renderTags(app$grid_ui))
    cat("SHINY CHECK OK: final snapshot, artifacts and UI loaded\n")
  } else runApp(app,host=Sys.getenv("GRID_SHINY_HOST","127.0.0.1"),port=as.integer(Sys.getenv("GRID_SHINY_PORT","3839")),launch.browser=FALSE)
}
