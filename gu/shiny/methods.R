# Shared method pages: a single click inspects a record without navigation.
gu_method_ui <- function(id) {
  ns <- NS(id)
  tagList(
    layout_columns(
      selectInput(ns("chr"),"Chromosome",c("ALL",chr_order),"ALL"),
      numericInput(ns("start"),"Start (0-based)",0,min=0),
      numericInput(ns("end"),"End",250000000,min=1),
      selectInput(ns("population"),"Super_pop",c("ALL","AFR","EAS","EUR","SAS","AMR")),
      col_widths=c(3,3,3,3)),
    tags$p(class="text-muted","按当前方法与区间展示最多 2,000 条片段；单击记录更新本页明细。坐标为 0-based half-open。"),
    gu_card(card_header("Regional segment landscape"),plotlyOutput(ns("plot"),height="400px")),
    gu_card(card_header("Segment calls"),DTOutput(ns("calls"),fill=FALSE)),
    gu_card(card_header("Selected record"),DTOutput(ns("selected"),fill=FALSE)),
    accordion(open=FALSE,
      accordion_panel("Per-sample burden",DTOutput(ns("burden"),fill=FALSE)),
      accordion_panel("Published callset overlay",DTOutput(ns("reference"),fill=FALSE)),
      accordion_panel("Exports",downloadButton(ns("download"),"Download filtered segments")))
  )
}
gu_method_server <- function(id,Q,dataset,build) {
  force(id)
  moduleServer(id,function(input,output,session) {
    calls <- reactive({
      req(input$chr,input$population,is.finite(input$start),is.finite(input$end))
      validate(need(input$end>input$start,"End 必须大于 Start。"))
      sql <- "SELECT s.sample_id,p.population,p.super_population,s.method,s.source_class,s.chr,s.start,s.end,s.length_bp,s.haplotype,s.score,s.posterior,s.locus_id FROM segments s LEFT JOIN sample_populations p ON p.dataset_id=s.dataset_id AND p.sample_id=s.sample_id WHERE s.dataset_id=? AND s.genome_build=? AND s.method=? AND s.end>? AND s.start<?"
      params <- list(dataset(),build(),id,input$start,input$end)
      if(input$chr!="ALL") { sql<-paste0(sql," AND s.chr=?");params<-c(params,input$chr) }
      if(input$population!="ALL") { sql<-paste0(sql," AND p.super_population=?");params<-c(params,input$population) }
      Q(paste0(sql," ORDER BY s.source_class,s.chr,s.start LIMIT 2000"),params)
    })
    output$calls <- renderDT({
      d<-calls()
      if(!nrow(d)) return(datatable(data.frame(状态="当前方法、区间或人群暂无结果"),rownames=FALSE,options=list(dom="t")))
      datatable(d,rownames=FALSE,selection=list(mode="single",selected=1),options=list(pageLength=12,scrollX=TRUE))
    })
    selected <- reactive({
      d<-calls();i<-input$calls_rows_selected
      if(!length(i) || i<1 || i>nrow(d)) i<-1L
      if(nrow(d)) d[i,,drop=FALSE] else d
    })
    output$selected <- renderDT({
      d<-selected()
      datatable(d,rownames=FALSE,class="compact stripe nowrap",selection="none",
        options=list(dom="t",scrollX=TRUE,ordering=FALSE,paging=FALSE))
    })
    output$plot <- renderPlotly({
      d<-calls();validate(need(nrow(d)>0,"当前筛选下暂无片段结果。"))
      d$row<-seq_len(nrow(d))
      d$label<-paste0(d$sample_id," · ",d$super_population," · chr",d$chr,":",d$start,"–",d$end)
      p<-plot_ly() %>% add_segments(data=d,x=~start,xend=~end,y=~row,yend=~row,color=~chr,text=~label,hoverinfo="text",line=list(width=3))
      layout(p,xaxis=list(title="Position (bp)"),yaxis=list(title="Segment calls",showticklabels=FALSE)) %>% config(displaylogo=FALSE)
    })
    output$burden <- renderDT({
      sql<-"SELECT b.sample_id,p.population,p.super_population,b.source_class,b.n_input_segments,b.n_merged_intervals,b.total_bp FROM sample_burden b LEFT JOIN sample_populations p ON p.dataset_id=b.dataset_id AND p.sample_id=b.sample_id WHERE b.dataset_id=? AND b.genome_build=? AND b.method=? AND b.burden_type='nonredundant_union'"
      params<-list(dataset(),build(),id)
      if(input$population!="ALL") {sql<-paste0(sql," AND p.super_population=?");params<-c(params,input$population)}
      datatable(Q(paste0(sql," ORDER BY b.total_bp DESC LIMIT 10000"),params),rownames=FALSE,options=list(pageLength=12,scrollX=TRUE))
    })
    output$reference <- renderDT({
      req(input$chr,input$end>input$start)
      sql<-"SELECT dataset_id,population,source_class,reference_role,chr,start,end FROM reference_callsets WHERE genome_build=? AND end>? AND start<?"
      params<-list(build(),input$start,input$end)
      if(input$chr!="ALL") {sql<-paste0(sql," AND chr=?");params<-c(params,input$chr)}
      datatable(Q(paste0(sql," ORDER BY chr,start LIMIT 2000"),params),rownames=FALSE,options=list(pageLength=12,scrollX=TRUE))
    })
    output$download <- downloadHandler(filename=function()paste0("gu-",id,"-",dataset(),"-",build(),".tsv"),content=function(file)data.table::fwrite(calls(),file,sep="\t"))
  })
}

# Overview panels share an explicit focus region; clicking a bin selects, never navigates.
gu_matching_ui <- function(id,label) {
  ns<-NS(id)
  tagList(uiOutput(ns("context")),plotlyOutput(ns("curve"),height="260px"))
}
gu_matching_server <- function(id,method,Q,dataset,build,region,on_select,on_open) {
  force(id);force(method)
  moduleServer(id,function(input,output,session) {
    available<-reactive({
      Q("SELECT 1 FROM segments WHERE dataset_id=? AND genome_build=? AND method=? LIMIT 1",list(dataset(),build(),method))
    })
    empty_message<-reactive({
      if(!nrow(available()))return("该方法暂无分析结果。")
      r<-region();req(!is.null(r))
      runs<-Q("SELECT 1 FROM method_runs WHERE dataset_id=? AND genome_build=? AND method=? AND chr=? AND status='complete' LIMIT 1",list(dataset(),build(),method,as.character(r$chr[[1]])))
      if(!nrow(runs)) "当前染色体尚无已完成的分析结果。" else "该方法在当前区间没有检出的片段。"
    })
    regional_calls <- reactive({
      if(!nrow(available()))return(data.frame())
      r<-region();req(!is.null(r),nrow(r)>0)
      Q("SELECT s.sample_id,p.super_population,s.start,s.end FROM segments s INDEXED BY idx_segments_region LEFT JOIN sample_populations p ON p.dataset_id=s.dataset_id AND p.sample_id=s.sample_id WHERE s.dataset_id=? AND s.genome_build=? AND s.chr=? AND s.start<? AND s.end>? AND s.method=?",
        list(dataset(),build(),as.character(r$chr[[1]]),r$end[[1]],r$start[[1]],method))
    })
    bins <- reactive({
      r<-region();req(!is.null(r))
      d<-regional_calls();if(!nrow(d)) return(data.frame())
      edges<-unique(round(seq(r$start[[1]],r$end[[1]],length.out=61)))
      populations<-c("AFR","EAS","EUR","SAS","AMR")
      out<-lapply(seq_len(length(edges)-1L),function(i) {
        x<-d[d$start<edges[i+1L] & d$end>edges[i],,drop=FALSE]
        data.frame(chr=as.character(r$chr[[1]]),start=edges[i],end=edges[i+1L],super_population=populations,
          carriers=vapply(populations,function(pop)length(unique(x$sample_id[!is.na(x$super_population)&x$super_population==pop])),integer(1)))
      })
      do.call(rbind,out)
    })
    output$context<-renderUI({
      r<-region();req(!is.null(r))
      tags$p(class="text-muted",paste0("chr",r$chr[[1]],":",format(r$start[[1]]+1,big.mark=","),"–",format(r$end[[1]],big.mark=","),
        " · 每个区间的非重复携带者人数；按 Super_pop 分类。单击联动 IGV 与 PhyML；双击查看方法明细。"))
    })
    output$curve<-renderPlotly({
      d<-bins()
      if(!nrow(d))return(plotly_empty() %>% layout(xaxis=list(visible=FALSE),yaxis=list(visible=FALSE),annotations=list(list(text=empty_message(),x=.5,y=.5,xref="paper",yref="paper",showarrow=FALSE))) %>% config(displayModeBar=FALSE))
      d$key<-seq_len(nrow(d));d$mid<-(d$start+d$end)/2
      p<-plot_ly(d,x=~mid,y=~carriers,color=~super_population,colors=c("#C47C25","#3585BB","#8556A5","#23977D","#D76478"),
        customdata=~key,type="scatter",mode="lines+markers",marker=list(size=4),
        text=~paste0(super_population," · chr",chr,":",start+1,"–",end),hovertemplate="%{text}<br>Carriers: %{y}<extra></extra>") %>%
        layout(xaxis=list(title="Position (bp)"),yaxis=list(title="Carriers",rangemode="tozero"),margin=list(l=60,r=20,t=10,b=45)) %>% config(displaylogo=FALSE,doubleClick=FALSE)
      htmlwidgets::onRender(p,sprintf("function(el){el.on('plotly_click',function(e){var k=e.points[0].customdata;if(k==null)return;Shiny.setInputValue('%s',k,{priority:'event'});var t=Date.now();if(el._last&&el._last.k===k&&t-el._last.t<450){Shiny.setInputValue('%s',k,{priority:'event'});el._last=null;}else{el._last={k:k,t:t};}});el.on('plotly_relayout',function(e){var a=e['xaxis.range']||[e['xaxis.range[0]'],e['xaxis.range[1]']];if(a[0]!=null&&a[1]!=null)Shiny.setInputValue('%s',{start:a[0],end:a[1]},{priority:'event'});});}",session$ns("select_bin"),session$ns("open_bin"),session$ns("view_range")))
    })
    observeEvent(input$view_range, {
      r<-region();req(!is.null(r))
      lo<-as.numeric(input$view_range$start);hi<-as.numeric(input$view_range$end)
      req(is.finite(lo),is.finite(hi),hi>lo)
      on_select(data.frame(chr=r$chr[[1]],start=max(0,floor(lo)),end=ceiling(hi),super_population="ALL"),method)
    })
    # Save the clicked row before shared context redraws the table.
    clicked<-reactiveVal(NULL)
    get_row<-function(i) {d<-isolate(bins());i<-as.integer(i);if(length(i)==1L&&!is.na(i)&&i>=1&&i<=nrow(d)) d[i,,drop=FALSE] else NULL}
    observeEvent(input$select_bin, {r<-get_row(input$select_bin);req(!is.null(r));clicked(r);on_select(r,method)})
    observeEvent(input$open_bin, {r<-isolate(clicked());if(is.null(r))r<-get_row(input$open_bin);req(!is.null(r));on_open(r,method)})
  })
}
