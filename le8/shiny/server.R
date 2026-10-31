# Stamp DT's rows_all events with the source generation. Late events from a
# replaced table cannot filter a new table/feature with the old row numbers.
filter_callback <- function(id,generation) {
  DT::JS(sprintf("var key = %s, generation = %d;
    var publish = function(rows) {
      Shiny.setInputValue(key + '_filter_state',
        {generation: generation, rows: rows, initialized: true}, {priority: 'event'});
    };
    $(document).off('shiny:inputchanged.le8_' + key)
      .on('shiny:inputchanged.le8_' + key, function(e) {
        if (e.name === key + '_rows_all' && e.value !== null) publish(e.value);
      });
    // A redraw may have identical indices, so DT may suppress rows_all's event.
    table.on('draw.dt.le8', function() {
      var json = table.ajax.json();
      if (json && json.DT_rows_all !== undefined) publish(json.DT_rows_all);
    });",jsonlite::toJSON(id,auto_unbox=TRUE),generation))
}

server <- function(input,output,session) {
  con<-open_store();session$onSessionEnded(function()DBI::dbDisconnect(con))
  catalogue<-read_catalogue(con)
  output$metadata<-renderPrint(DBI::dbReadTable(con,'metadata'))
  updateSelectInput(session,'trait',choices=c('all',sort(setdiff(unique(catalogue$trait),'all'))))
  modules<-c('all',sort(setdiff(unique(catalogue$module),'all')))
  updateSelectInput(session,'module',choices=modules)
  observeEvent(input$nav,{
    target<-switch(input$nav,twin='ABM twin',deviation='Molecular deviation',NULL)
    if(!is.null(target)&&target%in%modules)updateSelectInput(session,'module',selected=target)
  })
  scoped<-reactive({
    d<-catalogue
    for(pair in list(c('trait',input$trait),c('layer',input$layer),c('module',input$module)))
      if(length(pair)==2&&!is.null(pair[2])&&pair[2]!='all')d<-d[d[[pair[1]]]%in%c(pair[2],'all'),,drop=FALSE]
    d
  })
  observeEvent(scoped(),{
    d<-scoped();labels<-paste(d$name,'·',d$run_id,'·',d$status,'·',dirname(d$source))
    updateSelectizeInput(session,'table',choices=setNames(d$table_id,labels),server=TRUE,
       selected=if(!is.null(input$table)&&input$table%in%d$table_id)input$table else if(nrow(d))d$table_id[1] else character())
  })
  record<-reactive({req(input$table);d<-scoped();d<-d[d$table_id==input$table,,drop=FALSE];validate(need(nrow(d)==1,'Choose an available table.'));d})
  observeEvent(record(),{
    f<-features_for(con,record()$table_id)
    updateSelectizeInput(session,'feature',choices=c('ALL',f),selected=if(!is.null(input$feature)&&input$feature%in%f)input$feature else 'ALL',server=TRUE)
  })
  values<-reactive({req(record(),input$feature);query_table(con,record(),input$feature,trait=input$trait,layer=input$layer)})
  generation<-reactiveVal(0L)
  observeEvent(values(),{generation(isolate(generation())+1L)},priority=100)
  filtered_values<-function(id) reactive({
    d<-values();state<-input[[paste0(id,'_table_filter_state')]]
    # NULL before initialization is different from an initialized empty result.
    if(is.null(state)||!identical(as.integer(state$generation),generation())||!isTRUE(state$initialized))return(d)
    rows<-as.integer(unlist(state$rows,use.names=FALSE))
    rows<-unique(rows[!is.na(rows)&rows>=1L&rows<=nrow(d)])
    d[rows,,drop=FALSE]
  })
  panel_values<-setNames(lapply(c('explore','twin','deviation'),filtered_values),c('explore','twin','deviation'))
  observeEvent(values(),{
    n<-numeric_cols(values());updateSelectInput(session,'x',choices=n,selected=if(length(n))n[1] else character())
    updateSelectInput(session,'y',choices=n,selected=if(length(n)>1)n[2] else n)
    updateSelectInput(session,'group',choices=c('',names(values())),selected='')
  })
  panel_plot<-function(id)make_plot(panel_values[[id]](),record()$plot_kind,paste(record()$name,input$feature,sep=' · '),input$x,input$y,input$group)
  current_plot<-reactive({id<-input$nav;if(is.null(id)||!id%in%names(panel_values))id<-'explore';panel_plot(id)})
  write_filtered<-function(id,file)write.table(panel_values[[id]](),file,sep='\t',row.names=FALSE,quote=TRUE,na='NA')
  output$source_note<-renderText({req(record());paste('Run:',record()$run_id,'| Source:',record()$source,'| SHA:',substr(record()$source_sha256,1,12),'|',record()$note)})
  output$overview_table<-renderDT({d<-scoped();d<-d[,intersect(c('trait','layer','module','name','rows','status','run_id','table_id'),names(d)),drop=FALSE];datatable(d,rownames=FALSE,filter='top',selection='single',escape=TRUE,options=list(pageLength=15,scrollX=TRUE))},server=TRUE)
  observeEvent(input$overview_table_rows_selected,{
    i<-input$overview_table_rows_selected;if(length(i)&&i<=nrow(scoped())){updateSelectizeInput(session,'table',selected=scoped()$table_id[i]);updateNavbarPage(session,'nav',selected='explore')}
  })
  output$audit_table<-renderDT(datatable(DBI::dbReadTable(con,'audit'),rownames=FALSE,filter='top',escape=TRUE,options=list(pageLength=15,scrollX=TRUE)),server=TRUE)
  for(prefix in c('explore','twin','deviation'))local({
    id<-prefix
    output[[paste0(id,'_title')]]<-renderText({req(record());paste(record()$name,'—',input$feature)})
    output[[paste0(id,'_plot')]]<-renderPlotly({ggplotly(panel_plot(id),tooltip='all')})
    output[[paste0(id,'_note')]]<-renderText({d<-values();paste('筛选后',nrow(panel_values[[id]]()),'行（全部分页）；当前分子范围',nrow(d),'行；源表',record()$rows,'行。图与 TSV 使用相同筛选行。',if(isTRUE(attr(d,'truncated')))'当前选择超过100000行，只显示/下载前100000行；请进一步筛选。'else '')})
    output[[paste0(id,'_table')]]<-renderDT({d<-values();datatable(d,rownames=FALSE,filter='top',selection='single',escape=TRUE,callback=filter_callback(paste0(id,'_table'),generation()),options=list(pageLength=12,scrollX=TRUE,deferRender=TRUE,stateSave=FALSE))},server=TRUE)
    observeEvent(input[[paste0(id,'_table_rows_selected')]],{
      i<-input[[paste0(id,'_table_rows_selected')]];d<-values();col<-record()$feature_column
      if(length(i)&&nzchar(col)&&i<=nrow(d))updateSelectizeInput(session,'feature',selected=as.character(d[[col]][i]))
    })
    output[[paste0(id,'_data')]]<-downloadHandler(filename=function()paste0('LE8_filtered_',gsub('[^A-Za-z0-9_.-]','_',input$feature),'.tsv'),content=function(file){write_filtered(id,file)})
    output[[paste0(id,'_png')]]<-downloadHandler(filename=function()paste0('LE8_',gsub('[^A-Za-z0-9_.-]','_',input$feature),'.png'),content=function(file){ggsave(file,panel_plot(id),width=11,height=7,dpi=200)})
  })
}
