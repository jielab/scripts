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
  observeEvent(values(),{
    n<-numeric_cols(values());updateSelectInput(session,'x',choices=n,selected=if(length(n))n[1] else character())
    updateSelectInput(session,'y',choices=n,selected=if(length(n)>1)n[2] else n)
    updateSelectInput(session,'group',choices=c('',names(values())),selected='')
  })
  current_plot<-reactive(make_plot(values(),record()$plot_kind,paste(record()$name,input$feature,sep=' · '),input$x,input$y,input$group))
  output$source_note<-renderText({req(record());paste('Run:',record()$run_id,'| Source:',record()$source,'| SHA:',substr(record()$source_sha256,1,12),'|',record()$note)})
  output$overview_table<-renderDT({d<-scoped();d<-d[,intersect(c('trait','layer','module','name','rows','status','run_id','table_id'),names(d)),drop=FALSE];datatable(d,rownames=FALSE,filter='top',selection='single',escape=TRUE,options=list(pageLength=15,scrollX=TRUE))},server=TRUE)
  observeEvent(input$overview_table_rows_selected,{
    i<-input$overview_table_rows_selected;if(length(i)&&i<=nrow(scoped())){updateSelectizeInput(session,'table',selected=scoped()$table_id[i]);updateNavbarPage(session,'nav',selected='explore')}
  })
  output$audit_table<-renderDT(datatable(DBI::dbReadTable(con,'audit'),rownames=FALSE,filter='top',escape=TRUE,options=list(pageLength=15,scrollX=TRUE)),server=TRUE)
  for(prefix in c('explore','twin','deviation'))local({
    id<-prefix
    output[[paste0(id,'_title')]]<-renderText({req(record());paste(record()$name,'—',input$feature)})
    output[[paste0(id,'_plot')]]<-renderPlotly({ggplotly(current_plot(),tooltip='all')})
    output[[paste0(id,'_note')]]<-renderText({d<-values();paste('显示',nrow(d),'行；源表',record()$rows,'行。',if(isTRUE(attr(d,'truncated')))'当前选择超过100000行，只显示/下载前100000行；请进一步筛选。'else '')})
    output[[paste0(id,'_table')]]<-renderDT({d<-values();datatable(d,rownames=FALSE,filter='top',selection='single',escape=TRUE,options=list(pageLength=12,scrollX=TRUE,deferRender=TRUE))},server=TRUE)
    observeEvent(input[[paste0(id,'_table_rows_selected')]],{
      i<-input[[paste0(id,'_table_rows_selected')]];d<-values();col<-record()$feature_column
      if(length(i)&&nzchar(col)&&i<=nrow(d))updateSelectizeInput(session,'feature',selected=as.character(d[[col]][i]))
    })
    output[[paste0(id,'_data')]]<-downloadHandler(filename=function()paste0('LE8_',gsub('[^A-Za-z0-9_.-]','_',input$feature),'.tsv'),content=function(file){write.table(values(),file,sep='\t',row.names=FALSE,quote=TRUE,na='NA')})
    output[[paste0(id,'_png')]]<-downloadHandler(filename=function()paste0('LE8_',gsub('[^A-Za-z0-9_.-]','_',input$feature),'.png'),content=function(file){ggsave(file,current_plot(),width=11,height=7,dpi=200)})
  })
}
