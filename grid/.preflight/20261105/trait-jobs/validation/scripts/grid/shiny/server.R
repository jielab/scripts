filter_callback<-function(generation,input_id="table_filter_state",table_id="table_view"){
  DT::JS(sprintf("var generation=%d; var channel='%s'; var publish=function(rows){Shiny.setInputValue('%s',{generation:generation,rows:rows,initialized:true},{priority:'event'});};
    $(document).off('shiny:inputchanged.gridfilter_'+channel).on('shiny:inputchanged.gridfilter_'+channel,function(e){if(e.name==='%s_rows_all'&&e.value!==null)publish(e.value);});
    table.on('draw.dt.gridfilter',function(){var json=table.ajax.json();if(json&&json.DT_rows_all!==undefined)publish(json.DT_rows_all);});",generation,table_id,input_id,table_id))
}
server<-function(input,output,session){
  con<-open_store(DB_PATH);session$onSessionEnded(function()DBI::dbDisconnect(con))
  cat<-DBI::dbReadTable(con,'catalogue')
  updateSelectInput(session,'trait',choices=c('all',sort(setdiff(unique(cat$trait),'all'))))
  updateSelectInput(session,'module',choices=c('all',sort(unique(cat$module))))
  scoped<-reactive({d<-cat;if(!is.null(input$trait)&&input$trait!='all')d<-d[d$trait%in%c('all',input$trait),,drop=FALSE];if(!is.null(input$module)&&input$module!='all')d<-d[d$module==input$module,,drop=FALSE];d})
  observeEvent(scoped(),{
    d<-scoped();selected<-if(!is.null(input$table)&&input$table%in%d$table_id)input$table else if(nrow(d))d$table_id[1]else character()
    updateSelectizeInput(session,'table',choices=setNames(d$table_id,paste(d$name,d$trait,d$module,substr(d$run_id,1,12),sep=' · ')),selected=selected,server=TRUE)
  })
  record<-reactive({req(input$table);d<-scoped();d<-d[d$table_id==input$table,,drop=FALSE];shiny::validate(shiny::need(nrow(d)==1,'No available table in this scope.'));d})
  raw<-reactive({req(record());read_values(con,record())})
  choose<-function(values,preferred){values<-sort(unique(as.character(values[!is.na(values)])));if(preferred%in%values)preferred else 'ALL'}
  observeEvent(raw(),{
    d<-raw();ac<-first_col(d,c('ancestry','target'));sc<-first_col(d,c('subset'));mc<-first_col(d,c('method','model'))
    av<-if(!is.null(ac))d[[ac]]else character();sv<-if(!is.null(sc))d[[sc]]else character();mv<-if(!is.null(mc))sort(unique(as.character(d[[mc]])))else character()
    updateSelectInput(session,'ancestry',choices=unique(c('ALL',sort(unique(av)))),selected=choose(av,'ALL'))
    updateSelectInput(session,'subset',choices=unique(c('ALL',sort(unique(sv)))),selected=choose(sv,'all'))
    met<-if('metric'%in%names(d))sort(unique(as.character(d$metric)))else character()
    pref<-if('MSE'%in%met)'MSE'else if('Brier'%in%met)'Brier'else if(length(met))met[1]else'ALL'
    updateSelectInput(session,'source_metric',choices=unique(c('ALL',met)),selected=pref)
    updateSelectizeInput(session,'methods',choices=c('All methods'='__ALL__',setNames(mv,mv)),selected='__ALL__',server=TRUE)
    updateSelectInput(session,'metric',choices=metric_columns(d),selected=choose_metric(d))
  })
  source_values<-reactive({
    d<-raw();ac<-first_col(d,c('ancestry','target'));sc<-first_col(d,c('subset'));mc<-first_col(d,c('method','model'))
    if(!is.null(input$trait)&&input$trait!='all'&&'trait'%in%names(d))d<-d[d$trait%in%c(input$trait,'all'),,drop=FALSE]
    if(!is.null(ac)&&!is.null(input$ancestry)&&input$ancestry!='ALL')d<-d[!is.na(d[[ac]])&d[[ac]]==input$ancestry,,drop=FALSE]
    # ALL is a stored stratum as well as a UI sentinel: avoid pooling both when
    # an aggregate ALL stratum exists. Individual ancestries remain selectable.
    else if(!is.null(ac)&&any(d[[ac]]=='ALL',na.rm=TRUE))d<-d[!is.na(d[[ac]])&d[[ac]]=='ALL',,drop=FALSE]
    if(!is.null(sc)&&!is.null(input$subset)&&input$subset!='ALL')d<-d[!is.na(d[[sc]])&d[[sc]]==input$subset,,drop=FALSE]
    if('metric'%in%names(d)&&!is.null(input$source_metric)&&input$source_metric!='ALL')d<-d[!is.na(d$metric)&d$metric==input$source_metric,,drop=FALSE]
    if(!is.null(mc)&&!'__ALL__'%in%input$methods)d<-d[d[[mc]]%in%input$methods,,drop=FALSE]
    d
  })
  generation<-reactiveVal(0L)
  observeEvent(source_values(),{generation(isolate(generation())+1L)},priority=100)
  filtered<-reactive(filtered_rows(source_values(),input$table_filter_state,generation()))
  current_plot<-reactive({req(record());d<-filtered();metric<-input$metric;if(is.null(metric)||!nzchar(metric))metric<-choose_metric(d);make_scoped_plot(d,record()$plot_kind,metric,record()$name)})
  output$title<-renderText({req(record());paste(record()$trait,record()$name,sep=' · ')})
  output$scope_note<-renderText({req(record());'Frozen run identity and source TSV define the interpretation. Missing estimates and intervals remain unavailable.'})
  output$source<-renderText({req(record());paste('Run:',record()$run_id,'|',record()$source,'| SHA256:',record()$source_sha256,'|',record()$note)})
  output$plot<-renderPlot(current_plot())
  output$table_view<-renderDT({d<-source_values();datatable(d,rownames=FALSE,filter='top',selection='none',escape=TRUE,callback=filter_callback(generation()),options=list(pageLength=12,scrollX=TRUE,deferRender=TRUE,stateSave=FALSE))},server=TRUE)
  output$row_note<-renderText({paste(nrow(filtered()),'filtered rows across all pages;',nrow(source_values()),'rows in current scope.',if(isTRUE(attr(raw(),'truncated')))'Source exceeds 100,000 rows: current load is truncated; use a smaller aggregate table.'else'Figure and download use the same filtered rows.')})
  write_filtered<-function(file)write.table(filtered(),file,sep='\t',row.names=FALSE,quote=TRUE,na='NA')
  output$download_data<-downloadHandler(filename=function()'GRID_filtered.tsv',content=write_filtered)
  output$download_plot<-downloadHandler(filename=function()'GRID_current.png',content=function(file)ggsave(file,current_plot(),width=11,height=7,dpi=180))
  output$catalogue<-renderDT(datatable(scoped(),rownames=FALSE,filter='top',selection='single',escape=TRUE,options=list(pageLength=15,scrollX=TRUE)),server=TRUE)
  observeEvent(input$catalogue_rows_selected,{i<-input$catalogue_rows_selected;if(length(i)&&i>=1&&i<=nrow(scoped())){updateSelectizeInput(session,'table',selected=scoped()$table_id[i]);updateNavbarPage(session,'nav',selected='results')}})
  output$audit_table<-renderDT(datatable(DBI::dbReadTable(con,'audit'),rownames=FALSE,filter='top',escape=TRUE,options=list(pageLength=15,scrollX=TRUE)),server=TRUE)
  output$metadata<-renderPrint(DBI::dbReadTable(con,'metadata'))
  abm_catalogue<-cat[cat$module=='abm',,drop=FALSE]
  updateSelectInput(session,'abm_trait',choices=sort(unique(abm_catalogue$trait)))
  observeEvent(input$abm_trait,{
    runs<-sort(unique(abm_catalogue$run_id[abm_catalogue$trait==input$abm_trait]))
    updateSelectInput(session,'abm_run',choices=runs)
  })
  abm_raw<-reactive({req(input$abm_view);abm_read(con,cat,input$abm_trait,input$abm_run,input$abm_view)})
  abm_generation<-reactiveVal(0L)
  observeEvent(abm_raw(),{abm_generation(isolate(abm_generation())+1L)},priority=100)
  abm_filtered<-reactive(filtered_rows(abm_raw(),input$abm_filter_state,abm_generation()))
  abm_current_plot<-reactive(make_abm_plot(abm_filtered(),input$abm_view,paste(input$abm_trait,input$abm_run)))
  output$abm_table<-renderDT(datatable(abm_raw(),rownames=FALSE,filter='top',escape=TRUE,
    callback=filter_callback(abm_generation(),'abm_filter_state','abm_table'),
    options=list(pageLength=12,scrollX=TRUE)),server=TRUE)
  output$abm_plot<-renderPlot(abm_current_plot())
  output$abm_download_data<-downloadHandler(filename=function()paste0(input$abm_view,'.tsv'),
    content=function(file)write.table(abm_filtered(),file,sep='\t',row.names=FALSE,quote=TRUE,na='NA'))
  output$abm_download_plot<-downloadHandler(filename=function()paste0(input$abm_view,'.png'),
    content=function(file)ggsave(file,abm_current_plot(),width=12,height=7,dpi=180))
  output$abm_source<-renderText({r<-abm_records(cat,input$abm_trait,input$abm_run);r<-r[r$name==input$abm_view,,drop=FALSE];
    if(!nrow(r))'Unavailable: no completed table for this run.' else paste(r$source,r$source_sha256,collapse=' | ')})
}
