filter_callback<-function(generation){
  DT::JS(sprintf("var generation=%d; var publish=function(rows){Shiny.setInputValue('table_filter_state',{generation:generation,rows:rows,initialized:true},{priority:'event'});};
    $(document).off('shiny:inputchanged.gridfilter').on('shiny:inputchanged.gridfilter',function(e){if(e.name==='table_view_rows_all'&&e.value!==null)publish(e.value);});
    table.on('draw.dt.gridfilter',function(){var json=table.ajax.json();if(json&&json.DT_rows_all!==undefined)publish(json.DT_rows_all);});",generation))
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
  output$scope_note<-renderText({req(record());if(isTRUE(as.logical(record()$legacy)))'Historical saved result. Its metric and original evaluation population are unchanged.'else'Current twin run. Selection gain and same-subset training gain are evaluated separately.'})
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
  abm_values<-function(name){
    tr<-if(is.null(input$trait)||input$trait=='all')sort(setdiff(unique(cat$trait),'all'))[1]else input$trait
    if(length(tr)!=1L||is.na(tr))return(data.frame(message='Run the FINAL ABM, then reindex final.'))
    r<-cat[cat$name==name & cat$trait==tr & cat$module=='Twin test' & !as.logical(cat$legacy),,drop=FALSE]
    if(!nrow(r))return(data.frame(message='Run the FINAL ABM, then reindex final.'))
    read_values(con,r[nrow(r),,drop=FALSE])
  }
  abm_curve<-reactive(abm_values('selector_curve'))
  abm_panel<-function(purpose){d<-abm_curve();if(!'purpose'%in%names(d))return(d);d[d$purpose==purpose,,drop=FALSE]}
  abm_plot<-function(purpose){
    d<-abm_panel(purpose);if(!nrow(d)||!'actual_coverage'%in%names(d))return(empty_plot('No ABM selector curve in the current run.'))
    value<-if(purpose=='predictability')d$selector_delta_loss_vs_control else d$candidate_gain_same_subset
    label<-if(purpose=='predictability')'Fallback loss: selected - matched control (lower is better)' else 'Fallback loss - candidate loss (higher is better)'
    z<-data.frame(coverage=d$actual_coverage,value=value)
    ggplot(z,aes(coverage,value))+geom_hline(yintercept=0,linetype=2)+geom_line()+geom_point()+theme_minimal(base_size=12)+labs(x='Actual coverage',y=label)
  }
  output$abm_predictability_plot<-renderPlot(abm_plot('predictability'))
  output$abm_benefit_plot<-renderPlot(abm_plot('benefit'))
  output$abm_predictability_table<-renderDT(datatable(abm_panel('predictability'),rownames=FALSE,options=list(pageLength=6,scrollX=TRUE)),server=FALSE)
  output$abm_benefit_table<-renderDT(datatable(abm_panel('benefit'),rownames=FALSE,options=list(pageLength=6,scrollX=TRUE)),server=FALSE)
  output$abm_cleaning_plot<-renderPlot({d<-abm_values('reference_cleaning');if(!'delta_loss'%in%names(d))return(empty_plot('No reference-cleaning table.'));ggplot(d,aes(delta_loss,contrast))+geom_vline(xintercept=0,linetype=2)+geom_point(size=3)+geom_segment(aes(x=lower,xend=upper,yend=contrast),na.rm=TRUE)+theme_minimal(base_size=12)+labs(x='Paired loss difference (negative favours first model)',y=NULL)})
  output$abm_cleaning_table<-renderDT(datatable(abm_values('reference_cleaning'),rownames=FALSE,options=list(pageLength=8,scrollX=TRUE)),server=FALSE)
  output$abm_status<-renderDT(datatable(abm_values('abm_status'),rownames=FALSE,options=list(dom='t',scrollX=TRUE)),server=FALSE)
}
