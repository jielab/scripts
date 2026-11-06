args<-commandArgs(TRUE);DB_PATH<-args[1];root<-args[2]
for(p in c('shiny','DT','DBI','RSQLite','ggplot2','jsonlite'))suppressPackageStartupMessages(library(p,character.only=TRUE))
source(file.path(root,'data.R'));source(file.path(root,'plots.R'));source(file.path(root,'ui.R'));source(file.path(root,'server.R'))
con<-open_store(DB_PATH);catalogue<-dbReadTable(con,'catalogue')
r<-catalogue[catalogue$module=='abm'&catalogue$name=='benefit_curve',][1,]
stopifnot(!'legacy'%in%names(catalogue),nrow(r)==1)
for(view in c('predictability_curve','benefit_curve','policy_performance','count_incremental_contrasts')){
 d<-abm_read(con,catalogue,r$trait,r$run_id,view)
 z<-abm_plot_data(d,view)
 if(view=='benefit_curve'){
   stopifnot(isTRUE(all.equal(z$value[z$effect=='Gain enrichment'],d[['matched.gain_enrichment.estimate']])))
   stopifnot(isTRUE(all.equal(z$low[z$effect=='Gain enrichment'],d[['matched.gain_enrichment.ci_low']])))
 }
 if(view=='predictability_curve')stopifnot(isTRUE(all.equal(z$value[z$effect=='Selected - control'],d[['matched.predictability_loss_difference.estimate']])))
 if(view=='policy_performance'){
  policy<-d[grepl('^ABM\\(',d$method),]
  stopifnot(isTRUE(all.equal(z$value[z$effect=='Policy gain vs always-base'],policy[['policy_gain_vs_base.estimate']])))
 }
 png(tempfile(fileext='.png'));print(make_abm_plot(d,view));dev.off()
 png(tempfile(fileext='.png'));print(make_abm_plot(d[FALSE,,drop=FALSE],view));dev.off()
 if('matched.gain_enrichment.estimate'%in%names(d)){
   d[['matched.gain_enrichment.estimate']]<-NA_real_;d[['matched.gain_enrichment.ci_low']]<-NA_real_
   z<-abm_plot_data(d,view);stopifnot(all(is.na(z$value[z$effect=='Gain enrichment'])))
 }
}
dbDisconnect(con)
shiny::testServer(server,{
 session$setInputs(abm_trait=r$trait,abm_run=r$run_id,abm_view='benefit_curve')
 session$flushReact();stopifnot(nrow(abm_raw())==6, nrow(abm_filtered())==6)
 session$setInputs(abm_filter_state=list(generation=abm_generation(),initialized=TRUE,rows=c(1L,3L)));stopifnot(nrow(abm_filtered())==2)
 stopifnot(inherits(abm_current_plot(),'ggplot'))
 session$setInputs(abm_run='missing');stopifnot(abm_raw()$status=='unavailable_missing_table')
})
cat('Actual Shiny reactive page and PNG render: PASS\n')
