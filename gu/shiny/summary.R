# Static population map. All statistics use the cached exact Neanderthal union,
# independently of the Overview viewport, with equal weight per tested person.
gu_ibdmix_summary <- function(samples, burden) {
  d <- merge(samples, burden, by="sample_id", all.x=TRUE, sort=FALSE)
  good <- is.finite(d$coverage_pct) & d$tested_bp>0 & d$n_chromosomes>0
  measured <- d[good,,drop=FALSE]
  pops <- lapply(split(d, d$population), function(p) {
    x <- p[is.finite(p$coverage_pct) & p$tested_bp>0,,drop=FALSE]
    data.frame(population=p$population[[1]],super_population=p$super_population[[1]],
      n=nrow(x),n_total=nrow(p),pct=if(nrow(x))mean(x$coverage_pct) else NA_real_,
      mb=if(nrow(x))mean(x$neanderthal_bp)/1e6 else NA_real_)
  })
  list(samples=d, measured=measured, populations=if(length(pops))do.call(rbind,pops) else data.frame(),
       complete=nrow(d)>0 && all(good) && all(d$n_chromosomes==22),
       comparable_scope=nrow(measured)>0 && length(unique(measured$chromosomes))==1L,
       overall_pct=if(nrow(measured))mean(measured$coverage_pct) else NA_real_,
       overall_mb=if(nrow(measured))mean(measured$neanderthal_bp)/1e6 else NA_real_)
}

gu_ibdmix_map <- function(summary, land, locations) {
  old <- par(mar=c(1,0,2,0),family="sans",bg="white");on.exit(par(old))
  plot.new();plot.window(xlim=c(-168,181),ylim=c(-63,93),asp=1,xaxs="i",yaxs="i")
  for(d in split(land,land$polygon)) polygon(d$longitude,d$latitude,col="#f4f5f5",border="#bbc3c8",lwd=.65)
  title(main="Neanderthal sequence by population",adj=0,cex.main=1.05,col.main="#253746",line=.3)
  p <- merge(locations,summary$populations,by=c("population","super_population"),sort=FALSE)
  if(!nrow(p)) {text(0,0,"No mapped 1KG populations in this dataset",col="#607080");return(invisible(NULL))}
  segments(p$longitude,p$latitude,p$label_x,p$label_y,col="#aab3ba",lwd=.7)
  points(p$longitude,p$latitude,pch=20,col="#7f8a92",cex=.45)
  for(i in seq_len(nrow(p))) {
    r<-p[i,];radius<-4.5;theta<-seq(0,2*pi,length.out=100)
    polygon(r$label_x+radius*cos(theta),r$label_y+radius*sin(theta),col=if(is.finite(r$pct))"#e0e3e5" else "white",border="#89959e",lwd=.8)
    if(is.finite(r$pct) && r$pct>0) {
      theta<-seq(pi/2,pi/2-2*pi*min(r$pct,100)/100,length.out=40)
      polygon(c(r$label_x,r$label_x+radius*cos(theta)),c(r$label_y,r$label_y+radius*sin(theta)),col="#b3222b",border=NA)
    }
    label<-paste0(r$population,"\n",if(is.finite(r$pct))sprintf("%.2f%%",r$pct) else "N/A")
    text(r$label_x,r$label_y-6,label,adj=c(.5,1),cex=.78,col="#253746",font=2)
  }
  legend("bottomleft",inset=c(.015,0),legend=c("Neanderthal / diploid span","Remaining span","Not evaluated"),
    fill=c("#b3222b","#e0e3e5","white"),border="#89959e",bty="n",horiz=TRUE,cex=.75,text.col="#52616d")
}

gu_ibdmix_summary_text <- function(s, manifest) {
  pct <- function(x)if(is.finite(x))sprintf("%.2f%%",x) else "N/A"
  mb <- function(x)if(is.finite(x))sprintf("%.1f",x) else "—"
  d<-s$measured
  groups<-c("AFR (5)","EUR","EAS","SAS","AMR")
  rows<-lapply(groups,function(group) {
    p<-if(group=="AFR (5)") d[d$population %in% c("ESN","GWD","LWK","MSL","YRI"),,drop=FALSE] else d[d$super_population %in% group,,drop=FALSE]
    reference<-switch(group,`AFR (5)`="~17",EUR="~51",EAS="~55",SAS="~55","—")
    shiny::tags$tr(shiny::tags$td(group),shiny::tags$td(nrow(p)),
      shiny::tags$td(pct(if(nrow(p))mean(p$coverage_pct) else NA_real_)),
      shiny::tags$td(mb(if(nrow(p))mean(p$neanderthal_bp)/1e6 else NA_real_)),shiny::tags$td(reference))
  })
  scope<-if(s$complete)"22 / 22 autosomes · X excluded" else if(nrow(d))paste0("Partial coverage · ",paste(sort(unique(d$n_chromosomes)),collapse=", ")," / 22 autosomes per person") else "No certified whole-autosome results"
  shiny::tagList(
    shiny::tags$div(class="gu-summary-kicker","NEANDERTHAL / DIPLOID GENOME"),
    shiny::tags$div(class="gu-summary-value",pct(s$overall_pct)),
    shiny::tags$p(class="gu-summary-subvalue",paste0(mb(s$overall_mb)," Mb / person")),
    shiny::tags$p(paste0(format(nrow(d),big.mark=",")," / ",format(nrow(s$samples),big.mark=",")," individuals · ",sum(s$populations$n>0)," populations")),
    shiny::tags$p(class="gu-summary-scope",scope),
    if(length(manifest$diploid_autosome_bp))shiny::tags$p(class="text-muted",sprintf("完整常染色体二倍体分母：%.3f Gb（2 × 参考序列长度）。",manifest$diploid_autosome_bp/1e9)),
    if(!s$complete)shiny::tags$p(class="text-muted","仅统计已确认完成的整条常染色体；未分析的染色体、人群或缺少目标样本记录的运行记为 N/A，不外推全基因组比例。"),
    if(!s$comparable_scope && nrow(d))shiny::tags$p(class="text-muted","个体的已分析染色体范围不同，人群比例不能直接比较。"),
    shiny::tags$h6("Population summary"),
    shiny::tags$table(class="table table-sm gu-summary-table",
      shiny::tags$thead(shiny::tags$tr(lapply(c("Group","n","Coverage","Mb / ind.","Cell 2020"),shiny::tags$th))),shiny::tags$tbody(rows))
  )
}
