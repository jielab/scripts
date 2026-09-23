# Static Altai map and separate per-reference population burdens use cached exact
# unions, independently of the Overview viewport, with equal weight per tested person.
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
  title(main="Altai Neanderthal sequence by population",adj=0,cex.main=1.05,col.main="#253746",line=.3)
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
  scope<-if(s$complete)"22 / 22 autosomes · X excluded" else if(nrow(d))paste0("Partial coverage · ",paste(sort(unique(d$n_chromosomes)),collapse=", ")," / 22 autosomes per person") else "No certified whole-autosome results"
  filters<-manifest$ibdmix_filters
  daf_complete<-is.data.frame(filters) && nrow(filters)>0 &&
    any(filters$chrom %in% as.character(1:22)) &&
    all(filters$daf_status[filters$chrom %in% as.character(1:22)]=="applied")
  shiny::tagList(
    shiny::tags$div(class="gu-summary-kicker","ALTAI NEANDERTHAL / DIPLOID GENOME"),
    shiny::tags$div(class="gu-summary-value",pct(s$overall_pct)),
    shiny::tags$p(class="gu-summary-subvalue",paste0(mb(s$overall_mb)," Mb / person")),
    shiny::tags$p(paste0(format(nrow(d),big.mark=",")," / ",format(nrow(s$samples),big.mark=",")," individuals · ",sum(s$populations$n>0)," populations")),
    shiny::tags$p(class="gu-summary-scope",scope),
    shiny::tags$p(class="text-muted","Cell 2020 对照仅使用 Altai 和常染色体，不合并其他尼安德特人参考。"),
    shiny::tags$p(class="text-muted",if(daf_complete)"已应用 Altai 衍生等位基因比例最高 0.1% 窗口过滤。" else "尚未应用 Altai 衍生等位基因比例最高 0.1% 窗口过滤，不能视为论文严格过滤结果。"),
    if(length(manifest$diploid_autosome_bp))shiny::tags$p(class="text-muted",sprintf("完整常染色体二倍体分母：%.3f Gb（2 × 参考序列长度）。",manifest$diploid_autosome_bp/1e9)),
    if(!s$complete)shiny::tags$p(class="text-muted","仅统计已确认完成的整条常染色体；未分析的染色体、人群或缺少目标样本记录的运行记为 N/A，不外推全基因组比例。"),
    if(!s$comparable_scope && nrow(d))shiny::tags$p(class="text-muted","个体的已分析染色体范围不同，人群比例不能直接比较。")
  )
}

gu_ibdmix_summary_table <- function(s, burden) {
  refs<-c("Altai","Chagyr","Vindija","Denisova","Denisova25","Altai + Denisova","All Five")
  d<-merge(s$samples[,c("sample_id","population","super_population")],burden,by="sample_id",sort=FALSE)
  d<-d[is.finite(d$archaic_bp) & d$tested_bp>0 & d$n_chromosomes>0,,drop=FALSE]
  group_rows<-function(d,group) {
    if(group=="AFR (5)") d[d$population %in% c("ESN","GWD","LWK","MSL","YRI"),,drop=FALSE] else d[d$super_population %in% group,,drop=FALSE]
  }
  rows<-lapply(c("AFR (5)","EUR","EAS","SAS","AMR"),function(group) {
    p<-group_rows(d,group)
    altai<-p[p$reference=="Altai",,drop=FALSE]
    cells<-lapply(refs,function(ref) {
      x<-p[p$reference==ref,,drop=FALSE]
      if(!nrow(x))return(shiny::tags$td("N/A"))
      scope<-paste(sort(unique(x$n_chromosomes)),collapse=",")
      partial<-any(x$n_chromosomes!=22)
      detail<-paste0("n = ",nrow(x)," · ",scope," / 22 autosomes; chromosomes: ",paste(unique(x$chromosomes),collapse="; "))
      shiny::tags$td(title=detail,sprintf("%.1f",mean(x$archaic_bp)/1e6),
        if(partial || nrow(x)!=nrow(altai))shiny::tags$small(class="gu-summary-cell-scope",
          paste0("n = ",nrow(x),if(partial)paste0(" · partial ",scope,"/22 chr"))))
    })
    shiny::tags$tr(shiny::tags$td(group),shiny::tags$td(nrow(altai)),
      shiny::tags$td(if(nrow(altai))sprintf("%.2f%%",mean(altai$coverage_pct)) else "N/A"),cells,
      shiny::tags$td(switch(group,`AFR (5)`="~17",EUR="~51",EAS="~55",SAS="~55","—")))
  })
  headers<-c(list(shiny::tags$th("Group"),shiny::tags$th("n (Altai)"),shiny::tags$th("Coverage (Altai)")),
    lapply(refs,function(ref)shiny::tags$th(ref,shiny::tags$br(),"Mb / ind.")),
    list(shiny::tags$th("Cell 2020 (Altai)",shiny::tags$br(),"Mb / ind.")))
  shiny::tagList(
    shiny::tags$h6("Population summary · 5 archaic references + 2 unions"),
    shiny::tags$div(class="gu-summary-table-scroll",shiny::tags$table(class="table table-sm gu-summary-table",
      shiny::tags$thead(shiny::tags$tr(headers)),shiny::tags$tbody(rows))),
    shiny::tags$p(class="gu-summary-caption",
      "前五列分别统计各参考的片段并集长度；Altai + Denisova 仅合并这两个参考（不含 Denisova25），All Five 合并全部五参考。两种合并均先按个体取跨参考片段并集，重叠位置只计一次，再按已检测个体等权平均（Mb / ind.），不能直接相加前五列。排除 X，合并列仅统计所有成员参考均已确认完成的常染色体。"),
    shiny::tags$p(class="gu-summary-caption",
      "地图、Coverage 和右侧比例均仅基于 Altai。n 为 Altai 已检测人数；其他列人数不同或染色体不完整时在单元格标注，N/A 表示缺少该参考或组合的共同检测结果。鼠标悬停可查看各列人数与染色体范围。"),
    shiny::tags$p(class="gu-summary-caption",
      "Altai、Chagyr、Vindija 为 Neanderthal 参考；Denisova、Denisova25 为 Denisovan 参考。Denisovan 列为原生参考匹配片段，不能直接等同于已确认的渗入来源。AFR (5) 仅含 ESN、GWD、LWK、MSL、YRI。Cell 2020 为 Altai 文献近似值。")
  )
}
