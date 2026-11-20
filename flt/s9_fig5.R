#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig. 5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(ggplot2)
dir0 <- "D:/"
project <- file.path(dir0, "projects/flt")
out <- file.path(project,"analysis/main_plot")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fuma <- file.path(project,"analysis/03fuma")
csv <- function(p) read.csv(file.path(fuma,p),check.names=FALSE)
tsv <- function(p) read.delim(file.path(fuma,p),check.names=FALSE)
specific <- csv("MAGMA_GTEx54_all.csv"); general <- csv("MAGMA_GTEx_general_all.csv")
stopifnot(nrow(specific)==54,nrow(general)==30)
specific$q <- p.adjust(specific$P,"bonferroni"); general$q <- p.adjust(general$P,"bonferroni")
ev <- tsv("gene_mapping/G732__mdd_2025_Fig6b_gene_evidence.tsv")
ev <- ev %>% filter(Evidence_count >= 4)
sets <- csv("MAGMA_pathways_all.csv")
meta <- head(csv("metascape/Metascape_significant_cluster_representatives.csv"),20)
audit <- readRDS(file.path(project,"analysis/03phewas/cache/plot_audit.rds"))
cloud <- audit[["G732__mdd_2025_trait_n_loci_wordcloud_counts.tsv"]]
if(is.null(cloud)) cloud <- audit[["trait_display_counts.tsv"]]
stopifnot(nrow(cloud)>0, all(is.finite(cloud$n_loci)), all(cloud$n_loci > 0))
# A/B 顶部，C 横跨中间，D/E/F 底部。
boxes <- list(a=c(.005,.665,.755,.995), b=c(.680,.995,.755,.995),
              c=c(.005,.995,.525,.745), d=c(.005,.330,.005,.515),
              e=c(.335,.665,.005,.515), f=c(.670,.995,.005,.515))
ink <- "#243342"
muted <- "#65717E"
cols <- colorRampPalette(c("#77ACCB", "#CF5B67"))(256)
ss <- sets %>% filter(is.finite(P_BH), P_BH < .05) %>% arrange(P)
stopifnot(nrow(ss)>0, nrow(meta)>0, nrow(ev)>0)
score_limit <- max(1, ceiling(max(-log10(pmax(c(ss$P_BH, 10^meta[["Log(q-value)"]]), 1e-300)))))
methods <- names(ev)[3:8]
method_labels <- c("Positional mapping", "MAGMA gene", "eQTLGen cis-eQTL",
                   "GTEx v8 whole blood", "SMR - GTEx whole blood", "SMR - eQTLGen")
method_cols <- c("#298273", "#456A9B", "#6D76B2", "#9465A2", "#B9783F", "#A84E5A")
check_p <- function(d, label) {
  if (any(!is.finite(d$P) | d$P < 0 | d$P > 1))
    stop("Invalid P-values in ", label)
}
check_p(specific,"GTEx54"); check_p(general,"GTEx30")
panel <- function(box, tag, title) {
  par(fig=box, mar=c(0,0,0,0), new=TRUE, xpd=FALSE, cex=1)
  plot.new()
  plot.window(c(0,1), c(0,1), xaxs="i", yaxs="i")
  title_y <- if(tag=="f") .991 else .970
  text(0,title_y,toupper(tag),adj=c(0,1),font=2,cex=1.12,col=ink)
  text(.045,title_y,title,adj=c(0,1),font=2,cex=.92,col=ink)
}
# B：按位点数绘制词云。
draw_cloud <- function() {
  rect(.020,.055,.980,.875,border="#9BA7B3",lwd=.9)
 z <- cloud[order(-cloud$n_loci),]; occupied <- list(); placed <- list()
 for(j in seq_len(nrow(z))) {
  size <- .10+.20*sqrt(z$n_loci[j]/max(z$n_loci))
  word <- z$word[j]; ww <- strwidth(word,cex=size); hh <- strheight(word,cex=size)*1.35
  found <- FALSE
  for(k in 0:18000) {
   ang <- k*.23; radius <- .003*sqrt(k)
   xx <- .5+radius*cos(ang); yy <- .52+radius*sin(ang)
   r <- c(xx-ww/2,xx+ww/2,yy-hh/2,yy+hh/2)
   if(r[1]<0||r[2]>1||r[3]<.07||r[4]>.95) next
   if(any(vapply(occupied,function(b) r[1]<b[2]&&r[2]>b[1]&&r[3]<b[4]&&r[4]>b[3],logical(1)))) next
   placed[[j]] <- list(x=xx,y=yy,word=word,size=size,col=if("color"%in%names(z)) z$color[j] else "#377D70")
   occupied[[length(occupied)+1]] <- r; found <- TRUE; break
  }
  if(!found) stop("Word cloud did not fit: ",word)
 }
 bounds <- do.call(rbind,occupied)
 scale <- min(.82/(max(bounds[,2])-min(bounds[,1])),.66/(max(bounds[,4])-min(bounds[,3])))
 cx <- mean(range(bounds[,1:2])); cy <- mean(range(bounds[,3:4]))
 for(v in placed) text(.5+(v$x-cx)*scale,.52+(v$y-cy)*scale,v$word,cex=v$size*scale,col=v$col)
 text(.5,.015,"Size / colour: distinct loci; records >= 20",cex=.65)
}
# D/E：通路气泡图。
pathways <- function(box,tag,title,labels,x,score,count,xlabel,size_label,notes) {
  panel(box,tag,title)
  text(.01,.905,notes,adj=0,cex=.59,col=muted)
  labels <- vapply(labels,function(z)
    paste(strwrap(z,width=27),collapse="\n"),character(1))
  label_cex <- .67
  while (max(strwidth(labels,cex=label_cex)) > .405 && label_cex > .42)
    label_cex <- label_cex-.02
  yy <- seq(.845,.175,length.out=length(labels))
  if (length(yy)>1 && max(strheight(labels,cex=label_cex)) > .90*min(abs(diff(yy))))
    label_cex <- label_cex*.90*min(abs(diff(yy)))/max(strheight(labels,cex=label_cex))
  xmax <- max(pretty(c(0,x),n=3))
  if (!is.finite(xmax) || xmax<=0) xmax <- 1
  ticks <- pretty(c(0,xmax),n=3)
  ticks <- ticks[ticks>=0 & ticks<=xmax]
  x0 <- .455
  width <- .265
  for (tick in ticks)
    segments(x0+width*tick/xmax,.155,x0+width*tick/xmax,.865, col="#EEF0F2",lwd=.6)
  for (j in seq_along(labels)) {
    text(.01,yy[j],labels[j],adj=0,cex=label_cex,col=ink)
    segments(x0,yy[j],x0+width,yy[j],col="#EDF0F3",lwd=.6)
    ci <- 1+floor(255*min(max(score[j]/score_limit,0),1))
    points(x0+width*x[j]/xmax,yy[j],pch=21,bg=cols[ci],col="white",lwd=.4, cex=2.5*sqrt(count[j]/max(count)))
  }
  axis(1,at=x0+width*ticks/xmax,labels=ticks,pos=.115,
       cex.axis=.65,col=ink,col.axis=ink,tck=-.012,mgp=c(0,.45,0))
  text(x0+width/2,.065,xlabel,cex=.67,col=ink)
  # 右侧色条和点大小图例。
  text(.805,.855,"-log10(BH q)",adj=0,cex=.59,col=ink)
  rect(.810,.615+(0:255)*.205/256,.850,.615+(1:256)*.205/256,col=cols,border=NA)
  text(.865,c(.615,.820),c(0,score_limit),adj=0,cex=.59,col=muted)
  text(.805,.540,paste(strwrap(size_label,width=12),collapse="\n"), adj=c(0,1),cex=.59,col=ink)
  count_ticks <- unique(round(as.numeric(quantile(count,c(.25,.65,1)))))
  legend_y <- seq(.435,.255,length.out=length(count_ticks))
  points(rep(.835,length(count_ticks)),legend_y,pch=21,bg="#AAB8C2",col="white",
         cex=2.5*sqrt(count_ticks/max(count)))
  text(.900,legend_y,count_ticks,adj=0,cex=.59,col=muted)
}
# F：先在全部组织内校正，再筛选显著组织。
magma_tissue_selected <- rbind(
  transform(specific[specific$q < .05,c("name","P","q")],family="GTEx specific"),
  transform(general[general$q < .05,c("name","P","q")],family="General"))
draw_magma_tissues <- function() {
  tables <- list(general[general$q<.05,],specific[specific$q<.05,])
  titles <- c("General tissues","GTEx specific tissues")
  n <- sum(vapply(tables,nrow,integer(1)))
  if(!n) {
    text(.5,.5,"No tissues passed Bonferroni P < 0.05",cex=.8)
    return(invisible(NULL))
  }
  scores <- -log10(pmax(magma_tissue_selected$q,1e-300))
  xmax <- ceiling(max(scores))
  ticks <- pretty(c(0,xmax),n=4); ticks <- ticks[ticks>=0 & ticks<=xmax]
  x0 <- .53; width <- .45
  dy <- min(.040,.64/n)
  top <- .865
  plot_top <- top
  plot_bottom <- plot_top-n*dy-.09*(length(tables)-1)
  for(t in ticks) segments(x0+width*t/xmax,plot_bottom, x0+width*t/xmax,plot_top,col="#EEF0F2",lwd=.5)
  text(.005,.94,"Bonferroni P < 0.05; significant tissues only",adj=0,cex=.65,col=muted)
  for(i in seq_along(tables)) {
    d <- tables[[i]][order(tables[[i]]$q,tables[[i]]$P),]
    text(.005,top+.035,paste0(titles[i]," (",nrow(d),")"), adj=0,font=2,cex=.75,col=ink)
    if(!nrow(d)) { top <- top-.07; next }
    yy <- top-(seq_len(nrow(d))-.5)*dy
    bottom <- top-nrow(d)*dy
    labels <- vapply(gsub("_"," ",d$name),function(z)
      paste(strwrap(z,width=28),collapse="\n"),character(1))
    cex_label <- .70
    while(max(strwidth(labels,cex=cex_label))>.505 ||
          max(strheight(labels,cex=cex_label))>dy*.85) cex_label <- cex_label-.01
    text(x0-.018,yy,labels,adj=1,cex=cex_label,col=ink)
    score <- -log10(pmax(d$q,1e-300))
    rect(x0,yy-dy*.44,x0+width*score/xmax,yy+dy*.44, col="#75A5C0",border=NA)
    if(i<length(tables)) segments(x0,bottom-.045,x0+width,bottom-.045, col="#CBD2D8",lty=2,lwd=.6)
    top <- bottom-.09
  }
  rect(x0,plot_bottom,x0+width,plot_top,border="#90979E",lwd=.6)
  axis(1,at=x0+width*ticks/xmax,labels=ticks,pos=plot_bottom,
       cex.axis=.65,col.axis=ink,col="#90979E",tck=-.008,mgp=c(0,.4,0))
  text(.755,.055,"-log10(Bonferroni P)",cex=.72,col=ink)
  text(.005,.018,"Correction across all 54 / 30 tissues, separately",adj=0,cex=.62,col=muted)
}
# A：读取 PLACO 坐标与 P 值；FUMA 仅提供高亮 SNP 名单。
placo <- file.path(project,"analysis/03placo")
snp <- data.table::fread(file.path(placo,"G732__mdd_2025_PLACO_all.tsv.gz"),
                        select=c("SNP","CHR1","POS1","P_PLACO")) %>%
  as_tibble() %>% rename(CHR=CHR1,BP=POS1,P=P_PLACO) %>%
  mutate(CHR=as.integer(sub("^chr","",CHR))) %>%
  filter(CHR %in% 1:22,is.finite(BP),BP>0,is.finite(P),between(P,0,1))
stopifnot(!anyDuplicated(snp$SNP),nrow(snp)>0)
lead <- read.delim(file.path(placo,"FUMA_job787285/leadSNPs.txt"))$rsID
ind <- read.delim(file.path(placo,"FUMA_job787285/IndSigSNPs.txt"))$rsID
chr_axis <- snp %>% group_by(CHR) %>% summarise(length=max(BP),.groups="drop") %>%
  arrange(CHR) %>% mutate(offset=lag(cumsum(length+sum(length)*.003),default=0),center=offset+length/2)
snp <- snp %>% left_join(chr_axis,by="CHR") %>%
  mutate(x=BP+offset,logp=-log10(pmax(P,.Machine$double.xmin)),
    type=case_when(SNP %in% lead ~ "Lead SNPs",SNP %in% ind ~ "IndSig SNPs",.default="Background"))
# 仅对非显著背景点按显示分辨率抽稀，保留全部显著及高亮点。
highlights <- snp %>% filter(type != "Background")
background <- snp %>% filter(type == "Background") %>%
  mutate(bin_x=floor(x/(max(snp$x)/1800)),bin_y=floor(logp/.04)) %>%
  group_by(CHR,bin_x,bin_y) %>% filter(P<=1e-5 | P==min(P)) %>% ungroup()
chr_cols <- rep(c("#4594BD","#E5B63C","#51468D","#626E33","#BD287A", "#D86B32","#C581B5","#79BCA1","#782B44","#307B68"),length.out=22)
manhattan <- ggplot(background,aes(x,logp)) +
  geom_point(aes(color=factor(CHR)),size=.12,alpha=.65) +
  scale_color_manual(values=chr_cols,guide="none") +
  geom_hline(yintercept=c(5,-log10(5e-8)),linetype=c(3,2),color="grey50",linewidth=.2) +
  geom_point(data=highlights,aes(fill=type),shape=21,size=.7,stroke=.1) +
  scale_fill_manual(values=c("Lead SNPs"="#781B31","IndSig SNPs"="#397660"),name=NULL) +
  ggrepel::geom_text_repel(data=filter(highlights,type=="Lead SNPs"),aes(label=SNP),
    seed=732,size=1.3,max.overlaps=Inf,max.iter=30000,max.time=5,
    min.segment.length=0,segment.size=.1,box.padding=.12,ylim=c(6,12)) +
  scale_x_continuous(breaks=chr_axis$center,labels=chr_axis$CHR,expand=expansion(mult=.01)) +
  coord_cartesian(ylim=c(0,max(12,highlights$logp))) +
  labs(title="A  PLACO Manhattan plot",x="Chromosome",y=expression(-log[10](P))) + theme_classic(base_size=8) +
  theme(plot.title=element_text(size=8,face="bold"),axis.text=element_text(size=5),legend.position="inside",legend.position.inside=c(.85,.9),
    legend.text=element_text(size=5),legend.key.size=grid::unit(2,"mm"),
    plot.margin=margin(8,2,2,2))
draw <- function() {
  par(mar=rep(0,4),family="sans",bg="white",cex=1)
  plot.new()
  panel(boxes$a,"a","PLACO Manhattan plot")
  v <- boxes$a
  print(manhattan,newpage=FALSE,vp=grid::viewport(x=mean(v[1:2]),y=mean(v[3:4]),
    width=diff(v[1:2]),height=diff(v[3:4])))
  panel(boxes$b,"b","PheWAS")
  draw_cloud()
  # C：沿用 s5_2_gene_plot.R 的 Fig6b 格式：交替底色、完整方法名、大小点及竖排基因名。
  panel(boxes$c,"c","Genes meeting at least 4 of six criteria")
  text(.60,.89,"Large dot: meets criterion    Small grey dot: does not meet criterion",
       cex=1,col=muted)
  edges <- c(.255,.985)
  step <- diff(edges)/nrow(ev)
  xx <- edges[1]+(seq_len(nrow(ev))-.5)*step
  yy <- seq(.79,.35,length.out=6)
  for (j in seq_along(methods)) {
    if(j %% 2L == 1L) rect(edges[1],yy[j]-.043,edges[2],yy[j]+.043,col="#F1F4F8",border=NA)
    present <- as.logical(ev[[methods[j]]])
    points(xx,rep(yy[j],length(xx)),pch=16,cex=ifelse(present,1.1,.38),
           col=ifelse(present,method_cols[j],"#BEC9D6"))
  }
  text(edges[1]-.012,yy,method_labels,adj=1,cex=1,col=ink)
  text(xx,.29,ev$Gene,srt=90,adj=1,cex=1,col=ink)
  text(.985,.025,"MAGMA BH-FDR < 0.05; SMR: in Sheet8 for the tissue = significant; absent = not significant.",
       adj=1,cex=.8,col=muted)
  pathway_names <- gsub("_"," ",sub("^(GOCC_|GOBP_|GOMF_|GO_)","",ss$name))
  pathway_names <- paste0(toupper(substr(pathway_names,1,1)),tolower(substring(pathway_names,2)))
  pathways(boxes$d,"d","MAGMA gene-set enrichment",pathway_names,
           -log10(pmax(ss$P_BH,1e-300)),-log10(pmax(ss$P_BH,1e-300)),ss$NGENES,
           "-log10(BH q)","Genes in set",paste0("BH-FDR < 0.05; ",format(nrow(sets),big.mark=",")," tested gene sets"))
  pathways(boxes$e,"e","Metascape functional enrichment",meta$Description,
           meta$Enrichment,-meta[["Log(q-value)"]],meta[["#GeneInGOAndHitList"]],
           "Enrichment factor","Hit genes","Native BH q < 0.05; up to 20 cluster representatives")
  panel(boxes$f,"f","MAGMA tissue expression analysis")
  draw_magma_tissues()
}
# 仅保存合图 PDF；宽 210 mm，基础字号 8 pt。
pdf_path <- file.path(out,"Fig5.pdf")
cairo_pdf(pdf_path,width=210/25.4,height=240/25.4,pointsize=8)
tryCatch(draw(),finally=dev.off())
message("Saved: ",pdf_path)


