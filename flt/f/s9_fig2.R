source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig. 2: MR Manhattan plot from the current UVMR and MVMR runs.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(dplyr)

project <- Sys.getenv("FLT_OUT")
flt_file <- Sys.getenv("FLT_DICTIONARY")

output_dir <- file.path(Sys.getenv("FLT_PAPER"))
mr_method <- "IVW"
alpha <- 0.05

# 1. 读取 IVW 结果，按 GWAS ID 匹配食物名称和分类。
d <- flt_uv("forward") %>% filter(Method == mr_method) %>% mutate(FDR=suppressWarnings(as.numeric(`P-FDR`)))
meta <- readxl::read_excel(flt_file, sheet="item182") %>% select(gwas_id, name, item_raw, category)
# 分类配色与图 1 一致，各结局共用横坐标。
groups <- unique(meta$category)
palette <- c("#C54C3B", "#4BA7C1", "#80634F", "#415A82", "#DE897C", "#818DB6", "#80BFB3", "#C52C36", "#309B89", "#AA998C")
cols <- setNames(palette[seq_along(groups)], groups)
meta <- meta %>% filter(gwas_id %in% d$Exposure) %>% mutate(type=ifelse(startsWith(name, "F-"), "Derived", "Raw")) %>% arrange(match(category, groups), type == "Derived", tolower(item_raw)) %>% mutate(x=row_number())
outcomes <- c("All-cause dementia", "Stroke", "Sleep disorders", "Depression", "Anxiety disorders")
d <- d %>%
  mutate(label=meta$item_raw[match(Exposure, meta$gwas_id)],
         category=factor(meta$category[match(Exposure, meta$gwas_id)], levels=groups),
         x=meta$x[match(Exposure, meta$gwas_id)], disease=factor(Outcome, levels=outcomes),
         significant=!is.na(FDR) & FDR < alpha, IVW_P=suppressWarnings(as.numeric(`P-value`)),
         P_zero=!is.na(IVW_P) & IVW_P == 0,
         neglog10_P=-log10(pmax(IVW_P, .Machine$double.xmin)),
         neglog10_FDR=-log10(pmax(FDR, .Machine$double.xmin))) %>%
  arrange(disease, x) %>% mutate(pair_id=if_else(significant, cumsum(significant), NA_integer_))
groups <- groups[groups %in% meta$category]
g <- meta %>% group_by(category) %>% summarise(mid=mean(range(x)), .groups="drop") %>% arrange(match(category, groups))
counts <- d %>% group_by(disease, .drop=FALSE) %>% summarise(n=n(), n_FDR=sum(!is.na(FDR)), significant=sum(significant), .groups="drop")
# 2. Manhattan 图：按分类着色，标注全部 FDR 显著条目。
advanced_pairs <- data.frame(Exposure=c("G717","G824","G732"), Outcome=c("Sleep disorders","Sleep disorders","Depression"))
pair_key <- function(x) paste(x$Exposure,x$Outcome,sep="|")
d$advanced <- pair_key(d) %in% pair_key(advanced_pairs)
# 各结局独立设置纵轴范围。
if(any(d$FDR==0,na.rm=TRUE)) warning("Zero FDR floored for plotting; inspect input precision")
p <- ggplot(d,aes(x,neglog10_FDR,color=category)) +
  geom_hline(yintercept=0,color="#D6DCE1",linewidth=.3) +
  geom_hline(yintercept=-log10(alpha),linetype="dashed",color="#708898",linewidth=.35) +
  geom_point(data=d[!d$significant, ],size=1.15,alpha=.8,na.rm=TRUE) +
  geom_point(data=d[d$significant, ],size=2.8,na.rm=TRUE) +
  ggrepel::geom_label_repel(data=d[d$significant, ],aes(label=label),
    seed=42,size=3.1,fill="white",label.size=0,
    label.padding=grid::unit(.10,"lines"),box.padding=.35,point.padding=.2,
    min.segment.length=0,segment.size=.3,max.overlaps=Inf,
    max.iter=50000,max.time=5,show.legend=FALSE) +
  facet_wrap(~disease,ncol=1,scales="free_y") +
  scale_color_manual(values=cols,drop=FALSE) +
  scale_x_continuous(breaks=g$mid,labels=g$category, limits=c(-3,nrow(meta)+4),expand=expansion(mult=0)) +
  scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(.025,.18))) +
  labs(x=NULL,y=expression(-log[10](P[FDR])), caption="Labels: FDR < 0.05; colours: food categories.\nY ranges vary by outcome.") +
  theme_classic(base_size=13) +
  theme(legend.position="none",panel.background=element_rect(fill="white",color=NA),
        panel.spacing.y=grid::unit(.15,"lines"),
        strip.background=element_blank(),strip.placement="outside",
        strip.text.x=element_text(face="bold",size=12,hjust=0),
        axis.text.x=element_text(angle=35,hjust=1,size=11),
        axis.text.y=element_text(color="#444444"),axis.ticks.x=element_blank(),
        axis.line=element_line(color="#CCD2D7",linewidth=.3),
        plot.title=element_text(face="bold",size=14),
        plot.caption=element_text(size=10.5,hjust=0),plot.margin=margin(8,4,8,6))
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
# 3. 森林图：按 Exposure + Outcome 匹配 Sheet7 的 MV-IVW。
# 保留全部 IVW 显著条目，置信区间取自工作簿中的舍入值。
forest <- d %>% filter(significant)
if (!nrow(forest)) {
  empty <- ggplot() + annotate("text", x = 0, y = 0, label = "No forward IVW results met BH FDR < 0.05") +
    xlim(-1, 1) + ylim(-1, 1) + theme_void()
  ggsave(file.path(output_dir, "Fig2.png"), patchwork::wrap_plots(p, empty, widths = c(1, 1)),
         width = 15, height = 10, dpi = 300, bg = "white")
  quit(status = 0)
}
parse_or_ci <- function(value) {
  pattern <- "[0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?"
  values <- as.numeric(regmatches(value, gregexpr(pattern,value,perl=TRUE))[[1]])
  if(length(values)!=3L || any(!is.finite(values)) || any(values<=0) || values[2]>values[1] || values[3]<values[1]) stop("Invalid OR (95% CI): ",value)
  values
}
estimates <- t(vapply(forest[["OR (95% CI)"]],parse_or_ci,numeric(3)))
forest <- forest %>% mutate(OR=estimates[,1], CI_low=estimates[,2], CI_high=estimates[,3], OR_CI_label=sprintf("%.3f (%.3f - %.3f)", OR, CI_low, CI_high)) %>% arrange(disease, x)
mv <- flt_mv()
stopifnot(all(c("Exposure","Outcome","Method","OR (95%CI)","P-value") %in% names(mv)))
mv <- mv %>% filter(Method == "MV-IVW")
if(anyDuplicated(mv[c("Exposure","Outcome")])) stop("Duplicate MV-IVW exposure/outcome rows")
mi <- match(pair_key(forest),pair_key(mv))
forest <- forest %>% mutate(MVMR_status=ifelse(is.na(mi), "missing", "matched"), MVMR_method=mv$Method[mi], MVMR_P=as.numeric(mv[["P-value"]][mi]))
if(any(!is.na(forest$MVMR_P) & (!is.finite(forest$MVMR_P) | forest$MVMR_P<0 | forest$MVMR_P>1))) stop("Invalid MVMR P")
mv_estimates <- matrix(NA_real_,nrow(forest),3)
for(i in which(!is.na(mi))) mv_estimates[i, ] <- parse_or_ci(mv[["OR (95%CI)"]][mi[i]])
forest <- forest %>% mutate(MVMR_OR=mv_estimates[,1], MVMR_CI_low=mv_estimates[,2], MVMR_CI_high=mv_estimates[,3], MVMR_OR_CI_label=ifelse(is.na(mi), "NA", sprintf("%.3f (%.3f - %.3f)", MVMR_OR, MVMR_CI_low, MVMR_CI_high)))
message("MV-IVW matches: ",sum(!is.na(mi))," / ",nrow(forest))
if(anyNA(mi)) warning("Missing MVMR matches retained as NA")
shown_outcomes <- outcomes[outcomes %in% forest$Outcome]
# 为每个结局分配标题和间隔，文字与置信区间共用行坐标。
cursor <- 0
heading_rows <- numeric(length(shown_outcomes))
forest$row <- NA_real_
for (j in seq_along(shown_outcomes)) {
  cursor <- cursor+1
  heading_rows[j] <- cursor
  ii <- which(forest$Outcome==shown_outcomes[j])
  forest$row[ii] <- cursor+seq_along(ii)
  cursor <- cursor+length(ii)+.7
}
total_rows <- cursor+1.8
forest$y <- total_rows-forest$row
heading_y <- total_rows-heading_rows
# 两种方法共用 OR 轴，行内上下错开。
forest$UVMR_P <- forest$IVW_P
if(anyNA(forest$UVMR_P) || any(forest$UVMR_P<0 | forest$UVMR_P>1)) stop("Invalid UVMR P")
fmt_p <- function(x) vapply(x,function(v) if(is.na(v)) "NA" else
  if(v<.001) formatC(v,format="e",digits=2) else formatC(v,format="f",digits=3),character(1))
forest <- forest %>% mutate(UVMR_P_label=fmt_p(UVMR_P), FDR_label=fmt_p(FDR), MVMR_P_label=fmt_p(MVMR_P))
log_range <- range(log(c(forest$CI_low,forest$CI_high,forest$MVMR_CI_low,forest$MVMR_CI_high,1)),na.rm=TRUE)
log_range <- log_range+c(-1,1)*diff(log_range)*.06
# 保留 OR 数值标签，便于检查数据。
forest$UVMR_OR_label <- sprintf("%.3f", forest$OR)
forest$MVMR_OR_label <- ifelse(is.na(forest$MVMR_OR), "NA", sprintf("%.3f", forest$MVMR_OR))
# 数值列右对齐，OR 轴采用对数刻度。
table_x <- c(uv_ci=.80, uv_fdr=.95, mv_ci=1.19, mv_p=1.35)
right_edge <- 1.36
to_x <- function(or) .30+(log(or)-log_range[1])/diff(log_range)*.28
ticks <- c(.25,.5,1,2,4,8)
ticks <- ticks[log(ticks)>=log_range[1] & log(ticks)<=log_range[2]]
uv <- data.frame(y=forest$y+.16,OR=forest$OR,low=forest$CI_low,high=forest$CI_high, category=forest$category,method="UVMR-IVW")
mv_points <- data.frame(y=forest$y-.16,OR=forest$MVMR_OR,low=forest$MVMR_CI_low, high=forest$MVMR_CI_high,category=forest$category,method="MVMR-IVW")
points_data <- rbind(uv,mv_points)
points_data$method <- factor(points_data$method,levels=c("UVMR-IVW","MVMR-IVW"))
heads <- data.frame(y=heading_y,label=paste0(shown_outcomes," (n = ", vapply(shown_outcomes,function(z) sum(forest$Outcome==z),integer(1)),")"))
axes <- data.frame(x=to_x(ticks),label=format(ticks,trim=TRUE))
q <- ggplot() +
  geom_rect(data=heads,aes(xmin=0,xmax=right_edge,ymin=y-.4,ymax=y+.4),fill="#EDF2F5",color=NA) +
  geom_text(data=heads,aes(x=.005,y=y,label=label),hjust=0,fontface="bold",size=4.1) +
  annotate("segment",x=to_x(1),xend=to_x(1),y=min(forest$y)-.45,yend=max(forest$y)+.45, linetype=2,color="#708898",linewidth=.35) +
  geom_segment(data=points_data,aes(x=to_x(low),xend=to_x(high),y=y,yend=y,color=method,linetype=method),linewidth=.45,na.rm=TRUE) +
  geom_point(data=points_data,aes(x=to_x(OR),y=y,color=method,shape=method),size=2,na.rm=TRUE) +
  geom_rect(data=forest,aes(xmin=.007,xmax=.015,ymin=y-.30,ymax=y+.30,fill=category),color=NA) +
  geom_text(data=forest,aes(x=.027,y=y,label=label),hjust=0,size=3.8,color="#333333") +
  lapply(seq_along(table_x), function(i)
    geom_text(data=forest, aes(x=table_x[i], y=y, label=.data[[c("OR_CI_label", "FDR_label", "MVMR_OR_CI_label", "MVMR_P_label")[i]]]), hjust=1, size=3.6)) +
  annotate("text",x=c(.027,.44,table_x),y=total_rows+.4, label=c("Exposure","OR (95% CI)","UVMR OR (95% CI)","UVMR P-FDR","MVMR OR (95% CI)","MVMR P*"), hjust=c(0,.5,1,1,1,1),fontface="bold",size=3.8) +
  annotate("segment",x=0,xend=right_edge,y=total_rows-.05,yend=total_rows-.05,linewidth=.3) +
  annotate("segment",x=.30,xend=.58,y=1.1,yend=1.1,linewidth=.3) +
  geom_segment(data=axes,aes(x=x,xend=x,y=1.1,yend=.85),linewidth=.3) +
  geom_text(data=axes,aes(x=x,y=.4,label=label),size=3.6) +
  scale_fill_manual(values=cols,guide="none") +
  scale_color_manual(values=c("UVMR-IVW"="#2C5275","MVMR-IVW"="#D87936"),name=NULL) +
  scale_shape_manual(values=c("UVMR-IVW"=16,"MVMR-IVW"=15),name=NULL) +
  scale_linetype_manual(values=c("UVMR-IVW"="solid","MVMR-IVW"="dashed"),name=NULL) +
  scale_x_continuous(limits=c(0,right_edge),expand=expansion(mult=c(.005,.005))) +
  scale_y_continuous(limits=c(0,total_rows+1),expand=expansion(mult=0)) +
  labs(caption="* MVMR P is unadjusted. OR axis: log scale.\nColour strips: Fig. 1 food categories.") +
  theme_void(base_size=13) +
  theme(legend.position="top",legend.justification="left",legend.key.width=grid::unit(1.1,"cm"), plot.title=element_text(face="bold",size=14),plot.subtitle=element_text(size=10), plot.caption=element_text(size=10.5,hjust=0),plot.margin=margin(8,6,8,4))
# 4. 按 1:2 等高拼图并导出。
left_width <- 7
right_width <- 14
combined_width <- left_width+right_width
draw_combined <- function() {
  grid::grid.newpage()
  print(p,vp=grid::viewport(x=left_width/(2*combined_width),y=.5,width=left_width/combined_width,height=1))
  print(q,vp=grid::viewport(x=(left_width+right_width/2)/combined_width,y=.5,width=right_width/combined_width,height=1))
}
flt_png(file.path(output_dir, "Fig2.png"), combined_width, 10)
tryCatch(draw_combined(),finally=dev.off())
ggsave(file.path(output_dir,"Fig2_A.png"),p,width=left_width,height=10,dpi=180,bg="white")
ggsave(file.path(output_dir,"Fig2_B.png"),q,width=right_width,height=10,dpi=180,bg="white")
message("Combined figure saved: ",nrow(forest)," pairs")
print(counts)
