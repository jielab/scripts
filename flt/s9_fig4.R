#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig. 4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(ggplot2)
library(grid)
dir0 <- "D:/"
project <- file.path(dir0, "projects/flt")
output_dir <- file.path(project, "analysis/main_plot")
mri_dir <- file.path(project, "analysis/01lm/mri")
protein_dir <- file.path(project, "analysis/01protein")
alpha <- .05

read_result <- function(path) read.csv(path, check.names=FALSE)
map <- read_result(file.path(protein_dir, "limma/food_pair_mapping.csv"))
foods <- map$food
food_cols <- setNames(c("#39699E", "#328450", "#B96B35"), foods)
effect_cols <- c(Positive="#C83E35", Negative="#39699E", Other="#999999")
fmt <- function(x) ifelse(is.finite(x), formatC(x, digits=3, format="g"), "NA")

# A：复用 MRI 原始缓存；色阶按原渲染脚本计算，各食物分别显示。
mri <- read_result(file.path(mri_dir, "mri_lm_results.csv")) %>% filter(exposure %in% map$liking, family %in% c("ThicknessDK", "FA27", "MD27"))
stopifnot(!anyDuplicated(mri[c("exposure", "field_id")]), all(mri$status == "ok"))
images <- lapply(map$liking, function(id) {
  root <- file.path(mri_dir, "figures/anatomy", id, "cache")
  list(th=png::readPNG(file.path(root, "surface_renders/a_thickness_FDR05.png")),
       FA=png::readPNG(file.path(root, "tract_renders/FA27.png")),
       MD=png::readPNG(file.path(root, "tract_renders/MD27.png")))
})
text_at <- function(label, x, y, size=8, ...) grid.text(label, x, y, gp=gpar(fontsize=size, ...))
box <- function(x, y, w, h, fill, border=NA)
  grid.rect(x, y, w, h, gp=gpar(fill=fill, col=border, lwd=.5))
bar <- function(x, y, w, limit) {
  pal <- colorRampPalette(c("#053061", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#67001F"))(201)
  grid.raster(as.raster(matrix(pal, nrow=1)), x, y, width=w, height=.018)
  text_at(formatC(c(-limit, 0, limit), digits=3, format="f"), x+c(-.5,0,.5)*w, y-.025, 6)
}
draw_A <- function() {
  text_at("A", .008, .975, 10, fontface="bold")
  box(.515, .95, .95, .07, "#E3EAF0")
  text_at(c("Cortical thickness", "FA", "MD"), c(.36,.77,.92), .95)
  for(i in seq_along(foods)) {
    y <- .78-(i-1)*.27
    box(.026,y,.028,.27,food_cols[i])
    grid.text(gsub("Total salad vegetable", "Salad vegetables", foods[i]),.026,y,rot=90,
              gp=gpar(fontsize=8,col="white",fontface="bold"))
    # 原图左侧为四个脑视图；去除自带色条和空白，按原值重画色条。
    im <- images[[i]]$th[120:380,1:1400,,drop=FALSE]
    for(j in 1:4) grid.raster(im[,((j-1)*350+1):(j*350),,drop=FALSE],
                              .04+(j-.5)*.16,y+.025,height=.22)
    for(j in 1:2) grid.raster(images[[i]][[j+1]],c(.77,.92)[j],y+.025,height=.21)
    z <- mri %>% filter(exposure == map$liking[i])
    th_limit <- max(abs(z$beta[z$family == "ThicknessDK"]),na.rm=TRUE)
    wm <- z %>% filter(family != "ThicknessDK", FDR_family < alpha)
    wm_limit <- max(.02,abs(wm$beta),na.rm=TRUE)
    bar(.36,y-.09,.16,th_limit)
    bar(.845,y-.09,.14,wm_limit)
    box(.515,y,.95,.27,NA,"#777777")
  }
  text_at("Beta: MRI SD per liking-score point; row-specific scales. Coloured regions: family FDR < 0.05.",.52,.035,7)
}

# B：limma 高低组差异；仅标注偏好与摄入均显著且同向的蛋白。
shared <- read_result(file.path(protein_dir,"limma/shared_same_direction_proteins.csv")) %>% filter(adj.P.Val_liking < alpha, adj.P.Val_intake < alpha, effect_liking*effect_intake > 0)
set.seed(42)
strip <- read_result(file.path(protein_dir,"limma/all_limma_results.csv")) %>%
  filter(is.finite(effect), is.finite(adj.P.Val)) %>%
  mutate(food=factor(food,levels=foods), panel=factor(analysis,c("liking","intake"),c("Liking","Intake")),
         x=as.integer(food)+runif(n(),-.40,.40),
         direction=if_else(adj.P.Val < alpha,if_else(effect > 0,"Positive","Negative"),"Other"))
plots_B <- lapply(levels(strip$panel),function(panel_name) {
  z <- strip %>% filter(panel == panel_name)
  gap <- max(abs(z$effect))*.10
  z <- z %>% mutate(y=effect+sign(effect)*gap)
  labels <- z %>% semi_join(shared,by=c("food","protein")) %>%
    group_by(food,sign=sign(effect)) %>% arrange(x,.by_group=TRUE) %>%
    mutate(lx=as.integer(food)+if_else(row_number()<=ceiling(n()/2),-.23,.23),
           ly=sign*(max(abs(y))+gap*(1+2.1*((row_number()-1) %% ceiling(n()/2))))) %>% ungroup()
  bands <- z %>% group_by(food) %>% summarise(x=first(as.integer(food)),low=min(y),high=max(y),.groups="drop")
  ggplot(z,aes(x,y)) +
    geom_rect(data=bands,aes(xmin=x-.45,xmax=x+.45,ymin=low,ymax=high),inherit.aes=FALSE,fill="#EBECEF") +
    geom_point(aes(color=direction),size=.55,alpha=.75) +
    geom_rect(data=bands,aes(xmin=x-.5,xmax=x+.5,fill=food),ymin=-gap,ymax=gap,inherit.aes=FALSE) +
    geom_text(data=bands,aes(x,0,label=c("Butter on\nbread","Salad\nvegetables","Chicken")[x]),
              inherit.aes=FALSE,color="white",size=2.4,fontface="bold",lineheight=.9) +
    geom_segment(data=labels,aes(xend=lx,yend=ly),color="#7A8591",linewidth=.2) +
    geom_label(data=labels,aes(lx,ly,label=protein),size=2,
      label.padding=unit(.07,"lines"),linewidth=.15) +
    scale_color_manual(values=effect_cols) + scale_fill_manual(values=food_cols) +
    scale_x_continuous(limits=c(.4,3.6),breaks=NULL) +
    scale_y_continuous(breaks=function(lim) {v <- pretty(range(z$effect),4); v+sign(v)*gap},
      labels=function(v) format(sign(v)*pmax(0,abs(v)-gap),digits=2,trim=TRUE),expand=expansion(mult=.13)) +
    labs(title=panel_name,x=NULL,y="Adjusted protein difference") + theme_classic(base_size=8) +
    theme(legend.position="none",axis.line.x=element_blank(),plot.title=element_text(hjust=.5,size=8),
          plot.margin=margin(4,6,2,4))
})
draw_B <- function() {
  for(i in 1:2) print(plots_B[[i]],vp=viewport(x=.5,y=c(.76,.28)[i],width=1,height=.47))
  text_at("B",.015,.99,10,fontface="bold")
  text_at("Input scale; FDR < 0.05 coloured; gap at zero.",.52,.018,6)
}

# C：同向交集蛋白的连续偏好回归与后续 Cox；不重新筛选 Cox 显著性。
linear <- read_result(file.path(protein_dir,"linear/FLT_protein_linear_all.csv")) %>%
  semi_join(shared,by=c("food","protein")) %>% arrange(match(food,foods),protein)
cox <- list.files(file.path(protein_dir,"cox"),"^protein_cox_all.csv$",recursive=TRUE,full.names=TRUE) %>%
  lapply(read_result) %>% bind_rows() %>%
  inner_join(map %>% select(pair_id,food,outcome),by=c("pair_id","food","outcome"))
stopifnot(!anyDuplicated(linear[c("food","protein")]),!anyDuplicated(cox[c("food","protein")]))
forest <- linear %>% left_join(cox %>% select(food,protein,HR,lo=CI_low,hi=CI_high,FDR_outcome),by=c("food","protein")) %>%
  mutate(group=match(food,unique(food)),row=row_number()+group, y=1-row/(n()+4))
stopifnot(nrow(forest)>0,all(forest$status=="ok"),all(is.finite(forest$CI_low)),all(is.finite(forest$CI_high)))
draw_C <- function() {
  text_at("C",.01,.99,10,fontface="bold")
  text_at(c("Protein","Beta (95% CI)","P-FDR","HR (95% CI)","P-FDR"),c(.13,.36,.52,.73,.94),.985)
  beta_range <- range(c(0,forest$CI_low,forest$CI_high),na.rm=TRUE)*1.12
  hr_range <- range(log(c(1,forest$lo,forest$hi)),na.rm=TRUE)+c(-.03,.03)
  bx <- function(v) .25+.22*(v-beta_range[1])/diff(beta_range)
  hx <- function(v) .62+.25*(log(v)-hr_range[1])/diff(hr_range)
  grid.lines(c(bx(0),bx(0)),c(.07,.93),gp=gpar(col="#8797A5",lty=2))
  grid.lines(c(hx(1),hx(1)),c(.07,.93),gp=gpar(col="#8797A5",lty=2))
  step <- 1/(nrow(forest)+4)
  for(food in unique(forest$food)) {
    z <- forest %>% filter(.data$food == .env$food)
    top <- max(z$y)+step
    box(.365,top,.24,step*.8,"#E3EAF0"); box(.79,top,.37,step*.8,"#E3EAF0")
    text_at(food,.365,top,7)
    text_at(if(map$outcome[match(food,foods)]=="sleep_disorder") "Sleep disorders" else "Depression",.79,top,7)
    box(.025,mean(range(z$y)),.035,diff(range(z$y))+step*.9,food_cols[[food]])
    grid.text(food,.025,mean(range(z$y)),rot=90,gp=gpar(fontsize=8,col="white",fontface="bold"))
  }
  for(i in seq_len(nrow(forest))) {
    z <- forest[i,]; y <- z$y
    col <- if(z$FDR_food >= alpha) "#999999" else if(z$beta>0) effect_cols[1] else effect_cols[2]
    box(.14,y,.18,step*.9,"white","#777777")
    text_at(z$protein,.13,y,8)
    grid.segments(.215,y-sign(z$beta)*step*.16,.215,y+sign(z$beta)*step*.16,
                  arrow=arrow(length=unit(1.2,"mm"),type="closed"),gp=gpar(col=col,fill=col))
    grid.segments(bx(z$CI_low),y,bx(z$CI_high),y,gp=gpar(col=col))
    grid.points(bx(z$beta),y,pch=16,size=unit(1.7,"mm"),gp=gpar(col=col))
    text_at(fmt(z$FDR_food),.52,y)
    if(is.finite(z$HR)) {
      col <- if(z$FDR_outcome<alpha) "#AD412C" else "#4D5B6B"
      grid.segments(hx(z$lo),y,hx(z$hi),y,gp=gpar(col=col))
      grid.points(hx(z$HR),y,pch=16,size=unit(1.7,"mm"),gp=gpar(col=col))
    }
    text_at(if(is.finite(z$HR)) fmt(z$FDR_outcome) else "NA",.94,y)
  }
  for(kind in c("beta","hr")) {
    limits <- if(kind=="beta") beta_range else exp(hr_range)
    ticks <- pretty(limits,3); ticks <- ticks[ticks>=limits[1] & ticks<=limits[2]]
    fun <- if(kind=="beta") bx else hx
    grid.lines(fun(limits),c(.065,.065),gp=gpar(col="#657585"))
    text_at(format(ticks,trim=TRUE),fun(ticks),.045,7)
  }
  text_at("Beta: per liking-score point; HR: per protein SD. P-FDR: within food / outcome.",.52,.014,6)
}

# A 在上，B/C 在下；总高度为横版 A4 高度 210 mm，文字基础字号 8 pt。
save_pdf <- function(name,draw,width,height) {
  path <- file.path(output_dir,paste0(name,".pdf"))
  cairo_pdf(path,width=width/25.4,height=height/25.4,pointsize=8)
  tryCatch({grid.newpage(); draw()},finally=dev.off())
  message("Saved: ",path)
}
draw_all <- function() {
  regions <- list(c(.5,.80,1,.40),c(.17,.30,.34,.60),c(.67,.30,.66,.60))
  draws <- list(draw_A,draw_B,draw_C)
  for(i in 1:3) {
    v <- regions[[i]]
    pushViewport(viewport(x=v[1],y=v[2],width=v[3],height=v[4]))
    draws[[i]](); popViewport()
  }
}
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
save_pdf("Fig4_A",draw_A,270,84)
save_pdf("Fig4_B",draw_B,92,126)
save_pdf("Fig4_C",draw_C,178,126)
save_pdf("Fig4",draw_all,270,210)

