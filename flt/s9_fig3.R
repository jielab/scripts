#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig. 3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)

dir0 <- "D:"
project <- file.path(dir0, "projects/flt")
input_file <- file.path(project, "analysis/result.xlsx")
flt_file <- file.path(dir0, "files/flt.xlsx")
intake <- read.csv(file.path(project,"analysis/01intake/spearman_selected_actual_liking.csv"))
profile <- read.csv(file.path(project,"analysis/01intake/radar_project_intake.csv"))
lm_results <- read.csv(file.path(project,"analysis/01lm/lm_results.csv"))
output_dir <- file.path(project, "analysis/main_plot")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Fig. 3A
# 读取主要估计方法，并清理列名。
read_sheet <- function(sheet, method) {readxl::read_excel(input_file, sheet=sheet) %>% rename_with(trimws) %>% filter(trimws(Method) == method)}
pvalue <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & trimws(as.character(x)) != "" & is.na(v)) || any(!is.na(v) & (!is.finite(v) | v < 0 | v > 1))) stop("Invalid P value")
  v
}
uv <- read_sheet("Sheet2", "IVW") %>% mutate(`P-FDR`=pvalue(`P-FDR`)) %>% filter(`P-FDR` < alpha)
mv <- read_sheet("Sheet7", "MV-IVW")
rv <- read_sheet("Sheet3", "IVW")
meta <- as.data.frame(readxl::read_excel(flt_file, sheet = "item182"))
mi <- match(uv$Exposure, meta$gwas_id)
h <- readxl::read_excel(input_file, sheet = "Sheet1", col_names = FALSE, n_max = 2, .name_repair = "minimal")
stopifnot(h[[12]][1] == "M2", h[[12]][2] == "HR", h[[16]][2] == "P-FDR")
cx <- as.data.frame(readxl::read_excel(input_file, sheet = "Sheet1", skip = 2, col_names = FALSE, .name_repair = "minimal"))
cx <- cx[, c(1, 2, 4, 12, 16)]
names(cx) <- c("name", "label", "Outcome", "HR", "FDR")
key <- function(a, b) paste(trimws(a), trimws(b), sep = "||")
uk <- key(uv$Exposure, uv$Outcome)
ck <- key(cx$name, cx$Outcome)
mk <- key(mv$Exposure, mv$Outcome)
# Reverse MR swaps disease exposure and food outcome.
rk <- key(rv$Outcome, rv$Exposure)
stopifnot(!anyDuplicated(uk), !anyDuplicated(ck), !anyDuplicated(mk), !anyDuplicated(rk))
ci <- match(key(meta$name[mi], uv$Outcome), ck)
vi <- match(uk, mk)
ri <- match(uk, rk)
direction_or <- function(x) {
  v <- suppressWarnings(as.numeric(sub("^\\s*([0-9.eE+-]+).*", "\\1", as.character(x))))
  ifelse(is.na(v) | v <= 0 | v == 1, NA_character_, ifelse(v < 1, "protective", "harmful"))
}
d <- data.frame(Outcome = uv$Outcome, Exposure = uv$Exposure,
                Food_liking = meta$item_raw[mi], FDR_cox = pvalue(cx$FDR[ci]),
                Direction_cox = direction_or(cx$HR[ci]),
                P_FDR_IVW_mr = uv[["P-FDR"]],
                Direction_mr = direction_or(uv[["OR (95% CI)"]]),
                pval_mvmr = pvalue(mv[["P-value"]][vi]),
                Direction_mvmr = direction_or(mv[["OR (95%CI)"]][vi]),
                P_value_remr = pvalue(rv[["P-value"]][ri]), check.names = FALSE)
d <- d %>% rename(`Food liking`=Food_liking, `P-value_remr`=P_value_remr) %>% arrange(match(Outcome, unique(Outcome)))
colors <- c(protective = "#465D89", harmful = "#CC2529", nonsignificant = "#F5F5E9", missing = "#DADADA")
alpha <- 0.05
classify <- function(p, direction) {
  direction <- tolower(trimws(as.character(direction)))
  sig <- !is.na(p) & p < alpha
  if (any(sig & (is.na(direction) | !direction %in% c("protective", "harmful")))) stop("Significant result has missing/unknown direction.")
  case_when(is.na(p) ~ "missing", sig ~ direction, .default="nonsignificant")
}
status <- cbind(classify(d$FDR_cox, d$Direction_cox), classify(d$P_FDR_IVW_mr, d$Direction_mr), classify(d$pval_mvmr, d$Direction_mvmr))
checks <- !is.na(d[["P-value_remr"]]) & d[["P-value_remr"]] > alpha
draw_heatmap <- function() {
  n <- nrow(d)
  panel_width <- par("fin")[1]
  par(mai=c(.02,.98,.06,.04)*panel_width/(210/25.4),family="sans",lwd=.4,xaxs="i",yaxs="i",cex=1)
  plot.new()
  plot.window(xlim=c(0,n),ylim=c(-.65,7.7))
  legend("bottomright",legend=c("Positive","Negative","Non-significant"),col=colors[c("harmful","protective","nonsignificant")],pch=16,horiz=TRUE,bty="n",cex=1,inset=0,x.intersp=.6,y.intersp=.8)
  runs <- rle(as.character(d$Outcome))
  ends <- cumsum(runs$lengths); starts <- c(0,head(ends,-1))
  outcome_labels <- c(dementia_r13="Dementia",anystroke_eur="Stroke", sleep_disorder_r13="Sleep disorders",mdd_2025="Depression",anxiety_2026="Anxiety disorders")
  for(g in seq_along(ends)) {
    rect(starts[g],6.9,ends[g],7.7,col="#E3EAF0",border="#333333")
    lab <- outcome_labels[runs$values[g]]
    if(is.na(lab)) lab <- runs$values[g]
    text((starts[g]+ends[g])/2,7.3,lab)
    segments(c(starts[g],ends[g]),2.8,c(starts[g],ends[g]),6.9,col="#333333")
  }
  text(seq_len(n)-.5,3.02,d[["Food liking"]],srt=90,adj=c(0,.5),xpd=NA)
  titles <- c("Cox regression", "Univariable MR", "Multivariable MR", "Reverse MR")
  for(k in 1:4) {
    top <- 2.8-(k-1)*.7
    bottom <- top-.7
    text(-.16,(top+bottom)/2,titles[k],adj=c(1,.5),xpd=NA)
    fill <- if(k<=3) status[,k] else ifelse(is.na(d[["P-value_remr"]]), "missing", "nonsignificant")
    rect(seq_len(n)-1,bottom,seq_len(n),top,col=colors[fill],border="#333333")
    if(k==4) for(i in which(checks))
      lines(c(i-.76,i-.56,i-.22),bottom+c(.34,.16,.54),lwd=.8)
  }
  text(-.16,5.05,"Food-liking traits",adj=c(1,.5),xpd=NA)
}
save_pdf <- function(name, draw, height) {
  target <- file.path(output_dir, paste0(name, ".pdf"))
  cairo_pdf(target, width=210/25.4, height=height/25.4, pointsize=8)
  tryCatch(draw(), finally=dev.off())
  message("Saved: ", target)
}
save_pdf("Fig3_A", draw_heatmap, 62)
audit <- data.frame(Outcome = d$Outcome, Food_liking = d[["Food liking"]], Cox = status[,1], MR = status[,2], MVMR = status[,3], pval_mvmr = d$pval_mvmr, Cox_M2_FDR = d$FDR_cox, UVMR_FDR = d$P_FDR_IVW_mr, Exposure = d$Exposure, Reverse_MR_P = d[["P-value_remr"]], Reverse_MR_tick = checks)
write.csv(audit, file.path(output_dir, "Fig3A_audit.csv"), row.names = FALSE, na = "")

## Fig. 3B-D
liking <- c(`Butter on bread`="p20628", `Total salad vegetable`="F-Salad_vegetables", Chicken="p20638")
groups <- unique(intake$intake_field)
brain_blocks <- list(
	Cognitive=c("fluid_intelligence.z","reaction_time.z","pairs_matching.z","prospective_memory.b","numeric_memory.z"),
	`Mental health`=c("general_mental_health.b","life_satisfaction.z","subjective_wellbeing.z","phq4.z","phq9.z","cidi.z","gad7.z")
)
brain_labels <- c(
	"Fluid intelligence","Reaction time","Pairs matching","Prospective memory","Numeric memory",
	"General mental health","Life satisfaction","Subjective wellbeing","PHQ-4","PHQ-9","CIDI","GAD-7"
)
palette <- colorRampPalette(c("#29365F","#80ACD3","#FFFFFF","#EA8582","#852121"))(201)
wrap <- function(x,n=18) vapply(x,function(z) paste(strwrap(z,n),collapse="\n"),"")
glyph <- function(x,y,value,limit,nonsig) {
 rect(x-.5,y-.5,x+.5,y+.5,border="grey75",lwd=.4)
 if (!is.finite(value)) {text(x,y,"NA"); return(invisible(NULL))}
 v <- max(-1,min(1,value/limit))
 # Convert physical square size to each axis's units, even in wide panels.
 units <- c(diff(par("usr")[1:2]),diff(par("usr")[3:4]))/par("pin")
 half <- .46*min(1/units)*sqrt(abs(v))*units
 rect(x-half[1],y-half[2],x+half[1],y+half[2],col=palette[101+round(100*v)],border=NA)
 if (isTRUE(nonsig)) points(x,y,pch=4,cex=.8,lwd=.7)
}
scale_bar <- function(x,y,w,h,limit,title,clip=FALSE) {
 edges <- seq(x,x+w,length.out=202)
 rect(head(edges,-1),y,tail(edges,-1),y+h,col=palette,border=NA)
 labs <- sprintf("%.2f",c(-limit,0,limit))
 if(clip) labs[c(1,3)] <- paste0(c("<=",">="),labs[c(1,3)])
 text(c(x,x+w/2,x+w),y-.25,labs,adj=c(.5,1))
 text(x+w/2,y+h+.25,title,adj=c(.5,0))
}
draw_B <- function() {
 par(mai=c(1.25,1.25,.4,.06),family="sans",xaxs="i",yaxs="i",cex=1)
 plot(NA,xlim=c(.5,length(groups)+.5),ylim=c(-2,4.5),axes=FALSE,xlab="",ylab="")
 limit <- max(.01,abs(intake$cor),na.rm=TRUE)
 for(i in seq_len(nrow(intake))) glyph(match(intake$intake_field[i],groups),4-match(intake$food_liking[i],liking),intake$cor[i],limit,intake$p[i]>=.05)
 text(.3,3:1,names(liking),adj=1,xpd=NA)
 text(seq_along(groups),.3,groups,srt=90,adj=c(1,.5),xpd=NA)
 par(xpd=NA)
 scale_bar(6,3.9,4,.22,limit,"Spearman coefficient")
 mtext("Square area: |coefficient|; x: nominal P >= 0.05.",1,line=5)
}
draw_C <- function() {
 par(mar=c(3,1,2,1),family="sans",cex=1)
 theta <- pi/2-2*pi*(seq_along(groups)-1)/length(groups)
 x <- cos(theta); y <- sin(theta); closed <- c(seq_along(groups),1)
 plot(NA,xlim=c(-1.85,1.85),ylim=c(-1.55,1.4),asp=1,axes=FALSE,xlab="",ylab="")
 circle <- seq(0,2*pi,length.out=241)
 for(ring in seq(.2,1,.2)) lines(ring*cos(circle),ring*sin(circle),col="grey70",lty=2,lwd=.5)
 segments(0,0,x,y,col="grey70",lwd=.5)
 cols <- c("#3366A5","#16854B","#C66A27")
 for(i in seq_along(liking)) {
 rows <- profile %>% filter(food_liking == names(liking)[i])
 r <- rows$min_max[match(groups,rows$intake_group)]
 lines(r[closed]*x[closed],r[closed]*y[closed],col=cols[i],lwd=1)
 points(r*x,r*y,col=cols[i],pch=16,cex=.5)
 }
 for(i in seq_along(groups)) text(1.08*x[i],1.09*y[i],wrap(groups)[i],adj=c(if(abs(x[i])<.15) .5 else if(x[i]>0) 0 else 1,.5))
 legend("bottom",legend=names(liking),col=cols,lwd=1,bty="n",cex=1)
 mtext("Liking score = 9; within-axis min-max scaled mean intake.",1,line=1)
}
draw_D <- function() {
 par(mai=c(1.2,1.25,.6,.06),family="sans",cex=1)
 traits <- unlist(brain_blocks,use.names=FALSE)
 xs <- c(1:5,7:13)
 plot(NA,xlim=c(.5,13.5),ylim=c(-3,5.5),asp=1,axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
 rr <- lm_results %>% filter(trait %in% traits)
 stopifnot(nrow(rr)==36,!anyDuplicated(rr[c("trait","variable")]))
 for(i in seq_len(nrow(rr))) {
 valid <- rr$status[i]=="ok" && is.finite(rr$se[i]) && rr$se[i]>0 && is.finite(rr$FDR[i])
 glyph(xs[match(rr$trait[i],traits)],4-match(rr$variable[i],liking),if(valid) rr$beta[i] else NA_real_,.2,rr$FDR[i]>=.05)
 }
 text(.25,3:1,names(liking),adj=1,xpd=NA)
 text(xs,.25,brain_labels,srt=90,adj=c(1,.5),xpd=NA)
 text(c(3,10),4,names(brain_blocks),xpd=NA)
 par(xpd=NA)
 scale_bar(5,6,4,.25,.2,"Beta",TRUE)
 mtext("x: FDR >= 0.05 (BH within phenotype).",1,line=2)
 mtext("Colour / area capped at |beta| = 0.20.",1,line=3)
}
draw_figure3 <- function() {
 regions <- list(c(0,1,.70,1),c(0,1,.40,.70),c(0,.50,0,.40),c(.50,1,0,.40))
 draws <- list(draw_heatmap,draw_B,draw_C,draw_D)
 for(i in 1:4) {
 par(fig=regions[[i]],new=i>1,cex=1,ps=8,font=1,font.axis=1,font.lab=1,oma=rep(0,4))
 draws[[i]]()
 par(fig=regions[[i]],new=TRUE,mar=rep(0,4))
 plot.new(); plot.window(c(0,1),c(0,1))
 text(.01,.98,LETTERS[i],adj=c(0,1),cex=1,font=1)
 }
}
save_pdf("Fig3", draw_figure3, 297)
