dir0 <- "D:"
project <- file.path(dir0, "projects/flt")
outdir <- file.path(project, "analysis/01intake")
figdir <- file.path(outdir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# 直接读取统计脚本保存的 CSV，保留表中的显示顺序。
read_result <- function(name) read.csv(file.path(outdir, name), stringsAsFactors = FALSE)
group_results <- read_result("spearman_selected_actual_liking.csv")
selected_items <- read_result("spearman_selected_intake_items.csv")
radar_profile <- read_result("radar_project_intake.csv")

group_names <- unique(group_results$intake_field)
liking <- unique(group_results[c("food_liking", "food_liking_label")])
liking_selected <- setNames(liking$food_liking, liking$food_liking_label)
radar_order <- unique(radar_profile$intake_group)
all_fields <- unique(selected_items$intake_field)

# 食物名称取自字典；未匹配的字段使用原始 ID。
item_labels <- setNames(all_fields, all_fields)
dictionary <- file.path(dir0, "data/phe/ukb/phe/common/diet.lst")
if (file.exists(dictionary)) {
  lines <- readLines(dictionary, warn = FALSE)
  lines <- lines[grepl("^p26[0-9]+\\t.*#", lines)]
  names_by_id <- setNames(trimws(sub("^.*#", "", lines)), paste0("w.", sub("\\t.*$", "", lines)))
  matched <- intersect(all_fields, names(names_by_id))
  item_labels[matched] <- names_by_id[matched]
}
item_labels <- vapply(item_labels, function(x) paste(strwrap(x, width = 28), collapse = "\n"), "")
radar_colors <- c("#3366A5", "#16854B", "#C66A27")

save_plot <- function(draw, name, width, height, device = pdf) {
  device(file.path(figdir, paste0(name, ".pdf")), width = width, height = height)
  tryCatch(draw(), finally = dev.off())
}
draw_heatmap <- function(results, scores, fields=group_names, labels=fields,
                         limit=max(.01, abs(results$cor), na.rm=TRUE)) {
  # Keep item heatmaps readable at the same page height as the radar plot.
  if (length(fields)>24L) {
    for (idx in split(seq_along(fields), ceiling(seq_along(fields)/24))) {
      draw_heatmap(results[results$intake_field %in% fields[idx], ], scores, fields[idx], labels[idx], limit)
    }
    return(invisible(NULL))
  }
  nx <- length(scores); ny <- length(fields)
  x <- match(results$food_liking, scores)
  y <- ny + 1 - match(results$intake_field, fields)
  stopifnot(!anyNA(x), !anyNA(y), !anyDuplicated(paste(x, y)), nrow(results)==nx*ny)
  colors <- colorRampPalette(c("#29365F", "#80ACD3", "#FFFFFF", "#EA8582", "#852121"))(201)
  layout(matrix(1:2, nrow=1), widths=c(par("din")[1]-1.2, 1.2))
  par(oma=c(2,0,0,0), mar=c(1,13,8,1), family="sans")
  plot(NA, xlim=c(.5,nx+.5), ylim=c(.5,ny+.5), asp=1, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
  yrange <- diff(par("usr")[3:4])
  rect(x-.5, y-.5, x+.5, y+.5, border="grey75")
  for (i in seq_len(nrow(results))) {
    if (!is.finite(results$cor[i])) { text(x[i],y[i],"NA",cex=.65); next }
    value <- results$cor[i]/limit
    h <- .46*sqrt(abs(value))
    rect(x[i]-h,y[i]-h,x[i]+h,y[i]+h, col=colors[101+round(100*value)],border=NA)
    if (is.finite(results$p[i]) && results$p[i]>=.05) points(x[i],y[i],pch=4,cex=.8,lwd=1.4)
  }
  text(.3,ny:1,labels,adj=1,cex=.8,xpd=NA)
  text(seq_len(nx),ny+.7,names(scores),srt=55,adj=0,cex=.85,xpd=NA)
  par(mar=c(1,.5,8,3.5))
  plot(NA,xlim=c(0,1),ylim=c(-1,1)*limit*yrange/ny, axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
  edges <- seq(-limit,limit,length.out=202)
  rect(0,edges[-202],.5,edges[-1],col=colors,border=NA)
  ticks <- seq(-limit,limit,length.out=5)
  axis(4,at=ticks,labels=sprintf("%.2f",ticks),pos=.5,las=1,cex.axis=.75)
  mtext("coefficient",side=3,line=.5,cex=.7)
  mtext("Cells: Spearman coefficients; square area: |coefficient|.", outer=TRUE,side=1,line=.9,cex=.8)
  mtext("x: nominal P >= 0.05; NA: unavailable estimate.",outer=TRUE,side=1,line=0,cex=.8)
}

draw_radar <- function(profile) {
  theta <- pi / 2 - 2 * pi * (seq_along(radar_order) - 1) / length(radar_order)
  x <- cos(theta)
  y <- sin(theta)
  closed <- c(seq_along(theta), 1)
  circle <- seq(0, 2 * pi, length.out = 241)
  axis_labels <- vapply(radar_order, function(label) {paste(strwrap(label, width = 17), collapse = "\n")}, "")
  labels <- names(liking_selected)
  center <- .23

  old_par <- par(mar = rep(0, 4), family = "sans", xaxs = "i", yaxs = "i")
  on.exit(par(old_par))
  plot(NA, xlim = c(-1.84, 1.84), ylim = c(-1.24, 1.24), asp = 1, axes = FALSE, xlab = "", ylab = "")
  for (ring in seq(.2, 1, by = .2)) {lines(center + ring * cos(circle), ring * sin(circle), col = "#A7A7A7", lty = if (ring < 1) 2 else 1, lwd = .8)}
  segments(center, 0, center + x, y, col = "#A7A7A7", lwd = .72)

  for (i in seq_along(labels)) {
    rows <- profile[profile$food_liking == labels[i], ]
    r <- rows$min_max[match(radar_order, rows$intake_group)]
    lines(center + r[closed] * x[closed], r[closed] * y[closed], col = radar_colors[i], lwd = 2.3, lend = "round", ljoin = "round")
    points(center + r * x, r * y, pch = 16, col = radar_colors[i], cex = .75)
  }

  for (i in seq_along(axis_labels)) {
    align <- if (abs(x[i]) < .15) .5 else if (x[i] > 0) 0 else 1
    text(center + 1.055 * x[i], 1.075 * y[i], axis_labels[i], adj = c(align, .5), cex = 1.03, col = "#222222")
  }
  ns <- profile$n_liking[match(labels, profile$food_liking)]
  legend(x = -1.76, y = 1.15,
         legend = paste0(labels, " (N = ", format(ns, big.mark = ",", trim = TRUE), ")"),
         col = radar_colors, text.col = "#222222", lwd = 2.3, pch = 16, pt.cex = .75,
         seg.len = 1.7, x.intersp = .7, y.intersp = 1.22, bty = "n", cex = .88)
}
plot_height <- 7
save_plot(function() draw_heatmap(group_results, liking_selected), "spearman_actual_liking_heatmap", 7, plot_height)
save_plot(function() draw_heatmap(selected_items, liking_selected, all_fields, item_labels), "spearman_individual_intake_heatmap", 7, plot_height)
radar_device <- if (capabilities("cairo")) cairo_pdf else pdf
save_plot(function() draw_radar(radar_profile), "radar_intake", 10, plot_height, radar_device)
message("PDFs saved: ", figdir)
