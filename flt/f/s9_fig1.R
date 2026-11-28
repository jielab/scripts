source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fig. 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(readxl)

project <- Sys.getenv("FLT_OUT")
flt_file <-  Sys.getenv("FLT_DICTIONARY")

output_dir <- file.path(Sys.getenv("FLT_PAPER"))
alpha <- 0.05

hr_plot_max <- 1.15

# 读取并校验 M2 数据（L:P 列）。
d <- flt_cox()
meta <- as.data.frame(readxl::read_excel(flt_file, sheet = "item182"))
meta <- data.frame(name = meta$name, label = meta$item_raw, category = meta$category)
groups <- unique(meta$category)
stopifnot(length(groups) == 10L)
cols <- setNames(c("#C54C3B", "#4BA7C1", "#80634F", "#415A82", "#DE897C", "#818DB6", "#80BFB3", "#C52C36", "#309B89", "#AA998C"), groups)
meta$type <- ifelse(grepl("^p[0-9]+$", meta$name), "Raw", ifelse(startsWith(meta$name, "F-"), "Derived", "Coffee/Tea"))
meta <- meta[order(match(meta$category, groups), meta$type == "Derived", tolower(meta$label)), ]
meta$x <- seq_len(nrow(meta))
outcomes <- c("All-cause dementia", "Stroke", "Sleep disorders", "Depression", "Anxiety disorders")
stopifnot(nrow(d) == 910L, !anyDuplicated(d[c("name", "disease")]), setequal(d$disease, outcomes), setequal(d$name, meta$name))
idx <- match(d$name, meta$name)
stopifnot(all(d$category == meta$category[idx]), all(table(d$disease) == 182L))
d$x <- meta$x[idx]
stopifnot(all(is.na(d$HR) | (is.finite(d$HR) & d$HR > 0)), all(is.na(d$FDR) | (d$FDR >= 0 & d$FDR <= 1)),
          all(is.na(d$CI_low) | is.na(d$HR) | d$CI_low <= d$HR),
          all(is.na(d$CI_high) | is.na(d$HR) | d$CI_high >= d$HR))
d$significant <- !is.na(d$FDR) & d$FDR < alpha & !is.na(d$HR)
d <- d[order(match(d$disease, outcomes), d$x), ]
# 名称来自结果表，排序来自元数据。
stopifnot(!anyNA(d$label), all(nzchar(trimws(d$label))))
for (id in meta$name) stopifnot(length(unique(d$label[d$name == id])) == 1L)
meta$label <- d$label[match(meta$name, d$name)]
# 保留全部条目，FDR 仅决定点的颜色。
n_items <- nrow(meta)
stopifnot(n_items == 182L, nrow(d) == 910L, all(table(d$name) == length(outcomes)))
group_gap_slots <- 0.55
step <- 2*pi / (n_items + group_gap_slots*length(groups))
group_index <- match(meta$category, groups)
meta$theta <- pi/2 - step*(meta$x - 0.5 + group_gap_slots*(group_index - 0.5))
d$theta <- meta$theta[match(d$name, meta$name)]
categories <- do.call(rbind, lapply(groups, function(group) {
  theta <- meta$theta[meta$category == group]
  data.frame(category = group, low = min(theta) - step/2, high = max(theta) + step/2, middle = mean(range(theta)))
}))
stopifnot(is.na(hr_plot_max) || (is.finite(hr_plot_max) && hr_plot_max > 1))
d$log_hr <- log(d$HR)
if (is.na(hr_plot_max)) hr_plot_max <- exp(max(abs(d$log_hr), log(1.10), na.rm = TRUE))
log_limit <- log(hr_plot_max)
d$shown_log_hr <- pmax(-log_limit, pmin(log_limit, d$log_hr))
d$beyond <- !is.na(d$log_hr) & abs(d$log_hr) > log_limit + 1e-12
ring_width <- 0.37
ring_center <- 1.77 + (seq_along(outcomes) - 1)*0.46
track_outer <- max(ring_center) + ring_width/2
category_band_r <- track_outer + 0.115
effect_col <- c(lower = "#176E9B", higher = "#B64332", other = "#A8ADB1", neutral = "#424C53")
radial_hr <- function(value, k) ring_center[k] + ring_width/2 * value/log_limit

arc <- function(r, a, b, ...) {
  theta <- seq(a, b, length.out = 721)
  lines(r*cos(theta), r*sin(theta), ...)
}

draw_track_labels <- function() {
  # 疾病名称沿各自圆环外侧排列。
  outer_edges <- ring_center + ring_width/2
  next_edges <- c((ring_center - ring_width/2)[-1], category_band_r - 0.02)
  for (k in seq_along(outcomes)) {
    gap <- next_edges[k] - outer_edges[k]
    radius <- outer_edges[k] + gap/2
    chars <- strsplit(outcomes[k], "", fixed = TRUE)[[1]]
    label_cex <- 1.05
    label_cex <- label_cex * min(1, 0.95*gap / strheight("Ag", font = 2, cex = label_cex))
    widths <- strwidth(chars, font = 2, cex = label_cex)
    offsets <- cumsum(widths) - widths/2 - sum(widths)/2
    angles <- pi/2 - offsets/radius

    # 短引线指向对应圆环。
    half_span <- (sum(widths)/2 + 0.025)/radius
    ends <- pi/2 + c(-half_span, half_span)
    segments(radius*cos(ends), radius*sin(ends), outer_edges[k]*cos(ends), outer_edges[k]*sin(ends), col = "#58656C", lwd = 0.65)
    for (j in seq_along(chars)) {
      a <- angles[j]
      text(radius*cos(a), radius*sin(a), chars[j], srt = a*180/pi - 90, adj = c(0.5, 0.5), cex = label_cex, font = 2, col = "#23313A")
    }
  }
  # 中心分类图例。
  legend(0,0,xjust=0.5,yjust=0.5,legend=groups,col=cols[groups],pch=16,title="Groups",bty="o",box.col=NA,bg="white", cex=0.82,pt.cex=0.9,text.font=2,y.intersp=1.05,x.intersp=0.65)
}

draw_circular <- function(combined = FALSE) {
  par(mar = if (combined) c(3.2, 0.2, 1.5, 0.2) else c(4.0, 0.7, 2.1, 0.7), family = "sans", xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(c(-5.65, 5.65), c(-5.48, 5.48), asp = 1)
  if (!combined) title("Food-liking traits and incident brain disorders", cex.main = 1.75)

  # 各圆环采用相同的对称 log(HR) 刻度。
  for (k in rev(seq_along(outcomes))) {
    outer <- ring_center[k] + ring_width/2
    inner <- ring_center[k] - ring_width/2
    symbols(0, 0, circles = outer, inches = FALSE, add = TRUE, fg = "#E6E9EA", bg = if (k %% 2L == 1L) "#F4F6F7" else "#FAFBFC")
    symbols(0, 0, circles = inner, inches = FALSE, add = TRUE, fg = "#E6E9EA", bg = "white")
  }

  # 分类分隔线、效应点与超界三角标记。
  boundary <- categories$high + group_gap_slots*step/2
  for (k in seq_along(outcomes)) {
    inner <- ring_center[k] - ring_width/2
    outer <- ring_center[k] + ring_width/2
    segments(inner*cos(boundary), inner*sin(boundary), outer*cos(boundary), outer*sin(boundary), col = "#D9DEE0", lwd = 0.7)
    arc(ring_center[k], 0, 2*pi, col = "#737E85", lty = 2, lwd = 0.8)
    z <- d[d$disease == outcomes[k] & is.finite(d$log_hr), ]
    r <- radial_hr(z$shown_log_hr, k)
    x <- r*cos(z$theta)
    y <- r*sin(z$theta)
    mark_col <- ifelse(z$log_hr > 0, effect_col[["higher"]], ifelse(z$log_hr < 0, effect_col[["lower"]], effect_col[["neutral"]]))
    nonsig <- !z$significant & !z$beyond
    significant <- z$significant & !z$beyond
    points(x[nonsig], y[nonsig], pch = 16, cex = 0.46, col = effect_col[["other"]])
    if (any(z$significant)) {
      segments(ring_center[k]*cos(z$theta[z$significant]), ring_center[k]*sin(z$theta[z$significant]), x[z$significant], y[z$significant], col = adjustcolor(mark_col[z$significant], 0.7), lwd = 1.05)
    }
    points(x[significant], y[significant], pch = 16, cex = 0.74, col = mark_col[significant])
    for (i in which(z$beyond)) {
      angle <- z$theta[i]
      radial <- c(cos(angle), sin(angle))
      tangential <- c(-sin(angle), cos(angle))
      tip <- r[i]*radial
      base <- tip - sign(z$log_hr[i])*0.053*radial
      vertices <- rbind(tip, base + 0.022*tangential, base - 0.022*tangential)
      color <- if (z$significant[i]) mark_col[i] else effect_col[["other"]]
      polygon(vertices, col = color, border = NA)
    }
  }

  # 分类树的 182 个叶节点与圆环对齐。
  t <- seq(0, 1, length.out = 45)
  basis <- cbind((1-t)^3, 3*(1-t)^2*t, 3*(1-t)*t^2, t^3)
  a <- categories$middle
  segments(0.80*cos(a), 0.80*sin(a), 0.98*cos(a), 0.98*sin(a), col = cols[categories$category], lwd = 1.25)
  for (i in seq_len(n_items)) {
    a <- meta$theta[i]
    b <- categories$middle[match(meta$category[i], categories$category)]
    color <- cols[meta$category[i]]
    xy <- basis %*% rbind(0.98*c(cos(b), sin(b)), 1.13*c(cos(b), sin(b)), 1.28*c(cos(a), sin(a)), 1.49*c(cos(a), sin(a)))
    lines(xy, col = adjustcolor(color, 0.58), lwd = 0.72)
    points(1.49*cos(a), 1.49*sin(a), pch = 16, cex = 0.43, col = color)
  }

  # 食物名称及外侧分类弧线。
  label_r <- track_outer + 0.18
  label_cex <- if (combined) 0.75 else 0.85
  item_cex <- label_cex * pmin(1, (4.87 - label_r)/strwidth(meta$label, cex = label_cex))
  for (i in seq_len(n_items)) {
    a <- meta$theta[i]
    text(label_r*cos(a), label_r*sin(a), meta$label[i], srt = a*180/pi + if (cos(a) < 0) 180 else 0, adj = if (cos(a) < 0) 1 else 0, cex = item_cex[i], col = "#333E45")
  }
  for (i in seq_len(nrow(categories))) {
    g <- categories[i, ]
    color <- cols[g$category]
    arc(category_band_r, g$low, g$high, col = color, lwd = 4)
    # 分类名逐字沿弧排列，下半圆保持正向。
    outer_r <- 5.12
    group_cex <- if (combined) 1 else 1.35
    chars <- strsplit(as.character(g$category), "")[[1]]
    widths <- strwidth(chars, cex = group_cex) + 0.003
    available <- outer_r * (g$high - g$low - 0.06)
    group_cex <- group_cex * min(1, available/sum(widths))
    widths <- strwidth(chars, cex = group_cex) + 0.003
    flip <- sin(g$middle) < 0
    phi <- g$middle + ifelse(flip, 1, -1) * (cumsum(widths)-widths/2-sum(widths)/2)/outer_r
    for (j in seq_along(chars)) {text(outer_r*cos(phi[j]), outer_r*sin(phi[j]), chars[j], srt = phi[j]*180/pi + ifelse(flip, 90, -90), cex = group_cex, col = color)}
    gap <- sum(widths)/(2*outer_r) + 0.018
    ends <- list(c(g$low, min(g$high, g$middle-gap)), c(max(g$low, g$middle+gap), g$high))
    for (end in ends) if (diff(end) > 0)
      arc(outer_r, end[1], end[2], col = color, lwd = 1.2, lty = "longdash")
  }
  draw_track_labels()

  legend("bottom", inset = c(0, -0.028), xpd = NA, bty = "n", horiz = TRUE, legend = c("Positive", "Negative", "Non-significant"), col = effect_col[c("higher", "lower", "other")], pch = 16, cex = 1.05)
  limit_note <- sprintf("Shared log scale  |  displayed HR range %.3f to %.3f", 1/hr_plot_max, hr_plot_max)
  if (any(d$beyond)) {limit_note <- paste0(limit_note, "  |  boundary triangles: ", sum(d$beyond), " estimates beyond this range")}
  mtext(limit_note, side = 1, line = if (combined) 1.6 else 2.65, cex = if (combined) 0.7 else 0.95, col = "#45525A")
}

# A：人口学统计来自 s1_1_cox_plot.R 导出的 CSV。
draw_demographic <- function() {
  tab <- read.csv(file.path(project, "01cox/demographic_characteristics.csv"), stringsAsFactors = FALSE)
  vars <- unique(tab$characteristic)
  totals <- tapply(tab$n, tab$characteristic, sum)
  stopifnot(length(vars) == 9L, length(unique(totals)) == 1L, all(is.finite(tab$n)), all(tab$n >= 0))
  total <- unname(totals[1])
  palettes <- list(c("#9FC392", "#D1E4C2"), c("#79B3A8", "#A5D8C6"),
    c("#437F91", "#69A3B4", "#ACDCDE", "#D2EAEC"),
    c("#947EB5", "#AC98C5", "#C8B6D7"), c("#A78370", "#CAAA96", "#E4CBB7"),
    c("#835FA5", "#916EAF", "#B493C2", "#D7BEDC"),
    c("#937064", "#CCA08F"), c("#A57087", "#EAA6BE"),
    c("#DF702B", "#FBA043", "#F5C365"))
  par(mar = c(3, 0.2, 1.5, 0.3), family = "sans", xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(c(-total * 0.34, total * 1.01), c(0.35, 9.8))
  for (i in seq_along(vars)) {
    z <- tab[tab$characteristic == vars[i], ]
    y <- 10 - i
    colors <- c(palettes[[i]], if ("Missing" %in% z$category) "#CFD1D2")
    stopifnot(length(colors) == nrow(z))
    ends <- cumsum(z$n)
    rect(c(0, head(ends, -1)), y - 0.26, ends, y + 0.26, col = colors, border = NA)
    text(-total * 0.015, y, vars[i], adj = 1, cex = 1)
    labels <- z$category
    title_case <- labels %in% c("low", "middle", "high")
    labels[title_case] <- tools::toTitleCase(labels[title_case])
    # 较长图例分两行显示。
    legend(0, y + 0.72, legend = labels, fill = colors, border = NA, bty = "n", ncol = if (nrow(z) > 3) 3 else nrow(z), cex = 0.8, x.intersp = 0.4, y.intersp = 0.9, xjust = 0, yjust = 1, seg.len = 0.7)
  }
  ticks <- seq(0, total, by = 50000)
  axis(1, at = ticks, labels = format(ticks, big.mark = ",", trim = TRUE), pos = 0.74, cex.axis = 0.85, tck = -0.015)
  text(total/2, 0.15, "Number of participants", xpd = NA, cex = 1)
}

draw_fig1_AB <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 2))
  draw_demographic()
  mtext("A", side = 3, adj = 0, line = 0.4, font = 2, cex = 1.5)
  draw_circular(combined = TRUE)
  mtext("B", side = 3, adj = 0, line = -0.1, font = 2, cex = 1.5)
}
# 导出：单独 B，以及 A4 横版宽度的 A+B（基础字号 8 pt）。
save_png <- function(filename, draw, width, height, pointsize = 12) {
  path <- file.path(output_dir, filename)
  flt_png(path, width, height, pointsize = pointsize)
  tryCatch(draw(), finally = dev.off())
  message("Saved: ", path)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_png("Fig1_B.png", draw_circular, 16, 16)
save_png("Fig1.png", draw_fig1_AB, 297/25.4, 198/25.4, 8)
message("Plotted ", n_items, " items and ", sum(is.finite(d$HR)), " Cox estimates; ", sum(d$beyond), " estimates marked by boundary triangles.")
