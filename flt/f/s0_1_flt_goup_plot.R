source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 圆形分类树 :使用 item182 中的 182 个 FLT 及其 10 个分类
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(dplyr)

# 只需调整路径、图幅和字号
project <- Sys.getenv("FLT_OUT")
flt_file <-  Sys.getenv("FLT_DICTIONARY")
outdir <- file.path(project, "00data")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fig_size <- 16
leaf_size <- 3.2
group_size <- 5.2

d <- readxl::read_excel(flt_file, sheet = "item182") %>% transmute(variable = name, label = item_raw, category)
groups <- unique(d$category)
cols <- setNames(c("#C54C3B", "#4BA7C1", "#80634F", "#415A82", "#DE897C","#818DB6", "#80BFB3", "#C52C36", "#309B89", "#AA998C"), groups)
slots <- nrow(d) + 0.6 * length(groups)
d <- d %>%
  mutate(category = factor(category, levels = groups)) %>% arrange(category, startsWith(variable, "F-"), tolower(label)) %>%
  mutate(display_label = paste0(toupper(substr(label, 1, 1)), substring(label, 2)), theta = pi/2 - 2*pi*(row_number() + 0.6*(as.integer(category)-1) - 0.5)/slots, x = cos(theta), y = sin(theta), angle = theta*180/pi + ifelse(x < 0, 180, 0))
g <- d %>% group_by(category) %>%
  summarise(n = n(), lo = min(theta)-pi/slots, hi = max(theta)+pi/slots, theta = mean(range(theta)), .groups = "drop") %>%
  mutate(category = as.character(category), x = 0.43*cos(theta), y = 0.43*sin(theta))
text_width <- function(s, size, bold = FALSE) {
  grDevices::pdf(NULL, family = "sans")
  on.exit(grDevices::dev.off())
  graphics::strwidth(s, "inches", cex = size*72.27/25.4/12, font = if (bold) 2 else 1) * 25.4
}
w <- max(text_width(d$display_label, leaf_size, TRUE))
q <- 2*w/(fig_size*25.4 - 14)
stopifnot(q < 1)
lim <- max(1.6, 1.245/(1-q))
mm <- 2*lim/(fig_size*25.4 - 14)
radius <- 1.125 + w*mm
t <- seq(0, 1, length.out = 60)
B <- cbind((1-t)^3, 3*(1-t)^2*t, 3*(1-t)*t^2, t^3)
branches <- bind_rows(lapply(seq_len(nrow(d)), function(i) {
  a <- d[i, ]; b <- g[match(a$category, g$category), ]
  xy <- B %*% rbind(c(b$x, b$y), 0.67*c(cos(b$theta), sin(b$theta)), 0.84*c(a$x, a$y), c(a$x, a$y))
  data.frame(x = xy[,1], y = xy[,2], id = i, category = a$category)
}))
arcs <- labels <- list()
for (i in seq_len(nrow(g))) {
  a <- g[i, ]; chars <- strsplit(a$category, "")[[1]]
  widths <- text_width(chars, group_size) + 0.06
  flip <- sin(a$theta) < 0
  phi <- a$theta + ifelse(flip, 1, -1) * (cumsum(widths)-widths/2-sum(widths)/2)*mm/radius
  labels[[i]] <- data.frame(x = radius*cos(phi), y = radius*sin(phi), label = chars, angle = phi*180/pi + ifelse(flip, 90, -90), category = a$category)
  gap <- sum(widths)*mm/(2*radius) + 0.018
  ends <- list(c(a$lo, min(a$hi, a$theta-gap)), c(max(a$lo, a$theta+gap), a$hi))
  for (j in seq_along(ends)) if (diff(ends[[j]]) > 0) {
    phi <- seq(ends[[j]][1], ends[[j]][2], length.out = 80)
    arcs[[length(arcs)+1L]] <- data.frame(x = radius*cos(phi), y = radius*sin(phi), category = a$category, id = paste(i, j))
  }
}
p <- ggplot() +
  geom_segment(data = g, aes(x = 0, y = 0, xend = x, yend = y, colour = category), linewidth = 0.65, alpha = 0.62, lineend = "round", show.legend = FALSE) +
  geom_path(data = branches, aes(x, y, group = id, colour = category), linewidth = 0.46, alpha = 0.32, lineend = "round", show.legend = FALSE) +
  geom_point(data = d, aes(x, y, colour = category), shape = 16, size = 1.45) +
  geom_text(data = d, aes(1.025*x, 1.025*y, label = display_label, angle = angle, hjust = as.numeric(x < 0), colour = category), size = leaf_size, family = "sans", fontface = "bold", show.legend = FALSE) +
  geom_path(data = bind_rows(arcs), aes(x, y, group = id, colour = category), linewidth = 0.7, linetype = "longdash", show.legend = FALSE) +
  geom_text(data = bind_rows(labels), aes(x, y, label = label, angle = angle, colour = category), size = group_size, family = "sans", show.legend = FALSE) +
  scale_colour_manual("Groups", values = cols, breaks = groups) +
  guides(colour = guide_legend(override.aes = list(size = 3.6))) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE, clip = "off") +
  theme_void(base_family = "sans") +
  theme(legend.position = c(0.5, 0.5), legend.text = element_text(size = 10, face = "bold"), legend.title = element_text(size = 11, face = "bold"), legend.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(7, 7, 7, 7, "mm"))
pdf_file <- file.path(outdir, "FLT182_tree.pdf")
ggsave(pdf_file, p, device = "pdf", width = fig_size, height = fig_size, bg = "white", limitsize = FALSE)
print(g[c("category", "n")])
message("PDF saved: ", pdf_file)
