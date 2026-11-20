library(ggplot2)
library(patchwork)

project <- "D:/projects/flt"
outdir <- file.path(project, "analysis/01cox")
figdir <- file.path(outdir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 热图
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
result_file <- file.path(outdir, "cox_results.rds")
result <- readRDS(result_file)
cox_table <- result$table
meta <- result$meta
Ys <- result$Ys
out_map <- result$out_map
groups <- unique(meta$category)
meta <- meta[order(match(meta$category, groups), seq_len(nrow(meta))), ]

for (analysis_type in c("continuous", "medium_vs_low", "high_vs_low")) {
  z <- cox_table[cox_table$Analysis == analysis_type, ]
  res <- data.frame(variable = factor(z$name, levels = meta$name), category = factor(z$Category, levels = groups), outcome = factor(z$Disease, levels = rev(Ys)), beta = z$Beta, HR = z$HR, FDR = z$FDR, status = z$Status)
  res$valid <- res$status == "ok" & is.finite(res$beta) & is.finite(res$HR) & res$HR > 0 & is.finite(res$FDR)
## 画热图-----------------------------
# 警告模型已记录在结果表中，不纳入 FDR 或气泡图的有效估计
if (!any(res$valid)) { warning(paste("No valid estimates:", analysis_type)); next }
# 设置颜色：方向由 HR 决定，深浅由 FDR 决定 
bins <- c("<0.0001", "<0.001", "<0.01", "<0.05")
bin <- cut(res$FDR, c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), right = FALSE, labels = c(bins, ">=0.05"))
res$key <- ifelse(res$FDR >= 0.05, "Not significant", paste(ifelse(res$HR >= 1, "Higher risk", "Lower risk"), bin))
colors <- c(setNames(c("#8B1420", "#CB3543", "#E66F73", "#F3B4AF"), paste("Higher risk", bins)),setNames(c("#25345E", "#386AA4", "#77A5CD", "#BED3E5"), paste("Lower risk", bins)), "Not significant" = "#F4EEE6")
panel <- function(gs) {
  z <- res[res$category %in% gs, ]
  ggplot(z, aes(variable, outcome)) +
    geom_point(data = z[z$valid, ], aes(size = abs(beta), fill = key), shape = 21, color = "#8B8B8B", stroke = 0.3) +
    geom_point(data = z[!z$valid, ], shape = 4, size = 2.2, color = "grey45") +
    facet_grid(. ~ category, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = setNames(meta$Food_liking, meta$name)) +
    scale_y_discrete(drop = FALSE, labels = out_map) +
    scale_fill_manual(values = colors, limits = names(colors), drop = FALSE, name = "Direction / BH-FDR") +
    scale_size_continuous(range = c(2.2, 6.2), limits = c(0, max(abs(res$beta[res$valid]))), name = "|log(HR)|") +
    guides(fill = guide_legend(override.aes = list(size = 4.5)), size = guide_legend(order = 1)) +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1, size = 10, color = "black"), axis.text.y = element_text(face = "bold", color = "black"), strip.background = element_rect(fill = "#34495E"), strip.text = element_text(color = "white", face = "bold"))
}
group_rows <- split(groups, ceiling(seq_along(groups)/5))
p <- wrap_plots(lapply(group_rows, panel), ncol=1, guides="collect") +
	plot_annotation(
		title = "Food-liking traits and incident brain diseases",
		subtitle = paste(analysis_type, "| item182 | BH-FDR within each outcome (182 tests)"),
		caption = "Circle size indicates absolute log(HR), with a minimum size for visibility; color indicates direction and BH-FDR. Crosses indicate unavailable or flagged models."
	)
p <- p & theme(legend.position = "bottom", legend.key.height = grid::unit(0.7, "cm"))
ggsave(file.path(figdir, paste0("Fig_FLT_brain_", analysis_type, ".pdf")), p, width = 27, height = max(10, length(group_rows) * (3 + 0.32 * length(Ys))), bg = "white")
}
