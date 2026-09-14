# Render the same SNP axis for every ancestry; omitted LD stays transparent.
ld_plot <- function(result, race) {
  item <- Filter(function(x) identical(x$race, race), result$populations)[[1]]
  if (is.null(item$r2) || !length(result$pos)) return(NULL)
  z <- do.call(rbind, lapply(item$r2, function(row) vapply(row, function(v)
    if (is.null(v)) NA_real_ else as.numeric(v), numeric(1))))
  pos <- unlist(result$pos); snps <- unlist(result$snps)
  text <- outer(snps, snps, function(a, b) paste(a, "×", b))
  fig <- plotly::plot_ly(height = 310, x = pos, y = pos, z = z, text = text, type = "heatmap", zmin = 0, zmax = 1,
    colorscale = list(list(0, "#fff7ed"), list(.25, "#fed7aa"), list(.5, "#fb923c"), list(1, "#b91c1c")),
    colorbar = list(title = "r²", thickness = 12), hoverongaps = FALSE,
    hovertemplate = "%{text}<br>%{x} × %{y}<br>r² = %{z:.3f}<extra></extra>")
  fig <- plotly::layout(fig, margin = list(l = 68, r = 65, t = 16, b = 50),
    xaxis = list(title = paste0("chr", result$chr, " · GRCh", result$build, " (bp)"),
      range = c(result$start, result$end), tickformat = "~s", fixedrange = TRUE),
    yaxis = list(title = "Position (bp)", range = c(result$end, result$start), tickformat = "~s", fixedrange = TRUE),
    paper_bgcolor = "white", plot_bgcolor = "white")
  plotly::config(fig, displaylogo = FALSE, modeBarButtonsToRemove = c("zoom2d", "pan2d", "select2d", "lasso2d", "autoScale2d", "resetScale2d"))
}
