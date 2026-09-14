# One shared x axis keeps every population aligned while each has its own y axis.
bplot_locus <- function(chr, region) {
  if (chr == "All") return("all")
  paste0("chr", chr, ":", format(region[1], scientific = FALSE, trim = TRUE), "-",
    format(region[2], scientific = FALSE, trim = TRUE))
}

bplot_plot <- function(views, build, chr, region, threshold, page = "overview", selected = NULL) {
  n <- length(views)
  height <- max(350L, n * 175L + 65L)
  fig <- plotly::plot_ly(height = height)
  axes <- list(); annotations <- list(); shapes <- list()
  colors <- c("#2563a6", "#0f8a83")
  for (i in seq_len(n)) {
    view <- views[[i]]; data <- view$data
    axis <- if (i == 1L) "y" else paste0("y", i)
    key <- if (i == 1L) "yaxis" else paste0("yaxis", i)
    top <- 1 - (i - 1) / n; bottom <- top - 0.79 / n
    ymax <- max(10, if (nrow(data)) max(data$LOGP) * 1.08 else 0, -log10(threshold) * 1.1)
    axes[[key]] <- list(domain = c(bottom, top), range = c(0, ymax), fixedrange = TRUE,
      title = list(text = "−log₁₀(P)", font = list(size = 11)), tickfont = list(size = 10),
      gridcolor = "#eaf0f6", zeroline = FALSE, anchor = "x")
    if (length(view$blocks$edges)) {
      edges <- view$blocks$edges
      fig <- plotly::add_trace(fig, x = as.vector(rbind(edges, edges, NA_real_)),
        y = rep(c(0, ymax, NA_real_), length(edges)), type = "scatter", mode = "lines",
        yaxis = axis, xaxis = "x", line = list(color = "rgba(164,111,45,0.35)", width = 1),
        hoverinfo = "skip", showlegend = FALSE, name = paste0(view$track$race, " boundaries"))
    }
    if (nrow(data)) {
      # Explicit arrays avoid R formula/data inheritance across the stacked traces.
      fig <- plotly::add_trace(fig, x = data$X, y = data$LOGP, type = "scattergl", mode = "markers",
        xaxis = "x", yaxis = axis, name = view$track$label, showlegend = FALSE,
        marker = list(size = 3, opacity = 0.72,
          color = if (chr == "All") colors[(match(data$CHR, c(as.character(1:22), "X")) %% 2) + 1L] else colors[1]),
        customdata = as.matrix(data[, .(SNP, CHR, POS, P, BLOCK, ID)]),
        hovertemplate = "%{customdata[0]}<br>chr%{customdata[1]}:%{customdata[2]}<br>P = %{customdata[3]}<br>%{customdata[4]}<extra>%{fullData.name}</extra>")
    } else {
      fig <- plotly::add_trace(fig, x = numeric(), y = numeric(), type = "scatter", mode = "markers",
        xaxis = "x", yaxis = axis, showlegend = FALSE, hoverinfo = "skip", name = view$track$label)
    }
    shapes[[i]] <- list(type = "line", xref = "x", yref = axis,
      x0 = region[1], x1 = region[2], y0 = -log10(threshold), y1 = -log10(threshold),
      line = list(color = "#d86569", width = 1, dash = "dot"))
    label <- htmltools::htmlEscape(view$track$label)
    note <- htmltools::htmlEscape(view$blocks$note)
    annotations[[i]] <- list(xref = "paper", yref = "paper", x = 0, y = top,
      xanchor = "left", yanchor = "bottom", yshift = 5, showarrow = FALSE,
      text = paste0("<b>", label, "</b>  ·  ", note, "  ·  ",
        format(nrow(data), big.mark = ",", scientific = FALSE), " shown / ",
        format(view$total, big.mark = ",", scientific = FALSE), " in input region"),
      font = list(size = 12, color = "#334155"))
  }
  if (!is.null(selected) && chr != "All") shapes[[length(shapes) + 1L]] <- list(
    type = "rect", xref = "x", yref = "paper", x0 = selected[1], x1 = selected[2], y0 = 0, y1 = 1,
    fillcolor = "rgba(234,179,8,0.08)", line = list(width = 0), layer = "below")
  lengths <- bplot_lengths(build)
  offsets <- c(0, head(cumsum(lengths), -1))
  xaxis <- list(range = region, anchor = if (n == 1L) "y" else paste0("y", n),
    showgrid = FALSE, zeroline = FALSE, title = if (chr == "All") "Chromosome" else paste0("Chromosome ", chr, " · position (bp)"),
    tickfont = list(size = 11), fixedrange = FALSE)
  if (chr == "All") {
    xaxis$tickmode <- "array"; xaxis$tickvals <- offsets + lengths / 2; xaxis$ticktext <- names(lengths)
  } else xaxis$tickformat <- "~s"
  fig <- do.call(plotly::layout, c(list(p = fig, xaxis = xaxis, annotations = annotations,
    shapes = shapes, margin = list(l = 72, r = 25, t = 36, b = 60),
    paper_bgcolor = "#ffffff", plot_bgcolor = "#ffffff", dragmode = "zoom",
    hovermode = "closest", showlegend = FALSE), axes))
  fig <- plotly::config(fig, scrollZoom = TRUE, displaylogo = FALSE, doubleClick = FALSE,
    modeBarButtonsToRemove = c("select2d", "lasso2d", "autoScale2d", "resetScale2d"),
    toImageButtonOptions = list(format = "png", filename = paste0("bplot.GRCh", build, ".", chr), height = height, width = 1600))
  htmlwidgets::onRender(fig, "function(el, x, data) { window.bplot.bind(el, data); }",
    data = list(build = build, chr = chr, page = page,
      offsets = as.list(stats::setNames(offsets, names(lengths))),
      tracks = lapply(views, function(v) list(race = v$track$race,
        blocks = lapply(seq_len(nrow(v$blocks$data)), function(i) unname(as.list(v$blocks$data[i, .(CHR, START, END, ID)])))))))
}
