#!/usr/bin/env Rscript
# Run from the project root: Rscript --vanilla tests/bplot_hover.R
library(data.table)
source("shiny/data.R")
source("shiny/plot.R")

blocks <- data.table(CHR = "1", START = 100, END = 200, ID = 123L)
data <- data.table(SNP = c("rs1", "rs2"), CHR = "1", POS = c(150, 250),
  P = c(1e-9, 0.1), LOGP = c(9, 1), X = c(150, 250))
data[, ID := bplot_block_ids(data, blocks)]
data[, BLOCK := ifelse(is.na(ID), "", paste("block: AFR", ID))]

# Check the serialized trace consumed by the browser, not just plot_ly attrs.
# Include single-point regions and multiple tracks to catch array simplification
# and cross-track inheritance; the second SNP has no BED block.
for (page in c("overview", "blockview")) for (rows in list(1L, 1:2)) {
  view <- list(data = data[rows], total = length(rows),
    track = list(label = "height.AFR", race = "AFR"),
    blocks = list(data = blocks, edges = c(100, 200), note = "AFR blocks"))
  fig <- bplot_plot(list(view, view), "37", "1", c(100, 300), 5e-8, page)
  payload <- jsonlite::fromJSON(plotly::plotly_json(fig, jsonedit = FALSE),
    simplifyVector = FALSE)
  traces <- Filter(function(x) identical(x$type, "scattergl"), payload$data)
  stopifnot(length(traces) == 2L)
  for (trace in traces) {
    stopifnot(length(trace$customdata) == length(rows),
      identical(trace$customdata[[1]][[1]], "rs1"),
      identical(trace$customdata[[1]][[2]], "1"),
      trace$customdata[[1]][[3]] == 150,
      trace$customdata[[1]][[4]] == 1e-9,
      identical(trace$customdata[[1]][[5]], "block: AFR 123"),
      trace$customdata[[1]][[6]] == 123)
    if (length(rows) == 2L) stopifnot(
      identical(trace$customdata[[2]][[1]], "rs2"),
      identical(trace$customdata[[2]][[5]], ""),
      is.null(trace$customdata[[2]][[6]]))
  }
}
cat("Block hover and click payload checks passed.\n")

# Exercise six full tracks. A tiny fixture cannot expose the recursive-list
# slowdown that previously kept Shiny busy and left the panel blank.
large <- data[rep(1:2, 6000L)]
views <- lapply(seq_len(6L), function(i) list(data = large, total = nrow(large),
  track = list(label = paste0("track", i), race = "AFR"),
  blocks = list(data = blocks, edges = numeric(), note = "AFR blocks")))
elapsed <- system.time({
  fig <- bplot_plot(views, "37", "All", c(1, sum(bplot_lengths("37"))), 5e-8)
  payload <- jsonlite::fromJSON(plotly::plotly_json(fig, jsonedit = FALSE), simplifyVector = FALSE)
})[["elapsed"]]
stopifnot(length(payload$data) == 6L,
  all(vapply(payload$data, function(trace) length(trace$customdata) == 12000L, logical(1))))
cat(sprintf("Six-track / 72,000-point serialization completed in %.1f seconds.\n", elapsed))
