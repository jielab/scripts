#!/usr/bin/env Rscript
# Same pairwise heatmap idea as 0f/plot.f.R::plot_rg, using installed ggplot2.
# Keep raw estimates (including rg > 1) in labels/tables; saturate only the colors.
suppressPackageStartupMessages(library(ggplot2))
out <- commandArgs(trailingOnly = TRUE)[1]
rg <- read.delim(file.path(out, 'rg.tsv'), check.names = FALSE, na.strings = c('NA', 'nan'))
h2 <- read.delim(file.path(out, 'h2.tsv'), check.names = FALSE)
traits <- h2$trait
m <- expand.grid(p1 = traits, p2 = traits, stringsAsFactors = FALSE)
m$rg <- NA_real_; m$p <- NA_real_
for (i in seq_len(nrow(rg))) {
  hit <- (m$p1 == rg$p1[i] & m$p2 == rg$p2[i]) | (m$p1 == rg$p2[i] & m$p2 == rg$p1[i])
  if (rg$status[i] == 'ESTIMATED') {
    m$rg[hit] <- rg$rg[i]; m$p[hit] <- rg$p[i]
  }
}
available <- unique(c(rg$p1[rg$status == 'ESTIMATED'], rg$p2[rg$status == 'ESTIMATED'],
                      h2$trait[h2$status == 'ESTIMATED']))
diag <- m$p1 == m$p2 & m$p1 %in% available
m$rg[diag] <- 1; m$p[diag] <- 0
m$label <- ifelse(is.finite(m$rg), sprintf('%.2f', m$rg), 'NA')
m$label[is.finite(m$p) & m$p < 0.05 & m$p1 != m$p2] <- paste0(m$label[is.finite(m$p) & m$p < 0.05 & m$p1 != m$p2], '*')
m$p1 <- factor(m$p1, levels = traits)
m$p2 <- factor(m$p2, levels = rev(traits))
p <- ggplot(m, aes(p1, p2, fill = rg)) +
  geom_tile(color = 'white', linewidth = 0.4) + geom_text(aes(label = label), size = 3.3) +
  scale_fill_gradient2(low = '#2166ac', mid = 'white', high = '#b2182b', midpoint = 0,
                       limits = c(-1, 1), oob = scales::squish, na.value = 'grey90', name = 'rg') +
  coord_fixed() + labs(x = NULL, y = NULL, title = 'LDSC genetic correlation',
    subtitle = 'Chromosomes 1–22 · * nominal P < 0.05 · NA: unavailable',
    caption = 'Raw estimates are shown; color scale is capped at ±1. Diagonal = 1 by definition for available traits.') +
  theme_minimal(base_size = 12) + theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = 'bold'))
tmp <- file.path(out, 'rg.tmp.png')
ggsave(tmp, p, width = max(8, 3 + 0.65 * length(traits)), height = max(7, 2 + 0.65 * length(traits)), dpi = 180, bg = 'white')
if (!file.rename(tmp, file.path(out, 'rg.png'))) stop('Could not publish rg.png')
