# Base R manuscript summaries. Run within the staged ranked report directory.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
setwd(args[1L])
options(scipen = 999, digits = 12)
read_tsv <- function(stem) {
    read.delim(paste0(stem, ".tsv"), colClasses = "character", na.strings = NULL,
               comment.char = "", check.names = FALSE, fileEncoding = "UTF-8")
}
write_tsv <- function(x, stem) {
    write.table(x, paste0(stem, ".tsv"), sep = "\t", quote = FALSE,
                row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
ratio <- function(count, total, scale = 1) {
    ifelse(total > 0, scale * count / total, NA_real_)
}
s <- read_tsv("samples")
g <- read_tsv("transitions")
ev <- read_tsv("target_species_evidence")
stopifnot(nrow(s) > 0, !anyDuplicated(s$sample))
s$dataset <- sub("/[^/]+$", "", s$sample)
datasets <- unique(c(args[-1L], sort(unique(s$dataset))))
for (k in c("pairs", "species_pairs", "higher_rank_pairs", "unclassified_pairs",
            "human_pairs", "plant_pairs", "virus_pairs")) s[[k]] <- as.numeric(s[[k]])
g$pairs <- as.numeric(g$pairs)
ev$pairs <- as.numeric(ev$pairs)
stopifnot(all(s$pairs == s$species_pairs + s$higher_rank_pairs + s$unclassified_pairs))
s$species_pct <- ratio(s$species_pairs, s$pairs, 100)
groups <- c("human", "plant", "virus")

# Fill unobserved sample/group combinations with zero before computing denominators.
sample_targets <- do.call(rbind, lapply(groups, function(group) {
    selected <- g[g$kraken_group == group, , drop = FALSE]
    initial <- tapply(selected$pairs, selected$sample, sum)[s$sample]
    e <- ev[ev$group == group, , drop = FALSE]
    exact <- tapply(e$pairs, e$sample, sum)[s$sample]
    data.frame(dataset = s$dataset, sample = s$sample, group = group, input_pairs = s$pairs,
               initial_pairs = as.numeric(replace(initial, is.na(initial), 0)),
               final_pairs = s[[paste0(group, "_pairs")]],
               exact_pairs = as.numeric(replace(exact, is.na(exact), 0)),
               stringsAsFactors = FALSE)
}))
stopifnot(all(sample_targets$exact_pairs <= sample_targets$final_pairs))
write_tsv(sample_targets, "sample_targets")

# A separate scope column permits a dataset literally named "All".
overview <- targets <- plant_rows <- list()
plants <- ev[ev$group == "plant" & ev$pairs > 0, , drop = FALSE]
plant_ids <- unique(plants$taxid)
for (i in seq_len(length(datasets) + 1L)) {
    pooled <- i > length(datasets)
    dataset <- if (pooled) "" else datasets[i]
    key <- data.frame(scope = if (pooled) "all" else "dataset", dataset = dataset)
    d <- if (pooled) s else s[s$dataset == dataset, , drop = FALSE]
    usable <- d$sample[d$pairs > 0]
    n <- sum(d$pairs)
    sp <- sum(d$species_pairs)
    hi <- sum(d$higher_rank_pairs)
    un <- sum(d$unclassified_pairs)
    qq <- as.numeric(quantile(d$species_pct, c(.25, .5, .75), na.rm = TRUE, names = FALSE))
    overview[[i]] <- data.frame(key, samples = nrow(d), nonempty_samples = length(usable),
        input_pairs = n, species_pairs = sp, higher_rank_pairs = hi, unclassified_pairs = un,
        species_pct = ratio(sp, n, 100), higher_rank_pct = ratio(hi, n, 100),
        unclassified_pct = ratio(un, n, 100), sample_species_q25_pct = qq[1],
        sample_species_median_pct = qq[2], sample_species_q75_pct = qq[3])
    for (group in groups) {
        t <- sample_targets[sample_targets$sample %in% d$sample & sample_targets$group == group, ]
        init <- sum(t$initial_pairs)
        final <- sum(t$final_pairs)
        exact <- sum(t$exact_pairs)
        targets[[length(targets) + 1L]] <- data.frame(key, group = group, input_pairs = n,
            nonempty_samples = length(usable), initial_pairs = init, final_pairs = final,
            exact_pairs = exact, higher_rank_pairs = final - exact,
            initial_ppm = ratio(init, n, 1e6), final_ppm = ratio(final, n, 1e6),
            exact_ppm = ratio(exact, n, 1e6), exact_pct_of_final = ratio(exact, final, 100),
            samples_with_final = sum(t$final_pairs > 0), samples_with_exact = sum(t$exact_pairs > 0))
    }
    if (length(plant_ids)) {
        # One sample can supply several evidence types for the same taxid.
        dplants <- plants[plants$sample %in% usable, , drop = FALSE]
        for (taxid in plant_ids) {
            hits <- dplants[dplants$taxid == taxid, , drop = FALSE]
            positive <- length(unique(hits$sample))
            count <- sum(hits$pairs)
            plant_rows[[length(plant_rows) + 1L]] <- data.frame(key, taxid = taxid,
                name = plants$name[match(taxid, plants$taxid)], exact_pairs = count,
                input_pairs = n, exact_ppm = ratio(count, n, 1e6),
                positive_samples = positive, nonempty_samples = length(usable),
                detected_pct = ratio(positive, length(usable), 100))
        }
    }
}
overview <- do.call(rbind, overview)
targets <- do.call(rbind, targets)
plant_table <- if (length(plant_rows)) do.call(rbind, plant_rows) else
    data.frame(scope = character(), dataset = character(), taxid = character(),
               name = character(), exact_pairs = numeric(), input_pairs = numeric(), exact_ppm = numeric(),
               positive_samples = integer(), nonempty_samples = integer(), detected_pct = numeric())
write_tsv(overview, "table1_classification")
write_tsv(targets, "table2_target_groups")
write_tsv(plant_table, "table3_plant_species")

ink <- "#1C3342"
palette_rank <- c("#127C83", "#E5AB54", "#DCE3E6")
palette_target <- c("#516D87", "#168274", "#A24938")
co <- overview[overview$scope == "dataset", , drop = FALSE]
nd <- nrow(co)
yy <- rev(seq_len(nd))
ylabels <- paste0(co$dataset, "  (n=", co$samples, ")")
labels_y <- function() axis(2, at = yy, labels = ylabels, las = 1, tick = FALSE,
                             cex.axis = .85, line = -.2)
subtitle <- function(text) mtext(text, side = 3, line = .25, adj = 0, cex = .82, col = "#435967")
panel_title <- function(text) title(main = text, adj = 0, line = 1.5, cex.main = 1.1)
canvas <- function(xlim, xlab) {
    plot(NA_real_, NA_real_, xlim = xlim, ylim = c(.45, nd + .55), axes = FALSE,
         xlab = xlab, ylab = "", xaxs = "i", yaxs = "i")
}
save_figure <- function(stem, width, height, draw) {
    for (ext in c("pdf", "png")) {
        file <- paste0(stem, ".", ext)
        if (ext == "pdf") pdf(file, width = width, height = height, useDingbats = FALSE)
        else png(file, width = width, height = height, units = "in", res = 240, bg = "white")
        tryCatch({
            par(family = "sans", col = ink, col.axis = ink, col.lab = ink, col.main = ink,
                bty = "n", mgp = c(2.2, .6, 0), tcl = -.25)
            draw()
        }, finally = dev.off())
    }
}
label_margin <- max(7, min(22, max(nchar(ylabels)) * .52 + 1))
save_figure("fig1_classification", 8.2, max(5.4, 3.7 + .53 * nd), function() {
    par(mfrow = c(2, 1), mar = c(3.3, label_margin, 3.2, 1.0))
    canvas(c(0, 100), "Percent of all input read pairs")
    panel_title("A   Taxonomic resolution by dataset")
    for (i in seq_len(nd)) {
        if (!co$input_pairs[i]) next
        edges <- c(0, cumsum(as.numeric(co[i, c("species_pct", "higher_rank_pct", "unclassified_pct")])))
        for (j in 1:3) rect(edges[j], yy[i]-.21, edges[j+1], yy[i]+.21, col=palette_rank[j], border=NA)
    }
    axis(1, at = seq(0, 100, 20))
    labels_y()
    legend("topright", inset = c(0, -.30), xpd=NA, bty="n", horiz=TRUE, cex=.78,
           legend=c("Species", "Higher rank", "Unclassified"), fill=palette_rank, border=NA)
    values <- s$species_pct[is.finite(s$species_pct)]
    xmax <- max(c(1, values))
    canvas(c(-.015*xmax, 1.04*xmax), "Percent of each sample's input read pairs assigned to species")
    panel_title("B   Species-level assignment across samples")
    subtitle("Each point is one nonempty sample; black ticks mark dataset medians.")
    axis(1); labels_y()
    for (i in seq_len(nd)) {
        vals <- sort(na.omit(s$species_pct[s$dataset == co$dataset[i]]))
        if (!length(vals)) next
        offsets <- if (length(vals) == 1L) 0 else (((seq_along(vals)*.618034) %% 1)-.5)*.50
        points(vals, yy[i]+offsets, pch=16, cex=.78, col=adjustcolor(palette_rank[1], .65))
        segments(median(vals), yy[i]-.24, median(vals), yy[i]+.24, lwd=2.5, col=ink)
    }
})
save_figure("fig2_target_groups", 8.2, max(6.0, 3.6 + .9*nd), function() {
    par(mfrow=c(3,1), mar=c(2.5, label_margin, 2.4, 1), oma=c(1.4,0,1.2,0))
    cols <- c("initial_ppm", "final_ppm", "exact_ppm")
    actual <- targets[targets$scope == "dataset", ]
    top <- max(c(1, as.matrix(actual[cols])), na.rm=TRUE)
    ticks <- c(0, 10^seq(0, ceiling(log10(top))))
    xmax <- max(log10(1+top), .32)
    ticks <- ticks[log10(1+ticks) <= 1.04*xmax]
    for (j in seq_along(groups)) {
        gname <- groups[j]
        dat <- actual[actual$group==gname, ]
        canvas(c(-.02*xmax,1.05*xmax), "")
        title(main=c("A   Human","B   Plant","C   Virus")[j], adj=0, line=.65, cex.main=1.1)
        labels_y()
        axis(1, at=log10(1+ticks), labels=format(ticks, trim=TRUE, scientific=FALSE, big.mark=","))
        for (k in seq_along(cols)) points(log10(1+dat[[cols[k]]]), yy+(2-k)*.20,
            col=palette_target[k], pch=c(1,16,17)[k], cex=.92)
    }
    mtext("Read pairs per million input pairs (log10[1 + rate] axis)", side=1, outer=TRUE, line=.2, cex=.8)
    par(fig=c(0,1,0,1),new=TRUE,mar=rep(0,4)); plot.new()
    legend("topright", legend=c("Initial group","Final group","Exact species"), pch=c(1,16,17),
           col=palette_target, horiz=TRUE, bty="n", cex=.8)
})
top_plants <- plant_table[plant_table$scope == "all", , drop=FALSE]
top_plants <- head(top_plants[order(-top_plants$positive_samples, -top_plants$exact_pairs,
                                    as.numeric(top_plants$taxid)), , drop=FALSE], 12)
save_figure("fig3_plant_detection", max(7.5, 3.2 + .95*nd), max(3.2, 2.2 + .37*nrow(top_plants)), function() {
    if (!nrow(top_plants)) {
        par(mar=c(1,1,2.8,1)); plot.new()
        panel_title("Plant species detection")
        text(.5,.48,"No read pairs met the species-level plant criteria.", col=ink)
        return(invisible(NULL))
    }
    np <- nrow(top_plants)
    par(mar=c(max(4.2, min(12, max(nchar(co$dataset))*.45+1.2)),
              max(9, min(24, max(nchar(top_plants$name))*.50+1)), 3.8, 1.2))
    plot(NA, NA, xlim=c(.5,nd+.5), ylim=c(.5,np+.5), axes=FALSE, xlab="",ylab="",xaxs="i",yaxs="i")
    panel_title("Plant species detection by dataset")
    subtitle("Positive / nonempty samples; darker cells indicate higher frequency.")
    colors <- colorRampPalette(c("#EEF4F3", "#79B7AF", "#006A66"))(101)
    details <- plant_table[plant_table$scope == "dataset", ]
    for (i in seq_len(np)) {
        taxon <- details[details$taxid == top_plants$taxid[i], ]
        taxon <- taxon[match(co$dataset, taxon$dataset), ]
        y <- np-i+1
        for (j in seq_len(nd)) {
            pct <- taxon$detected_pct[j]
            shade <- if (is.na(pct)) "#DDDDDD" else colors[1+round(pct)]
            rect(j-.49,y-.48,j+.49,y+.48, col=shade,border=NA)
            text(j, y, if (is.na(pct)) "NA" else paste0(taxon$positive_samples[j], "/",taxon$nonempty_samples[j]),
                 cex=.84, col=if (is.finite(pct) && pct>65) "white" else ink)
        }
    }
    axis(2,at=rev(seq_len(np)),labels=top_plants$name,las=1,font=3,tick=FALSE,cex.axis=.85)
    axis(1,at=seq_len(nd),labels=co$dataset,las=if(nd>5) 2 else 1,tick=FALSE,cex.axis=.85)
})
