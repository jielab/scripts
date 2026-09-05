#!/usr/bin/env Rscript

# Render completed PhyML Newick trees as static, publication-friendly overview
# PNGs.  Large trees retain every branch but omit thousands of overlapping tip
# labels; archaic tips remain marked and named.

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) return(NULL)
    args[[i + 1L]]
}

out_arg <- value_after("--out")
if (is.null(out_arg)) stop("Usage: Rscript tree_plot.R --out PHYML_OUTPUT_DIR", call. = FALSE)
out <- normalizePath(out_arg, winslash = "/", mustWork = FALSE)
loci_root <- file.path(out, "loci")

flat <- file.path(loci_root, "haplotypes.phy")
phy_inputs <- if (file.exists(flat)) {
    flat
} else if (dir.exists(loci_root)) {
    list.files(loci_root, pattern = "^haplotypes\\.phy$", recursive = TRUE, full.names = TRUE)
} else {
    character()
}
phy_inputs <- sort(normalizePath(phy_inputs, winslash = "/", mustWork = FALSE))

positive_int <- function(value, default, minimum, maximum) {
    parsed <- suppressWarnings(as.integer(value))
    if (is.na(parsed) || parsed < minimum || parsed > maximum) default else parsed
}

plot_px <- positive_int(Sys.getenv("PHYML_TREE_PLOT_PX", "4800"), 4800L, 1600L, 12000L)
plot_dpi <- positive_int(Sys.getenv("PHYML_TREE_PLOT_DPI", "300"), 300L, 72L, 600L)
label_limit <- positive_int(Sys.getenv("PHYML_TREE_LABEL_LIMIT", "500"), 500L, 0L, 5000L)

read_label_map <- function(path) {
    if (!file.exists(path) || file.info(path)$size <= 0) return(NULL)
    x <- read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("phy_label", "label") %in% names(x))) return(NULL)
    x
}

truthy <- function(x) toupper(trimws(as.character(x))) %in% c("1", "TRUE", "T", "YES", "Y")

descendant_tips <- function(tree, node) {
    tips <- integer()
    frontier <- as.integer(node)
    while (length(frontier)) {
        children <- tree$edge[tree$edge[, 1] %in% frontier, 2]
        tips <- c(tips, children[children <= Ntip(tree)])
        frontier <- unique(children[children > Ntip(tree)])
    }
    unique(tips)
}

archaic_lineage <- function(labels) {
    out <- rep("", length(labels))
    out[grepl("Altai|Chagyr|Vindija|Neander", labels, ignore.case = TRUE)] <- "Neanderthal"
    out[grepl("Denis", labels, ignore.case = TRUE)] <- "Denisovan"
    out
}

candidate_edge_for_tree <- function(tree) {
    n_tip <- Ntip(tree)
    if (n_tip < 3L) return(NULL)
    internal_children <- unique(tree$edge[tree$edge[, 2] > n_tip, 2])
    best <- NULL
    key_is_better <- function(left, right) {
        if (is.null(right)) return(TRUE)
        for (i in seq_along(left)) {
            if (left[[i]] < right[[i]]) return(TRUE)
            if (left[[i]] > right[[i]]) return(FALSE)
        }
        FALSE
    }
    for (node in internal_children) {
        down <- descendant_tips(tree, node)
        support <- if (length(tree$node.label)) {
            suppressWarnings(as.numeric(tree$node.label[[node - n_tip]]))
        } else NA_real_
        support_rank <- if (is.finite(support)) support else -1
        sides <- list(descendant = down, complement = setdiff(seq_len(n_tip), down))
        for (side_name in names(sides)) {
            side <- sides[[side_name]]
            side_labels <- tree$tip.label[side]
            if ("Ancestral" %in% tree$tip.label && "Ancestral" %in% side_labels) next
            lineages <- archaic_lineage(side_labels)
            archaic <- side_labels[nzchar(lineages)]
            modern <- setdiff(side_labels, archaic)
            lineage <- unique(lineages[nzchar(lineages)])
            if (length(lineage) != 1L || length(modern) < 2L || length(side) >= n_tip) next
            # Match tree_summary.py: highest support, fewer modern tips,
            # then more references from the same archaic lineage.
            key <- c(-support_rank, length(modern), -length(archaic),
                     if (identical(side_name, "complement")) 0 else 1, node)
            if (key_is_better(key, if (is.null(best)) NULL else best$key)) {
                best <- list(key = key, tips = side_labels, support = support, side = side_name,
                             lineage = lineage[[1]], pass = is.finite(support) && support >= 70)
            }
        }
    }
    best
}

redraw_phylogram_edges <- function(tree, plot_state, color, width) {
    edge <- tree$edge
    for (i in seq_len(nrow(edge))) {
        parent <- edge[i, 1]; child <- edge[i, 2]
        segments(plot_state$xx[parent], plot_state$yy[child],
                 plot_state$xx[child], plot_state$yy[child], col = color, lwd = width)
    }
    for (parent in unique(edge[, 1])) {
        children <- edge[edge[, 1] == parent, 2]
        segments(plot_state$xx[parent], min(plot_state$yy[children]),
                 plot_state$xx[parent], max(plot_state$yy[children]), col = color, lwd = width)
    }
}

locus_id_for_phy <- function(phy) {
    if (dirname(phy) != loci_root) return(basename(dirname(phy)))
    loci_file <- file.path(out, "final", "loci.tsv")
    if (!file.exists(loci_file)) return(basename(out))
    loci <- read.delim(loci_file, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
    ids <- unique(as.character(loci$locus_id[!is.na(loci$locus_id) & nzchar(loci$locus_id)]))
    if (length(ids) == 1L) ids[[1]] else basename(out)
}

render_one <- function(phy) {
    tree_file <- paste0(phy, "_phyml_tree.txt")
    meta_file <- paste0(phy, ".meta.tsv")
    if (!file.exists(tree_file) || file.info(tree_file)$size <= 0 ||
        !file.exists(meta_file) || file.info(meta_file)$size <= 0) return(NULL)

    png_file <- sub("\\.txt$", ".png", tree_file)
    if (file.exists(png_file) && file.info(png_file)$size > 0) {
        message("PHYML TREE PLOT: status=exists output=", png_file)
        return(png_file)
    }

    tr <- read.tree(tree_file)
    if (is.null(tr) || Ntip(tr) < 2L) stop("No drawable tree in ", tree_file)
    labels <- read_label_map(meta_file)
    if (!is.null(labels)) {
        idx <- match(tr$tip.label, labels$phy_label)
        replace <- !is.na(idx) & !is.na(labels$label[idx]) & nzchar(labels$label[idx])
        tr$tip.label[replace] <- labels$label[idx[replace]]
    }
    tr$tip.label <- make.unique(trimws(tr$tip.label))
    candidate_edge <- candidate_edge_for_tree(tr)
    candidate_side <- if (is.null(candidate_edge) || !candidate_edge$pass) character() else candidate_edge$tips
    outside <- setdiff(tr$tip.label, candidate_side)
    if (length(candidate_side) >= 2L && length(outside)) {
        outgroup <- if ("Ancestral" %in% outside) "Ancestral" else outside[[1]]
        tr <- tryCatch(root(unroot(tr), outgroup = outgroup, resolve.root = FALSE), error = function(e) tr)
    }
    tr <- ladderize(tr)

    n_tip <- Ntip(tr)
    show_all_labels <- n_tip <= label_limit
    archaic <- grep("Altai|Chagyr|Vindija|Denisova|Neander|archaic", tr$tip.label,
                    ignore.case = TRUE)
    candidate_idx <- match(candidate_side, tr$tip.label)
    candidate_idx <- candidate_idx[!is.na(candidate_idx)]
    candidate_side_idx <- match(candidate_side, tr$tip.label)
    candidate_side_idx <- candidate_side_idx[!is.na(candidate_side_idx)]
    candidate_node <- if (length(candidate_side_idx) >= 2L) getMRCA(tr, candidate_side_idx) else NA_integer_
    candidate_desc <- if (is.finite(candidate_node)) descendant_tips(tr, candidate_node) else integer()
    locus_id <- locus_id_for_phy(phy)
    tmp <- tempfile(pattern = paste0(".", basename(png_file), "."),
                    tmpdir = dirname(png_file), fileext = ".png")
    device_open <- FALSE
    tryCatch({
        bitmap_type <- if (capabilities("cairo")) "cairo" else getOption("bitmapType", "Xlib")
        png(tmp, width = plot_px, height = plot_px, res = plot_dpi, type = bitmap_type)
        device_open <- TRUE
        par(mar = c(1.1, 0.8, 2.6, 8.0), xpd = NA, fg = "black", col = "black")
        tip_cex <- if (show_all_labels) max(0.10, min(0.65, 35 / sqrt(n_tip))) else 0.1
        tip_col <- rep("#5f6b73", n_tip)
        tip_col[candidate_idx] <- "#b8322a"
        tip_col[archaic] <- "#111111"
        plot.phylo(
            tr, type = "phylogram", direction = "rightwards", use.edge.length = TRUE,
            show.tip.label = show_all_labels, cex = tip_cex,
            tip.color = tip_col, edge.color = "#42484d", edge.width = if (n_tip > 1000L) 0.25 else 0.55,
            no.margin = FALSE
        )

        # Nature 2020 Fig. 2-style emphasis: shade the smallest displayed
        # clade containing the passing modern haplotypes and their best archaic
        # reference, then redraw the topology above the translucent band.
        if (length(candidate_desc) && is.finite(candidate_node)) {
            pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
            y0 <- min(pp$yy[candidate_desc]) - 0.42
            y1 <- max(pp$yy[candidate_desc]) + 0.42
            rect(pp$xx[candidate_node], y0, max(pp$xx, na.rm = TRUE) * 1.015, y1,
                 col = adjustcolor("#ef9a9a", alpha.f = 0.28), border = NA)
            redraw_phylogram_edges(tr, pp, "#42484d", if (n_tip > 1000L) 0.25 else 0.55)
        }

        tiplabels(pch = 16, cex = if (n_tip > 1000L) 0.10 else 0.22,
                  col = tip_col, frame = "none")
        if (length(archaic)) {
            tiplabels(tip = archaic, pch = 21, bg = "white", col = "black",
                      cex = 0.75, lwd = 0.7, frame = "none")
            if (!show_all_labels) {
                legend("bottomleft", legend = paste("Archaic:", paste(tr$tip.label[archaic], collapse = ", ")),
                       bty = "n", cex = 0.65, text.col = "black")
            }
        }
        if (!show_all_labels && length(candidate_idx)) {
            tiplabels(text = tr$tip.label[candidate_idx], tip = candidate_idx,
                      frame = "none", adj = c(-0.05, 0.5), cex = 0.55, col = "#b8322a")
        }

        # Bootstrap labels are useful on small trees but become illegible on a
        # thousands-tip overview.  Show only support >=70 when space permits.
        if (n_tip <= 250L && length(tr$node.label)) {
            support <- suppressWarnings(as.numeric(tr$node.label))
            keep <- which(is.finite(support) & support >= 70)
            if (length(keep)) {
                nodelabels(support[keep], node = n_tip + keep, frame = "none",
                           cex = 0.45, col = "black")
            }
        }
        if (is.finite(candidate_node) && !is.null(candidate_edge) && is.finite(candidate_edge$support)) {
            if (!is.na(candidate_edge$support)) {
                nodelabels(candidate_edge$support, node = candidate_node, frame = "rect",
                           bg = "white", cex = 0.60, col = "#8e1b17")
            }
        }
        if (!is.null(tr$edge.length) && any(is.finite(tr$edge.length))) {
            try(add.scale.bar(cex = 0.55, lwd = 0.8), silent = TRUE)
        }
        title(main = sprintf("%s: PhyML candidate-clade phylogeny (%s tips)", locus_id,
                             format(n_tip, big.mark = ",")), cex.main = 0.85)
        note <- if ("Ancestral" %in% tr$tip.label) {
            "ML topology rooted on the 1KG INFO/AA inferred ancestral sequence"
        } else {
            "Unrooted ML topology; display-rooted outside the candidate clade"
        }
        if (!show_all_labels) note <- paste0(note, "; non-candidate labels omitted")
        mtext(note, side = 1, line = 0.1, cex = 0.55)
        dev.off()
        device_open <- FALSE
        if (!file.exists(tmp) || file.info(tmp)$size <= 0) stop("PNG device produced an empty file")
        if (file.exists(png_file)) unlink(png_file)
        if (!file.rename(tmp, png_file)) stop("Could not atomically install plot: ", png_file)
    }, error = function(e) {
        if (device_open) try(dev.off(), silent = TRUE)
        if (file.exists(tmp)) unlink(tmp)
        stop("Failed to render ", tree_file, ": ", conditionMessage(e), call. = FALSE)
    })
    message("PHYML TREE PLOT: status=created tips=", n_tip, " output=", png_file)
    png_file
}

plots <- Filter(Negate(is.null), lapply(phy_inputs, render_one))
message("PHYML TREE PLOTS: created_or_current=", length(plots), " output=", out)
