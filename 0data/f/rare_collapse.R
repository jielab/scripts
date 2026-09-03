#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: rare_collapse.R INPUT_TSV GENE_BED NCBI_GENE_LOC OUTPUT_TSV", call. = FALSE)
}

input <- args[[1L]]
bed_file <- args[[2L]]
gene_loc_file <- args[[3L]]
output <- args[[4L]]

read_tsv <- function(path, header = TRUE) {
  utils::read.delim(path, header = header, quote = "\"", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE)
}

x <- read_tsv(input)
required <- c("Gene", "Protein", "p-value")
missing_required <- setdiff(required, names(x))
if (length(missing_required)) {
  stop("Input is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

p <- suppressWarnings(as.numeric(x[["p-value"]]))
if (any(!is.finite(p))) {
  stop("Input contains ", sum(!is.finite(p)), " non-numeric or missing p-value entries", call. = FALSE)
}

# Sorting by p-value and then original row number gives deterministic first-row
# tie handling. duplicated() then retains exactly one row per Protein/Gene pair.
pair_key <- paste(x$Protein, x$Gene, sep = "\034")
ord <- order(p, seq_len(nrow(x)))
keep <- ord[!duplicated(pair_key[ord])]
collapsed <- x[keep, , drop = FALSE]

bed <- read_tsv(bed_file, header = FALSE)
if (ncol(bed) < 4L) stop("BED must contain at least four columns", call. = FALSE)
bed <- bed[, 1:4, drop = FALSE]
names(bed) <- c("Gene_chr", "Gene_Start", "Gene_end", "Gene")
if (anyDuplicated(bed$Gene)) {
  duplicate_genes <- unique(bed$Gene[duplicated(bed$Gene)])
  stop("BED column 4 contains duplicate gene names (for example: ",
       paste(utils::head(duplicate_genes, 5L), collapse = ", "), ")", call. = FALSE)
}

bed_index <- match(collapsed$Gene, bed$Gene)
bed_matched <- !is.na(bed_index)

gene_loc <- read_tsv(gene_loc_file, header = FALSE)
if (ncol(gene_loc) < 6L) stop("NCBI gene.loc must contain at least six columns", call. = FALSE)
gene_loc <- gene_loc[, c(2L, 3L, 4L, 6L), drop = FALSE]
names(gene_loc) <- c("Gene_chr", "Gene_Start", "Gene_end", "Gene")
target_loc <- gene_loc$Gene %in% collapsed$Gene[!bed_matched]
if (anyDuplicated(gene_loc$Gene[target_loc])) {
  duplicate_genes <- unique(gene_loc$Gene[target_loc][duplicated(gene_loc$Gene[target_loc])])
  stop("NCBI gene.loc contains duplicate target gene names (for example: ",
       paste(utils::head(duplicate_genes, 5L), collapse = ", "), ")", call. = FALSE)
}
loc_index <- match(collapsed$Gene, gene_loc$Gene)
loc_matched <- !bed_matched & !is.na(loc_index)
matched <- bed_matched | loc_matched

coordinates <- data.frame(Gene_chr = rep(NA_character_, nrow(collapsed)),
                          Gene_Start = rep(NA_real_, nrow(collapsed)),
                          Gene_end = rep(NA_real_, nrow(collapsed)),
                          stringsAsFactors = FALSE)
coordinates[bed_matched, ] <- bed[bed_index[bed_matched],
                                  c("Gene_chr", "Gene_Start", "Gene_end"), drop = FALSE]
coordinates[loc_matched, ] <- gene_loc[loc_index[loc_matched],
                                       c("Gene_chr", "Gene_Start", "Gene_end"), drop = FALSE]
out <- cbind(collapsed, coordinates)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(output, ".tmp.", Sys.getpid())
on.exit(unlink(tmp, force = TRUE), add = TRUE)
utils::write.table(out, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
if (file.exists(output)) unlink(output, force = TRUE)
if (!file.rename(tmp, output)) stop("Could not install output: ", output, call. = FALSE)

message("rare collapse: input_rows=", nrow(x),
        " unique_Protein_Gene=", nrow(collapsed),
        " BED_matched=", sum(bed_matched),
        " NCBI_fallback_matched=", sum(loc_matched),
        " unmatched=", sum(!matched),
        " output=", output)
if (any(!matched)) {
  message("rare collapse: unmatched genes retained with NA coordinates: ",
          paste(utils::head(sort(unique(collapsed$Gene[!matched])), 20L), collapse = ", "),
          if (length(unique(collapsed$Gene[!matched])) > 20L) ", ..." else "")
}
