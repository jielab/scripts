#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value for ", name, call. = FALSE)
  args[i + 1L]
}

model_dir <- get_arg("--dir")
pops <- strsplit(get_arg("--pops", "AFR,EAS,EUR,SAS"), ",", fixed = TRUE)[[1L]]
model_id <- get_arg("--model-id", basename(model_dir))
phi <- get_arg("--phi", NA_character_)
manifest <- get_arg("--manifest")
qc_file <- get_arg("--qc", file.path(model_dir, "score_qc.tsv"))
if (is.null(model_dir) || !dir.exists(model_dir)) stop("Missing --dir directory", call. = FALSE)

manifest_rows <- list()
qc_rows <- list()
for (pop in pops) {
  pop_dir <- file.path(model_dir, pop)
  files <- file.path(pop_dir, paste0(pop, ".chr", 1:22, ".sscore"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing score files: ", paste(missing, collapse = ", "), call. = FALSE)

  scores <- rbindlist(lapply(seq_along(files), function(chr) {
    x <- fread(files[chr], showProgress = FALSE)
    id_col <- intersect(c("IID", "#IID"), names(x))[1L]
    allele_col <- grep("ALLELE_CT$", names(x), value = TRUE)[1L]
    score_col <- grep("_SUM$", names(x), value = TRUE)[1L]
    if (anyNA(c(id_col, allele_col, score_col))) stop("Unexpected PLINK2 columns in ", files[chr], call. = FALSE)
    if (anyDuplicated(x[[id_col]])) stop("Duplicate IDs in ", files[chr], call. = FALSE)
    x[, .(eid = as.character(get(id_col)), allele = as.numeric(get(allele_col)),
          score = as.numeric(get(score_col)), chr = chr)]
  }))

  coverage <- scores[, .(n_chr = uniqueN(chr)), by = eid]
  incomplete <- coverage[n_chr != 22L]
  if (nrow(incomplete)) stop(nrow(incomplete), " IDs are absent from one or more chromosome score files", call. = FALSE)

  result <- scores[, .(allele_cnt = sum(allele), PRS = sum(score)), by = eid]
  if (!all(is.finite(result$PRS))) stop("Non-finite PRS values for ", pop, call. = FALSE)
  setnames(result, c("allele_cnt", "PRS"), c(paste0("allele_cnt.", pop), paste0("score_sum.", pop)))
  outfile <- file.path(model_dir, paste0(pop, ".prs.gz"))
  fwrite(result, outfile, sep = "\t", quote = FALSE)

  qc_rows[[pop]] <- data.table(model_id = model_id, phi = phi, population = pop,
    n = nrow(result), mean_prs = mean(result[[paste0("score_sum.", pop)]]),
    sd_prs = sd(result[[paste0("score_sum.", pop)]]),
    min_prs = min(result[[paste0("score_sum.", pop)]]),
    max_prs = max(result[[paste0("score_sum.", pop)]]),
    mean_allele_count = mean(result[[paste0("allele_cnt.", pop)]]), file = normalizePath(outfile))
  manifest_rows[[pop]] <- data.table(model_id = model_id, phi = phi, pop = pop,
    file = normalizePath(outfile), score_col = paste0("score_sum.", pop))
}

qc <- rbindlist(qc_rows)
dir.create(dirname(qc_file), recursive = TRUE, showWarnings = FALSE)
fwrite(qc, qc_file, sep = "\t", quote = FALSE)
if (!is.null(manifest)) {
  m <- rbindlist(manifest_rows)
  dir.create(dirname(manifest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(manifest)) {
    old <- fread(manifest)
    mid <- model_id
    old <- old[model_id != mid]
    m <- rbind(old, m, fill = TRUE)
  }
  setorder(m, model_id, pop)
  fwrite(m, manifest, sep = "\t", quote = FALSE)
}
cat("Combined PRS for", paste(pops, collapse = ","), "under", model_id, "\n")
