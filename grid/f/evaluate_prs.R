#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value after ", key, call. = FALSE)
  args[[i + 1L]]
}
as_bool <- function(x) toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")

trait <- get_arg("--trait")
phe_file <- get_arg("--phe")
pca_file <- get_arg("--pca")
med_file <- get_arg("--med")
csx_manifest_file <- get_arg("--csx-manifest")
outdir <- get_arg("--out")
ancestry_file <- get_arg("--ancestry-file", "")
ancestry_id_col <- get_arg("--ancestry-id-col", "eid")
ancestry_col <- get_arg("--ancestry-col", "")
ancestry_prob_col <- get_arg("--ancestry-prob-col", "")
self_report_col <- get_arg("--self-report-col", "")
target_ids_file <- get_arg("--target-ids", "")
exclude_ids_file <- get_arg("--exclude-ids", "")
prscsx_repeats <- as.integer(get_arg("--prscsx-repeats", get_arg("--repeats", "100")))
disco_repeats <- as.integer(get_arg("--disco-repeats", "20"))
validation_frac <- as.numeric(get_arg("--validation-frac", "0.5"))
fine_tune_n <- as.integer(get_arg("--fine-tune-n", "1300"))
prscsx_prob_min <- as.numeric(get_arg("--prscsx-prob-min", "0.9"))
fine_prob_min <- as.numeric(get_arg("--fine-tune-prob-min", "0.999999"))
require_self_match <- as_bool(get_arg("--require-self-match", "TRUE"))
cov_n_pc <- as.integer(get_arg("--cov-pcs", get_arg("--pcs", "20")))
distance_n_pc <- as.integer(get_arg("--distance-pcs", "10"))
min_test_group <- as.integer(get_arg("--min-test-group", get_arg("--min-group", "200")))
seed <- as.integer(get_arg("--seed", "20260802"))
run_prscsx <- as_bool(get_arg("--run-prscsx", "TRUE"))
run_disco <- as_bool(get_arg("--run-disco", "TRUE"))
write_predictions <- as_bool(get_arg("--write-predictions", "FALSE"))
disco_a_arg <- get_arg("--disco-a", "1,1,1,1")

required <- list(trait, phe_file, pca_file, med_file, csx_manifest_file, outdir)
if (any(vapply(required, is.null, logical(1)))) stop("Required: --trait --phe --pca --med --csx-manifest --out", call. = FALSE)
if (!run_prscsx && !run_disco) stop("At least one of --run-prscsx or --run-disco must be TRUE", call. = FALSE)
if (!is.finite(prscsx_repeats) || prscsx_repeats < 2L) stop("--prscsx-repeats must be >=2", call. = FALSE)
if (!is.finite(disco_repeats) || disco_repeats < 2L) stop("--disco-repeats must be >=2", call. = FALSE)
if (!is.finite(validation_frac) || validation_frac <= .1 || validation_frac >= .9) stop("--validation-frac must be between 0.1 and 0.9", call. = FALSE)
if (!is.finite(fine_tune_n) || fine_tune_n < 100L) stop("--fine-tune-n must be >=100", call. = FALSE)
if (!is.finite(prscsx_prob_min) || prscsx_prob_min < 0 || prscsx_prob_min > 1) stop("--prscsx-prob-min must be within [0,1]", call. = FALSE)
if (!is.finite(fine_prob_min) || fine_prob_min < 0 || fine_prob_min > 1) stop("--fine-tune-prob-min must be within [0,1]", call. = FALSE)
if (!is.finite(cov_n_pc) || cov_n_pc < 5L) stop("--cov-pcs must be >=5", call. = FALSE)
if (!is.finite(distance_n_pc) || distance_n_pc < 5L || distance_n_pc > cov_n_pc) stop("Invalid --distance-pcs", call. = FALSE)
if (!is.finite(min_test_group) || min_test_group < 50L) stop("--min-test-group must be >=50", call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

POP4 <- c("AFR", "EAS", "EUR", "SAS")
ALL_GROUPS <- c("AFR", "EAS", "EUR", "SAS", "AMR", "OTH")
M_CSX <- "PRS-CSx tuned"
M_BEST <- "Best single component"
M_CONV <- "Conventional matched/similar"
M_DISCO <- "DiscoDivas actual centers"
M_DISCO_1KG <- "DiscoDivas 1KG centers"
M_NEAREST <- "Nearest single-ancestry model"

quantile_safe <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(rep(NA_real_, length(p)))
  unname(quantile(x, p, names = FALSE, type = 7))
}
fit_fast <- function(X, y) {
  X <- as.matrix(X); y <- as.numeric(y)
  f <- lm.fit(X, y)
  b <- f$coefficients; b[!is.finite(b)] <- 0
  fitted <- as.numeric(X %*% b)
  residuals <- y - fitted
  list(coefficients = b, fitted = fitted, residuals = residuals, sse = sum(residuals^2), rank = f$rank)
}
paper_r2 <- function(y, Xcov, score) {
  y <- as.numeric(y); Xcov <- as.matrix(Xcov); score <- as.numeric(score)
  if (length(y) != length(score) || nrow(Xcov) != length(y)) stop("Metric input length mismatch", call. = FALSE)
  if (length(y) <= ncol(Xcov) + 10L) {
    return(list(n = length(y), partial_r2 = NA_real_, partial_r = NA_real_, delta_total_r2 = NA_real_,
                covariate_r2 = NA_real_, full_model_r2 = NA_real_, score_beta = NA_real_, score_beta_std = NA_real_))
  }
  f0 <- fit_fast(Xcov, y)
  f1 <- fit_fast(cbind(Xcov, PRS = score), y)
  fs <- fit_fast(Xcov, score)
  r <- suppressWarnings(cor(f0$residuals, fs$residuals)); if (!is.finite(r)) r <- NA_real_
  sst <- sum((y - mean(y))^2)
  beta <- tail(f1$coefficients, 1L)
  list(n = length(y), partial_r2 = (f0$sse - f1$sse) / f0$sse, partial_r = r,
       delta_total_r2 = (f0$sse - f1$sse) / sst, covariate_r2 = 1 - f0$sse / sst,
       full_model_r2 = 1 - f1$sse / sst, score_beta = beta,
       score_beta_std = beta * sd(score) / sd(y))
}
normalise_group <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z <- gsub("[_-]", " ", z)
  out <- rep(NA_character_, length(z))
  # Preserve explicit canonical labels first, then map common descriptions.
  canonical <- z %chin% ALL_GROUPS
  out[canonical] <- z[canonical]
  set_if_missing <- function(pattern, value) {
    hit <- is.na(out) & grepl(pattern, z)
    out[hit] <<- value
  }
  set_if_missing("^(AFR|AFRICAN|BLACK)", "AFR")
  set_if_missing("^(EAS|EAST ASIAN|CHINESE|JAPANESE|KOREAN)", "EAS")
  set_if_missing("^(EUR|EUROPEAN|WHITE)", "EUR")
  set_if_missing("^(SAS|SOUTH ASIAN|INDIAN|PAKISTANI|BANGLADESHI)", "SAS")
  set_if_missing("^(AMR|ADMIXED AMERICAN|HISPANIC|LATINO|LATIN AMERICAN)", "AMR")
  set_if_missing("^(OTH|OTHER|MIXED|ADMIX|UNCLASS|UNKNOWN)", "OTH")
  out
}
row_euclidean <- function(X, centers) {
  X <- as.matrix(X); centers <- as.matrix(centers)
  out <- vapply(seq_len(nrow(centers)), function(j) sqrt(rowSums(sweep(X, 2L, centers[j, ], FUN = "-")^2)), numeric(nrow(X)))
  if (is.null(dim(out))) out <- matrix(out, ncol = nrow(centers))
  colnames(out) <- rownames(centers)
  out
}
shrinkage_vector <- function(centers) {
  D <- as.matrix(dist(as.matrix(centers), upper = TRUE, diag = TRUE))
  d <- tryCatch(solve(D, rep(1, nrow(D))), error = function(e) NULL)
  if (is.null(d) || any(!is.finite(d))) d <- tryCatch(qr.solve(D, rep(1, nrow(D)), tol = 1e-12), error = function(e) NULL)
  if (is.null(d) || any(!is.finite(d))) stop("Fine-tuning-center distance matrix is singular", call. = FALSE)
  list(d = as.numeric(d), D = D, condition_number = kappa(D))
}
calculate_disco_weights <- function(target_pc, centers, A) {
  centers <- as.matrix(centers); target_pc <- as.matrix(target_pc)
  sv <- shrinkage_vector(centers)
  dist_mat <- row_euclidean(target_pc, centers)
  dist_safe <- pmax(dist_mat, sqrt(.Machine$double.eps))
  w <- sweep(1 / dist_safe, 2L, sv$d * A, FUN = "*")
  den <- rowSums(w); scale0 <- rowSums(abs(w))
  bad <- !is.finite(den) | !is.finite(scale0) | scale0 == 0 | abs(den) <= .Machine$double.eps * scale0
  if (any(bad)) stop("Invalid DiscoDivas denominator in ", sum(bad), " rows", call. = FALSE)
  a <- sweep(w, 1L, den, FUN = "/")
  if (any(!is.finite(a))) stop("Non-finite DiscoDivas weights", call. = FALSE)
  colnames(a) <- rownames(centers)
  list(weights = a, distances = dist_mat, d = sv$d, center_distance = sv$D,
       condition_number = sv$condition_number, row_sum_error = max(abs(rowSums(a) - 1)))
}
residualize_and_scale <- function(score_matrix, pc_matrix) {
  score_matrix <- as.matrix(score_matrix)
  X <- cbind(`(Intercept)` = 1, as.matrix(pc_matrix))
  z <- matrix(NA_real_, nrow(score_matrix), ncol(score_matrix), dimnames = dimnames(score_matrix))
  coef <- matrix(NA_real_, ncol(X), ncol(score_matrix), dimnames = list(colnames(X), colnames(score_matrix)))
  mu <- ss <- numeric(ncol(score_matrix))
  for (j in seq_len(ncol(score_matrix))) {
    f <- fit_fast(X, score_matrix[, j]); r <- f$residuals
    mu[j] <- mean(r); ss[j] <- sd(r)
    if (!is.finite(ss[j]) || ss[j] <= 0) stop("Zero residual SD for fine-tuned PRS", call. = FALSE)
    z[, j] <- (r - mu[j]) / ss[j]; coef[, j] <- f$coefficients
  }
  names(mu) <- names(ss) <- colnames(score_matrix)
  list(score = z, coefficients = coef, residual_mean = mu, residual_sd = ss)
}
fit_combination <- function(idx, raw_csx, Xcov, y, component_names) {
  raw <- raw_csx[idx, , drop = FALSE]
  mu <- colMeans(raw); ss <- apply(raw, 2L, sd)
  if (any(!is.finite(ss) | ss <= 0)) stop("Zero/non-finite PRS SD in fine-tuning data", call. = FALSE)
  z <- sweep(sweep(raw, 2L, mu, FUN = "-"), 2L, ss, FUN = "/")
  fit <- fit_fast(cbind(Xcov[idx, , drop = FALSE], z), y[idx])
  beta_std <- tail(fit$coefficients, ncol(z)); beta_std[!is.finite(beta_std)] <- 0
  names(beta_std) <- component_names
  raw_weights <- beta_std / ss; names(raw_weights) <- component_names
  score <- as.numeric(raw %*% raw_weights)
  val <- paper_r2(y[idx], Xcov[idx, , drop = FALSE], score)
  single <- vapply(seq_len(ncol(raw)), function(j) paper_r2(y[idx], Xcov[idx, , drop = FALSE], raw[, j])$partial_r2, numeric(1))
  best_j <- which.max(replace(single, !is.finite(single), -Inf))
  list(n = length(idx), component_mean = mu, component_sd = ss,
       standardized_coefficients = beta_std, raw_weights = raw_weights,
       validation_partial_r2 = val$partial_r2, single_validation_r2 = single,
       best_component = component_names[best_j])
}
apply_combination <- function(model, raw) as.numeric(as.matrix(raw) %*% model$raw_weights)
make_decile <- function(x, n = 10L) {
  if (length(x) < n) return(rep(NA_integer_, length(x)))
  pmin(n, pmax(1L, ceiling(n * rank(x, ties.method = "average", na.last = "keep") / sum(is.finite(x)))))
}
metric_summary <- function(z, group_cols) {
  z[, {
    q <- quantile_safe(partial_r2, c(.025, .975))
    .(mean_partial_r2 = mean(partial_r2, na.rm = TRUE), sd_partial_r2 = sd(partial_r2, na.rm = TRUE),
      q025_partial_r2 = q[1L], q975_partial_r2 = q[2L],
      mean_delta_total_r2 = mean(delta_total_r2, na.rm = TRUE), mean_covariate_r2 = mean(covariate_r2, na.rm = TRUE),
      mean_full_model_r2 = mean(full_model_r2, na.rm = TRUE), mean_n = mean(n, na.rm = TRUE), repeats = .N)
  }, by = group_cols]
}
paired_summary <- function(metrics, method_a, method_b, design_label) {
  x <- metrics[method %chin% c(method_a, method_b)]
  w <- dcast(x, design + eval_group + `repeat` ~ method, value.var = "partial_r2")
  if (!all(c(method_a, method_b) %in% names(w))) return(data.table())
  w[, difference := get(method_b) - get(method_a)]
  w[, relative_change := fifelse(abs(get(method_a)) > 1e-12, difference / get(method_a), NA_real_)]
  w[, {
    qd <- quantile_safe(difference, c(.025, .975)); qr <- quantile_safe(relative_change, c(.025, .975))
    tp <- tryCatch(t.test(difference)$p.value, error = function(e) NA_real_)
    wp <- tryCatch(wilcox.test(difference, exact = FALSE)$p.value, error = function(e) NA_real_)
    .(comparator = method_a, method = method_b, mean_difference = mean(difference, na.rm = TRUE),
      sd_difference = sd(difference, na.rm = TRUE), q025_difference = qd[1L], q975_difference = qd[2L],
      mean_relative_change = mean(relative_change, na.rm = TRUE), median_relative_change = median(relative_change, na.rm = TRUE),
      q025_relative_change = qr[1L], q975_relative_change = qr[2L],
      fraction_repeats_favoring_method = mean(difference > 0, na.rm = TRUE), paired_t_p = tp, wilcoxon_p = wp, repeats = .N)
  }, by = .(design, eval_group)][, requested_design := design_label][]
}
file_info <- function(label, path) {
  ok <- length(path) == 1L && nzchar(path) && file.exists(path)
  fi <- if (ok) file.info(path) else NULL
  data.table(label = label, path = as.character(path), exists = ok,
             bytes = if (ok) as.numeric(fi$size) else NA_real_,
             modified = if (ok) format(fi$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_)
}

message("Reading phenotype, PCA, ancestry information, and PRS-CSx scores ...")
phe0 <- readRDS(phe_file)
base_phe_cols <- unique(c("eid", trait, "age", "sex",
                          if (!nzchar(ancestry_file)) c(ancestry_col, ancestry_prob_col, self_report_col) else character()))
base_phe_cols <- base_phe_cols[nzchar(base_phe_cols)]
miss <- setdiff(base_phe_cols, names(phe0))
if (length(miss)) stop("Phenotype RDS is missing: ", paste(miss, collapse = ", "), call. = FALSE)
phe <- as.data.table(phe0[, base_phe_cols, drop = FALSE])
setnames(phe, trait, "phenotype")
phe[, eid := as.character(eid)]

if (nzchar(ancestry_file)) {
  if (!file.exists(ancestry_file)) stop("Missing ancestry file: ", ancestry_file, call. = FALSE)
  anc <- fread(ancestry_file, showProgress = FALSE)
  needed <- unique(c(ancestry_id_col, ancestry_col, ancestry_prob_col, self_report_col))
  needed <- needed[nzchar(needed)]
  if (!all(needed %in% names(anc))) stop("Ancestry file lacks: ", paste(setdiff(needed, names(anc)), collapse = ", "), call. = FALSE)
  anc <- anc[, ..needed]
  setnames(anc, ancestry_id_col, "eid")
  anc[, eid := as.character(eid)]
  if (anyDuplicated(anc$eid)) stop("Duplicate IDs in ancestry file", call. = FALSE)
  phe <- merge(phe, anc, by = "eid", all.x = TRUE, sort = FALSE)
}

if (nzchar(target_ids_file)) {
  ids <- as.character(fread(target_ids_file, header = FALSE, showProgress = FALSE)[[1L]])
  phe <- phe[eid %chin% ids]
}
if (nzchar(exclude_ids_file)) {
  ids <- as.character(fread(exclude_ids_file, header = FALSE, showProgress = FALSE)[[1L]])
  phe <- phe[!eid %chin% ids]
}
if (!nrow(phe)) stop("No phenotype participants remain after ID filters", call. = FALSE)

pca <- fread(pca_file, showProgress = FALSE)
id_col <- intersect(c("IID", "#IID", "eid"), names(pca))[1L]
if (length(id_col) != 1L || is.na(id_col)) stop("PCA file needs IID/#IID/eid", call. = FALSE)
setnames(pca, id_col, "eid"); pca[, eid := as.character(eid)]
cov_pc_cols <- paste0("PC", seq_len(cov_n_pc))
dist_pc_cols <- paste0("PC", seq_len(distance_n_pc))
if (!all(cov_pc_cols %in% names(pca))) stop("PCA file lacks covariate PCs: ", paste(setdiff(cov_pc_cols, names(pca)), collapse = ", "), call. = FALSE)

med <- fread(med_file, showProgress = FALSE)
pop_col <- names(med)[1L]
if (!all(dist_pc_cols %in% names(med))) stop("Reference median file lacks: ", paste(setdiff(dist_pc_cols, names(med)), collapse = ", "), call. = FALSE)
reference_names <- toupper(as.character(med[[pop_col]]))
reference_centers <- as.matrix(med[, ..dist_pc_cols]); storage.mode(reference_centers) <- "double"; rownames(reference_centers) <- reference_names
if (!all(POP4 %in% reference_names)) stop("Reference medians must include AFR,EAS,EUR,SAS", call. = FALSE)
reference_centers <- reference_centers[POP4, , drop = FALSE]

csx_manifest <- fread(csx_manifest_file, showProgress = FALSE)
need_csx <- c("model_id", "phi", "pop", "file", "score_col")
if (!all(need_csx %in% names(csx_manifest))) stop("CSx manifest needs: ", paste(need_csx, collapse = ", "), call. = FALSE)
csx_manifest[, pop := toupper(as.character(pop))]
if (uniqueN(csx_manifest$model_id) != 1L || !setequal(csx_manifest$pop, POP4)) stop("Expected exactly one four-population PRS-CSx model", call. = FALSE)
model_id <- as.character(unique(csx_manifest$model_id)); phi <- as.character(unique(csx_manifest$phi))
csx <- NULL
for (g in POP4) {
  z <- csx_manifest[pop == g]
  x <- fread(z$file, select = c("eid", z$score_col), showProgress = FALSE)
  x[, eid := as.character(eid)]; setnames(x, z$score_col, paste0("PRS_", g))
  if (anyDuplicated(x$eid)) stop("Duplicate IDs in PRS score for ", g, call. = FALSE)
  csx <- if (is.null(csx)) x else merge(csx, x, by = "eid", all = FALSE, sort = FALSE)
}
csx_cols <- paste0("PRS_", POP4)

message("Merging complete participants ...")
dat <- Reduce(function(x, y) merge(x, y, by = "eid", all = FALSE, sort = FALSE), list(phe, pca, csx))
num_cols <- c("phenotype", "age", cov_pc_cols, csx_cols, if (nzchar(ancestry_prob_col)) ancestry_prob_col else character())
num_cols <- intersect(num_cols, names(dat))
dat[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
dat <- dat[complete.cases(dat[, c("phenotype", "age", cov_pc_cols, csx_cols), with = FALSE]) & !is.na(sex)]
if (nrow(dat) < 1000L) stop("Only ", nrow(dat), " complete participants remain", call. = FALSE)
if (uniqueN(dat$phenotype) <= 10L) stop("This evaluator is for quantitative traits", call. = FALSE)
dat[, sex := factor(sex)]

pc_dist <- as.matrix(dat[, ..dist_pc_cols]); storage.mode(pc_dist) <- "double"
assign_dist <- row_euclidean(pc_dist, reference_centers)
assign_i <- max.col(-assign_dist, ties.method = "first")
dat[, `:=`(closest_reference = POP4[assign_i], nearest_reference_distance = assign_dist[cbind(seq_len(.N), assign_i)])]
sorted_assign <- t(apply(assign_dist, 1L, sort, partial = 1:2))
dat[, reference_margin := sorted_assign[, 2L] - sorted_assign[, 1L]]

formal_ancestry <- nzchar(ancestry_col) && ancestry_col %in% names(dat)
if (formal_ancestry) {
  dat[, eval_group := normalise_group(get(ancestry_col))]
  grouping_basis <- paste0("genetic ancestry column '", ancestry_col, "'")
} else {
  dat[, eval_group := closest_reference]
  grouping_basis <- "nearest 1KG center fallback; exploratory and not a paper-level ancestry assignment"
  warning("No valid --ancestry-col supplied. OTH/AMR cannot be identified; results are exploratory.", call. = FALSE)
}
if (nzchar(ancestry_prob_col) && ancestry_prob_col %in% names(dat)) dat[, ancestry_probability := as.numeric(get(ancestry_prob_col))] else dat[, ancestry_probability := NA_real_]
if (nzchar(self_report_col) && self_report_col %in% names(dat)) dat[, self_report_group := normalise_group(get(self_report_col))] else dat[, self_report_group := NA_character_]
dat <- dat[!is.na(eval_group) & eval_group %chin% ALL_GROUPS]
dat[, row_id := .I]

has_prob <- any(is.finite(dat$ancestry_probability))
has_self <- any(!is.na(dat$self_report_group))
dat[, prscsx_eligible := TRUE]
if (has_prob) dat[, prscsx_eligible := ancestry_probability >= prscsx_prob_min]
dat[, fine_eligible_single := eval_group %chin% POP4]
if (has_prob) dat[, fine_eligible_single := fine_eligible_single & ancestry_probability >= fine_prob_min]
if (require_self_match && has_self) dat[, fine_eligible_single := fine_eligible_single & self_report_group == eval_group]
dat[, fine_eligible_oth := eval_group == "OTH"]

Xcov <- model.matrix(reformulate(c("age", "sex", cov_pc_cols)), data = dat)
y <- as.numeric(dat$phenotype)
raw_csx <- as.matrix(dat[, ..csx_cols]); storage.mode(raw_csx) <- "double"; colnames(raw_csx) <- POP4
pc_cov <- as.matrix(dat[, ..cov_pc_cols]); storage.mode(pc_cov) <- "double"
pc_dist <- as.matrix(dat[, ..dist_pc_cols]); storage.mode(pc_dist) <- "double"

A <- as.numeric(strsplit(disco_a_arg, ",", fixed = TRUE)[[1L]])
if (length(A) != 4L || any(!is.finite(A)) || any(A < 0) || sum(A) <= 0) stop("--disco-a must contain four non-negative values", call. = FALSE)
names(A) <- POP4

group_counts <- dat[, .(N = .N, eligible_single = sum(fine_eligible_single), eligible_oth = sum(fine_eligible_oth),
                        mean_age = mean(age), phenotype_mean = mean(phenotype), phenotype_sd = sd(phenotype)), by = eval_group]
group_counts[, group_order := match(eval_group, ALL_GROUPS)]; setorder(group_counts, group_order); group_counts[, group_order := NULL]
report_groups <- group_counts[N >= min_test_group]$eval_group
prscsx_groups <- report_groups[report_groups %chin% c(POP4, "AMR")]
if (!length(report_groups)) stop("No test group has at least --min-test-group participants", call. = FALSE)

ancestry_qc <- dat[, {
  pq <- if (any(is.finite(ancestry_probability))) quantile_safe(ancestry_probability, c(.025, .5, .975)) else rep(NA_real_, 3L)
  dq <- quantile_safe(nearest_reference_distance, c(.025, .5, .975))
  mq <- quantile_safe(reference_margin, c(.025, .5, .975))
  .(N = .N, nearest_center_agreement = mean(eval_group == closest_reference),
    ancestry_probability_mean = mean(ancestry_probability, na.rm = TRUE),
    ancestry_probability_q025 = pq[1L], ancestry_probability_median = pq[2L], ancestry_probability_q975 = pq[3L],
    self_report_match = mean(eval_group == self_report_group, na.rm = TRUE),
    distance_mean = mean(nearest_reference_distance), distance_q025 = dq[1L], distance_median = dq[2L], distance_q975 = dq[3L],
    margin_mean = mean(reference_margin), margin_q025 = mq[1L], margin_median = mq[2L], margin_q975 = mq[3L])
}, by = eval_group]
ancestry_confusion <- dat[, .(N = .N), by = .(eval_group, closest_reference)]
ancestry_confusion[, percent_within_group := 100 * N / sum(N), by = eval_group]
component_correlations <- rbindlist(lapply(report_groups, function(g) {
  idx <- which(dat$eval_group == g)
  C <- suppressWarnings(cor(raw_csx[idx, , drop = FALSE], use = "pairwise.complete.obs"))
  z <- as.data.table(C, keep.rownames = "component1")
  melt(z, id.vars = "component1", variable.name = "component2", value.name = "correlation")[, eval_group := g][]
}), fill = TRUE)

run_info <- data.table(
  item = c("trait", "complete_n", "phi", "grouping_basis", "covariate_pcs", "distance_pcs",
           "PRS-CSx_repeats", "PRS-CSx_validation_fraction", "DiscoDivas_repeats", "fine_tune_n",
           "PRS-CSx_probability_min", "fine_probability_min", "require_self_report_match", "has_probability_column", "has_self_report_column",
           "DiscoDivas_A", "target_ids", "exclude_ids"),
  value = c(trait, nrow(dat), phi, grouping_basis, cov_n_pc, distance_n_pc,
            prscsx_repeats, validation_frac, disco_repeats, fine_tune_n, prscsx_prob_min, fine_prob_min,
            require_self_match, has_prob, has_self, paste(A, collapse = ","), target_ids_file, exclude_ids_file)
)


# 🚩 Analysis 1: PRS-CSx paper-style repeated 50:50 validation/testing splits.
prscsx_metrics <- data.table(); prscsx_weights <- data.table(); prscsx_components <- data.table()
if (run_prscsx) {
  message("Running PRS-CSx repeated split evaluation: ", prscsx_repeats, " repeats ...")
  met_rows <- list(); wt_rows <- list(); comp_rows <- list(); mc <- wc <- sc <- 0L
  for (r in seq_len(prscsx_repeats)) {
    set.seed(seed + r)
    for (g in prscsx_groups) {
      idx <- which(dat$eval_group == g & dat$prscsx_eligible)
      n_val <- floor(length(idx) * validation_frac)
      if (n_val <= ncol(Xcov) + ncol(raw_csx) + 20L || length(idx) - n_val < min_test_group) next
      val <- sample(idx, n_val, replace = FALSE); test <- setdiff(idx, val)
      mdl <- fit_combination(val, raw_csx, Xcov, y, POP4)
      tuned <- apply_combination(mdl, raw_csx[test, , drop = FALSE])
      best <- raw_csx[test, match(mdl$best_component, POP4)]
      for (method in c(M_CSX, M_BEST)) {
        score <- if (method == M_CSX) tuned else best
        met <- paper_r2(y[test], Xcov[test, , drop = FALSE], score)
        mc <- mc + 1L
        met_rows[[mc]] <- as.data.table(c(list(design = "PRS-CSx_50_50", eval_group = g, `repeat` = r,
                                                    method = method, n_validation = length(val), n_testing = length(test)), met))
      }
      for (j in seq_along(POP4)) {
        wc <- wc + 1L
        wt_rows[[wc]] <- data.table(design = "PRS-CSx_50_50", eval_group = g, `repeat` = r,
                                    component = POP4[j], standardized_coefficient = mdl$standardized_coefficients[j],
                                    raw_weight = mdl$raw_weights[j], validation_partial_r2 = mdl$validation_partial_r2,
                                    best_component = mdl$best_component)
        sc <- sc + 1L
        comp_rows[[sc]] <- data.table(design = "PRS-CSx_50_50", eval_group = g, `repeat` = r,
                                      component = POP4[j], validation_partial_r2 = mdl$single_validation_r2[j])
      }
    }
    if (r %% 10L == 0L || r == prscsx_repeats) message("  PRS-CSx repeat ", r, "/", prscsx_repeats)
  }
  prscsx_metrics <- rbindlist(met_rows, fill = TRUE)
  prscsx_weights <- rbindlist(wt_rows, fill = TRUE)
  prscsx_components <- rbindlist(comp_rows, fill = TRUE)
}


# 🚩 Analysis 2: DiscoDivas paper design
# Fixed 1,300-person fine-tuning cohorts; remaining participants are used for
# testing, and the OTH conventional model is used for OTH/AMR.
disco_metrics <- data.table(); disco_sensitivity <- data.table(); disco_weights <- data.table()
transfer_metrics <- data.table(); continuum_metrics <- data.table(); calibration_metrics <- data.table()
fine_centers_long <- data.table(); disco_model_weights <- data.table(); disco_qc <- data.table()
first_predictions <- NULL
if (run_disco) {
  eligible_counts <- setNames(vapply(POP4, function(g) sum(dat$eval_group == g & dat$fine_eligible_single), integer(1)), POP4)
  if (any(eligible_counts < fine_tune_n)) {
    stop("Not enough eligible samples for paper-style fine-tuning: ",
         paste0(names(eligible_counts), "=", eligible_counts, collapse = ", "),
         "; requested fine_tune_n=", fine_tune_n, call. = FALSE)
  }
  oth_available <- sum(dat$fine_eligible_oth) >= fine_tune_n
  if (!oth_available) warning("Fewer than ", fine_tune_n, " OTH samples: matched OTH/AMR conventional comparator is unavailable", call. = FALSE)
  message("Running DiscoDivas paper-style evaluation: ", disco_repeats, " repeats; fine-tuning N=", fine_tune_n, " ...")
  met_rows <- list(); sens_rows <- list(); weight_rows <- list(); transfer_rows <- list()
  continuum_rows <- list(); calib_rows <- list(); center_rows <- list(); model_w_rows <- list(); qc_rows <- list()
  mc <- sc <- wc <- tc <- dc <- cc <- cec <- mwc <- qc_cur <- 0L

  for (r in seq_len(disco_repeats)) {
    set.seed(seed + 100000L + r)
    fine_idx <- vector("list", length(POP4)); names(fine_idx) <- POP4
    for (g in POP4) {
      pool <- which(dat$eval_group == g & dat$fine_eligible_single)
      fine_idx[[g]] <- sample(pool, fine_tune_n, replace = FALSE)
    }
    oth_idx <- integer()
    if (oth_available) oth_idx <- sample(which(dat$fine_eligible_oth), fine_tune_n, replace = FALSE)
    excluded <- sort(unique(c(unlist(fine_idx, use.names = FALSE), oth_idx)))
    test_idx <- setdiff(seq_len(nrow(dat)), excluded)
    test_group <- dat$eval_group[test_idx]

    models <- vector("list", length(POP4)); names(models) <- POP4
    centers <- matrix(NA_real_, 4L, distance_n_pc, dimnames = list(POP4, dist_pc_cols))
    fine_scores <- matrix(NA_real_, length(test_idx), 4L, dimnames = list(NULL, POP4))
    for (g in POP4) {
      mdl <- fit_combination(fine_idx[[g]], raw_csx, Xcov, y, POP4)
      models[[g]] <- mdl
      centers[g, ] <- apply(pc_dist[fine_idx[[g]], , drop = FALSE], 2L, median)
      fine_scores[, g] <- apply_combination(mdl, raw_csx[test_idx, , drop = FALSE])
      for (j in seq_along(POP4)) {
        mwc <- mwc + 1L
        model_w_rows[[mwc]] <- data.table(`repeat` = r, fine_tuning_group = g, component = POP4[j],
                                          standardized_coefficient = mdl$standardized_coefficients[j],
                                          raw_weight = mdl$raw_weights[j], validation_partial_r2 = mdl$validation_partial_r2)
      }
      for (j in seq_along(dist_pc_cols)) {
        cec <- cec + 1L
        center_rows[[cec]] <- data.table(`repeat` = r, fine_tuning_group = g, PC = dist_pc_cols[j], median = centers[g, j])
      }
    }
    oth_model <- NULL; oth_score <- rep(NA_real_, length(test_idx))
    if (oth_available) {
      oth_model <- fit_combination(oth_idx, raw_csx, Xcov, y, POP4)
      oth_score <- apply_combination(oth_model, raw_csx[test_idx, , drop = FALSE])
      for (j in seq_along(POP4)) {
        mwc <- mwc + 1L
        model_w_rows[[mwc]] <- data.table(`repeat` = r, fine_tuning_group = "OTH", component = POP4[j],
                                          standardized_coefficient = oth_model$standardized_coefficients[j],
                                          raw_weight = oth_model$raw_weights[j], validation_partial_r2 = oth_model$validation_partial_r2)
      }
    }

    transform <- residualize_and_scale(fine_scores, pc_dist[test_idx, , drop = FALSE])
    actual <- calculate_disco_weights(pc_dist[test_idx, , drop = FALSE], centers, A)
    fixed <- calculate_disco_weights(pc_dist[test_idx, , drop = FALSE], reference_centers, A)
    disco_actual <- rowSums(transform$score * actual$weights)
    disco_fixed <- rowSums(transform$score * fixed$weights)
    nearest_i <- max.col(-actual$distances, ties.method = "first")
    nearest_model <- POP4[nearest_i]
    nearest_score <- fine_scores[cbind(seq_along(test_idx), nearest_i)]
    nearest_distance <- actual$distances[cbind(seq_along(test_idx), nearest_i)]
    sorted_distance <- t(apply(actual$distances, 1L, sort, partial = 1:2))
    distance_margin <- sorted_distance[, 2L] - sorted_distance[, 1L]
    max_weight <- apply(actual$weights, 1L, max)
    weight_hhi <- rowSums(actual$weights^2)
    nonnegative <- apply(actual$weights, 1L, function(z) all(z >= 0))
    weight_entropy <- rep(NA_real_, length(test_idx))
    if (any(nonnegative)) {
      p <- actual$weights[nonnegative, , drop = FALSE]
      p[p <= 0] <- NA_real_
      weight_entropy[nonnegative] <- -rowSums(p * log(p), na.rm = TRUE) / log(ncol(p))
    }

    conventional_score <- nearest_score
    conventional_model <- nearest_model
    for (g in POP4) {
      pos <- which(test_group == g)
      if (length(pos)) {
        conventional_score[pos] <- fine_scores[pos, g]
        conventional_model[pos] <- g
      }
    }
    if (oth_available) {
      pos <- which(test_group %chin% c("OTH", "AMR"))
      conventional_score[pos] <- oth_score[pos]
      conventional_model[pos] <- "OTH"
    }

    for (g in report_groups) {
      pos <- which(test_group == g)
      if (length(pos) < min_test_group) next
      score_list <- list(); score_list[[M_CONV]] <- conventional_score[pos]; score_list[[M_DISCO]] <- disco_actual[pos]
      for (method in names(score_list)) {
        met <- paper_r2(y[test_idx[pos]], Xcov[test_idx[pos], , drop = FALSE], score_list[[method]])
        mc <- mc + 1L
        met_rows[[mc]] <- as.data.table(c(list(design = "DiscoDivas_1300", eval_group = g, `repeat` = r,
                                                   method = method, n_fine_tuning = length(excluded), n_testing = length(pos)), met))
      }
      for (method in c(M_DISCO_1KG, M_NEAREST)) {
        score <- if (method == M_DISCO_1KG) disco_fixed[pos] else nearest_score[pos]
        met <- paper_r2(y[test_idx[pos]], Xcov[test_idx[pos], , drop = FALSE], score)
        sc <- sc + 1L
        sens_rows[[sc]] <- as.data.table(c(list(design = "DiscoDivas_1300", eval_group = g, `repeat` = r,
                                                     method = method, n_testing = length(pos)), met))
      }
      for (j in seq_along(POP4)) {
        w0 <- actual$weights[pos, j]; q <- quantile_safe(w0, c(.025, .5, .975))
        wc <- wc + 1L
        weight_rows[[wc]] <- data.table(`repeat` = r, eval_group = g, center_source = "actual fine-tuning cohorts",
                                        component = POP4[j], n = length(pos), mean_weight = mean(w0), sd_weight = sd(w0),
                                        q025 = q[1L], median = q[2L], q975 = q[3L], negative_fraction = mean(w0 < 0))
        w1 <- fixed$weights[pos, j]; q1 <- quantile_safe(w1, c(.025, .5, .975))
        wc <- wc + 1L
        weight_rows[[wc]] <- data.table(`repeat` = r, eval_group = g, center_source = "1KG approximate centers",
                                        component = POP4[j], n = length(pos), mean_weight = mean(w1), sd_weight = sd(w1),
                                        q025 = q1[1L], median = q1[2L], q975 = q1[3L], negative_fraction = mean(w1 < 0))
      }

      transfer_list <- setNames(lapply(POP4, function(m) fine_scores[pos, m]), paste0("Fine-tuned in ", POP4))
      if (oth_available) transfer_list[["Fine-tuned in OTH"]] <- oth_score[pos]
      transfer_list[[M_DISCO]] <- disco_actual[pos]
      for (method in names(transfer_list)) {
        met <- paper_r2(y[test_idx[pos]], Xcov[test_idx[pos], , drop = FALSE], transfer_list[[method]])
        tc <- tc + 1L
        transfer_rows[[tc]] <- data.table(`repeat` = r, eval_group = g, model = method, n = length(pos), partial_r2 = met$partial_r2)
      }

      d_near <- make_decile(nearest_distance[pos]); d_entropy <- make_decile(weight_entropy[pos]); d_max <- make_decile(-max_weight[pos])
      for (axis_name in c("nearest_distance", "weight_entropy", "low_max_weight")) {
        dd <- switch(axis_name, nearest_distance = d_near, weight_entropy = d_entropy, low_max_weight = d_max)
        for (d0 in sort(unique(dd[is.finite(dd)]))) {
          pp <- which(dd == d0)
          if (length(pp) <= ncol(Xcov) + 20L) next
          for (method in c(M_CONV, M_DISCO)) {
            score <- if (method == M_CONV) conventional_score[pos[pp]] else disco_actual[pos[pp]]
            met <- paper_r2(y[test_idx[pos[pp]]], Xcov[test_idx[pos[pp]], , drop = FALSE], score)
            dc <- dc + 1L
            continuum_rows[[dc]] <- data.table(`repeat` = r, eval_group = g, axis = axis_name, decile = d0,
                                               method = method, n = length(pp), partial_r2 = met$partial_r2,
                                               mean_nearest_distance = mean(nearest_distance[pos[pp]]),
                                               mean_entropy = mean(weight_entropy[pos[pp]], na.rm = TRUE),
                                               mean_max_weight = mean(max_weight[pos[pp]]))
          }
        }
      }

      for (method in c(M_CONV, M_DISCO)) {
        score <- if (method == M_CONV) conventional_score[pos] else disco_actual[pos]
        yy <- fit_fast(Xcov[test_idx[pos], , drop = FALSE], y[test_idx[pos]])$residuals
        ss <- fit_fast(Xcov[test_idx[pos], , drop = FALSE], score)$residuals
        dec <- make_decile(ss)
        for (d0 in sort(unique(dec[is.finite(dec)]))) {
          pp <- which(dec == d0)
          cc <- cc + 1L
          calib_rows[[cc]] <- data.table(`repeat` = r, eval_group = g, method = method, score_decile = d0,
                                         n = length(pp), mean_score_residual = mean(ss[pp]),
                                         mean_phenotype_residual = mean(yy[pp]), sd_phenotype_residual = sd(yy[pp]))
        }
      }
    }

    qc_cur <- qc_cur + 1L
    qc_rows[[qc_cur]] <- data.table(`repeat` = r, n_fine_tuning = length(excluded), n_testing = length(test_idx),
                                    actual_center_condition_number = actual$condition_number,
                                    fixed_center_condition_number = fixed$condition_number,
                                    actual_weight_row_sum_error = actual$row_sum_error,
                                    fixed_weight_row_sum_error = fixed$row_sum_error,
                                    min_actual_weight = min(actual$weights), max_actual_weight = max(actual$weights),
                                    fraction_negative_actual_weights = mean(actual$weights < 0),
                                    mean_weight_entropy = mean(weight_entropy, na.rm = TRUE),
                                    mean_weight_hhi = mean(weight_hhi), oth_model_available = oth_available)

    if (r == 1L && write_predictions) {
      first_predictions <- data.table(eid = dat$eid[test_idx], eval_group = test_group,
                                      conventional_model = conventional_model, nearest_single_ancestry_model = nearest_model,
                                      nearest_center_distance = nearest_distance, distance_margin = distance_margin,
                                      weight_entropy = weight_entropy, max_weight = max_weight,
                                      PRS_conventional = conventional_score, PRS_DiscoDivas = disco_actual,
                                      PRS_DiscoDivas_1KG = disco_fixed)
      first_predictions <- cbind(first_predictions, as.data.table(actual$weights))
      setnames(first_predictions, POP4, paste0("Disco_weight_", POP4))
    }
    if (r %% 5L == 0L || r == disco_repeats) message("  DiscoDivas repeat ", r, "/", disco_repeats)
  }
  disco_metrics <- rbindlist(met_rows, fill = TRUE)
  disco_sensitivity <- rbindlist(sens_rows, fill = TRUE)
  disco_weights <- rbindlist(weight_rows, fill = TRUE)
  transfer_metrics <- rbindlist(transfer_rows, fill = TRUE)
  continuum_metrics <- rbindlist(continuum_rows, fill = TRUE)
  calibration_metrics <- rbindlist(calib_rows, fill = TRUE)
  fine_centers_long <- rbindlist(center_rows, fill = TRUE)
  disco_model_weights <- rbindlist(model_w_rows, fill = TRUE)
  disco_qc <- rbindlist(qc_rows, fill = TRUE)
}


# 🚩 Summaries
prscsx_main <- if (nrow(prscsx_metrics)) metric_summary(prscsx_metrics, c("design", "eval_group", "method")) else data.table()
prscsx_paired <- if (nrow(prscsx_metrics)) paired_summary(prscsx_metrics, M_BEST, M_CSX, "PRS-CSx_50_50") else data.table()
prscsx_weight_summary <- if (nrow(prscsx_weights)) prscsx_weights[, {
  q <- quantile_safe(standardized_coefficient, c(.025, .975))
  .(mean_coefficient = mean(standardized_coefficient), sd_coefficient = sd(standardized_coefficient),
    q025 = q[1L], q975 = q[2L], mean_validation_partial_r2 = mean(validation_partial_r2),
    selected_as_best_fraction = mean(best_component == component), repeats = .N)
}, by = .(eval_group, component)] else data.table()

disco_main <- if (nrow(disco_metrics)) metric_summary(disco_metrics, c("design", "eval_group", "method")) else data.table()
disco_paired <- if (nrow(disco_metrics)) paired_summary(disco_metrics, M_CONV, M_DISCO, "DiscoDivas_1300") else data.table()
sensitivity_main <- if (nrow(disco_sensitivity)) metric_summary(disco_sensitivity, c("design", "eval_group", "method")) else data.table()
transfer_summary <- if (nrow(transfer_metrics)) transfer_metrics[, {
  q <- quantile_safe(partial_r2, c(.025, .975))
  .(mean_partial_r2 = mean(partial_r2), sd_partial_r2 = sd(partial_r2), q025 = q[1L], q975 = q[2L], mean_n = mean(n), repeats = .N)
}, by = .(eval_group, model)] else data.table()
disco_weight_summary <- if (nrow(disco_weights)) disco_weights[, {
  q <- quantile_safe(mean_weight, c(.025, .975))
  .(mean_weight = mean(mean_weight), sd_across_repeats = sd(mean_weight), q025 = q[1L], q975 = q[2L],
    mean_negative_fraction = mean(negative_fraction), repeats = .N)
}, by = .(eval_group, center_source, component)] else data.table()
continuum_summary <- if (nrow(continuum_metrics)) continuum_metrics[, {
  q <- quantile_safe(partial_r2, c(.025, .975))
  .(mean_partial_r2 = mean(partial_r2), sd_partial_r2 = sd(partial_r2), q025 = q[1L], q975 = q[2L],
    mean_n = mean(n), mean_nearest_distance = mean(mean_nearest_distance),
    mean_entropy = mean(mean_entropy, na.rm = TRUE), mean_max_weight = mean(mean_max_weight), repeats = .N)
}, by = .(eval_group, axis, decile, method)] else data.table()
calibration_summary <- if (nrow(calibration_metrics)) calibration_metrics[, .(
  mean_score_residual = mean(mean_score_residual), mean_phenotype_residual = mean(mean_phenotype_residual),
  sd_across_repeats = sd(mean_phenotype_residual), mean_n = mean(n), repeats = .N
), by = .(eval_group, method, score_decile)] else data.table()
center_summary <- if (nrow(fine_centers_long)) fine_centers_long[, {
  q <- quantile_safe(median, c(.025, .975))
  .(mean_median = mean(median), sd_median = sd(median), q025 = q[1L], q975 = q[2L])
}, by = .(fine_tuning_group, PC)] else data.table()
disco_model_weight_summary <- if (nrow(disco_model_weights)) disco_model_weights[, {
  q <- quantile_safe(standardized_coefficient, c(.025, .975))
  .(mean_coefficient = mean(standardized_coefficient), sd_coefficient = sd(standardized_coefficient),
    q025 = q[1L], q975 = q[2L], mean_validation_partial_r2 = mean(validation_partial_r2), repeats = .N)
}, by = .(fine_tuning_group, component)] else data.table()

all_metric_rows <- rbindlist(list(prscsx_metrics, disco_metrics, disco_sensitivity), fill = TRUE)
r2_decomposition <- if (nrow(all_metric_rows)) all_metric_rows[, .(
  mean_partial_r2 = mean(partial_r2, na.rm = TRUE), mean_delta_total_r2 = mean(delta_total_r2, na.rm = TRUE),
  mean_covariate_r2 = mean(covariate_r2, na.rm = TRUE), mean_full_model_r2 = mean(full_model_r2, na.rm = TRUE),
  mean_score_beta_std = mean(score_beta_std, na.rm = TRUE), mean_n = mean(n, na.rm = TRUE), repeats = .N
), by = .(design, eval_group, method)] else data.table()


# 🚩 Full-data deployment model
# Performance claims remain based on held-out analyses.
message("Fitting full-data deployment bundle ...")
full_models <- vector("list", 4L); names(full_models) <- POP4
full_centers <- matrix(NA_real_, 4L, distance_n_pc, dimnames = list(POP4, dist_pc_cols))
full_scores <- matrix(NA_real_, nrow(dat), 4L, dimnames = list(NULL, POP4))
for (g in POP4) {
  idx <- which(dat$eval_group == g & dat$fine_eligible_single)
  if (length(idx) < 100L) idx <- which(dat$eval_group == g)
  mdl <- fit_combination(idx, raw_csx, Xcov, y, POP4)
  full_models[[g]] <- mdl
  full_centers[g, ] <- apply(pc_dist[idx, , drop = FALSE], 2L, median)
  full_scores[, g] <- apply_combination(mdl, raw_csx)
}
full_oth_model <- NULL; full_oth_score <- rep(NA_real_, nrow(dat))
if (sum(dat$eval_group == "OTH") >= 100L) {
  idx <- which(dat$eval_group == "OTH")
  full_oth_model <- fit_combination(idx, raw_csx, Xcov, y, POP4)
  full_oth_score <- apply_combination(full_oth_model, raw_csx)
}
full_transform <- residualize_and_scale(full_scores, pc_dist)
full_parts <- calculate_disco_weights(pc_dist, full_centers, A)
full_disco_score <- rowSums(full_transform$score * full_parts$weights)
full_nearest_i <- max.col(-full_parts$distances, ties.method = "first")
full_conventional <- full_scores[cbind(seq_len(nrow(dat)), full_nearest_i)]
for (g in POP4) {
  pos <- which(dat$eval_group == g); full_conventional[pos] <- full_scores[pos, g]
}
if (!is.null(full_oth_model)) {
  pos <- which(dat$eval_group %chin% c("OTH", "AMR")); full_conventional[pos] <- full_oth_score[pos]
}
probs <- seq(0, 1, length.out = 1001L)
calibration <- list()
for (g in c(report_groups, "ALL_GLOBAL")) {
  idx <- if (g == "ALL_GLOBAL") seq_len(nrow(dat)) else which(dat$eval_group == g)
  if (length(idx) < min_test_group) next
  cal <- list(n = length(idx), cov_colnames = colnames(Xcov), quantile_probs = probs, methods = list())
  for (method in c(M_CONV, M_DISCO)) {
    score <- if (method == M_CONV) full_conventional[idx] else full_disco_score[idx]
    fit <- fit_fast(cbind(Xcov[idx, , drop = FALSE], score = score), y[idx])
    cal$methods[[method]] <- list(coefficients = fit$coefficients,
                                  score_quantiles = quantile_safe(score, probs),
                                  prediction_quantiles = quantile_safe(fit$fitted, probs))
  }
  calibration[[g]] <- cal
}
deployment <- list(
  trait = trait, phi = phi, model_id = model_id,
  grouping_basis = grouping_basis, fine_groups = POP4, report_groups = report_groups,
  cov_pc_cols = cov_pc_cols, distance_pc_cols = dist_pc_cols, cov_terms = c("age", "sex", cov_pc_cols),
  csx_score_cols = csx_cols, sex_levels = levels(dat$sex), assignment_reference_centers = reference_centers,
  fine_tuning_centers = full_centers, disco_A = A, disco_d = full_parts$d,
  fine_models = full_models, oth_model = full_oth_model,
  disco_transform = list(coefficients = full_transform$coefficients,
                         residual_mean = full_transform$residual_mean, residual_sd = full_transform$residual_sd),
  calibration = calibration,
  note = "Full-data deployment fit; predictive performance must be taken from held-out repeated analyses"
)
saveRDS(deployment, file.path(outdir, "deployment_models.rds"), compress = "xz")
template_cols <- c("eid", "eval_group", "age", "sex", cov_pc_cols, csx_cols)
template <- as.data.table(setNames(replicate(length(template_cols), character(), simplify = FALSE), template_cols))
fwrite(template, file.path(outdir, "individual_input_template.tsv"), sep = "\t", quote = FALSE)
if (!is.null(first_predictions)) fwrite(first_predictions, file.path(outdir, "first_repeat_predictions.tsv.gz"), sep = "\t", quote = FALSE)


# 🚩 Workbook
literature_design <- data.table(
  analysis = c("PRS-CSx repeated validation/testing", "DiscoDivas empirical UKBB-style", "DiscoDivas center sensitivity"),
  implementation = c(
    paste0(prscsx_repeats, " random 50:50 splits within each ancestry among participants meeting probability >=", prscsx_prob_min, " when supplied; four PRS-CSx components standardized and combined in validation; tested with age, sex and PC1-PC", cov_n_pc, " covariates"),
    paste0(disco_repeats, " repeats; ", fine_tune_n, " fine-tuning participants from AFR/EAS/EUR/SAS and OTH when available; remaining participants tested; distance uses PC1-PC", distance_n_pc),
    "Primary uses medians of actual fine-tuning cohorts; sensitivity uses released 1KG four-population medians"
  ),
  limitation = c(
    "Fixed phi is prespecified rather than selected across a grid; discovery-target overlap must be excluded by the user",
    if (formal_ancestry) "Accuracy depends on the supplied ancestry labels and fine-tuning eligibility columns" else "Nearest-center fallback cannot reproduce OTH/AMR or high-confidence random-forest ancestry assignment",
    "The four input PRS models are PRS-CSx-based rather than the LDpred2 pipeline used in the DiscoDivas empirical paper"
  )
)
input_files <- rbindlist(list(file_info("phenotype_rds", phe_file), file_info("pca", pca_file),
                              file_info("reference_medians", med_file), file_info("csx_manifest", csx_manifest_file),
                              file_info("ancestry_file", ancestry_file), file_info("target_ids", target_ids_file),
                              file_info("exclude_ids", exclude_ids_file)), fill = TRUE)
csx_model_dir <- dirname(as.character(csx_manifest$file[1L]))
csx_root_dir <- dirname(csx_model_dir)
score_qc_file <- file.path(csx_model_dir, "score_qc.tsv")
chromosome_qc_file <- file.path(csx_model_dir, "chromosome_qc.tsv")
score_qc <- if (file.exists(score_qc_file)) fread(score_qc_file, showProgress = FALSE) else data.table()
chromosome_qc <- if (file.exists(chromosome_qc_file)) fread(chromosome_qc_file, showProgress = FALSE) else data.table()
sumstats_qc <- rbindlist(lapply(POP4, function(g) {
  f <- file.path(csx_root_dir, "input", g, paste0(g, ".sumstats_qc.tsv"))
  if (file.exists(f)) fread(f, showProgress = FALSE) else data.table()
}), fill = TRUE)
qc_overview <- rbindlist(list(
  data.table(check = "complete_participants", value = nrow(dat), note = "Complete phenotype, age, sex, PC1-PC20 and all four PRS components"),
  data.table(check = "formal_ancestry_available", value = as.numeric(formal_ancestry), note = grouping_basis),
  data.table(check = "fine_tuning_probability_used", value = as.numeric(has_prob), note = ancestry_prob_col),
  data.table(check = "self_report_match_used", value = as.numeric(require_self_match && has_self), note = self_report_col),
  data.table(check = "discovery_target_overlap_exclusion_file", value = as.numeric(nzchar(exclude_ids_file)), note = exclude_ids_file),
  group_counts[, .(check = paste0("group_", eval_group), value = N,
                   note = paste0("eligible_single=", eligible_single, "; eligible_oth=", eligible_oth))]
), fill = TRUE)

wb <- createWorkbook()
header_style <- createStyle(fontColour = "white", fgFill = "#1F4E78", textDecoration = "bold", halign = "center")
write_sheet <- function(name, z) {
  name <- substr(name, 1L, 31L); addWorksheet(wb, name)
  if (is.null(z) || !nrow(z)) { writeData(wb, name, "No rows"); return(invisible(NULL)) }
  writeDataTable(wb, name, z, tableStyle = "TableStyleMedium2")
  addStyle(wb, name, header_style, rows = 1, cols = seq_len(ncol(z)), gridExpand = TRUE)
  freezePane(wb, name, firstRow = TRUE); setColWidths(wb, name, cols = seq_len(ncol(z)), widths = "auto")
}
write_sheet("Run_info", run_info)
write_sheet("Literature_design", literature_design)
write_sheet("Group_counts", group_counts)
write_sheet("Ancestry_QC", ancestry_qc)
write_sheet("Ancestry_vs_nearest", ancestry_confusion)
write_sheet("Component_correlations", component_correlations)
write_sheet("PRSCSx_main", prscsx_main)
write_sheet("PRSCSx_split", prscsx_metrics)
write_sheet("PRSCSx_paired", prscsx_paired)
write_sheet("PRSCSx_weights", prscsx_weight_summary)
write_sheet("PRSCSx_components", prscsx_components)
write_sheet("Disco_main", disco_main)
write_sheet("Disco_split", disco_metrics)
write_sheet("Disco_paired", disco_paired)
write_sheet("Center_sensitivity", sensitivity_main)
write_sheet("R2_decomposition", r2_decomposition)
write_sheet("Transfer_matrix", transfer_summary)
write_sheet("Transfer_raw", transfer_metrics)
write_sheet("Interpolation_weights", disco_weight_summary)
write_sheet("Interpolation_raw", disco_weights)
write_sheet("Continuum", continuum_summary)
write_sheet("Calibration", calibration_summary)
write_sheet("Fine_tuning_centers", center_summary)
write_sheet("Disco_model_weights", disco_model_weight_summary)
write_sheet("Disco_QC", disco_qc)
write_sheet("QC_overview", qc_overview)
write_sheet("CSx_score_QC", score_qc)
write_sheet("CSx_chromosome_QC", chromosome_qc)
write_sheet("Sumstats_QC", sumstats_qc)
write_sheet("Input_files", input_files)
saveWorkbook(wb, file.path(outdir, "evaluation_results.xlsx"), overwrite = TRUE)


# 🚩 Figures
ord <- ALL_GROUPS[ALL_GROUPS %chin% report_groups]
theme_eval <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")

if (nrow(prscsx_main)) {
  p1a <- ggplot(prscsx_main, aes(factor(eval_group, levels = ord), mean_partial_r2, fill = method)) +
    geom_col(position = position_dodge(.75), width = .68) +
    geom_errorbar(aes(ymin = q025_partial_r2, ymax = q975_partial_r2), position = position_dodge(.75), width = .16) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "A. Held-out partial R2", subtitle = paste0(prscsx_repeats, " repeated ancestry-specific 50:50 validation/testing splits"), x = NULL, y = "Partial R2", fill = NULL) + theme_eval
  p1b <- ggplot(prscsx_weights, aes(`repeat`, standardized_coefficient, group = component)) +
    geom_line(alpha = .35) + facet_wrap(~eval_group, scales = "free_y") +
    labs(title = "B. PRS-CSx combination-coefficient stability", x = "Repeat", y = "Coefficient") + theme_eval
  p1c <- ggplot(prscsx_weight_summary, aes(component, mean_coefficient, fill = component)) +
    geom_col(show.legend = FALSE) + geom_errorbar(aes(ymin = q025, ymax = q975), width = .16) +
    facet_wrap(~eval_group, scales = "free_y") + labs(title = "C. Mean validation-tuned coefficients", x = NULL, y = "Coefficient") + theme_eval
  p1d <- ggplot(prscsx_components, aes(component, validation_partial_r2, fill = component)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) + facet_wrap(~eval_group, scales = "free_y") +
    scale_y_continuous(labels = scales::percent) + labs(title = "D. Single-component validation performance", x = NULL, y = "Partial R2") + theme_eval
  fig1 <- (p1a | p1b) / (p1c | p1d) + plot_annotation(title = paste0(trait, ": PRS-CSx paper-style repeated validation/testing evaluation"))
  ggsave(file.path(outdir, "Fig1.PRSCSx_evaluation.png"), fig1, width = 15, height = 10, dpi = 220, bg = "white")
}

if (nrow(disco_main)) {
  paired_plot <- copy(disco_paired); paired_plot[, eval_group := factor(eval_group, levels = ord)]
  p2a <- ggplot(disco_main, aes(factor(eval_group, levels = ord), mean_partial_r2, fill = method)) +
    geom_col(position = position_dodge(.75), width = .68) +
    geom_errorbar(aes(ymin = q025_partial_r2, ymax = q975_partial_r2), position = position_dodge(.75), width = .16) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "A. Held-out partial R2", subtitle = paste0("", fine_tune_n, " participants per fine-tuning cohort; remaining participants tested"), x = NULL, y = "Partial R2", fill = NULL) + theme_eval
  p2b <- ggplot(paired_plot, aes(eval_group, mean_difference)) + geom_col() +
    geom_errorbar(aes(ymin = q025_difference, ymax = q975_difference), width = .16) + geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(labels = scales::percent) + labs(title = "B. Paired absolute R2 difference", subtitle = "DiscoDivas minus conventional", x = NULL, y = "Difference in partial R2") + theme_eval
  p2c <- ggplot(paired_plot, aes(eval_group, mean_relative_change)) + geom_col() +
    geom_errorbar(aes(ymin = q025_relative_change, ymax = q975_relative_change), width = .16) + geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(labels = scales::percent) + labs(title = "C. Relative R2 change", x = NULL, y = "Relative change") + theme_eval
  decomp_plot <- r2_decomposition[design == "DiscoDivas_1300" & method %chin% c(M_CONV, M_DISCO)]
  decomp_long <- melt(decomp_plot, id.vars = c("eval_group", "method"),
                      measure.vars = c("mean_covariate_r2", "mean_delta_total_r2", "mean_full_model_r2"),
                      variable.name = "metric", value.name = "R2")
  p2d <- ggplot(decomp_long, aes(factor(eval_group, levels = ord), R2, fill = metric)) +
    geom_col(position = "dodge") + facet_wrap(~method) + scale_y_continuous(labels = scales::percent) +
    labs(title = "D. R2 decomposition", x = NULL, y = "R2", fill = NULL) + theme_eval
  fig2 <- (p2a | p2b) / (p2c | p2d) + plot_annotation(title = paste0(trait, ": DiscoDivas paper-aligned primary comparison"))
  ggsave(file.path(outdir, "Fig2.DiscoDivas_primary.png"), fig2, width = 15, height = 10, dpi = 220, bg = "white")
}

if (nrow(transfer_summary)) {
  p3a <- ggplot(transfer_summary, aes(model, factor(eval_group, levels = rev(ord)), fill = mean_partial_r2)) +
    geom_tile() + geom_text(aes(label = sprintf("%.1f%%", 100 * mean_partial_r2)), size = 3) +
    labs(title = "A. Model-transfer matrix", subtitle = "Rows: testing groups; columns: fine-tuning model", x = NULL, y = NULL, fill = "Partial R2") +
    theme_eval + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  p3b <- ggplot(disco_model_weight_summary, aes(component, mean_coefficient, fill = component)) +
    geom_col(show.legend = FALSE) + geom_errorbar(aes(ymin = q025, ymax = q975), width = .16) +
    facet_wrap(~fine_tuning_group, scales = "free_y") +
    labs(title = "B. Fine-tuned model coefficients", x = NULL, y = "Standardized coefficient") + theme_eval
  p3c <- ggplot(disco_model_weights, aes(`repeat`, validation_partial_r2, group = fine_tuning_group)) +
    geom_line() + facet_wrap(~fine_tuning_group, scales = "free_y") + scale_y_continuous(labels = scales::percent) +
    labs(title = "C. Fine-tuning performance stability", x = "Repeat", y = "Validation partial R2") + theme_eval
  p3d <- ggplot(group_counts[eval_group %chin% ord], aes(factor(eval_group, levels = ord), N)) + geom_col() +
    geom_text(aes(label = format(N, big.mark = ",")), vjust = -.2, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, .15))) + labs(title = "D. Testing-group sample size", x = NULL, y = "Participants") + theme_eval
  fig3 <- (p3a | p3b) / (p3c | p3d) + plot_annotation(title = "Why performance differs across testing ancestries")
  ggsave(file.path(outdir, "Fig3.model_transfer.png"), fig3, width = 15, height = 10, dpi = 220, bg = "white")
}

if (nrow(disco_weight_summary)) {
  p4a <- ggplot(disco_weight_summary[center_source == "actual fine-tuning cohorts"],
                 aes(component, mean_weight, fill = component)) + geom_col(show.legend = FALSE) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = .16) + facet_wrap(~eval_group) +
    labs(title = "A. Participant-specific interpolation weights", subtitle = "Actual fine-tuning-cohort medians", x = NULL, y = "Mean weight") + theme_eval
  p4b <- ggplot(disco_weight_summary, aes(component, mean_weight, fill = center_source)) +
    geom_col(position = "dodge") + facet_wrap(~eval_group) +
    labs(title = "B. Actual versus 1KG approximate centers", x = NULL, y = "Mean weight", fill = NULL) + theme_eval
  sens_all <- rbindlist(list(disco_main[, .(eval_group, method, mean_partial_r2, q025_partial_r2, q975_partial_r2)],
                             sensitivity_main[, .(eval_group, method, mean_partial_r2, q025_partial_r2, q975_partial_r2)]), fill = TRUE)
  p4c <- ggplot(sens_all, aes(factor(eval_group, levels = ord), mean_partial_r2, fill = method)) +
    geom_col(position = position_dodge(.8), width = .72) +
    geom_errorbar(aes(ymin = q025_partial_r2, ymax = q975_partial_r2), position = position_dodge(.8), width = .14) +
    scale_y_continuous(labels = scales::percent) + labs(title = "C. Center and comparator sensitivity", x = NULL, y = "Partial R2", fill = NULL) + theme_eval
  p4d <- ggplot(disco_qc, aes(`repeat`, mean_weight_entropy)) + geom_line() + geom_point() +
    labs(title = "D. Interpolation entropy stability", x = "Repeat", y = "Mean normalized entropy") + theme_eval
  fig4 <- (p4a | p4b) / (p4c | p4d) + plot_annotation(title = "DiscoDivas interpolation diagnostics")
  ggsave(file.path(outdir, "Fig4.interpolation_diagnostics.png"), fig4, width = 15, height = 10, dpi = 220, bg = "white")
}

if (nrow(continuum_summary)) {
  diff_cont <- dcast(continuum_summary[method %chin% c(M_CONV, M_DISCO)],
                     eval_group + axis + decile ~ method, value.var = "mean_partial_r2")
  if (all(c(M_CONV, M_DISCO) %in% names(diff_cont))) diff_cont[, difference := get(M_DISCO) - get(M_CONV)]
  p5a <- ggplot(continuum_summary[axis == "nearest_distance"], aes(decile, mean_partial_r2, group = method, shape = method)) +
    geom_line() + geom_point() + facet_wrap(~eval_group, scales = "free_y") + scale_y_continuous(labels = scales::percent) +
    labs(title = "A. Performance by within-group distance decile", x = "Nearest-center distance decile", y = "Partial R2", shape = NULL) + theme_eval
  p5b <- ggplot(diff_cont[axis == "nearest_distance"], aes(decile, difference)) + geom_line() + geom_point() +
    geom_hline(yintercept = 0, linetype = 2) + facet_wrap(~eval_group, scales = "free_y") + scale_y_continuous(labels = scales::percent) +
    labs(title = "B. DiscoDivas advantage by distance", x = "Distance decile", y = "R2 difference") + theme_eval
  p5c <- ggplot(continuum_summary[axis == "weight_entropy"], aes(decile, mean_partial_r2, group = method, shape = method)) +
    geom_line() + geom_point() + facet_wrap(~eval_group, scales = "free_y") + scale_y_continuous(labels = scales::percent) +
    labs(title = "C. Performance by ancestry-weight entropy", x = "Entropy decile", y = "Partial R2", shape = NULL) + theme_eval
  p5d <- ggplot(diff_cont[axis == "weight_entropy"], aes(decile, difference)) + geom_line() + geom_point() +
    geom_hline(yintercept = 0, linetype = 2) + facet_wrap(~eval_group, scales = "free_y") + scale_y_continuous(labels = scales::percent) +
    labs(title = "D. DiscoDivas advantage by entropy", x = "Entropy decile", y = "R2 difference") + theme_eval
  fig5 <- (p5a | p5b) / (p5c | p5d) + plot_annotation(title = "Performance along the genetic-ancestry continuum")
  ggsave(file.path(outdir, "Fig5.continuum.png"), fig5, width = 15, height = 11, dpi = 220, bg = "white")
}

if (nrow(calibration_summary)) {
  p6a <- ggplot(calibration_summary, aes(mean_score_residual, mean_phenotype_residual, group = method, shape = method)) +
    geom_line() + geom_point() + facet_wrap(~eval_group, scales = "free") +
    labs(title = "A. Covariate-adjusted calibration by score decile", x = "Mean PRS residual", y = "Mean phenotype residual", shape = NULL) + theme_eval
  p6b <- ggplot(calibration_summary, aes(score_decile, mean_phenotype_residual, group = method, shape = method)) +
    geom_line() + geom_point() + facet_wrap(~eval_group, scales = "free_y") +
    labs(title = "B. Phenotype residual across PRS deciles", x = "PRS decile", y = "Mean phenotype residual", shape = NULL) + theme_eval
  stability <- disco_metrics[eval_group %chin% report_groups]
  p6c <- ggplot(stability, aes(`repeat`, partial_r2, group = method, shape = method)) + geom_line(alpha = .6) +
    facet_wrap(~eval_group, scales = "free_y") + scale_y_continuous(labels = scales::percent) +
    labs(title = "C. Testing-set performance stability", x = "Repeat", y = "Partial R2", shape = NULL) + theme_eval
  p6d <- ggplot(disco_qc, aes(`repeat`, actual_center_condition_number)) + geom_line() + geom_point() +
    labs(title = "D. Fine-tuning-center distance-matrix conditioning", x = "Repeat", y = "Condition number") + theme_eval
  fig6 <- (p6a | p6b) / (p6c | p6d) + plot_annotation(title = "Calibration and stability diagnostics")
  ggsave(file.path(outdir, "Fig6.calibration_stability.png"), fig6, width = 15, height = 10, dpi = 220, bg = "white")
}

if (nrow(ancestry_qc) && nrow(component_correlations)) {
  p7a <- ggplot(ancestry_confusion, aes(closest_reference, factor(eval_group, levels = rev(ord)), fill = percent_within_group)) +
    geom_tile() + geom_text(aes(label = sprintf("%.1f%%", percent_within_group)), size = 3) +
    labs(title = "A. Supplied group versus nearest 1KG center", x = "Nearest 1KG center", y = "Evaluation group", fill = "% within group") + theme_eval
  prob_plot <- copy(ancestry_qc); prob_plot[, eval_group := factor(eval_group, levels = ord)]
  if (has_prob) {
    p7b <- ggplot(prob_plot, aes(eval_group, ancestry_probability_median)) + geom_point(size = 2) +
      geom_errorbar(aes(ymin = ancestry_probability_q025, ymax = ancestry_probability_q975), width = .15) +
      labs(title = "B. Ancestry-assignment probability", x = NULL, y = "Probability: median and 2.5%-97.5%") + theme_eval
  } else {
    p7b <- ggplot(prob_plot, aes(eval_group, distance_median)) + geom_point(size = 2) +
      geom_errorbar(aes(ymin = distance_q025, ymax = distance_q975), width = .15) +
      labs(title = "B. Nearest-center distance", subtitle = "No ancestry-probability column supplied", x = NULL, y = "Distance: median and 2.5%-97.5%") + theme_eval
  }
  p7c <- ggplot(component_correlations, aes(component1, component2, fill = correlation)) +
    geom_tile() + geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.5) + facet_wrap(~eval_group) +
    labs(title = "C. Correlation among raw PRS-CSx components", x = NULL, y = NULL, fill = "r") +
    theme_eval + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p7d <- ggplot(group_counts[eval_group %chin% ord], aes(factor(eval_group, levels = ord), phenotype_sd)) + geom_col() +
    geom_text(aes(label = sprintf("%.2f", phenotype_sd)), vjust = -.2, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, .15))) +
    labs(title = "D. Phenotype standard deviation by group", x = NULL, y = paste0("SD of ", trait)) + theme_eval
  fig7 <- (p7a | p7b) / (p7c | p7d) + plot_annotation(title = "Ancestry assignment and input-data diagnostics")
  ggsave(file.path(outdir, "Fig7.ancestry_data_QC.png"), fig7, width = 15, height = 11, dpi = 220, bg = "white")
}

readme <- c(
  paste0("Trait: ", trait), paste0("Complete N: ", nrow(dat)),
  paste0("Grouping: ", grouping_basis), paste0("Covariates: age, sex, PC1-PC", cov_n_pc),
  paste0("DiscoDivas distance and PRS residualization PCs: PC1-PC", distance_n_pc), "",
  "Primary performance metric:",
  "  partial R2 = (SSE_covariates - SSE_covariates+PRS) / SSE_covariates in held-out testing data.",
  "  delta_total_R2 and full-model R2 are also reported to avoid misinterpreting partial R2.", "",
  "Analyses:",
  paste0("  1. PRS-CSx: ", prscsx_repeats, " random 50:50 validation/testing splits within each ancestry group."),
  paste0("  2. DiscoDivas: ", disco_repeats, " repeats with ", fine_tune_n, " fine-tuning participants from AFR/EAS/EUR/SAS; remaining participants tested."),
  "  3. If OTH is available, an OTH-fine-tuned conventional model is used for OTH and AMR, matching the paper design.",
  "  4. Primary DiscoDivas centers are medians of actual fine-tuning cohorts; released 1KG medians are a sensitivity analysis.", "",
  "Critical interpretation:",
  "  A supplied high-confidence genetic-ancestry file is required for a paper-level AFR/EAS/EUR/SAS/AMR/OTH analysis.",
  "  Without it, nearest-center groups are fallback exploratory groups and OTH/AMR cannot be evaluated.",
  "  Discovery GWAS participants and their close relatives must be excluded from target evaluation using --exclude-ids.",
  "  deployment_models.rds is a full-data application bundle and is not used for performance claims."
)
writeLines(readme, file.path(outdir, "README.txt"))
cat("Evaluation complete\n")
cat("  output:", normalizePath(outdir), "\n")
cat("  grouping:", grouping_basis, "\n")
cat("  complete participants:", nrow(dat), "\n")
