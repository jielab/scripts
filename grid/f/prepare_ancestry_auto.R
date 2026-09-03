#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value after ", key, call. = FALSE)
  args[[i + 1L]]
}
as_bool <- function(x) toupper(x) %in% c("TRUE", "T", "1", "YES", "Y")

phe_file <- get_arg("--phe")
pca_file <- get_arg("--pca")
med_file <- get_arg("--med")
out_file <- get_arg("--out")
outdir <- get_arg("--outdir", dirname(out_file))
ethnicity_requested <- get_arg("--ethnicity-col", "auto")
n_pc <- as.integer(get_arg("--n-pc", "20"))
distance_n_pc <- as.integer(get_arg("--distance-pcs", "10"))
prob_threshold <- as.numeric(get_arg("--prob-threshold", "0.90"))
fine_prob_threshold <- as.numeric(get_arg("--fine-prob-threshold", "0.999999"))
anchor_max <- as.integer(get_arg("--anchor-max-per-group", "10000"))
anchor_margin_quantile <- as.numeric(get_arg("--anchor-margin-quantile", "0.25"))
seed <- as.integer(get_arg("--seed", "20260803"))
allow_geometry_fallback <- as_bool(get_arg("--allow-geometry-fallback", "TRUE"))

if (any(vapply(list(phe_file, pca_file, med_file, out_file), is.null, logical(1)))) {
  stop("Required: --phe FILE --pca FILE --med FILE --out FILE", call. = FALSE)
}
for (f in c(phe_file, pca_file, med_file)) if (!file.exists(f)) stop("Missing input: ", f, call. = FALSE)
if (!is.finite(n_pc) || n_pc < 5L || !is.finite(distance_n_pc) || distance_n_pc < 5L || distance_n_pc > n_pc) {
  stop("Invalid --n-pc/--distance-pcs", call. = FALSE)
}
if (!is.finite(prob_threshold) || prob_threshold <= 0.5 || prob_threshold >= 1) stop("--prob-threshold must be in (0.5,1)", call. = FALSE)
if (!is.finite(fine_prob_threshold) || fine_prob_threshold < prob_threshold || fine_prob_threshold > 1) stop("Invalid --fine-prob-threshold", call. = FALSE)
if (!is.finite(anchor_max) || anchor_max < 100L) stop("--anchor-max-per-group must be >=100", call. = FALSE)
if (!is.finite(anchor_margin_quantile) || anchor_margin_quantile < 0 || anchor_margin_quantile >= 0.9) stop("Invalid --anchor-margin-quantile", call. = FALSE)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)
POP4 <- c("AFR", "EAS", "EUR", "SAS")

map_ukb_ethnicity <- function(x) {
  raw <- trimws(as.character(x))
  z <- toupper(raw)
  z <- gsub("[_-]", " ", z)
  code <- suppressWarnings(as.integer(sub("[.]0+$", "", z)))
  out <- rep(NA_character_, length(z))

  out[code %in% c(1L, 1001L, 1002L, 1003L)] <- "EUR"
  out[code %in% c(2L, 2001L, 2002L, 2003L, 2004L)] <- "OTH"
  out[code %in% c(3L, 3001L, 3002L, 3003L, 3004L)] <- "SAS"
  out[code %in% c(4L, 4001L, 4002L, 4003L)] <- "AFR"
  out[code %in% 5L] <- "EAS"
  out[code %in% 6L] <- "OTH"

  set_missing <- function(pattern, value) {
    hit <- is.na(out) & grepl(pattern, z)
    out[hit] <<- value
  }
  set_missing("^(EUR|EUROPEAN|WHITE|BRITISH|IRISH)", "EUR")
  set_missing("^(AFR|AFRICAN|BLACK|CARIBBEAN)", "AFR")
  set_missing("^(EAS|EAST ASIAN|CHINESE|JAPANESE|KOREAN)", "EAS")
  set_missing("^(SAS|SOUTH ASIAN|ASIAN|INDIAN|PAKISTANI|BANGLADESHI)", "SAS")
  set_missing("^(AMR|ADMIXED AMERICAN|LATIN|LATINO|HISPANIC|MIXED|OTHER|OTH|UNKNOWN|UNCLASS)", "OTH")
  out[z %in% c("", "NA", "NAN", "NULL", "PREFER NOT TO ANSWER", "DO NOT KNOW", "-1", "-3")] <- NA_character_
  out
}

find_ethnicity_column <- function(d, requested = "auto") {
  if (!identical(tolower(requested), "auto")) {
    if (!requested %in% names(d)) stop("Requested ethnicity column is absent: ", requested, call. = FALSE)
    return(requested)
  }
  preferred <- c("ethnicity", "ethnic.c", "ethnic", "self_report_ancestry", "self_report_ethnicity",
                 "p21000", "p21000_i0", "p21000_a1", "f21000_0_0", "field21000")
  regex_hits <- grep("(^|[._])ethnic|21000", names(d), value = TRUE, ignore.case = TRUE)
  candidates <- unique(c(preferred[preferred %in% names(d)], regex_hits))
  if (!length(candidates)) return(NA_character_)
  score <- vapply(candidates, function(v) sum(!is.na(map_ukb_ethnicity(d[[v]]))), numeric(1))
  candidates[which.max(score)]
}

row_euclidean <- function(X, centers) {
  X <- as.matrix(X); centers <- as.matrix(centers)
  out <- vapply(seq_len(nrow(centers)), function(j) sqrt(rowSums(sweep(X, 2L, centers[j, ], FUN = "-")^2)), numeric(nrow(X)))
  if (is.null(dim(out))) out <- matrix(out, ncol = nrow(centers))
  colnames(out) <- rownames(centers)
  out
}

geometry_posterior <- function(dist_mat, anchors = NULL) {
  powers <- seq(1, 20, by = 0.5)
  chosen <- 4
  if (!is.null(anchors) && nrow(anchors)) {
    quality <- vapply(powers, function(p) {
      inv <- 1 / pmax(dist_mat[anchors$row_id, , drop = FALSE], sqrt(.Machine$double.eps))^p
      pr <- inv / rowSums(inv)
      true_j <- match(anchors$group, colnames(pr))
      true_p <- pr[cbind(seq_len(nrow(pr)), true_j)]
      abs(unname(quantile(true_p, 0.05, na.rm = TRUE)) - 0.90)
    }, numeric(1))
    chosen <- powers[which.min(quality)]
  }
  inv <- 1 / pmax(dist_mat, sqrt(.Machine$double.eps))^chosen
  pr <- inv / rowSums(inv)
  list(posterior = pr, power = chosen)
}

message("Reading phenotype and automatically locating UKB self-reported ethnicity ...")
phe0 <- readRDS(phe_file)
if (!"eid" %in% names(phe0)) stop("Phenotype RDS lacks eid", call. = FALSE)
ethnicity_col <- find_ethnicity_column(phe0, ethnicity_requested)
if (is.na(ethnicity_col)) {
  message("No usable ethnicity column was found; ancestry will use a geometry-only fallback.")
  phe <- data.table(eid = as.character(phe0$eid), self_report_ancestry = NA_character_)
} else {
  message("Using self-reported ethnicity column: ", ethnicity_col)
  phe <- data.table(eid = as.character(phe0$eid), self_report_ancestry = map_ukb_ethnicity(phe0[[ethnicity_col]]))
}
if (anyDuplicated(phe$eid)) stop("Duplicate eid in phenotype RDS", call. = FALSE)
rm(phe0); invisible(gc())

message("Reading projected PCs and released 1KG centers ...")
pca <- fread(pca_file, showProgress = FALSE)
id_col <- intersect(c("IID", "#IID", "eid"), names(pca))[1L]
if (length(id_col) != 1L || is.na(id_col)) stop("PCA file needs IID/#IID/eid", call. = FALSE)
setnames(pca, id_col, "eid"); pca[, eid := as.character(eid)]
pc_cols <- paste0("PC", seq_len(n_pc)); dist_pc_cols <- paste0("PC", seq_len(distance_n_pc))
if (!all(pc_cols %in% names(pca))) stop("PCA file lacks: ", paste(setdiff(pc_cols, names(pca)), collapse = ", "), call. = FALSE)
if (anyDuplicated(pca$eid)) stop("Duplicate eid in PCA file", call. = FALSE)

med <- fread(med_file, showProgress = FALSE)
pop_col <- names(med)[1L]
if (!all(dist_pc_cols %in% names(med))) stop("Median file lacks distance PCs", call. = FALSE)
med[, (pop_col) := toupper(as.character(get(pop_col)))]
if (!all(POP4 %in% med[[pop_col]])) stop("Median file must contain AFR,EAS,EUR,SAS", call. = FALSE)
centers <- as.matrix(med[match(POP4, get(pop_col)), ..dist_pc_cols]); storage.mode(centers) <- "double"; rownames(centers) <- POP4

dat <- merge(pca[, c("eid", pc_cols), with = FALSE], phe, by = "eid", all.x = TRUE, sort = FALSE)
pc_dist <- as.matrix(dat[, ..dist_pc_cols]); storage.mode(pc_dist) <- "double"
dist_mat <- row_euclidean(pc_dist, centers)
nearest_j <- max.col(-dist_mat, ties.method = "first")
nearest <- POP4[nearest_j]
sorted_dist <- t(apply(dist_mat, 1L, sort, partial = 1:2))
nearest_distance <- dist_mat[cbind(seq_len(nrow(dist_mat)), nearest_j)]
second_distance <- sorted_dist[, 2L]
margin <- second_distance - nearest_distance
margin_ratio <- margin / pmax(second_distance, sqrt(.Machine$double.eps))

dat[, `:=`(nearest_1kg = nearest, nearest_1kg_distance = nearest_distance,
           second_minus_first_distance = margin, relative_distance_margin = margin_ratio, row_id = .I)]

anchor_candidates <- dat[self_report_ancestry %chin% POP4 & self_report_ancestry == nearest_1kg & is.finite(relative_distance_margin),
                         .(row_id, group = self_report_ancestry, relative_distance_margin)]
anchors <- rbindlist(lapply(POP4, function(g) {
  z <- anchor_candidates[group == g]
  if (!nrow(z)) return(z)
  cut <- unname(quantile(z$relative_distance_margin, anchor_margin_quantile, na.rm = TRUE))
  keep <- z[relative_distance_margin >= cut]
  if (nrow(keep) < 200L) keep <- z
  keep
}), use.names = TRUE)
anchor_counts_initial <- anchors[, .(available = .N), by = group]

classifier_method <- "geometry_only"
validation_confusion <- data.table()
validation_accuracy <- NA_real_
posterior <- NULL
model_note <- "No self-report anchors were available."

can_lda <- nrow(anchor_counts_initial) == 4L && min(anchor_counts_initial$available) >= 200L && requireNamespace("MASS", quietly = TRUE)
if (can_lda) {
  balanced_n <- min(anchor_max, min(anchor_counts_initial$available))
  balanced <- rbindlist(lapply(POP4, function(g) {
    z <- anchors[group == g]
    z[sample.int(nrow(z), balanced_n)]
  }))
  balanced[, fold_test := FALSE]
  balanced[, fold_test := seq_len(.N) <= max(50L, floor(.N * 0.2)), by = group]
  # Shuffle the test designation independently within group.
  balanced[, fold_test := sample(fold_test), by = group]
  # Use explicit predicates for compatibility with data.table >= 1.17, where a
  # lone symbol in i is resolved in calling scope instead of as a DT column.
  train_ids <- balanced[fold_test == FALSE, row_id]
  test_ids <- balanced[fold_test == TRUE, row_id]

  Xtrain <- as.matrix(dat[train_ids, ..pc_cols]); storage.mode(Xtrain) <- "double"
  mu <- colMeans(Xtrain); ss <- apply(Xtrain, 2L, sd)
  keep_pc <- is.finite(ss) & ss > 0
  Xtrain <- sweep(sweep(Xtrain[, keep_pc, drop = FALSE], 2L, mu[keep_pc], "-"), 2L, ss[keep_pc], "/")
  ytrain <- factor(dat$self_report_ancestry[train_ids], levels = POP4)
  fit0 <- tryCatch(MASS::lda(x = Xtrain, grouping = ytrain, prior = rep(0.25, 4), tol = 1e-8), error = function(e) e)

  if (!inherits(fit0, "error")) {
    Xtest <- as.matrix(dat[test_ids, ..pc_cols]); storage.mode(Xtest) <- "double"
    Xtest <- sweep(sweep(Xtest[, keep_pc, drop = FALSE], 2L, mu[keep_pc], "-"), 2L, ss[keep_pc], "/")
    pred_test <- predict(fit0, Xtest)
    validation_confusion <- data.table(truth = dat$self_report_ancestry[test_ids], predicted = as.character(pred_test$class))[, .N, by = .(truth, predicted)]
    validation_accuracy <- mean(as.character(pred_test$class) == dat$self_report_ancestry[test_ids])

    # Refit on all balanced anchors, then predict the entire cohort.
    all_ids <- balanced$row_id
    Xall_train <- as.matrix(dat[all_ids, ..pc_cols]); storage.mode(Xall_train) <- "double"
    mu <- colMeans(Xall_train); ss <- apply(Xall_train, 2L, sd); keep_pc <- is.finite(ss) & ss > 0
    Xall_train <- sweep(sweep(Xall_train[, keep_pc, drop = FALSE], 2L, mu[keep_pc], "-"), 2L, ss[keep_pc], "/")
    yall <- factor(dat$self_report_ancestry[all_ids], levels = POP4)
    fit <- tryCatch(MASS::lda(x = Xall_train, grouping = yall, prior = rep(0.25, 4), tol = 1e-8), error = function(e) e)
    if (!inherits(fit, "error")) {
      Xall <- as.matrix(dat[, ..pc_cols]); storage.mode(Xall) <- "double"
      Xall <- sweep(sweep(Xall[, keep_pc, drop = FALSE], 2L, mu[keep_pc], "-"), 2L, ss[keep_pc], "/")
      posterior <- predict(fit, Xall)$posterior
      posterior <- posterior[, POP4, drop = FALSE]
      classifier_method <- "balanced_LDA_from_self_report_and_1KG_concordant_anchors"
      model_note <- paste0("Balanced anchors per group=", balanced_n, "; holdout accuracy=", sprintf("%.4f", validation_accuracy),
                           "; PCs=", paste(pc_cols[keep_pc], collapse = ","))
    }
  } else {
    model_note <- paste0("LDA failed: ", conditionMessage(fit0))
  }
}

if (is.null(posterior)) {
  if (!allow_geometry_fallback) stop("Automatic ancestry model could not be trained and geometry fallback is disabled", call. = FALSE)
  gp <- geometry_posterior(dist_mat, if (nrow(anchors)) anchors else NULL)
  posterior <- gp$posterior[, POP4, drop = FALSE]
  classifier_method <- paste0("1KG_center_inverse_distance_power_", gp$power)
  model_note <- paste0(model_note, " Geometry-only fallback used.")
}

pred_j <- max.col(posterior, ties.method = "first")
pred_core <- POP4[pred_j]
max_prob <- posterior[cbind(seq_len(nrow(posterior)), pred_j)]
genetic_ancestry <- ifelse(max_prob >= prob_threshold, pred_core, "OTH")

out <- data.table(
  eid = dat$eid,
  genetic_ancestry = genetic_ancestry,
  ancestry_probability = max_prob,
  predicted_core_ancestry = pred_core,
  self_report_ancestry = dat$self_report_ancestry,
  nearest_1kg = dat$nearest_1kg,
  nearest_1kg_distance = dat$nearest_1kg_distance,
  second_minus_first_distance = dat$second_minus_first_distance,
  relative_distance_margin = dat$relative_distance_margin,
  classifier_method = classifier_method
)
for (g in POP4) out[, (paste0("prob_", g)) := posterior[, g]]
out[, ancestry_self_report_match := genetic_ancestry == self_report_ancestry]
out[, fine_tune_eligible := genetic_ancestry %chin% POP4 & ancestry_probability >= fine_prob_threshold &
                                  !is.na(self_report_ancestry) & genetic_ancestry == self_report_ancestry]
if (anyDuplicated(out$eid)) stop("Internal duplicate output IDs", call. = FALSE)
fwrite(out, out_file, sep = "\t", quote = FALSE)

counts <- out[, .(N = .N, percent = 100 * .N / nrow(out), mean_probability = mean(ancestry_probability),
                  median_probability = median(ancestry_probability), self_report_available = sum(!is.na(self_report_ancestry)),
                  self_report_match = sum(ancestry_self_report_match, na.rm = TRUE), fine_tune_eligible = sum(fine_tune_eligible)),
              by = genetic_ancestry]
counts[, order := match(genetic_ancestry, c(POP4, "OTH"))]; setorder(counts, order); counts[, order := NULL]
self_confusion <- out[!is.na(self_report_ancestry), .N, by = .(self_report_ancestry, genetic_ancestry)]
nearest_confusion <- out[, .N, by = .(nearest_1kg, genetic_ancestry)]
model_info <- data.table(
  item = c("method", "ethnicity_column", "probability_threshold", "fine_probability_threshold", "anchor_margin_quantile",
           "anchor_max_per_group", "holdout_accuracy", "n_pc", "distance_n_pc", "note"),
  value = c(classifier_method, ifelse(is.na(ethnicity_col), "NONE", ethnicity_col), prob_threshold, fine_prob_threshold,
            anchor_margin_quantile, anchor_max, validation_accuracy, n_pc, distance_n_pc, model_note)
)

wb <- createWorkbook()
write_sheet <- function(name, x) {
  addWorksheet(wb, name)
  if (!nrow(x)) { writeData(wb, name, "No rows"); return(invisible(NULL)) }
  writeDataTable(wb, name, x, tableStyle = "TableStyleMedium2")
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(x)), widths = "auto")
}
write_sheet("Model", model_info)
write_sheet("Assigned_counts", counts)
write_sheet("Anchor_counts", anchor_counts_initial)
write_sheet("Holdout_confusion", validation_confusion)
write_sheet("Selfreport_confusion", self_confusion)
write_sheet("Nearest1KG_confusion", nearest_confusion)
saveWorkbook(wb, file.path(outdir, "ancestry_auto_qc.xlsx"), overwrite = TRUE)

plot_theme <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
p1 <- ggplot(counts, aes(genetic_ancestry, N)) + geom_col() +
  geom_text(aes(label = format(N, big.mark = ",")), vjust = -0.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "A. Automatically assigned genetic ancestry", x = NULL, y = "Participants") + plot_theme

p2 <- ggplot(out, aes(genetic_ancestry, ancestry_probability)) + geom_boxplot(outlier.alpha = .08) +
  geom_hline(yintercept = prob_threshold, linetype = 2) +
  labs(title = "B. Assignment probability", subtitle = paste0("Dashed line: classification threshold ", prob_threshold), x = NULL, y = "Maximum posterior probability") + plot_theme

heat <- copy(self_confusion)
if (nrow(heat)) heat[, pct := 100 * N / sum(N), by = self_report_ancestry]
p3 <- if (nrow(heat)) {
  ggplot(heat, aes(self_report_ancestry, genetic_ancestry, fill = pct)) + geom_tile() +
    geom_text(aes(label = sprintf("%.1f", pct)), size = 3) +
    labs(title = "C. Self-report versus automatic assignment", x = "Self-reported ancestry", y = "Automatic genetic ancestry", fill = "% within self-report") + plot_theme
} else {
  ggplot() + annotate("text", x = 0, y = 0, label = "No self-reported ethnicity column found") +
    labs(title = "C. Self-report versus automatic assignment") + theme_void()
}

set.seed(seed + 1L)
idx <- if (nrow(dat) > 50000L) sample.int(nrow(dat), 50000L) else seq_len(nrow(dat))
plot_pc <- data.table(PC1 = dat$PC1[idx], PC2 = dat$PC2[idx], group = genetic_ancestry[idx])
center_dt <- data.table(group = POP4, PC1 = centers[, "PC1"], PC2 = centers[, "PC2"])
p4 <- ggplot(plot_pc, aes(PC1, PC2, shape = group)) + geom_point(alpha = .15, size = .45) +
  geom_point(data = center_dt, aes(PC1, PC2), inherit.aes = FALSE, shape = 21, fill = "white", size = 3, stroke = 1) +
  geom_text(data = center_dt, aes(PC1, PC2, label = group), inherit.aes = FALSE, vjust = -1, fontface = "bold") +
  labs(title = "D. PC projection and assigned groups", subtitle = "Random sample of up to 50,000 participants", shape = NULL) + plot_theme

fig <- (p1 | p2) / (p3 | p4) + plot_annotation(title = "Automatic UK Biobank ancestry construction")
ggsave(file.path(outdir, "Fig2.Ancestry_QC.png"), fig, width = 14, height = 10, dpi = 220, bg = "white")

cat("Automatic ancestry file:", out_file, "\n")
cat("Classifier:", classifier_method, "\n")
cat("Self-reported ethnicity column:", ifelse(is.na(ethnicity_col), "NONE", ethnicity_col), "\n")
cat("Holdout accuracy:", validation_accuracy, "\n")
print(counts)
