source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
library(dplyr)

project <- Sys.getenv("FLT_OUT")
input <- Sys.getenv("FLT_COHORT")
flt_file <- Sys.getenv("FLT_DICTIONARY")
outdir <- file.path(project, "01intake")

run_full_scan <- FALSE  # TRUE also scans item182. The three selected FLTs always run.
radar_score <- 9  # Exact FPQ score used to define each radar group.

liking_selected <- c("Butter on bread" = "p20628", "Total salad vegetable" = "F-Salad_vegetables", "Chicken" = "p20638")
# 1. WebQ and FPQ ----------------------------------------------------------
diet_fields <- list(
  `Alcoholic Beverages` = c(26067,26138,26151,26152,26153),
  `Non-alcoholic Beverages` = c(26081,26082,26095,26124,26126,26127,26141,26142,26148),
  Dairy = c(26084,26087,26096,26099,26102,26103,26131,26133,26150,26154),
  Eggs = 26088,
  `Fats and Oils` = c(26062,26063,26110,26111,26112),
  Fruits = 26089:26094,
  Vegetables = c(26065,26098,26123,26125,26143,26144,26146,26147,26101,26115,26118,26119,26120),
  `Soybean Products` = c(26086,26136,26137),
  `Nuts and Seeds` = 26106:26108,
   `Grains` = c(26071,26074,26076,26077,26105,26114,26072,26073,26075,26078,26079,26083,26097,26113,26128,26139),
  `Red and Processed Meat` = c(26066,26100,26104,26117,26122),
  Poultry = c(26069,26121),
  `Seafood and Fish` = c(26070,26109,26132,26149),
  `Meat Substitutes` = 26145,
  `Sweets and Desserts` = c(26064,26068,26080,26085,26140),
  `Savoury Foods and Condiments` = c(26116,26129,26130,26134,26135)
)
group_names <- names(diet_fields)
all_fields <- paste0("w.p", 26062:26154)
mapping <- data.frame(field_id = unlist(diet_fields, use.names = FALSE), food_consumption_group = rep(group_names, lengths(diet_fields)))
mapping$column <- paste0("w.p", mapping$field_id)

cfg <- readRDS(file.path(input, "analysis_config.rds"))
flt <- as.data.frame(readRDS(file.path(input, cfg$primary_data_file)))
diet <- as.data.frame(readRDS(file.path(Sys.getenv("FLT_YGB"), "Rdata/diet.fudan.rds")))
meta <- data.frame(name = unname(liking_selected), item_raw = names(liking_selected))

if (run_full_scan) {
  meta <- as.data.frame(readxl::read_excel(flt_file, sheet = "item182"))
  stopifnot(all(c("name", "item_raw") %in% names(meta)), !anyNA(meta$name), !anyDuplicated(meta$name), nrow(meta) > 0)
  extra <- !liking_selected %in% meta$name
  meta <- rbind(meta[c("name", "item_raw")], data.frame(name = unname(liking_selected[extra]), item_raw = names(liking_selected)[extra]))
}
require_columns <- function(x, columns, source) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop(source, " missing columns: ", paste(missing, collapse = ", "))
}
score_fields <- unique(c(unname(liking_selected), meta$name))
require_columns(flt, c("eid", score_fields), "FPQ")
require_columns(diet, c("eid", all_fields), "WebQ")
flt$eid <- as.character(flt$eid)
diet$eid <- as.character(diet$eid)
for (field in all_fields) {
  x <- diet[[field]]
  if (!is.numeric(x) || any(!is.na(x) & (!is.finite(x) | x < 0))) {stop("Invalid WebQ weight: ", field)}
}
for (field in score_fields) {
  x <- flt[[field]]
  if (!is.numeric(x)) stop("Non-numeric liking score: ", field)
  x[!is.finite(x)] <- NA_real_
  if (grepl("^p20[67][0-9]{2}$", field)) x[!x %in% 1:9] <- NA_real_
  flt[[field]] <- x
}
# Keep recorded zeros. A missing component makes its group total missing.
intake <- data.frame(eid = diet$eid)
for (group in group_names) {intake[[group]] <- rowSums(diet[paste0("w.p", diet_fields[[group]])])}
intake[all_fields] <- diet[all_fields]
dat <- flt %>% select(eid, all_of(score_fields)) %>% inner_join(intake, by = "eid", relationship = "one-to-one", na_matches = "never")
if (!nrow(dat)) stop("No matching FPQ and WebQ participant IDs.")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name) {write.csv(x, file.path(outdir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")}
write_table(mapping, "food_weight_mapping.csv")
write_table(meta, "food_liking_metadata.csv")
write_table(data.frame(stage = c("FPQ cohort", "WebQ file", "Matched IDs"), N = c(nrow(flt), nrow(diet), nrow(dat))), "sample_counts.csv")

# 2. Pairwise Spearman correlations -----------------------------------------
spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  rho <- p <- NA_real_
  if (n >= 3 && any(x != x[1]) && any(y != y[1])) {
    fit <- cor.test(x, y, method = "spearman", exact = FALSE)
    rho <- unname(fit$estimate)
    p <- fit$p.value
  }
  data.frame(n = n, cor = rho, p = p, zero_intake_pct = if (n) 100 * mean(y == 0) else NA_real_)
}
run_grid <- function(scores, fields) {
  pairs <- expand.grid(intake_field = fields, food_liking = unname(scores), stringsAsFactors = FALSE)
  pairs$food_liking_label <- names(scores)[match(pairs$food_liking, scores)]
  estimates <- Map(function(x, y) spearman(dat[[x]], dat[[y]]), pairs$food_liking, pairs$intake_field)
  results <- bind_cols(pairs, bind_rows(estimates))
  results$p_fdr <- p.adjust(results$p, method = "BH", n = nrow(results))
  results
}

# 2.1 Individual intake fields, including eggs and water.
item_liking <- if (run_full_scan) setNames(meta$name, meta$item_raw) else liking_selected
item_results <- run_grid(item_liking, all_fields)
item_results$food_consumption_group <- mapping$food_consumption_group[match(item_results$intake_field, mapping$column)]
item_results$food_consumption_group[item_results$intake_field == "w.p26088"] <- "Eggs (outside 17 groups)"
item_results$food_consumption_group[item_results$intake_field == "w.p26148"] <- "Water (outside 17 groups)"
selected_items <- item_results %>% filter(food_liking %in% liking_selected) %>% mutate(p_fdr_selected_items = p.adjust(p, method = "BH"))
write_table(selected_items, "spearman_selected_intake_items.csv")
# 2.2 Intake groups
group_results <- run_grid(liking_selected, group_names)
write_table(group_results, "spearman_selected_actual_liking.csv")
if (run_full_scan) {
  write_table(run_grid(item_liking, group_names), "spearman_all_actual_liking.csv")
  write_table(item_results, "spearman_all_intake_items.csv")
}
# 3. Descriptive intake profiles ------------------------------------------
# Use one exact liking score. Normalize the three group means within each axis.
# Follow the active classification automatically; no hard-coded category names.
radar_order <- names(diet_fields)
radar_profile <- bind_rows(lapply(names(liking_selected), function(label) {
  score <- dat[[liking_selected[[label]]]]
  grams <- dat[which(score == radar_score), radar_order, drop = FALSE]
  means <- colMeans(grams, na.rm = TRUE)
  means[is.nan(means)] <- NA_real_
  data.frame(food_liking = label, liking_score = radar_score, intake_group = radar_order, n_liking = nrow(grams), n_intake = colSums(!is.na(grams)), mean_g_day = means, row.names = NULL)
}))

min_max <- function(means) {
  if (anyNA(means) || diff(range(means)) == 0) return(rep(NA_real_, length(means)))
  (means - min(means)) / diff(range(means))
}
radar_profile$min_max <- ave(radar_profile$mean_g_day, radar_profile$intake_group, FUN = min_max)
write_table(radar_profile, "radar_project_intake.csv")

message("Results saved: ", outdir)
