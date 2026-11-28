source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
options(uvmr.skip_main = TRUE)
source(file.path(Sys.getenv("FLT_CODE"), "s2_1_uvmr.R"))
library(data.table)
source_dir <- file.path(Sys.getenv("FLT_TMP"), "uvmr", Sys.getenv("FLT_PANEL"))
dirs <- file.path(source_dir, c(paste0("forward-", flt_panel()$id), "reverse"))
planned <- lapply(dirs, function(d) fread(file.path(d, "planned_pairs.csv")))
offsets <- c(0L, head(cumsum(lengths(lapply(planned, function(d) d$pair_id))), -1L))
read_parts <- function(name) {
    rbindlist(Map(function(d, offset) {
        x <- fread(file.path(d, paste0(name, ".csv")))
        if ("pair_id" %in% names(x)) x[, pair_id := pair_id + offset]
        x
    }, dirs, offsets), use.names = TRUE, fill = TRUE)
}
pairs <- read_parts("planned_pairs")
stopifnot(!anyDuplicated(pairs[, .(direction, exposure_key, outcome_key)]))
out <- flt_run("uvmr")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (name in c("planned_pairs", "steiger_trait_info_used", "uvmr_results", "uvmr_hits",
               "uvmr_snp_counts", "uvmr_failures", "uvmr_steiger", "uvmr_main_table",
               "uvmr_sensitivity_table", "uvmr_heterogeneity", "uvmr_egger",
               "uvmr_issues", "uvmr_leave_one_out")) {
    d <- read_parts(name)
    if (name == "steiger_trait_info_used") d <- unique(d, by = "key")
    fwrite(d, file.path(out, paste0(name, ".csv")))
}
raw <- fread(file.path(out, "uvmr_results.csv"))
counts <- fread(file.path(out, "uvmr_snp_counts.csv"))
fwrite(compare_directions(raw, counts, pairs), file.path(out, "uvmr_bidirectional.csv"))
