# Use every member of the analysis cohort. AOD's patient file is case-enriched.
args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)
source(file.path(args[4], "phenotype.R"), local = TRUE)
d <- as.data.frame(readRDS(args[1]))
trait <- args[2]
stopifnot(all(c("eid", "birth_date", "date_attend") %in% names(d)), !anyNA(d$eid), !anyDuplicated(d$eid))
fields <- c("birth_date", "date_attend", "date_lost", "date_death")
for (src in c("icd10", "ref")) fields <- c(fields, paste0(c("fod_", "cnt_"), src, "_", trait))
x <- d[, c("eid", intersect(c(fields,"age","sex","drug.lipid","bmi","sbp","dbp","bb_TC","bb_HDL","bb_HBA1C"), names(d))), drop = FALSE]
launch <- commandArgs(FALSE)
entry <- sub("^--file=", "", grep("^--file=",launch,value=TRUE)[1])
if (is.na(entry)) entry <- launch[match("-f",launch)+1L]
source(file.path(dirname(normalizePath(entry)), "baseline_medication.R"), local=TRUE)
phe_root <- dirname(dirname(normalizePath(args[1])))
med <- le8_baseline_med_values(le8_attach_baseline_medication(d, root=phe_root))
x$baseline_lipid_i0 <- med$values[, "drug.lipid"]
for (nm in intersect(fields[!startsWith(fields, "cnt_")], names(x))) {
    raw <- x[[nm]]
    x[[nm]] <- .phenotype_date(raw)
    if (startsWith(nm, "fod_")) x[[paste0("invalid_", nm)]] <- !is.na(raw) & is.na(x[[nm]])
}
data.table::fwrite(x, file.path(args[3], "cohort.tsv.gz"), sep = "\t", na = "", compress = "gzip")
writeLines(as.character(date_follow_end), file.path(args[3], "followup.txt"))
