source("f/ukb_followup.R")
raw <- readRDS("D:/data/ukb/phe/Rdata/diet0.rds")
window <- maha_webq_window(raw)
cat("Raw date rows:", nrow(window), "; undated recalls:", sum(window$webq_undated), "\n")
rm(raw); gc()
d <- readRDS("D:/data/ukb/phe/Rdata/all.rds")
d <- d[!is.na(d$ethnic.c) & d$ethnic.c == "White", ]
scores <- paste0("diet.", c("maha", "dash", "mind", "medi24"), ".sum")
d <- maha_set_landmark(d, window, scores)
cat("Landmark set. Scored participants:", sum(!is.na(d$diet.maha.sum)), "\n")
for (y in c("cvd_cad", "cvd_stroke_i", "cvd_hfail", "t2dm", "ckd", "death")) {
  d <- maha_outcome_followup(d, y, as.Date("2023-04-01"))
  print(attr(d, "maha_outcome_audit"))
}
cv <- c("age.landmark", "sex.f", "center", "tdi", "PC1", "PC2", "bmi.pts", "bp.pts", "nonhdl.pts", "smoke.pts", "pa.pts", "sleep.pts")
cols <- c("eid", "diet_landmark", "date_attend", "age", cv, scores, "death.t2e", "death.Yt2e")
saveRDS(d[, unique(cols)], "reviewer_revision/ukb_preflight_data.rds")
cat("Input timing preflight passed. No full analysis run.\n")
