# Exploratory food-description proxy; not a validated NOVA classification.
# Word boundaries ensure "ham" does not classify "hamburger" as processed meat.
maha_processed_meat_flag <- function(description) {
  grepl("\\b(bacon|sausages?|ham|hot[ -]?dogs?|frankfurters?|pepperoni|salami|bologna|pastrami|corned beef|luncheon meat)\\b",
        tolower(description), perl = TRUE)
}

nhanes_meat_sensitivity <- function(dat, covariates) {
  score_map <- c(maha = "Original MAHA", maha_processed = "MAHA: processed-meat proxy only",
                 maha_nomeat = "MAHA: meat omitted from adverse proxy", dash = "DASH", mind = "MIND", medi = "MEDI")
  score_cols <- paste0("diet.", names(score_map), ".sum")
  # Freeze covariates and use the same participants for every compared score.
  cv <- good_covs(dat, covariates)
  need <- unique(c("death", "death_time_y", "wt", "strata", "psu", score_cols, cv))
  d <- dat[complete.cases(dat[, need]), need, drop = FALSE]
  d <- d[is.finite(d$death_time_y) & d$death_time_y > 0 & is.finite(d$wt) & d$wt > 0, ]
  if (nrow(d) < 200 || sum(d$death) < 20) stop("Insufficient common sample for meat sensitivity")
  fits <- lapply(names(score_map), function(nm) {
    dd <- d
    dd$score <- as.numeric(scale(dd[[paste0("diet.", nm, ".sum")]]))
    fm <- reformulate(c("score", cv), "survival::Surv(death_time_y, death)")
    fit <- survey::svycoxph(fm, design = make_design(dd))
    b <- unname(coef(fit)["score"]); se <- sqrt(vcov(fit)["score", "score"])
    data.frame(score = score_map[[nm]], N = nrow(dd), events = sum(dd$death),
      HR_per_SD = exp(b), lower95 = exp(b - 1.96 * se), upper95 = exp(b + 1.96 * se),
      P = 2 * pnorm(-abs(b / se)), adjustment = paste(cv, collapse = "; "))
  })
  changes <- lapply(c("maha_processed", "maha_nomeat"), function(nm) {
    a <- dat$diet.maha.sum; b <- dat[[paste0("diet.", nm, ".sum")]]
    ok <- is.finite(a) & is.finite(b)
    data.frame(score = score_map[[nm]], N_scored = sum(ok), N_changed = sum(abs(a[ok] - b[ok]) > 1e-10),
      Spearman_unweighted = cor(a[ok], b[ok], method = "spearman"))
  })
  list(mortality_common_sample = do.call(rbind, fits), score_changes = do.call(rbind, changes),
    metadata = data.frame(item = c("Purpose", "Original adverse proxy", "Processed-only proxy", "Meat-omitted proxy", "Other components", "Scaling", "Weighting", "Limitations"),
      value = c("Post-review sensitivity; original primary score is retained",
        "refined grain + sweets/pastries + SSB + fried/fast food + red/processed meat",
        "Same four non-meat terms + word-boundary processed-meat description proxy",
        "Same four non-meat terms, with no separate meat term",
        "Same eight components and minimum completeness rule; adverse component quintiles recomputed in original scoring population",
        "Per one unweighted SD within the identical complete-case mortality sample",
        "Survey-weighted Cox: original pooled dietary weights, strata and PSU",
        "Food groups may overlap; mixed dishes remain; proxy is not validated NOVA; dairy and fat-ratio limitations remain")))
}
