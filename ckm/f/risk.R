# PREVENT base 10-year total-CVD risk, with unrounded probabilities.
# Scientific source: Khan SS et al., Circulation 2024;149:430-449.
# DOI https://doi.org/10.1161/CIRCULATIONAHA.123.067626 (supplementary base model).
# Coefficients/transforms independently cross-checked against preventr 0.12.2:
# https://cran.r-project.org/src/contrib/preventr_0.12.2.tar.gz
# https://github.com/martingmayer/preventr (R/estimate_risk.R; base_10yr data).
# Reference package copyright 2024 preventr authors; MIT license below.
# This vectorized implementation computes ONLY the specified base/10yr/total-CVD
# equation. BMI has zero coefficients in this endpoint, but its validated input
# range is retained to match preventr's base-model eligibility contract.
# Total CVD here means ASCVD + heart failure, not broad CKM stage 4 (AF/PAD).
#
# MIT License (reference implementation notice)
# Copyright (c) 2024 preventr authors
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

ckm_prevent_coefficients <- function() {
    terms <- c("age", "non_hdl", "hdl", "sbp_low", "sbp_high", "diabetes", "smoking",
        "bmi_low", "bmi_high", "egfr_low", "egfr_high", "bp_tx", "statin",
        "bp_tx_sbp", "statin_non_hdl", "age_non_hdl", "age_hdl", "age_sbp",
        "age_diabetes", "age_smoking", "age_bmi", "age_egfr", "constant")
    female <- c(.7939329,.0305239,-.1606857,-.2394003,.3600781,.8667604,
        .5360739,0,0,.6045917,.0433769,.3151672,-.1477655,-.0663612,
        .1197879,-.0819715,.0306769,-.0946348,-.2705700,-.0787150,0,
        -.1637806,-3.3077280)
    male <- c(.7688528,.0736174,-.0954431,-.4347345,.3362658,.7692857,
        .4386871,0,0,.5378979,.0164827,.2888790,-.1337349,-.0475924,
        .1502730,-.0517874,.0191169,-.1049477,-.2251948,-.0895067,0,
        -.1543702,-3.0311680)
    data.frame(term=terms,female=female,male=male,stringsAsFactors=FALSE)
}

ckm_prevent_base <- function(d) {
    a <- (d$age-55)/10
    nh <- d$total_chol-d$hdl-3.5
    h <- (d$hdl-1.3)/.3
    slo <- (pmin(d$sbp,110)-110)/20
    shi <- (pmax(d$sbp,110)-130)/20
    blo <- (pmin(d$bmi,30)-25)/5
    bhi <- (pmax(d$bmi,30)-30)/5
    elo <- (pmin(d$egfr,60)-60)/(-15)
    ehi <- (pmax(d$egfr,60)-90)/(-15)
    dm <- as.numeric(d$diabetes); sm <- as.numeric(d$smoking)
    bp <- as.numeric(d$bp_tx); st <- as.numeric(d$statin)
    x <- cbind(a,nh,h,slo,shi,dm,sm,blo,bhi,elo,ehi,bp,st,bp*shi,st*nh,
        a*nh,a*h,a*shi,a*dm,a*sm,a*bhi,a*elo,1)
    cf <- ckm_prevent_coefficients()
    lp <- ifelse(as.character(d$sex)=="F",drop(x %*% cf$female),drop(x %*% cf$male))
    stats::plogis(lp)
}

ckm_prevent_risk <- function(d) {
    d <- as.data.frame(d); n <- nrow(d)
    out <- data.frame(prevent_10y=rep(NA_real_,n),
        prevent_reason=rep("missing_or_ineligible_input",n),
        prevent_source=rep(NA_character_,n),stringsAsFactors=FALSE)
    prior <- d$prevalent_cvd %in% TRUE
    out$prevent_reason[prior] <- "prevalent_cvd_not_primary_prevention"
    out$prevent_reason[is.na(d$prevalent_cvd)] <- "uncertain_baseline_cvd_timing"
    required <- c("age","sex","sbp","bp_tx","total_chol","hdl","statin","diabetes","smoking","egfr","bmi")
    if (length(setdiff(required,names(d)))) stop("PREVENT input columns are missing")
    eligible <- complete.cases(d[,required,drop=FALSE]) & d$prevalent_cvd %in% FALSE
    limits <- list(age=c(30,79),sbp=c(90,180),total_chol=c(3.36,8.28),
        hdl=c(.52,2.59),egfr=c(15,140),bmi=c(18.5,39.9))
    for (nm in names(limits)) {
        bounds <- limits[[nm]]
        invalid <- !is.na(d[[nm]]) & (!is.finite(d[[nm]]) | d[[nm]] < bounds[1L] | d[[nm]] > bounds[2L])
        out$prevent_reason[invalid & !prior] <- paste0(nm,"_outside_",bounds[1L],"_",bounds[2L])
        eligible[invalid] <- FALSE
    }
    eligible[!as.character(d$sex) %in% c("F","M")] <- FALSE
    for (nm in c("bp_tx","statin","diabetes","smoking"))
        eligible[!d[[nm]] %in% c(FALSE,TRUE,0,1)] <- FALSE
    i <- which(eligible)
    if (length(i)) {
        out$prevent_10y[i] <- ckm_prevent_base(d[i,,drop=FALSE])
        out$prevent_reason[i] <- "ok"
        out$prevent_source[i] <- "PREVENT base total-CVD 10yr; unrounded; verified against preventr 0.12.2"
    }
    # Explicit external values may be used for curated input or labelled demo.
    override <- if ("prevent_10y" %in% names(d)) is.finite(d$prevent_10y) else rep(FALSE,n)
    if (any(override)) {
        if (!all(c("risk_provenance","risk_index_date") %in% names(d)))
            stop("Precomputed prevent_10y requires risk_provenance and risk_index_date")
        valid <- d$prevent_10y>=0 & d$prevent_10y<=1 &
            !is.na(d$risk_provenance) & nzchar(trimws(d$risk_provenance)) &
            !is.na(d$risk_index_date) & d$risk_index_date==d$baseline_date
        if (any(override & !(valid %in% TRUE)))
            stop("Precomputed PREVENT probabilities must be 0..1 with a documented baseline input date")
        j <- which(override & d$prevalent_cvd %in% FALSE)
        out$prevent_10y[j] <- d$prevent_10y[j]
        out$prevent_source[j] <- d$risk_provenance[j]
        out$prevent_reason[j] <- "ok_precomputed"
    }
    out
}

ckm_score2_risk <- function(d) {
    d <- as.data.frame(d)
    # SCORE2 2021, low-risk European region (UK), original Table 1.
    # doi:10.1093/eurheartj/ehab309; coefficients also in authors' corrections.
    # SCORE2 is an explicitly adapted CKM sensitivity definition, not AHA PREVENT.
    # No clipping or age extrapolation. SCORE2-OP, SCORE2-Diabetes and SCORE2-CKD
    # are separate models and are not simulated by changing these covariates.
    n <- nrow(d)
    ans <- data.frame(score2_10y = rep(NA_real_, n),
        score2_reason = rep("not_validated_for_this_profile", n), stringsAsFactors = FALSE)
    inputs <- d[, c("age", "sex", "sbp", "total_chol", "hdl", "smoking", "diabetes", "egfr", "uacr")]
    eligible <- complete.cases(inputs) & d$age >= 40 & d$age < 70 &
        d$diabetes %in% FALSE & d$prevalent_cvd %in% FALSE &
        d$egfr >= 60 & d$uacr < 30 & d$sbp > 0 & d$total_chol > 0 & d$hdl > 0
    i <- which(eligible)
    if (length(i)) {
        a <- (d$age[i] - 60) / 5
        s <- (d$sbp[i] - 120) / 20
        tc <- d$total_chol[i] - 6
        h <- (d$hdl[i] - 1.3) / .5
        sm <- as.numeric(d$smoking[i])
        male <- as.character(d$sex[i]) == "M"
        lp_m <- .3742*a + .6012*sm + .2777*s + .1458*tc - .2698*h -
            .0755*a*sm - .0255*a*s - .0281*a*tc + .0426*a*h
        lp_f <- .4648*a + .7744*sm + .3131*s + .1002*tc - .2606*h -
            .1088*a*sm - .0277*a*s - .0226*a*tc + .0613*a*h
        lp <- ifelse(male, lp_m, lp_f)
        # Algebraically equivalent stable complementary-log-log recalibration.
        eta <- log(-log(ifelse(male, .9605, .9776))) + lp
        z <- ifelse(male, -.5699 + .7476 * eta, -.7380 + .7019 * eta)
        ans$score2_10y[i] <- -expm1(-exp(z))
        ans$score2_reason[i] <- "ok_low_risk_Europe_original_SCORE2"
    }
    override <- if ("score2_10y" %in% names(d)) is.finite(d$score2_10y) else rep(FALSE, n)
    if (any(override)) {
        if (!all(c("risk_provenance", "risk_index_date") %in% names(d)))
            stop("Precomputed score2_10y requires risk_provenance and risk_index_date")
        valid <- d$score2_10y >= 0 & d$score2_10y <= 1 &
            !is.na(d$risk_provenance) & nzchar(trimws(d$risk_provenance)) &
            !is.na(d$risk_index_date) & d$risk_index_date == d$baseline_date
        if (any(override & !(valid %in% TRUE))) stop("Invalid precomputed SCORE2 probabilities")
        # Precomputed alternative SCORE2-family models require an explicit flag
        # for their validated high-risk threshold; original SCORE2 uses age cutoffs.
        j <- which(override & d$prevalent_cvd %in% FALSE)
        outside_original <- j[!(eligible[j] %in% TRUE)]
        if (length(outside_original) && (!"score2_high_risk" %in% names(d) ||
            any(is.na(d$score2_high_risk[outside_original]))))
            stop("Precomputed SCORE2-family risk outside original eligibility requires its own score2_high_risk flag")
        ans$score2_10y[j] <- d$score2_10y[j]
        ans$score2_reason[j] <- paste0("precomputed: ", d$risk_provenance[j])
    }
    ans
}
