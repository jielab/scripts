#############################################################
####                    Calibration                      ####
#############################################################
# This script generates
## SFig 13: Calibration plot and metrics

rm(list=ls())
options(stringsAsFactors = F)
options(scipen = 1)
setwd("<path>")

## --> packages needed <-- ##
require(data.table)
require(survival)
require(tidyverse)
require(pmcalibration)
require(readxl)
require(dplyr)
require(magrittr)
require(ggsci)
require(patchwork)
require(ggridges)
require(survminer)
require(riskRegression)

## import data
ukb.phe <-  fread("<path>")
## import labels
lab.phe <- fread("<path>")

ukb.phe.original <- ukb.phe

###########################################
####            Calibration            ####
###########################################

EO_list_raw <- list()
EO_list_platt <- list()
calibrate <- function(outc, cdat, ten_year, nbins){
  
  # 1) remove early events
  dt <- ukb.phe.original %>%
    filter(!(((eval(as.name(cdat)) < 0.5)) & (eval(as.name(outc)) == 1))) %>% as.data.frame()
  
  # 2) splits
  set.seed(123)
  ##--  create data splits  --##
  # optimization
  o.idx <- ukb.phe$f.eid[ukb.phe$split=="o.idx"]
  # remaining samples for validation
  v.idx <- ukb.phe$f.eid[ukb.phe$split=="v.idx"]
  # convert to data frame
  ukb.phe  <- as.data.frame(ukb.phe)
  
  # 3) models
  ## model trained and tested in 25% opt set (o.idx)
  tr_name <- load(file = paste0("<path>", outc, ".RData"))
  train.model <- get(tr_name)
  ## model trained in 25% opt set and tested in validation set (v.idx)
  te_name <- load(file = paste0("<path>", outc, ".RData"))
  feature.com <- get(te_name)
  
  # 4) build train and test df
  ukb.train <- dt %>% dplyr::filter(f.eid %in% o.idx) %>% as.data.table()
  ukb.train[, lin.pred.train := train.model$linear.predictor]
  ukb.test  <- dt %>% dplyr::filter(f.eid %in% v.idx) %>% as.data.table()
  ukb.test[,  lin.pred.test  := feature.com$linear.predictor]
  
  # 5) keep relevant columns
  train_df <- ukb.train[, .(time = get(cdat), status = as.integer(get(outc)), lp = lin.pred.train)]
  test_df  <- ukb.test[,  .(time = get(cdat), status = as.integer(get(outc)), lp = lin.pred.test)]
  setDT(train_df); setDT(test_df)
  
  # 6) baseline risks
  fit_off <- coxph(Surv(time, status) ~ offset(lp), data = train_df) ## offset(lp) to get baseline hazards
  ## plain baseline hazard
  bh <- basehaz(fit_off, centered = FALSE)
  ## function approximate the baseline hazard step function and apply to LPs to get absolute risk
  predict_surv_h <- function(lp_vec, t_star) {
    H0t <- approx(bh$time, bh$hazard, xout = t_star, rule = 2)$y
    exp(-H0t * exp(lp_vec))
  }
  ## apply the function to training lp (to later train platt scaling here)
  train_df[, risk_raw := 1 - predict_surv_h(lp, ten_year)]
  ## apply the function to test lp (to plot non-platt raw risks in the test set)
  test_df[,  risk_raw := 1 - predict_surv_h(lp, ten_year)]
  ## rename events
  train_df[, event_h  := status == 1]
  test_df[,  event_h  := status == 1]
  
  # 7) Platt scaling
  ## small number to add/subtract, because logits of 0 and 1 are Inf
  clip <- 1e-6
  train_df[, p_clip := pmin(pmax(risk_raw, clip), 1 - clip)]
  
  ## fit logistic regression on logit of raw estimated risks in the training set
  ## https://pmc.ncbi.nlm.nih.gov/articles/PMC7075534/ and ChatGPT
  platt_fit <- glm(event_h ~ qlogis(p_clip), family = binomial(), data = train_df)
  
  ## function to apply the fit on the raw risks of test set
  cal_platt <- function(p, fit = platt_fit, clip = 1e-6){
    ## add clip to 0 or subtract from 1
    p <- pmin(pmax(p, clip), 1 - clip)
    ## we fit the model:
    ## coef(fit)[1] is the intercept from our training fit
    ## coef(fit)[2] is the slope from our training fit
    ## qlogis(p) is logit of out-of-the-box abs risks (in test set)
    ## plogis is inverse logit, so that the output is transformed back to absolute risk
    plogis(coef(fit)[1] + coef(fit)[2] * qlogis(p))
  }
  
  ## apply platt to raw risks
  test_df[,  risk_platt := cal_platt(risk_raw)]
  
  # 8) deciles of predictions (only on risk_platt, since ranking is the same it is the same as risk_raw for plot)
  qbreaks <- quantile(test_df$risk_platt, probs = seq(0,1,length.out = nbins+1), na.rm = TRUE)
  test_df[, bin := cut(risk_platt, breaks = qbreaks, include.lowest = TRUE, labels = FALSE)]
  
  ## we get the means of each bin (mean predicted risk)
  pred_by_bin <- test_df[, .(
    pred_raw   = mean(risk_raw,   na.rm = TRUE),
    pred_platt = mean(risk_platt, na.rm = TRUE)
  ), by = bin]
  
  # 9) We fit KM to get observed risk of each bin
  sf <- survfit(Surv(time, status) ~ bin, data = test_df, conf.type = "log-log")
  ## get KM results, the time we need is 10y, so ten_year = 10
  ## extend = T is to ensure that the latest risk at latest event is extrapolated horizontally to 10 y
  ## e.g. if there was no death on the exact 10y day, it takes the risk at last one before that
  ## otherwise NA if there was no exact event on 10y
  sm <- summary(sf, times = ten_year, extend = TRUE)
  ## write these to table
  obs_km <- data.table(
    bin      = as.integer(gsub("bin=", "", sm$strata)),
    S        = sm$surv,
    S_low    = sm$lower,
    S_high   = sm$upper
  )[, `:=`(obs      = 1 - S,
           obs_low  = 1 - S_high,
           obs_high = 1 - S_low)]
  
  calib_dt <- merge(obs_km, pred_by_bin, by = "bin")[order(bin)]
  
  # 10) metrics (EO / calibration-in-the-large)
  ## get observed risk
  sf_all        <- survfit(Surv(time, status) ~ 1, data = test_df, conf.type = "log-log")
  sm_all        <- summary(sf_all, times = ten_year, extend = TRUE)
  obs_all       <- 1 - sm_all$surv[1]
  obs_all_low   <- 1 - sm_all$upper[1]
  obs_all_high  <- 1 - sm_all$lower[1]
  
  ## estimated (mean predicted) risks
  E_raw         <- mean(test_df$risk_raw,   na.rm = TRUE)
  E_platt       <- mean(test_df$risk_platt, na.rm = TRUE)
  ## calculate estimated by observed ratio (EO)
  EO_raw        <- E_raw   / obs_all
  EO_platt      <- E_platt / obs_all
  EminusO_raw   <- E_raw   - obs_all
  EminusO_platt <- E_platt - obs_all
  
  caption_txt <- sprintf(
    "(E-O raw=%.2f) | (E-O platt=%.2f)\n(E/O raw=%.2f) | (E/O platt=%.2f)",
    EminusO_raw, EminusO_platt,EO_raw, EO_platt
  )
  
  # 11) calibration plot
  x_max <- quantile(test_df$risk_platt, probs = 0.995, na.rm = TRUE)
  
  p_cal <- ggplot() +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_errorbar(
      data = calib_dt,
      aes(x = pred_raw, ymin = obs_low, ymax = obs_high),
      width = 0,
      color = "darkgrey"
    ) +
    geom_errorbar(
      data = calib_dt,
      aes(x = pred_platt, ymin = obs_low, ymax = obs_high),
      width = 0,
      color = "skyblue"
    ) +
    geom_point(
      data = calib_dt,
      aes(x = pred_raw, y = obs),
      size = 1.5,
      color = "darkgrey"
    ) +
    geom_line(
      data = calib_dt,
      aes(x = pred_raw, y = obs),
      size = 0.75,
      color = "darkgrey"
    ) +
    geom_point(
      data = calib_dt,
      aes(x = pred_platt, y = obs),
      size = 1.5,
      color = "skyblue"
    ) +
    geom_line(
      data = calib_dt,
      aes(x = pred_platt, y = obs),
      size = 0.75,
      color = "skyblue"
    ) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, x_max)) +
    labs(x = paste("Predicted risk", outc),
         y = "Observed risk (95% CI)",
         subtitle = caption_txt) +
    theme(legend.position = "none") +
    theme_bw()
  
  p_cal
  
}

outc_list <- c("<outcomes>")

cdat_list <- c("<followup>")

## outcomes to plot
outcome_map <- setNames(cdat_list[1:2], outc_list[1:2])

## calibration plot for each outcome
plots <- lapply(names(outcome_map), function(outc) {
  cdat <- outcome_map[[outc]]
  calibrate(outc = outc, cdat = cdat, ten_year = 10, nbins = 10)
})
names(plots) <- names(outcome_map)
EO_list_platt

layout = "
ABC
DEF
GHI
JKL
MNO
PRS
"

pdf("<path>/calibrationplot.pdf", height=8.5, width=7)
wrap_plots(plots) + plot_layout(design = layout) & theme(text=element_text(size=6))
dev.off()
