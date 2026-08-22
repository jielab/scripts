################################################
####----  Julia Carrasco-Zanini 040625  ----####
################################################

## external validation of the OBScore model

rm(list = ls())

setwd("<path>")

library(tidyr)
library(dplyr)

epic <- readstata13::read.dta13("<path>")

## filter all individuals elegible for the drug 
obs <- epic %>% 
  filter(bmi>=27&antidep_1hc!=1)

req.variables <- read.csv("<path>")

obs <- merge(obs, req.variables, by.x = "ll2id", by.y = "LL2ID")

library(lubridate)
obs <- obs %>% 
  mutate(inc_mace = ifelse(INC_MI_DERIVED==1|inc_stroke==1,1,0),
         inc_ihd_date = decimal_date(as.Date(inc_ihd_date, "%Y-%m-%d")),
         inc_stroke_date = decimal_date(as.Date(inc_stroke_date, "%Y-%m-%d")),
         inc_diabetes_date = decimal_date(as.Date(inc_diabetes_date, "%Y-%m-%d")),
         inc_mace_date = ifelse(inc_mace==1&INC_MI_DERIVED_DATE_Y10TH<inc_stroke_date,INC_MI_DERIVED_DATE_Y10TH,
                                ifelse(inc_mace==1&INC_MI_DERIVED_DATE_Y10TH>inc_stroke_date,inc_stroke_date,
                                       ifelse(inc_mace==0,NA,NA))),
         whtr = waist/height) %>% 
  rename("inc.401.1" = INC_HYPERT,
         "inc.411.4" = INC_CORONARY_ATHERO2,
         "inc.411.3" = INC_ANGINA_PECTORIS2,
         "inc.411.8" = INC_ISCH_HEART2,
         "inc.272.11" = INC_HYPERCHOL2,
         "inc.274.1" = INC_GOUT2,
         "inc.530.11" = INC_GERD2,
         "inc.550.2" = INC_DIAPH_HERNIA,
         "inc.574.1" = INC_CHOLELITHIASIS2,
         "inc.571.5" = INC_NAFLD2,
         "inc.716.9" = INC_ARTHROPATHY2,
         "inc.585.3" = INC_RENAL_DIS2,
         "inc.327.3" = INC_SLEEP_APNOEA2,
         "t2d" = inc_diabetes,
         "MI" = INC_MI_DERIVED,
         "stroke" = inc_stroke,
         "CVdeath" = CV_DEATH_AMORT,
         "MACE" = inc_mace)

## define names for outcomes 
outc_list <- c(
  "inc.401.1", "inc.411.4", "inc.411.3", "inc.411.8",
  "inc.272.11", "inc.274.1", "inc.530.11", "inc.550.2", "inc.574.1",
  "inc.571.5", "inc.716.9", "inc.585.3","inc.327.3",
  "t2d","MI","stroke","CVdeath", "MACE")

new_labels <- c(
  "inc.411.3"  =  "Angina Pectoris",
  "inc.716.9"  =  "Arthropathy",
  "inc.574.1"  =  "Cholelithiasis",
  "inc.585.3"  =  "Chronic Renal Failure",
  "inc.411.4"  =  "Coronary Athero",
  "CVdeath"    =  "CV Death",
  "inc.550.2"  =  "Diaph Hernia",
  "inc.401.1"  =  "Hypertension",
  "inc.530.11" =  "GERD",
  "inc.274.1"  =  "Gout",
  "inc.272.11" =  "Hypercholesterolemia",
  "MACE"       =  "MACE",
  "MI"         =  "MI",
  "inc.411.8"  =  "Ischemic Heart",
  "inc.571.5"  =  "NAFLD",
  "Premdeath"  =  "Premature Death",
  "inc.327.3"  =  "Sleep Apnea",
  "stroke"     =  "Stroke",
  "t2d"        =  "Type 2 Diabetes"
)

## rename predictor variable 
pred.vars <- c("age", "urate","hba1c","sexFemale",
  "chol","hdl_chol","bin_401.1","family_heart",
  "log_crea_s","log_alt","Overall_health_ratingFair",
  "smokingCurrent","bin_418","Overall_health_ratingPoor","log_ggt","whtr",
  "Long_standing_illness__disability_or_infirmityYes","bin_785")

obs <- obs %>% 
  mutate("age" = age,
         "urate" = uricacid_n,
         "hba1c" = hba1c,
         "sexFemale" = ifelse(sex==2,1,0),
         "chol" = cholestl,
         "hdl_chol" = hdl,
         "bin_401.1" = ifelse(antihyp_1hc==1|systol>=140|diastol>=90,1,0),
         "family_heart" = FH_MI,
         "log_crea_s"= creat_n,
         "log_alt" = alt_n,
         "Overall_health_ratingFair" = HEALTH,
         "smokingCurrent" = ifelse(cigstat==1,1,0),
         "bin_418" = CHEST_PAIN,
         "Overall_health_ratingPoor" = ifelse(HEALTH==4,1,0),
         "log_ggt"=ggt_n,
         "whtr"=whtr,
         "Long_standing_illness__disability_or_infirmityYes" = ifelse(PSYCHIATRIC_ILLNESS==1|ARTHRITIS==1|CA==1|ARRHYTHMIA==1,1,0),
         "bin_785" = ifelse(INC_ABDO_PAIN==1&INC_ABDO_PAIN_DATE_Y10TH<as.numeric(format(as.Date(obs$ndate, format="%Y-%m-%d"),"%Y")),1,0))

## define continuos and categorical variables 
cont.vars <- c("age","urate","hba1c","chol","hdl_chol","log_crea_s","log_alt","log_ggt","whr")
cat.vars <- pred.vars[-which(pred.vars%in%cont.vars)]

for(i in cont.vars){
  print(i)
  print(summary(obs[,i]))
}
for(i in cat.vars){
  # print(i)
  print(summary(as.factor(obs[,i])))
}

## remove participants with missing data 
obs.na <- obs %>%
  filter(if_all(pred.vars, complete.cases))

for(i in outc_list){
  print(paste0(i, " : ", sum(obs.na[,i])))
}

## defined outcomes date variable names 
outc_date_list <- c("INC_HYPERT_DATE_Y10TH","INC_CORONARY_ATHERO2_DATE_Y10TH",
                    "INC_ANGINA_PECTORIS2_DATE_Y10TH","INC_ISCH_HEART2_DATE_Y10TH",
                    "INC_HYPERCHOL2_DATE_Y10TH","INC_GOUT2_DATE_Y10TH","INC_GERD2_DATE_Y10TH",
                    "INC_DIAPH_HERNIA_DATE_Y10TH","INC_CHOLELITHIASIS2_DATE_Y10TH","INC_NAFLD2_DATE_Y10TH",
                    "INC_ARTHROPATHY2_DATE_Y10TH","INC_RENAL_DIS2_DATE_Y10TH","INC_SLEEP_APNOEA2_DATE_Y10TH",
                    "inc_diabetes_date","INC_MI_DERIVED_DATE_Y10TH","inc_stroke_date","DOD_Y10TH","inc_mace_date")

prev_exclusion <- c("INC_HYPERT_DERIVED1","INC_CORONARY_ATHERO1","INC_ANGINA_PECTORIS1","INC_IHD","INC_HYPERCHOL1",
                    "INC_GOUT1","INC_GERD1","inc.550.2","INC_CHOLELITHIASIS1","INC_NAFLD1","INC_ARTHROPATHY1",
                    "INC_RENAL_DIS1","INC_SLEEP_APNOEA1","t2d","MI","stroke","CVdeath","MACE")

prev.date <- c("INC_HYPERT_DERIVED1_DATE_Y10TH","INC_CORONARY_ATHERO1_DATE_Y10TH",
               "INC_ANGINA_PECTORIS1_DATE_Y10TH","INC_IHD_DATE_Y10TH",
               "INC_HYPERCHOL1_DATE_Y10TH","INC_GOUT1_DATE_Y10TH","INC_GERD1_DATE_Y10TH",
               "INC_DIAPH_HERNIA_DATE_Y10TH","INC_CHOLELITHIASIS1_DATE_Y10TH","INC_NAFLD1_DATE_Y10TH",
               "INC_ARTHROPATHY1_DATE_Y10TH","INC_RENAL_DIS1_DATE_Y10TH","INC_SLEEP_APNOEA1_DATE_Y10TH",
               "inc_diabetes_date","INC_MI_DERIVED_DATE_Y10TH","inc_stroke_date","DOD_Y10TH","inc_mace_date")


## convert hba1c to mmol.mol
obs.na$hba1c <- 10.929 * (obs.na$hba1c - 2.15)
summary(obs.na$hba1c)

## convert uric acid to umol/L
obs.na$urate <- obs.na$urate*(1000/160)
## log transform variables that are needed
obs.na$log_crea_s <- log10(obs.na$log_crea_s)
obs.na$log_alt <- log10(obs.na$log_alt)
obs.na$log_ggt <- log10(obs.na$log_ggt)

## force UKB scaling 
obs.na <- obs.na %>% 
  mutate(age = (age - 56.92) / 7.96,
         whtr = (whtr- 0.59)/0.06,
         urate = (urate-336.45) /78.38,
         chol = (chol-5.67)/1.18,
         hdl_chol = (hdl_chol - 1.33 )/0.33,
         hba1c = (hba1c- 37.22)/7.69,
         log_crea_s = (log_crea_s- 1.86)/0.09,
         log_alt = (log_alt- 1.38 )/0.20 ,
         log_ggt = (log_ggt- 1.53 )/0.27
  )


## scale
# library(data.table)
# obs.na[, (cont.vars) := lapply(.SD, scale), .SDcols = cont.vars]
obs.na[,cont.vars] <- sapply(obs.na[,cont.vars], scale)

## perform validation 
require(survival)

obs.na$baseline_date <-  decimal_date(as.Date(obs.na$ndate, "%Y-%m-%d"))

## define function to load and rename RData 
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

val.res <- lapply(1:length(outc_list), function(x){
  print(x)
  ## exclude prevalent cases and incident cases within the first 6 months
  dat <- obs.na
  dat$fol_prev <- dat[,prev.date[x]] - dat$baseline_date
  ## generate follow-up variable for incident status
  dat$fol <- dat[,outc_date_list[x]] - dat$baseline_date
  ## exclude prevalent and incident cases within 6 months
  if(sum((dat[,prev_exclusion[x]]==1&dat$fol_prev<0)|(dat[,outc_list[x]]==1&dat$fol<0.5), na.rm = T)>0){
    dat <- dat[-which((dat[,prev_exclusion[x]]==1&dat$fol_prev<0)|(dat[,outc_list[x]]==1&dat$fol<0.5))]
  }
  
  ## censor follow up date to 31/03/2024
  dat$fol <- ifelse(is.na(dat$fol),2024.246-dat$baseline_date,dat$fol)
  
  ## generate survival object 
  surv.dat <- Surv(dat$fol,dat[, outc_list[x]])
  
  ## boot samples
  set.seed(100)
  boot.s <- lapply(1:1000, function(x) sample(nrow(dat), round(nrow(dat)*1),replace = T))
  
  ## load trained OBScore for the specific outcome
  OBSCORE.tmp <- loadRData(paste0("<path>",outc_list[x],".RData"))
  
  ## generate c-index bootstapped estimates
  obscore.pred <- list()
  obscore.cindex <- NULL
  
  for (i in 1:1000) {
    ## clinical models
    obscore.pred[[i]] <- as.numeric(predict(OBSCORE.tmp$glmnet.opt, 
                                         type = "response",
                                         newx = as.matrix(dat[boot.s[[i]],row.names(OBSCORE.tmp$opt.coefficients)]),
                                         s = OBSCORE.tmp$lambda.opt))
    obscore.cindex <- c(obscore.cindex,glmnet::Cindex(obscore.pred[[i]], as.matrix(surv.dat[boot.s[[i]]])))
  }
  ## summarise in data frame 
  res.list <- data.frame(
    "disease.code" = outc_list[x],
    "disease.label" = new_labels[outc_list[x]],
    "cohort" = "EPIC-N",
    "obscore.cindex.mean" = mean(obscore.cindex),
    "obscore.cindex.ci.lo" = quantile(obscore.cindex, 0.025),
    "obscore.cindex.ci.upp" = quantile(obscore.cindex, 0.975),
    "n.cases" = sum(dat[,outc_list[x]]==1, na.rm = T))
})
val.res <- do.call(rbind, val.res)
# write.table(val.res, file = "<path>", sep = "\t", row.names = F)
ukb.cindex <- read.delim("<path>")

ukb.cindex <- merge(val.res, ukb.cindex, by.x = "disease.label", by.y = "Outcome")
ggplot(ukb.cindex, aes(color = disease.label))+
  geom_point(aes(x = Cindex.mean, y = obscore.cindex.mean))+
  geom_errorbar(aes(x = Cindex.mean, y = obscore.cindex.mean, xmin = ci.low, xmax = ci.upp), width=0)+
  geom_errorbar(aes(x = Cindex.mean, y = obscore.cindex.mean,ymin = obscore.cindex.ci.lo, ymax = obscore.cindex.ci.upp), width=0)+
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  theme_bw()+
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.key.size = unit(0,"cm"),)+
  xlab("C-index in UK biobank (95% CI)")+
  ylab("C-index in EPIC-N (95% CI)")+
  ylim(0.5,0.9)+
  xlim(0.5,0.9)+
  guides(color = guide_legend(ncol = 3))

cor.test(ukb.cindex$obscore.cindex.mean, ukb.cindex$Cindex.mean, method = "pearson")

val.res.10 <- lapply(1:length(outc_list), function(x){
  print(x)
  ## exclude prevalent cases and incident cases within the first 6 months
  dat <- obs.na
  dat$fol_prev <- dat[,prev.date[x]] - dat$baseline_date
  ## generate follow-up variable for incident status
  dat$fol <- dat[,outc_date_list[x]] - dat$baseline_date
  ## exclude prevalent and incident cases within 6 months
  if(sum((dat[,prev_exclusion[x]]==1&dat$fol_prev<0)|(dat[,outc_list[x]]==1&dat$fol<0.5), na.rm = T)>0){
    dat <- dat[-which((dat[,prev_exclusion[x]]==1&dat$fol_prev<0)|(dat[,outc_list[x]]==1&dat$fol<0.5))]
  }
  
  if(sum(dat[,outc_list[x]]==1, na.rm = T)>50){
  ## censor follow up date to 31/03/2024
  dat[,outc_list[x]] <- ifelse(dat[,outc_list[x]]==1&dat$fol<10,1,0)
  dat$fol <- ifelse(dat$fol<=10,dat$fol,10)
  dat$fol <- ifelse(is.na(dat$fol),10,dat$fol)
  
  ## generate survival object 
  surv.dat <- Surv(dat$fol,dat[, outc_list[x]])
  
  ## boot samples
  set.seed(100)
  boot.s <- lapply(1:1000, function(x) sample(nrow(dat), round(nrow(dat)*1),replace = T))
  
  ## load trained OBScore for the specific outcome
  OBSCORE.tmp <- loadRData(paste0("<path>",outc_list[x],".RData"))
  # load(paste0("<path>",outc_list[x],".RData"))
  
  ## generate c-index bootstapped estimates
  obscore.pred <- list()
  obscore.cindex <- NULL
  
  for (i in 1:1000) {
    ## clinical models
    obscore.pred[[i]] <- as.numeric(predict(OBSCORE.tmp$glmnet.opt, 
                                            type = "response",
                                            newx = as.matrix(dat[boot.s[[i]],row.names(OBSCORE.tmp$opt.coefficients)]),
                                            s = OBSCORE.tmp$lambda.opt))
    obscore.cindex <- c(obscore.cindex,glmnet::Cindex(obscore.pred[[i]], as.matrix(surv.dat[boot.s[[i]]])))
  }
  ## summarise in data frame 
  res.list <- data.frame(
    "disease.code" = outc_list[x],
    "disease.label" = new_labels[outc_list[x]],
    "cohort" = "EPIC-N",
    "obscore.cindex.mean" = mean(obscore.cindex),
    "obscore.cindex.ci.lo" = quantile(na.omit(obscore.cindex), 0.025),
    "obscore.cindex.ci.upp" = quantile(na.omit(obscore.cindex), 0.975),
    "n.cases" = sum(dat[,outc_list[x]]==1, na.rm = T))
  return(res.list)
}
  
  
})
val.res.10 <- do.call(rbind, val.res.10)

ukb.cindex <- read.delim("<path>")
ukb.cindex.10 <- merge(val.res.10, ukb.cindex, by.x = "disease.label", by.y = "Outcome")
library(pals)
cols25(n=16)
pdf("<path>", width = 6,height = 4)
ukb.cindex.10 %>% 
  # filter(n.cases>50) %>%
  ggplot( aes(color = disease.label))+
  geom_point(aes(x = Cindex.mean, y = obscore.cindex.mean))+
  geom_errorbar(aes(x = Cindex.mean, y = obscore.cindex.mean, xmin = ci.low, xmax = ci.upp), width=0)+
  geom_errorbar(aes(x = Cindex.mean, y = obscore.cindex.mean,ymin = obscore.cindex.ci.lo, ymax = obscore.cindex.ci.upp), width=0)+
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  theme_bw()+
  theme(legend.position = c(0.78,0.12),
        text = element_text(size = 9),
        plot.title = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.key.size = unit(0,"cm"),)+
  xlab("C-index in UK biobank (95% CI)")+
  ylab("C-index in EPIC-N (95% CI)")+
  # ylim(0.49,0.95)+
  # xlim(0.49,0.95)+
  guides(color = guide_legend(ncol = 2))+
  scale_color_manual(values = c("Angina Pectoris"="#1F78C8",
                                "Arthropathy"="#ff0000",
                                "Cholelithiasis"="#33a02c",
                                "Chronic Renal Failure"="#6A33C2",
                                "Coronary Athero"="#ff7f00",
                                "CV Death"="#565656",
                                "Diaph Hernia"="#FFD700",
                                "GERD"="#a6cee3",
                                "Gout"="#FB6496",
                                "Hypercholesterolemia"="#b2df8a",
                                "Hypertension"="#CAB2D6",
                                "Ischemic Heart"="#FDBF6F",
                                "MACE"="#999999",
                                "MI"="#EEE685",
                                "Stroke"="#C8308C",
                                "Type 2 Diabetes"="#FF83FA"))
dev.off()

cor.test(ukb.cindex.10$obscore.cindex.mean[which(ukb.cindex.10$n.cases>20)], ukb.cindex.10$Cindex.mean[which(ukb.cindex.10$n.cases>20)], method = "pearson")
cor.test(ukb.cindex.10$obscore.cindex.mean, ukb.cindex.10$Cindex.mean, method = "pearson")

# val.res.10 <- val.res.10 %>% 
#   filter(n.cases>50)
write.table(val.res.10, file = "<path>", sep = "\t", row.names = F)

sample.size <-c(5,10,20,50,70,100,200,300,390)
## downsample hypertension 
hyper.sim <- lapply(sample.size, function(x){
  print(x)
  ## exclude prevalent cases and incident cases within the first 6 months
  dat <- obs.na
  dat$fol_prev <- dat[,prev.date[1]] - dat$baseline_date
  ## generate follow-up variable for incident status
  dat$fol <- dat[,outc_date_list[1]] - dat$baseline_date
  ## exclude prevalent and incident cases within 6 months
  if(sum((dat[,prev_exclusion[1]]==1&dat$fol_prev<0)|(dat[,outc_list[1]]==1&dat$fol<0.5), na.rm = T)>0){
    dat <- dat[-which((dat[,prev_exclusion[1]]==1&dat$fol_prev<0)|(dat[,outc_list[1]]==1&dat$fol<0.5))]
  }
  
  if(sum(dat[,outc_list[1]]==1, na.rm = T)>50){
    ## censor follow up date to 31/03/2024
    dat[,outc_list[1]] <- ifelse(dat[,outc_list[1]]==1&dat$fol<10,1,0)
    dat$fol <- ifelse(dat$fol<=10,dat$fol,10)
    dat$fol <- ifelse(is.na(dat$fol),10,dat$fol)
    
    ## downsample
    ex.tmp <- which(dat[,outc_list[1]]==1)
    ex.tmp <- sample(ex.tmp, x, replace = F)
    
    dat <- dat[c(ex.tmp,which(dat[,outc_list[1]]==0)),]
    
    ## generate survival object 
    surv.dat <- Surv(dat$fol,dat[, outc_list[1]])
    
    ## boot samples
    set.seed(100)
    boot.s <- lapply(1:1000, function(x) sample(nrow(dat), round(nrow(dat)*1),replace = T))
    
    ## load trained OBScore for the specific outcome
    OBSCORE.tmp <- loadRData(paste0("<path>",outc_list[1],".RData"))
    # load(paste0("<path>".RData"))
    
    ## generate c-index bootstapped estimates
    obscore.pred <- list()
    obscore.cindex <- NULL
    
    for (i in 1:1000) {
      ## clinical models
      obscore.pred[[i]] <- as.numeric(predict(OBSCORE.tmp$glmnet.opt, 
                                              type = "response",
                                              newx = as.matrix(dat[boot.s[[i]],row.names(OBSCORE.tmp$opt.coefficients)]),
                                              s = OBSCORE.tmp$lambda.opt))
      obscore.cindex <- c(obscore.cindex,glmnet::Cindex(obscore.pred[[i]], as.matrix(surv.dat[boot.s[[i]]])))
    }
    ## summarise in data frame 
    res.list <- data.frame(
      "disease.code" = outc_list[1],
      "disease.label" = new_labels[outc_list[1]],
      "cohort" = "EPIC-N",
      "obscore.cindex.mean" = mean(obscore.cindex),
      "obscore.cindex.ci.lo" = quantile(na.omit(obscore.cindex), 0.025),
      "obscore.cindex.ci.upp" = quantile(na.omit(obscore.cindex), 0.975),
      "n.cases" = sum(dat[,outc_list[1]]==1, na.rm = T))
    return(res.list)
  }
  
  
})
hyper.sim <- do.call(rbind, hyper.sim)
ggplot(hyper.sim, aes(x = n.cases, y = obscore.cindex.mean))+
  geom_point()
