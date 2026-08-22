######################################################
####                   G&H ext val                ####
######################################################

# This script generates:
## external validation for incident T2D in G&H

rm(list=ls())
setwd("<path>")
options(stringsAsFactors = F)

## --> packages needed <-- ##
require(data.table)
require(arrow)
require(readxl)
require(tidyverse)
require(doMC)
require(magrittr)
require(arrow)
require(survival)
require(ggsci)
require(glmnet)

################################
####       import data      ####
################################
#------------------------------#
##--      baseline covs     --##
#------------------------------#
#-- Age, Sex, Overall Health Poor, Overall Health Average --#
## get data for IID and OrageneID
elgh.genetic       <- fread("<path>")
## import phenotypic data
elgh.id      <- fread("<path>")
## import date of birth and sex, ethnic group to get age at sample was taken
elgh.dob     <- fread("<path>")
## add oragene ID
elgh.id      <- elgh.id[ OrageneID %in% elgh.genetic$OrageneID | OrageneID %in% elgh.dob$S1QST_Oragene_ID ]
## make sure that all samples with genetics are covered
elgh.id[, gwas.id := ifelse(OrageneID %in% elgh.genetic$OrageneID, 1, 0)]
elgh.id      <- elgh.id[ order(`pseudonhs_2024-07-10`, -gwas.id, OrageneID)]
elgh.id[, ind := 1:.N, by="pseudonhs_2024-07-10"]

## merge with dataset
elgh.genetic <- merge(elgh.genetic, unique(elgh.id))
elgh.genetic %<>% left_join(.,elgh.dob[,c( "S1QST_MM-YYYY_ofBirth","S1QST_Oragene_ID", "S1QST_HealthWellbeing", "S1QST_ParentsBloodRelated", "S1QST_ParentsBloodRelated2", "S1QST_ParentsBloodRelatedOther")],
                           by=c("OrageneID"="S1QST_Oragene_ID" ))
## add baseline date
elgh.doe     <- fread("<path>")
elgh.genetic %<>% left_join(., elgh.doe[,1:2], by=c("OrageneID" = "LABCO_OrageneID"))
elgh.genetic$`LABCO_YYYY-MM_DateInLab` <- as.Date(paste0(elgh.genetic$`LABCO_YYYY-MM_DateInLab`, "-01"))

#------------------------------#
##--     continuous vars    --##
#------------------------------#
#-- BMI, ALT, Crea, Waist, Height, HDL-Chol, HbA1c, GGT, Tot-Chol --#
## run in parallel
registerDoMC(4)

elgh.phenotypes    <- c("2025_04_BMI_all_readings_at_unique_timepoints.csv",
                        "2025_04_ALT_all_readings_at_unique_timepoints.csv", 
                        "2025_04_creatinine_all_readings_at_unique_timepoints.csv", 
                        "2025_04_Waist_circumference_all_readings_at_unique_timepoints.csv", 
                        "2025_04_Height_all_readings_at_unique_timepoints.csv", 
                        "2025_04_HDL-C_all_readings_at_unique_timepoints.csv", 
                        "2025_04_HbA1c_all_readings_at_unique_timepoints.csv", 
                        "2025_04_GGT_all_readings_at_unique_timepoints.csv", 
                        "2025_04_Total_cholesterol_all_readings_at_unique_timepoints.csv")

## import the relevant data
elgh.phenotypes    <- mclapply(elgh.phenotypes, function(x){
  ## import the relevant data
  tmp       <- fread(paste0("<path>", x))
  ## filter to dataset
  tmp       <- merge(tmp, elgh.genetic[, .(`pseudonhs_2024-07-10`, `LABCO_YYYY-MM_DateInLab`)], by.x="pseudo_nhs_number", by.y = "pseudonhs_2024-07-10")
  ## keep only the most frequent unit to avoid any doubt downstream
  ii        <- table(tmp$unit)
  tmp       <- tmp[ unit == names(ii[which.max(ii)])]
  ## compute date difference
  tmp[, date.diff := as.numeric(as.IDate(date) - as.IDate(`LABCO_YYYY-MM_DateInLab`))]
  ## only keep before baseline
  tmp <- tmp[date.diff < 0]
  ## order by differences and keep only the most recent one
  tmp <- tmp[ order(pseudo_nhs_number, -date.diff)]
  tmp[, ind := 1:.N, by="pseudo_nhs_number"]
  tmp        <- tmp[ ind == 1]
  ## edit names
  tmp        <- tmp[, .(pseudo_nhs_number, value, date.diff)]
  names(tmp) <- c("pseudo_nhs_number", gsub("2025_04_|_all_readings_at_unique_timepoints.csv", "", x), paste0(gsub("2025_04_|_all_readings_at_unique_timepoints.csv", "", x), ".date.diff"))
  ## return
  return(tmp)
}, mc.cores=4)

## do full join to avoid dropouts from sex-specific outcomes
elgh.phenotypes    <- elgh.phenotypes %>% reduce(full_join, by = "pseudo_nhs_number")

## missingness
missing_perc <- colMeans(is.na(elgh.phenotypes[,-1])) * 100
nonmissing_n <- colSums(!is.na(elgh.phenotypes[,-1]))

missing_data <- data.frame(
  variable=names(missing_perc),
  missing_percentage = missing_perc,
  available_count = nonmissing_n)

## merge with covariate data
elgh.dat         <- merge(elgh.phenotypes, elgh.genetic, by.y="pseudonhs_2024-07-10", by.x = "pseudo_nhs_number", all.x=T)

#------------------------------#
##--     categorical vars   --##
#------------------------------#
##---- Abd pain ----#
abd.pain.dat <- fread("<path>")
## merge with data for date at plasma sampling
abd.pain.dat <- merge(abd.pain.dat, elgh.dat[, c('pseudo_nhs_number', 'LABCO_YYYY-MM_DateInLab')], by.x = "nhs_number", by.y = "pseudo_nhs_number")
## compute number of days since/to diagnosis
abd.pain.dat[, date.diff := as.numeric(as.IDate(date) - as.IDate(`LABCO_YYYY-MM_DateInLab`))]
## only keep those before baseline
abd.pain.dat <- abd.pain.dat[date.diff<0]
## add column for event
abd.pain.dat[, event := 1]
## edit names
abd.pain.dat <- abd.pain.dat[, .(nhs_number, event, date.diff)]
names(abd.pain.dat) <- c("pseudo_nhs_number", "abd.pain", "abd.pain.date.diff")
## combine abd.pain with the whole dataset
elgh.dat %<>% left_join(., abd.pain.dat, by="pseudo_nhs_number")
## make 0 if abd.pain is NA
elgh.dat$abd.pain[is.na(elgh.dat$abd.pain)] <- 0

##---- Smoking ----#
smoking.dat <- 
  fread("<path>")
## merge with data for date at plasma sampling
smoking.dat <- merge(smoking.dat, elgh.dat[, c('pseudo_nhs_number', 'LABCO_YYYY-MM_DateInLab')], by.x = "nhs_number", by.y = "pseudo_nhs_number")
## compute number of days since/to diagnosis
smoking.dat[, date.diff := as.numeric(as.IDate(date) - as.IDate(`LABCO_YYYY-MM_DateInLab`))]
## only keep those before baseline
smoking.dat <- smoking.dat[date.diff<0]
## add column for event
smoking.dat[, event := 1]
## edit names
smoking.dat <- smoking.dat[, .(nhs_number, event, date.diff)]
names(smoking.dat) <- c("pseudo_nhs_number", "smoker", "smoker.date.diff")
## combine smoking with the whole dataset
elgh.dat %<>% left_join(., smoking.dat, by="pseudo_nhs_number")
## make 0 if smoking is NA
elgh.dat$smoker[is.na(elgh.dat$smoker)] <- 0

##---- Hypertension ----#
ht.dat <- fread("<path>")
## merge with data for date at plasma sampling
ht.dat <- merge(ht.dat, elgh.dat[, c('pseudo_nhs_number', 'LABCO_YYYY-MM_DateInLab')], by.x = "nhs_number", by.y = "pseudo_nhs_number")
## compute number of days since/to diagnosis
ht.dat[, date.diff := as.numeric(as.IDate(date) - as.IDate(`LABCO_YYYY-MM_DateInLab`))]
## only keep those before baseline
ht.dat <- ht.dat[date.diff<0]
## add column for event
ht.dat[, event := 1]
## edit names
ht.dat <- ht.dat[, .(nhs_number, event, date.diff)]
names(ht.dat) <- c("pseudo_nhs_number", "prev_hypertension", "hypertension.date.diff")
## combine ht with the whole dataset
elgh.dat %<>% left_join(., ht.dat, by="pseudo_nhs_number")
## make 0 if smoking is NA
elgh.dat$prev_hypertension[is.na(elgh.dat$prev_hypertension)] <- 0

##---- Long-illness ----#
## 20 of most common in G&H
binary.paths    <- c("Dermatitis_atopc_contact_other_unspecified_summary_report.csv",
                     "Hypertension_summary_report.csv",
                     "Type_2_Diabetes_summary_report.csv",
                     "Depression_summary_report.csv",
                     "Enthesopathies__synovial_disorders_summary_report.csv",
                     "Gastro-oesophageal_reflux_disease_summary_report.csv",
                     "Gastritis_and_duodenitis_summary_report.csv",
                     "Anxiety_and_phobia_summary_report.csv",
                     "Chronic_sinusitis_summary_report.csv",
                     "Iron_deficiency_with_and_without_anaemia_summary_report.csv",
                     "Asthma_summary_report.csv",
                     "Osteoarthritis_excl_spine_summary_report.csv",
                     "Migraine_summary_report.csv",
                     "Menorrhagia_and_polymenorrhoea_summary_report.csv",
                     "Diabetic_eye_disease_summary_report.csv",
                     "Post-traumatic_stress_and_stress-related_disorders_summary_report.csv",
                     "Non-alcoholic_fatty_liver_disease_and_steatohepatitis_summary_report.csv",
                     "Alcohol_dependence_and_related_disease_summary_report.csv",
                     "Coronary_heart_disease_summary_report.csv",
                     "Thyroid_disease_summary_report.csv")

registerDoMC(2)

## function to process each disease folder
binary.phenotypes <- mclapply(binary.paths, function(x){

  ## import the relevant data
  tmp       <- fread(paste0("<path>", x))
  ## restrict to dataset
  tmp       <- merge(tmp, elgh.genetic[, .(`pseudonhs_2024-07-10`, `LABCO_YYYY-MM_DateInLab`)], by.x="nhs_number", by.y = "pseudonhs_2024-07-10")
  ## compute number of days since/to diagnosis
  tmp[, date.diff := as.numeric(as.IDate(date) - as.IDate(`LABCO_YYYY-MM_DateInLab`))]
  ## order by differences and keep only the most recent one (first occurence data anyways)
  tmp <- tmp[order(nhs_number, date.diff)]
  tmp[, ind := 1:.N, by = "nhs_number"]
  ## only keep prevalent
  tmp <- tmp[date.diff<0]
  ## edit names
  tmp <- tmp[, .(nhs_number, date, date.diff)]
  names(tmp) <- c("pseudo_nhs_number",paste0(x, ".date"), paste0(x, ".date.diff"))
  
  return(tmp)
}, mc.cores=2)

## merge
binary.phenotypes    <- binary.phenotypes %>% reduce(full_join, by = "pseudo_nhs_number")

## feature based on most common diseases
elgh.dat$long_illness <- ifelse(
  elgh.dat$pseudo_nhs_number %in% binary.phenotypes$pseudo_nhs_number,
  1, 0)

elgh.dat %<>% select(pseudo_nhs_number, `LABCO_YYYY-MM_DateInLab`, BMI, AgeAtRecruitment, S1QST_Gender, S1QST_HealthWellbeing, smoker, abd.pain, 
                     prev_hypertension, ALT, creatinine, Waist_circumference, Height, `HDL-C`,HbA1c,
                     GGT, Total_cholesterol,long_illness) %>% filter(BMI >= 22) %>% na.omit()

#------------------------------#
##--      incident T2D      --##
#------------------------------#
t2d.dat <- fread("<path>")
## merge with data for date at plasma sampling
t2d.dat <- merge(t2d.dat, elgh.dat[, c('pseudo_nhs_number', 'LABCO_YYYY-MM_DateInLab')],
                by.x = "nhs_number", by.y = "pseudo_nhs_number")
## compute number of days since/to diagnosis
t2d.dat[, date.diff := as.numeric(as.IDate(date) - as.IDate(`LABCO_YYYY-MM_DateInLab`))]
## add column for event
t2d.dat[, event := 1]
## edit names
t2d.dat <- t2d.dat[, .(nhs_number, event, date.diff, date)]
names(t2d.dat) <- c("pseudo_nhs_number", "t2d", "t2d.date.diff", "t2d.date")
## join to dataset
elgh.dat <- elgh.dat %>% left_join(., t2d.dat, by="pseudo_nhs_number")
## make 0 if t2d is NA
elgh.dat$t2d[is.na(elgh.dat$t2d)] <- 0
## remove prevalent T2D
elgh.dat %<>% filter(HbA1c<48) %>% filter(!(t2d == 1 & t2d.date.diff <182.5))

## get death data
deaths <- read_ipc_file("<path>")
## merge
elgh.dat %<>% left_join(deaths, by="pseudo_nhs_number")
## end of fu
elgh.dat$t2d.date[is.na(elgh.dat$t2d.date)] <- "2024-11-20"
## if there is a death, and no t2d event, replace t2d.date with date of death
elgh.dat$t2d.date[!is.na(elgh.dat$date_ons_death) & elgh.dat$t2d!=1] <- elgh.dat$date_ons_death[!is.na(elgh.dat$date_ons_death) & elgh.dat$t2d!=1]
## calculate fu time
elgh.dat[, t2d.date.diff := as.numeric(as.Date(t2d.date) - as.Date(`LABCO_YYYY-MM_DateInLab`))]

#------------------------------#
##--        validation      --##
#------------------------------#
## prep predictors
## mean(SD) from derivation cohort
elgh.dat %<>% mutate(
  age = (AgeAtRecruitment - 56.92) / 7.96,
  hba1c = (HbA1c - 37.22) / 7.69,
  sexFemale = if_else(S1QST_Gender == 2, 1L, 0L),
  chol = (Total_cholesterol-5.67) / 1.18,
  hdl_chol = (`HDL-C`-1.33) / 0.33,
  bin_401.1 = prev_hypertension,
  log_alt = (log10(ALT)-1.38) / 0.20,
  log_crea_s = (log10(creatinine)-1.86) / 0.09,
  Overall_health_ratingFair = if_else(S1QST_HealthWellbeing == 2, 1L, 0L),
  Overall_health_ratingPoor = if_else(S1QST_HealthWellbeing == 1, 1L, 0L),
  Long_standing_illness__disability_or_infirmityYes = long_illness,
  smokingCurrent = smoker,
  bin_785 = abd.pain,
  log_ggt = (log10(GGT)-1.53) / 0.27,
  whtr = ((Waist_circumference/Height)-0.59)/0.06,
  bmi = (BMI-31.16)/3.91)

## create surv object
surv.dat <- Surv(elgh.dat$t2d.date.diff, elgh.dat$t2d)
## set seed
set.seed(100)
## bootstrap samples
boot.s <- lapply(1:1000, function(x) sample(nrow(elgh.dat), round(nrow(elgh.dat)*1),replace = T))


## -----> OBSCORE <----- ##
## define function to load and rename RData 
loadRData <- function(fileName){load(fileName)get(ls()[ls() != "fileName"])}
## load model
OBSCORE.tmp <- loadRData("<path>")

## empty list to store
obscore.pred <- list()
obscore.cindex <- NULL
## coef names
feat <- rownames(OBSCORE.tmp$opt.coefficients)
## design matrix
X <- as.matrix(elgh.dat[, ..feat])
## predict
for (i in 1:1000) {
  obscore.pred[[i]] <- as.numeric(predict(OBSCORE.tmp$glmnet.opt,
                                          type = "link",
                                          newx = X[boot.s[[i]], , drop = FALSE],
                                          s = OBSCORE.tmp$lambda.opt))
  obscore.cindex <- c(obscore.cindex,glmnet::Cindex(obscore.pred[[i]], as.matrix(surv.dat[boot.s[[i]]])))
}

## -----> BMI-based <----- ##
## load model
agesexbmi.tmp <- loadRData("<path>")

## empty list to store
agesexbmi.pred <- list()
agesexbmi.cindex <- NULL
## coef names
feat_agesexbmi <- rownames(agesexbmi.tmp$opt.coefficients)
## design matrix
X <- as.matrix(elgh.dat[, ..feat_agesexbmi])
## predict
for (i in 1:1000) {
  agesexbmi.pred[[i]] <- as.numeric(predict(agesexbmi.tmp$glmnet.opt,
                                          type = "link",
                                          newx = X[boot.s[[i]], , drop = FALSE],
                                          s = agesexbmi.tmp$lambda.opt))
  agesexbmi.cindex <- c(agesexbmi.cindex,glmnet::Cindex(agesexbmi.pred[[i]], as.matrix(surv.dat[boot.s[[i]]])))
}
