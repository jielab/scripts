# ==========================================================================
# Modifiable Risk Factors of PFS and FI
# Xueqing Jia, 2025
# ==========================================================================

###### Load data ######
load(file = "Analysis sample_split.Rdata")
rownames(temp_data_split)<-temp_data_split$n_eid

pred_FI<-read.csv(file = "pred_vs_true_allfolds.csv")
colnames(pred_FI)[4]<-"FI_sum_0_pre"

data<-merge(pred_FI[,c("n_eid","FI_sum_0_pre")],temp_data_split[,c(1:4,12,2927:2946)],by="n_eid")

data$FI_sum_0_pre<-scale(data$FI_sum_0_pre)
data$FIsum_0<-scale(data$FIsum_0)

###### Earlylife ######
load(file = "./Datasets/Ukb_modifiable_factors/Earlylife.Rdata")
Earlylife$n_eid<-rownames(Earlylife)

exposure_dat<-na.omit(Earlylife[,c("n_eid","n_1687_0_0")])
tem_dat<-merge(data,exposure_dat,by="n_eid")
tem_dat[,c("n_1687_0_0")]<-as.factor(tem_dat[,c("n_1687_0_0")])
FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+","n_1687_0_0"))
fit<-lm(FM,data = tem_dat)
GSum<-summary(fit)

lm_results<-matrix(ncol = 6,nrow = 2)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:2) {
  exposure<-c("n_1677_0_0","n_1787_0_0")[i]
  exposure_dat<-na.omit(Earlylife[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
  
}
results_1<-lm_results

lm_results<-matrix(ncol = 6,nrow = 4)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:2) {
  exposure<-c("n_1687_0_0","n_1697_0_0")[i]
  exposure_dat<-na.omit(Earlylife[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(2*i-1,2*i),1]<-rownames(GSum$coefficients)[5:6]
  lm_results[c(2*i-1,2*i),2]<-nrow(tem_dat)
  lm_results[c(2*i-1,2*i),3]<-round(coef(fit)[5:6],3)
  lm_results[c(2*i-1),4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[c(2*i),4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[c(2*i-1,2*i),5]<-round(GSum$coefficients[5:6,2],4)
  lm_results[c(2*i-1,2*i),6]<-GSum$coefficients[5:6,4]
}
results_2<-lm_results
lm_results<-rbind(results_1,results_2)

write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_earlylife.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_earlylife.csv")

###### Psychosocial factors ######
load(file = "./Datasets/Ukb_modifiable_factors/Psychosocial.Rdata")
Psychosocial$n_eid<-rownames(Psychosocial)

lm_results<-matrix(ncol = 6,nrow = 13)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:13) {
  exposure<-colnames(Psychosocial)[i]
  exposure_dat<-na.omit(Psychosocial[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
}
write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_Psychosocial.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_Psychosocial.csv")

###### Physical measures ######
load(file = "./Datasets/Ukb_modifiable_factors/Physical.Rdata")
Physical<-Physical[,c("n_23099_0_0","n_23105_0_0","leg_FP","arm_FP","n_23127_0_0","FVC","FEV1","PEF","lung_func","DBP","SBP","HGS")]
Physical$n_eid<-rownames(Physical)

lm_results<-matrix(ncol = 6,nrow = 12)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:12) {
  exposure<-colnames(Physical)[i]
  exposure_dat<-na.omit(Physical[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-scale(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
}
write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_Physical_measures.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_Physical_measures.csv")

###### Local environment ######
load(file = "./Datasets/Ukb_modifiable_factors/Localenvironment.Rdata")
Localenvironment<-Localenvironment[,c("mean_nitro","mean_pm10","n_24006_0_0","n_24008_0_0","n_24020_0_0","n_24021_0_0",
                                      "n_24022_0_0","n_24023_0_0","n_24024_0_0","n_24500_0_0","n_24501_0_0","n_24502_0_0",
                                      "n_24503_0_0","n_24504_0_0","n_24505_0_0","n_24506_0_0","n_24507_0_0")]
Localenvironment$n_eid<-rownames(Localenvironment)

temp_data<-merge(data,Localenvironment,by="n_eid")

lm_results<-matrix(ncol = 6,nrow = 17)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:17) {
  exposure<-colnames(Localenvironment)[c(1:17)][i]
  exposure_dat<-na.omit(Localenvironment[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-scale(log(tem_dat[,c(exposure)]+1))
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
}

write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_Local_environment.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_Local_environment.csv")

###### SES ######
XC_data<-read.csv(file = "./Datasets/Ukb_Covariates_XC.csv",header = T)
SES<-XC_data[,c("n_eid","TD_index","Edu_XC","income_0","occup_0_XC")]
SES$Edu_XC[SES$Edu_XC==1]<-0
SES$Edu_XC[SES$Edu_XC==2|SES$Edu_XC==3]<-1

lm_results<-matrix(ncol = 6,nrow = 1)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1) {
  exposure<-colnames(SES)[c(2:5)][i]
  exposure_dat<-na.omit(SES[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-scale(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
}
lm_results1<-lm_results

lm_results<-matrix(ncol = 6,nrow = 2)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:2) {
  exposure<-c("Edu_XC","occup_0_XC")[i]
  exposure_dat<-na.omit(SES[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
}
lm_results2<-rbind(lm_results,lm_results1)


lm_results<-matrix(ncol = 6,nrow = 4)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1) {
  exposure<-c("income_0")[i]
  exposure_dat<-na.omit(SES[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[1:4,1]<-rownames(GSum$coefficients)[5:8]
  lm_results[1:4,2]<-nrow(tem_dat)
  lm_results[1:4,3]<-round(coef(fit)[5:8],3)
  lm_results[1,4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[2,4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[3,4] <-paste0(round(coef(fit)[7],3)," ","(",round(coef(fit)[7]-1.96*GSum$coefficients[7,2],3),", ",round(coef(fit)[7]+1.96*GSum$coefficients[7,2],3),")")
  lm_results[4,4] <-paste0(round(coef(fit)[8],3)," ","(",round(coef(fit)[8]-1.96*GSum$coefficients[8,2],3),", ",round(coef(fit)[8]+1.96*GSum$coefficients[8,2],3),")")
  lm_results[1:4,5]<-round(GSum$coefficients[5:8,2],4)
  lm_results[1:4,6]<-GSum$coefficients[5:8,4]
}
lm_results<-rbind(lm_results,lm_results2)

write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_SES.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_SES.csv")

###### Lifestyle #####
load(file = "./Datasets/Ukb_modifiable_factors/lifestyle_data.Rdata")
lifestyle_data<-lifestyle_data[,c("n_1210_0_0","ever_smoke","alcohol_intake","n_2159_0_0","n_22036_0_0",
                                  "mod_sleep","n_1190_0_0","n_1200_0_0","n_1220_0_0","n_1269_0_0","n_1279_0_0","n_1120_0_0",
                                  "n_1130_0_0","fruit","fish","processed_meat","red_meat","whole_grain","white_grain",
                                  "vegetables","heal_diet","n_1050_0_0","n_1060_0_0","n_2139_0_0","n_2149_0_0","n_22037_0_0",
                                  "n_22038_0_0","n_22039_0_0","n_22040_0_0","n_1070_0_0","n_1080_0_0","n_1488_0_0",
                                  "n_1498_0_0","n_1528_0_0","diet_index","lei_activity","VA","VB","VC","VD","VE","VB9","V_Multi",
                                  "Fish_oil","Glucosamine","Calcium","Zinc","Iron","Selenium")]
lifestyle_data$n_eid<-rownames(lifestyle_data)

## Multiple-choice categorical variables
lm_results<-matrix(ncol = 6,nrow = 14)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:14) {
  exposure<-c("lei_activity","VA","VB","VC","VD","VE","VB9","V_Multi",
              "Fish_oil","Glucosamine","Calcium","Zinc","Iron","Selenium")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,lifestyle_data,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
}
lm_results1<-lm_results


## Single-choice categorical variables
lm_results<-matrix(ncol = 6,nrow = 15)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:15) {
  exposure<-c("n_1210_0_0","ever_smoke","n_2159_0_0","n_22036_0_0",
              "mod_sleep","n_1269_0_0","n_1279_0_0","fruit","fish","processed_meat",
              "red_meat","whole_grain","white_grain","vegetables","heal_diet")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,lifestyle_data,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
}
lm_results2<-lm_results

lm_results<-matrix(ncol = 6,nrow = 4)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:2) {
  exposure<-c("n_1190_0_0","n_1200_0_0")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,lifestyle_data,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(2*i-1,2*i),1]<-rownames(GSum$coefficients)[5:6]
  lm_results[c(2*i-1,2*i),2]<-nrow(tem_dat)
  lm_results[c(2*i-1,2*i),3]<-round(coef(fit)[5:6],3)
  lm_results[c(2*i-1),4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[c(2*i),4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[c(2*i-1,2*i),5]<-round(GSum$coefficients[5:6,2],4)
  lm_results[c(2*i-1,2*i),6]<-GSum$coefficients[5:6,4]
}
lm_results3<-lm_results


lm_results<-matrix(ncol = 6,nrow = 3)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1) {
  exposure<-c("n_1220_0_0")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,lifestyle_data,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(1:3),1]<-rownames(GSum$coefficients)[5:7]
  lm_results[c(1:3),2]<-nrow(tem_dat)
  lm_results[c(1:3),3]<-round(coef(fit)[5:7],3)
  lm_results[1,4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[2,4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[3,4] <-paste0(round(coef(fit)[7],3)," ","(",round(coef(fit)[7]-1.96*GSum$coefficients[7,2],3),", ",round(coef(fit)[7]+1.96*GSum$coefficients[7,2],3),")")
  lm_results[c(1:3),5]<-round(GSum$coefficients[5:7,2],4)
  lm_results[c(1:3),6]<-GSum$coefficients[5:7,4]
}
lm_results4<-lm_results


lm_results<-matrix(ncol = 6,nrow = 5)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1) {
  exposure<-c("n_1120_0_0")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,lifestyle_data,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(1:5),1]<-rownames(GSum$coefficients)[5:9]
  lm_results[c(1:5),2]<-nrow(tem_dat)
  lm_results[c(1:5),3]<-round(coef(fit)[5:9],3)
  lm_results[1,4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[2,4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[3,4] <-paste0(round(coef(fit)[7],3)," ","(",round(coef(fit)[7]-1.96*GSum$coefficients[7,2],3),", ",round(coef(fit)[7]+1.96*GSum$coefficients[7,2],3),")")
  lm_results[4,4] <-paste0(round(coef(fit)[8],3)," ","(",round(coef(fit)[8]-1.96*GSum$coefficients[8,2],3),", ",round(coef(fit)[8]+1.96*GSum$coefficients[8,2],3),")")
  lm_results[5,4] <-paste0(round(coef(fit)[9],3)," ","(",round(coef(fit)[9]-1.96*GSum$coefficients[9,2],3),", ",round(coef(fit)[9]+1.96*GSum$coefficients[9,2],3),")")
  lm_results[c(1:5),5]<-round(GSum$coefficients[5:9,2],4)
  lm_results[c(1:5),6]<-GSum$coefficients[5:9,4]
}
lm_results5<-lm_results


lm_results<-matrix(ncol = 6,nrow = 4)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1) {
  exposure<-c("n_1130_0_0")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,lifestyle_data,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(1:4),1]<-rownames(GSum$coefficients)[5:8]
  lm_results[c(1:4),2]<-nrow(tem_dat)
  lm_results[c(1:4),3]<-round(coef(fit)[5:8],3)
  lm_results[1,4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[2,4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[3,4] <-paste0(round(coef(fit)[7],3)," ","(",round(coef(fit)[7]-1.96*GSum$coefficients[7,2],3),", ",round(coef(fit)[7]+1.96*GSum$coefficients[7,2],3),")")
  lm_results[4,4] <-paste0(round(coef(fit)[8],3)," ","(",round(coef(fit)[8]-1.96*GSum$coefficients[8,2],3),", ",round(coef(fit)[8]+1.96*GSum$coefficients[8,2],3),")")
  lm_results[c(1:4),5]<-round(GSum$coefficients[5:8,2],4)
  lm_results[c(1:4),6]<-GSum$coefficients[5:8,4]
}
lm_results6<-lm_results

lm_results<-matrix(ncol = 6,nrow = 15)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:15) {
  exposure<-c("alcohol_intake","n_1050_0_0","n_1060_0_0","n_2139_0_0","n_2149_0_0","n_22037_0_0",
              "n_22038_0_0","n_22039_0_0","n_22040_0_0","n_1070_0_0","n_1080_0_0","n_1488_0_0",
              "n_1498_0_0","n_1528_0_0","diet_index")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-scale(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
}
lm_results7<-lm_results

lm_results<-rbind(lm_results1,lm_results2,lm_results3,lm_results4,lm_results5,lm_results6,lm_results7)

write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_lifestyle.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_lifestyle.csv")

###### Physical measures-tertiles######
load(file = "./Datasets/Ukb_modifiable_factors/Physical.Rdata")
Physical<-Physical[,c("n_23099_0_0","n_23105_0_0","leg_FP","arm_FP","n_23127_0_0","FVC","FEV1","PEF","lung_func","DBP","SBP","HGS")]
Physical$n_eid<-rownames(Physical)

lm_results<-matrix(ncol = 6,nrow = 36)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:12) {
  exposure<-colnames(Physical)[i]
  exposure_dat<-na.omit(Physical[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  per<-quantile(tem_dat[,c(exposure)],c(1/3,2/3))
  tem_dat[,c(exposure)]<-as.factor(ifelse(tem_dat[,c(exposure)]<=per[1],1,ifelse(tem_dat[,c(exposure)]<=per[2],2,3)))
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(3*i-2),1]<-exposure
  lm_results[c(3*i-1,3*i),1]<-rownames(GSum$coefficients)[5:6]
  lm_results[c(3*i-2,3*i-1,3*i),2]<-c(paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==1)),"/",nrow(tem_dat)),paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==2)),"/",nrow(tem_dat)),paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==3)),"/",nrow(tem_dat)))
  lm_results[c(3*i-1,3*i),3]<-round(coef(fit)[5:6],3)
  lm_results[c(3*i-1),4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[c(3*i),4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[c(3*i-1,3*i),5]<-round(GSum$coefficients[5:6,2],4)
  lm_results[c(3*i-1,3*i),6]<-GSum$coefficients[5:6,4]
}
write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_Physical_measures_tertiles.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_Physical_measures_tertiles.csv")

###### Local environment-tertiles ######
load(file = "./Datasets/Ukb_modifiable_factors/Localenvironment.Rdata")
Localenvironment<-Localenvironment[,c("mean_nitro","mean_pm10","n_24006_0_0","n_24008_0_0","n_24020_0_0","n_24021_0_0",
                                      "n_24022_0_0","n_24023_0_0","n_24024_0_0","n_24500_0_0","n_24501_0_0","n_24502_0_0",
                                      "n_24503_0_0","n_24504_0_0","n_24505_0_0","n_24506_0_0","n_24507_0_0")]
Localenvironment$n_eid<-rownames(Localenvironment)

temp_data<-merge(data,Localenvironment,by="n_eid")

lm_results<-matrix(ncol = 6,nrow = 51)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:17) {
  exposure<-colnames(Localenvironment)[c(1:17)][i]
  exposure_dat<-na.omit(Localenvironment[,c("n_eid",exposure)])  
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  per<-quantile(tem_dat[,c(exposure)],c(1/3,2/3))
  tem_dat[,c(exposure)]<-as.factor(ifelse(tem_dat[,c(exposure)]<=per[1],1,ifelse(tem_dat[,c(exposure)]<=per[2],2,3)))
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(3*i-2),1]<-exposure
  lm_results[c(3*i-1,3*i),1]<-rownames(GSum$coefficients)[5:6]
  lm_results[c(3*i-2,3*i-1,3*i),2]<-c(paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==1)),"/",nrow(tem_dat)),paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==2)),"/",nrow(tem_dat)),paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==3)),"/",nrow(tem_dat)))
  lm_results[c(3*i-1,3*i),3]<-round(coef(fit)[5:6],3)
  lm_results[c(3*i-1),4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[c(3*i),4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[c(3*i-1,3*i),5]<-round(GSum$coefficients[5:6,2],4)
  lm_results[c(3*i-1,3*i),6]<-GSum$coefficients[5:6,4]
}
write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_Local_environment_tertiles.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_Local_environment_tertiles.csv")

###### SES-tertiles ######
XC_data<-read.csv(file = "./Datasets/Ukb_Covariates_XC.csv",header = T)
SES<-XC_data[,c("n_eid","TD_index","Edu_XC","income_0","occup_0_XC")]
SES$Edu_XC[SES$Edu_XC==1]<-0
SES$Edu_XC[SES$Edu_XC==2|SES$Edu_XC==3]<-1

lm_results<-matrix(ncol = 6,nrow = 3)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1) {
  exposure<-colnames(SES)[c(2:5)][i]
  exposure_dat<-na.omit(SES[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  per<-quantile(tem_dat[,c(exposure)],c(1/3,2/3))
  tem_dat[,c(exposure)]<-as.factor(ifelse(tem_dat[,c(exposure)]<=per[1],1,ifelse(tem_dat[,c(exposure)]<=per[2],2,3)))
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(3*i-2),1]<-exposure
  lm_results[c(3*i-1,3*i),1]<-rownames(GSum$coefficients)[5:6]
  lm_results[c(3*i-2,3*i-1,3*i),2]<-c(paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==1)),"/",nrow(tem_dat)),paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==2)),"/",nrow(tem_dat)),paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==3)),"/",nrow(tem_dat)))
  lm_results[c(3*i-1,3*i),3]<-round(coef(fit)[5:6],3)
  lm_results[c(3*i-1),4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[c(3*i),4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[c(3*i-1,3*i),5]<-round(GSum$coefficients[5:6,2],4)
  lm_results[c(3*i-1,3*i),6]<-GSum$coefficients[5:6,4]
}
lm_results1<-lm_results

lm_results<-matrix(ncol = 6,nrow = 2)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:2) {
  exposure<-c("Edu_XC","occup_0_XC")[i]
  exposure_dat<-na.omit(SES[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[i,1]<-rownames(GSum$coefficients)[5]
  lm_results[i,2]<-nrow(tem_dat)
  lm_results[i,3]<-round(coef(fit)[5],3)
  lm_results[i,4]<-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[i,5]<-round(GSum$coefficients[5,2],4)
  lm_results[i,6]<-GSum$coefficients[5,4]
}
lm_results2<-rbind(lm_results,lm_results1)


lm_results<-matrix(ncol = 6,nrow = 4)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1) {
  exposure<-c("income_0")[i]
  exposure_dat<-na.omit(SES[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[1:4,1]<-rownames(GSum$coefficients)[5:8]
  lm_results[1:4,2]<-nrow(tem_dat)
  lm_results[1:4,3]<-round(coef(fit)[5:8],3)
  lm_results[1,4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[2,4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[3,4] <-paste0(round(coef(fit)[7],3)," ","(",round(coef(fit)[7]-1.96*GSum$coefficients[7,2],3),", ",round(coef(fit)[7]+1.96*GSum$coefficients[7,2],3),")")
  lm_results[4,4] <-paste0(round(coef(fit)[8],3)," ","(",round(coef(fit)[8]-1.96*GSum$coefficients[8,2],3),", ",round(coef(fit)[8]+1.96*GSum$coefficients[8,2],3),")")
  lm_results[1:4,5]<-round(GSum$coefficients[5:8,2],4)
  lm_results[1:4,6]<-GSum$coefficients[5:8,4]
}
lm_results<-rbind(lm_results,lm_results2)

write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_SES_tertiles.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_SES_tertiles.csv")

###### Lifestyle-tertiles ######
load(file = "./Datasets/Ukb_modifiable_factors/lifestyle_data.Rdata")
lifestyle_data<-lifestyle_data[,c("n_1210_0_0","ever_smoke","alcohol_intake","n_2159_0_0","n_22036_0_0",
                                  "mod_sleep","n_1190_0_0","n_1200_0_0","n_1220_0_0","n_1269_0_0","n_1279_0_0","n_1120_0_0",
                                  "n_1130_0_0","fruit","fish","processed_meat","red_meat","whole_grain","white_grain",
                                  "vegetables","heal_diet","n_1050_0_0","n_1060_0_0","n_2139_0_0","n_2149_0_0","n_22037_0_0",
                                  "n_22038_0_0","n_22039_0_0","n_22040_0_0","n_1070_0_0","n_1080_0_0","n_1488_0_0",
                                  "n_1498_0_0","n_1528_0_0","diet_index","lei_activity","VA","VB","VC","VD","VE","VB9","V_Multi",
                                  "Fish_oil","Glucosamine","Calcium","Zinc","Iron","Selenium")]
lifestyle_data$n_eid<-rownames(lifestyle_data)

lm_results<-matrix(ncol = 6,nrow = 45)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1:15) {
  exposure<-c("alcohol_intake","n_1050_0_0","n_1060_0_0","n_2139_0_0","n_2149_0_0","n_22037_0_0",
              "n_22038_0_0","n_22039_0_0","n_22040_0_0","n_1070_0_0","n_1080_0_0","n_1488_0_0",
              "n_1498_0_0","n_1528_0_0","diet_index")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,exposure_dat,by="n_eid")
  per<-quantile(tem_dat[,c(exposure)],c(1/3,2/3))
  tem_dat[,c(exposure)]<-as.factor(ifelse(tem_dat[,c(exposure)]<=per[1],1,ifelse(tem_dat[,c(exposure)]<=per[2],2,3)))
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(3*i-2),1]<-exposure
  lm_results[c(3*i-1,3*i),1]<-rownames(GSum$coefficients)[5:6]
  lm_results[c(3*i-2,3*i-1,3*i),2]<-c(paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==1))," / ",nrow(tem_dat)),paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==2)),"/",nrow(tem_dat)),paste0(nrow(subset(tem_dat,tem_dat[,c(exposure)]==3)),"/",nrow(tem_dat)))
  lm_results[c(3*i-1,3*i),3]<-round(coef(fit)[5:6],3)
  lm_results[c(3*i-1),4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[c(3*i),4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[c(3*i-1,3*i),5]<-round(GSum$coefficients[5:6,2],4)
  lm_results[c(3*i-1,3*i),6]<-GSum$coefficients[5:6,4]
}

write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_lifestyle_tertiles.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_lifestyle_tertiles.csv")

###### Lifestyle-supplementary ######
### smoke + alcohol comsumption
XC_data<-read.csv(file = "./Datasets/Ukb_Covariates_XC.csv",header = T)
lifestyle_data<-XC_data[,c("n_eid","smok_0","alco_0_XL")]

lifestyle_data$smok_0<-as.factor(lifestyle_data$smok_0)
lifestyle_data$alco_0_XL<-as.factor(lifestyle_data$alco_0_XL)

lm_results<-matrix(ncol = 6,nrow = 2)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1) {
  exposure<-c("smok_0")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,lifestyle_data,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(1:2),1]<-rownames(GSum$coefficients)[5:6]
  lm_results[c(1:2),2]<-nrow(tem_dat)
  lm_results[c(1:2),3]<-round(coef(fit)[5:6],3)
  lm_results[1,4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[2,4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[c(1:2),5]<-round(GSum$coefficients[5:6,2],4)
  lm_results[c(1:2),6]<-GSum$coefficients[5:6,4]
}
lm_results1<-lm_results


lm_results<-matrix(ncol = 6,nrow = 3)
colnames(lm_results)<-c("exposure","n_sample","Bvalue","Bvalue(95%CI)","SE","Pvalue")
for (i in 1) {
  exposure<-c("alco_0_XL")[i]
  exposure_dat<-na.omit(lifestyle_data[,c("n_eid",exposure)])
  tem_dat<-merge(data,lifestyle_data,by="n_eid")
  tem_dat[,c(exposure)]<-as.factor(tem_dat[,c(exposure)])
  FM<-as.formula(paste0("FIsum_0~Age_XC+sex+ethnic_v1+",exposure))
  fit<-lm(FM,data = tem_dat)
  GSum<-summary(fit)
  lm_results[c(1:3),1]<-rownames(GSum$coefficients)[5:7]
  lm_results[c(1:3),2]<-nrow(tem_dat)
  lm_results[c(1:3),3]<-round(coef(fit)[5:7],3)
  lm_results[1,4] <-paste0(round(coef(fit)[5],3)," ","(",round(coef(fit)[5]-1.96*GSum$coefficients[5,2],3),", ",round(coef(fit)[5]+1.96*GSum$coefficients[5,2],3),")")
  lm_results[2,4] <-paste0(round(coef(fit)[6],3)," ","(",round(coef(fit)[6]-1.96*GSum$coefficients[6,2],3),", ",round(coef(fit)[6]+1.96*GSum$coefficients[6,2],3),")")
  lm_results[3,4] <-paste0(round(coef(fit)[7],3)," ","(",round(coef(fit)[7]-1.96*GSum$coefficients[7,2],3),", ",round(coef(fit)[7]+1.96*GSum$coefficients[7,2],3),")")
  lm_results[c(1:3),5]<-round(GSum$coefficients[5:7,2],4)
  lm_results[c(1:3),6]<-GSum$coefficients[5:7,4]
}
lm_results2<-lm_results

lm_results<-rbind(lm_results1,lm_results2)

write.csv(lm_results,file = "./Results/Modi_factors/ProteoFI_lifestyle_smok&alco.csv")
write.csv(lm_results,file = "./Results/Modi_factors/FI_lifestyle_smok&alco.csv")
