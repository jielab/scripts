# ============================================================================
# Proteome-wide association analysis using linear regression
# Xueqing Jia, 2025
# ============================================================================


### Load R packages
library(plyr)
library(fdrci)

### Load data
load(file = "Analysis sample_part1.Rdata")
colnames(reg_data)<-gsub("-","_",colnames(reg_data))

# Standardize Frailty Index
reg_data$FIsum_0 <- scale(reg_data$FIsum_0)

# Convert covariates to factors
reg_data$Edu_XC<-as.factor(reg_data$Edu_XC)
reg_data$smok_0<-as.factor(reg_data$smok_0)
reg_data$alco_0_XL<-as.factor(reg_data$alco_0_XL)

# first 20 genetic PCs
PCs<-paste0(colnames(reg_data)[2927:2946],collapse="+")

### linear regression
## Simple model (Model 1) controlling for age, sex, and the first 20 genetic principal components
UniLm <- function(x){
  FML <- as.formula(paste0("FIsum_0~",x,'+Age_XC+sex+',PCs))
  fit <- lm(FML,data = reg_data)
  GSum <- summary(fit)
  res_conti = as.data.frame(t(GSum$coefficients[2,]))
  res = res_conti |>
    mutate(protein = x)
  return(res)
}

VarNames<-colnames(reg_data)[16:2926]
UniVar<-lapply(VarNames, UniLm)
lm_results<-ldply(UniVar, data.frame)

# Multiple testing correction
lm_results$FDR<-p.adjust(lm_results$Pr...t.., method = 'fdr',n=length(lm_results$Pr...t..))
lm_results$P_bon<-p.adjust(lm_results$Pr...t.., method = 'bonferroni',n=length(lm_results$Pr...t..))

# Summary of significant associations
nrow(lm_results[which(lm_results$FDR<=0.05),])
nrow(lm_results[which(lm_results$P_bon<=0.05),])

write.csv(lm_results,file = "lm_proteins&FI_Model1.csv")


## Full model (Model 2) further controlling for ethnicity, educational level, Townsend deprivation index, smoking status, alcohol intake frequency, regular exercise, health diet, and body mass index based on Model 1
UniLm <- function(x){
  FML <- as.formula(paste0("FIsum_0~",x,'+Age_XC+sex+ethnic_v1+Edu_XC+TD_index+BMI_0+smok_0+alco_0_XL+hediet_0+exerc_0_XL+',PCs))
  fit <- lm(FML,data = reg_data)
  GSum <- summary(fit)
  res_conti = as.data.frame(t(GSum$coefficients[2,]))
  res = res_conti |>
    mutate(protein = x)
  return(res)
}

VarNames<-colnames(reg_data)[16:2926]
UniVar<-lapply(VarNames, UniLm)
lm_results<-ldply(UniVar, data.frame)

# Multiple testing correction
lm_results$FDR<-p.adjust(lm_results$Pr...t.., method = 'fdr',n=length(lm_results$Pr...t..))
lm_results$P_bon<-p.adjust(lm_results$Pr...t.., method = 'bonferroni',n=length(lm_results$Pr...t..))

# Summary of significant associations
nrow(lm_results[which(lm_results$FDR<=0.05),])#2051
nrow(lm_results[which(lm_results$P_bon<=0.05),])#1339

write.csv(lm_results,file = "lm_proteins&FI_Model2.csv")

