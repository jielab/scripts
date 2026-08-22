# ==========================================================================
# Predictive Performance for Health-related Outcomes
# Xueqing Jia and Jingyun Zhang, 2025
# ==========================================================================

### Load R packages
rm(list = ls())
library(readxl)
library(survival)
library(compareC)
library(ggplot2)

############# Predictive performance for 655 incident diseases ################################
load("./UKB_dataset_FrailtyAndProtein.Rdata")
# df_used: all data used in the analyses,including 
# Dis_name: 655 disease names
# Dis_in_female: diseases only incident in females
# ConvFactorName_new: Conventional factors names and related data type

df_used$FIsum_0_z <- scale(df_used$FIsum_0)
df_used$FI_sum_0_pre_z <- scale(df_used$FI_sum_0_pre)

### defining Cox model function ###
Dis_in_female_index<-match(Dis_in_female,Dis_name)
Dis_in_female_index<-Dis_in_female_index[is.na(Dis_in_female_index)==F]

Cox_model<-function(df_used,Dis_name,ExpName1,ExpName2,ConvFactorName,save_dir){
  cov_type<-ifelse(ConvFactorName$type=="num",ConvFactorName$name,paste("as.factor(",ConvFactorName$name,")",sep = ""))
  temp_compare<-NULL
  temp_cox<-NULL
  for (i in setdiff(c(1:length(Dis_name)),Dis_in_female_index)) {
    var_in<-c("n_eid",paste(Dis_name[i],"_ba",sep = ""), paste(Dis_name[i],"_in",sep = ""),paste("Timeto",Dis_name[i],sep = ""),"Age_XC", "sex", ExpName1,ExpName2, ConvFactorName$name)
    df_in<-df_used[,var_in]
    df_in<-na.omit(df_in)
    if(sum(df_in[,3])<10){
      next
    } # excluding disease if events<10
    
    fml1<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~ Age_XC + as.factor(sex)",sep = ""))
    fml2<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~ ",paste(cov_type,collapse = " + "),sep = ""))
    fml3<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1,sep = ""))
    fml4<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1," + Age_XC + as.factor(sex)",sep = ""))
    fml5<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1,"+",paste(cov_type,collapse = " + "),sep = ""))
    fml6<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName2,sep = ""))
    fml7<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName2," + Age_XC + as.factor(sex)",sep = ""))
    fml8<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName2,"+",paste(cov_type,collapse = " + "),sep = ""))
    fml9<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1,"+",ExpName2," + Age_XC + as.factor(sex)",sep = ""))
    fml10<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1,"+",ExpName2,"+",paste(cov_type,collapse = " + "),sep = ""))
    fml_list<-list(fml1,fml2,fml3,fml4,fml5,fml6,fml7,fml8,fml9,fml10)
    
    temp<-NULL
    fit_list<-list()
    for(j in 1:(length(fml_list))){
      fml<-fml_list[[j]]
      fit<-coxph(fml,data=df_in)
      sfit<-summary(fit)
      res<-data.frame(summary=paste(fit$nevent,"/",fit$n,sep = ""),
                      AIC=AIC(fit),
                      BIC=BIC(fit),
                      C_index=sfit$concordance["C"],
                      se_C=sfit$concordance["se(C)"])
      temp<-rbind.data.frame(temp,res)
      fit_list[[j]]<-fit
    }
    
    temp$Dis_name<-Dis_name[i]
    temp$model<-c("Age+Sex","Conv","PredFI","PredFI+Age+Sex","PredFI+Conv","FI","FI+Age+Sex","FI+Conv","PredFI+FI+Age+Sex","PredFI+FI+Conv")
    temp_cox<-rbind.data.frame(temp_cox,temp)
    
    # comparing c-index
    time<-df_in[,paste("Timeto",Dis_name[i],sep = "")]
    event<-df_in[,paste(Dis_name[i],"_in",sep = "")]
    pred2<-predict(fit_list[[2]])
    pred3<-predict(fit_list[[3]])
    pred5<-predict(fit_list[[5]])
    pred6<-predict(fit_list[[6]])
    pred8<-predict(fit_list[[8]])
    pred10<-predict(fit_list[[10]])
    
    compare_result1 <- compareC(time,event,pred3,pred2) # PredFI vs Conv
    compare_result2 <- compareC(time,event,pred5,pred2) # PredFI+Conv vs Conv
    compare_result3 <- compareC(time,event,pred3,pred6) # PredFI vs FI
    compare_result4 <- compareC(time,event,pred6,pred2) # FI vs Conv
    compare_result5 <- compareC(time,event,pred5,pred8) # PredFI+Conv vs FI+Conv
    compare_result6 <- compareC(time,event,pred10,pred8) # PredFI+FI+Conv vs FI+Conv 
    compare_result7 <- compareC(time,event,pred10,pred5) # PredFI+FI+Conv vs PredFI+Conv
    
    res_compare<-rbind.data.frame( compare_result1[c(2,8)],
                                   compare_result2[c(2,8)],
                                   compare_result3[c(2,8)],
                                   compare_result4[c(2,8)],
                                   compare_result5[c(2,8)],
                                   compare_result6[c(2,8)],
                                   compare_result7[c(2,8)])
    res_compare$Dis_name<-Dis_name[i]
    res_compare$Pairs<-c("PredFI vs Conv","PredFI+Conv vs Conv","PredFI vs FI","FI vs Conv","PredFI+Conv vs FI+Conv","PredFI+FI+Conv vs FI+Conv","PredFI+FI+Conv vs PredFI+Conv")
    temp_compare<-rbind.data.frame(temp_compare,res_compare)
    
    print(i)
    
  }
  
  for (i in Dis_in_female_index) {
    var_in<-c("n_eid",paste(Dis_name[i],"_ba",sep = ""), paste(Dis_name[i],"_in",sep = ""),paste("Timeto",Dis_name[i],sep = ""),"Age_XC", "sex", ExpName1,ExpName2, ConvFactorName$name)
    df_in<-df_used[,var_in]
    df_in<-na.omit(df_in)
    if(sum(df_in[,3])<10){
      next
    }
    cov_type<-cov_type[-2]
    fml1<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~ Age_XC",sep = ""))
    fml2<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~ ",paste(cov_type,collapse = " + "),sep = ""))
    fml3<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1,sep = ""))
    fml4<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1," + Age_XC",sep = ""))
    fml5<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1,"+",paste(cov_type,collapse = " + "),sep = ""))
    fml6<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName2,sep = ""))
    fml7<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName2," + Age_XC",sep = ""))
    fml8<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName2,"+",paste(cov_type,collapse = " + "),sep = ""))
    fml9<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1,"+",ExpName2," + Age_XC ",sep = ""))
    fml10<-formula(paste("Surv(Timeto",Dis_name[i],",",Dis_name[i],"_in) ~",ExpName1,"+",ExpName2," +", paste(cov_type,collapse = " + "),sep = ""))
    fml_list<-list(fml1,fml2,fml3,fml4,fml5,fml6,fml7,fml8,fml9,fml10)
    
    temp<-NULL
    fit_list<-list()
    for(j in 1:(length(fml_list))){
      fml<-fml_list[[j]]
      fit<-coxph(fml,data=df_in)
      sfit<-summary(fit)
      res<-data.frame(summary=paste(fit$nevent,"/",fit$n,sep = ""),
                      AIC=AIC(fit),
                      BIC=BIC(fit),
                      C_index=sfit$concordance["C"],
                      se_C=sfit$concordance["se(C)"])
      temp<-rbind.data.frame(temp,res)
      fit_list[[j]]<-fit
    }
    temp$Dis_name<-Dis_name[i]
    temp$model<-c("Age+Sex","Conv","PredFI","PredFI+Age+Sex","PredFI+Conv","FI","FI+Age+Sex","FI+Conv","PredFI+FI+Age+Sex","PredFI+FI+Conv")
    temp_cox<-rbind.data.frame(temp_cox,temp)
    
    # comparing c-index
    time<-df_in[,paste("Timeto",Dis_name[i],sep = "")]
    event<-df_in[,paste(Dis_name[i],"_in",sep = "")]
    pred2<-predict(fit_list[[2]])
    pred3<-predict(fit_list[[3]])
    pred5<-predict(fit_list[[5]])
    pred6<-predict(fit_list[[6]])
    pred8<-predict(fit_list[[8]])
    pred10<-predict(fit_list[[10]])
    
    compare_result1 <- compareC(time,event,pred3,pred2) # PredFI vs Conv
    compare_result2 <- compareC(time,event,pred5,pred2) # PredFI+Conv vs Conv
    compare_result3 <- compareC(time,event,pred3,pred6) # PredFI vs FI
    compare_result4 <- compareC(time,event,pred6,pred2) # FI vs Conv
    compare_result5 <- compareC(time,event,pred5,pred8) # PredFI+Conv vs FI+Conv
    compare_result6 <- compareC(time,event,pred10,pred8) # PredFI+FI+Conv vs FI+Conv
    compare_result7 <- compareC(time,event,pred10,pred5) # PredFI+FI+Conv vs PredFI+Conv
    
    res_compare<-rbind.data.frame( compare_result1[c(2,8)],
                                   compare_result2[c(2,8)],
                                   compare_result3[c(2,8)],
                                   compare_result4[c(2,8)],
                                   compare_result5[c(2,8)],
                                   compare_result6[c(2,8)],
                                   compare_result7[c(2,8)])
    res_compare$Dis_name<-Dis_name[i]
    res_compare$Pairs<-c("PredFI vs Conv","PredFI+Conv vs Conv","PredFI vs FI","FI vs Conv","PredFI+Conv vs FI+Conv","PredFI+FI+Conv vs FI+Conv","PredFI+FI+Conv vs PredFI+Conv")
    temp_compare<-rbind.data.frame(temp_compare,res_compare)
    
    print(i)
  }
  
  write.csv(temp_cox,file = paste0(save_dir,"Association_Disease_FI.csv"))
  write.csv(temp_compare,file = paste0(save_dir,"Association_Disease_FI_CIndex_Compare.csv"))
}

### Associations between FI/PFS and diseases ###
save_dir<-"./diseases/"
Cox_model(df_used,Dis_name,"FI_sum_0_pre_z","FIsum_0_z",ConvFactorName_new)


############# Predictive performance for mortality ##############################
death_data<-read.csv(file = "UKB_Deathcause.csv")
death_data<-death_data[,c(1,7,9,11:16)]

load("UKB_dataset_FrailtyAndProtein.Rdata")
df_used<-merge(data,df_used[c(1:5,1971:1986)], by = "n_eid")
df_used<-merge(death_data,df_used, by = "n_eid")
df_used$FIsum_0_z <- scale(df_used$FIsum_0)
df_used$FI_sum_0_pre_z <- scale(df_used$FI_sum_0_pre)

Dis_name<-colnames(death_data)[c(2,4:9)]

### defining Cox model function ###
Cox_model<-function(df_used,Dis_name,ExpName1,ExpName2,ConvFactorName,save_dir){
  cov_type<-ifelse(ConvFactorName$type=="num",ConvFactorName$name,paste("as.factor(",ConvFactorName$name,")",sep = ""))
  temp_compare<-NULL
  temp_cox<-NULL
  for (i in 1:length(Dis_name)) {
    var_in<-c("n_eid",Dis_name[i],"Timetodeath_y","Age_XC", "sex", ExpName1,ExpName2, ConvFactorName$name)
    df_in<-df_used[,var_in]
    df_in<-na.omit(df_in)
    
    fml1<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~ Age_XC + as.factor(sex)",sep = ""))
    fml2<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~ ",paste(cov_type,collapse = " + "),sep = ""))
    fml3<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~",ExpName1,sep = ""))
    fml4<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~",ExpName1," + Age_XC + as.factor(sex)",sep = ""))
    fml5<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~",ExpName1,"+",paste(cov_type,collapse = " + "),sep = ""))
    fml6<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~",ExpName2,sep = ""))
    fml7<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~",ExpName2," + Age_XC + as.factor(sex)",sep = ""))
    fml8<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~",ExpName2,"+",paste(cov_type,collapse = " + "),sep = ""))
    fml9<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~",ExpName1,"+",ExpName2," + Age_XC + as.factor(sex)",sep = ""))
    fml10<-formula(paste("Surv(Timetodeath_y,",Dis_name[i],") ~",ExpName1,"+",ExpName2,"+",paste(cov_type,collapse = " + "),sep = ""))
    fml_list<-list(fml1,fml2,fml3,fml4,fml5,fml6,fml7,fml8,fml9,fml10)
    
    temp<-NULL
    fit_list<-list()
    for(j in 1:(length(fml_list))){
      fml<-fml_list[[j]]
      fit<-coxph(fml,data=df_in)
      sfit<-summary(fit)
      res<-data.frame(summary=paste(fit$nevent,"/",fit$n,sep = ""),
                      AIC=AIC(fit),
                      BIC=BIC(fit),
                      C_index=sfit$concordance["C"],
                      se_C=sfit$concordance["se(C)"])
      temp<-rbind.data.frame(temp,res)
      fit_list[[j]]<-fit
    }
    
    temp$Dis_name<-Dis_name[i]
    temp$model<-c("Age+Sex","Conv","PredFI","PredFI+Age+Sex","PredFI+Conv","FI","FI+Age+Sex","FI+Conv","PredFI+FI+Age+Sex","PredFI+FI+Conv")
    temp_cox<-rbind.data.frame(temp_cox,temp)
    
    # comparing c-index
    time<-df_in[,"Timetodeath_y"]
    event<-df_in[,Dis_name[i]]
    pred2<-predict(fit_list[[2]])
    pred3<-predict(fit_list[[3]])
    pred5<-predict(fit_list[[5]])
    pred6<-predict(fit_list[[6]])
    pred8<-predict(fit_list[[8]])
    pred10<-predict(fit_list[[10]])
    
    compare_result1 <- compareC(time,event,pred3,pred2) # PredFI vs Conv
    compare_result2 <- compareC(time,event,pred5,pred2) # PredFI+Conv vs Conv
    compare_result3 <- compareC(time,event,pred3,pred6) # PredFI vs FI
    compare_result4 <- compareC(time,event,pred6,pred2) # FI vs Conv
    compare_result5 <- compareC(time,event,pred5,pred8) # PredFI+Conv vs FI+Conv
    compare_result6 <- compareC(time,event,pred10,pred8) # PredFI+FI+Conv vs FI+Conv 
    compare_result7 <- compareC(time,event,pred10,pred5) # PredFI+FI+Conv vs PredFI+Conv
    
    res_compare<-rbind.data.frame( compare_result1[c(2,8)],
                                   compare_result2[c(2,8)],
                                   compare_result3[c(2,8)],
                                   compare_result4[c(2,8)],
                                   compare_result5[c(2,8)],
                                   compare_result6[c(2,8)],
                                   compare_result7[c(2,8)])
    res_compare$Dis_name<-Dis_name[i]
    res_compare$Pairs<-c("PredFI vs Conv","PredFI+Conv vs Conv","PredFI vs FI","FI vs Conv","PredFI+Conv vs FI+Conv","PredFI+FI+Conv vs FI+Conv","PredFI+FI+Conv vs PredFI+Conv")
    temp_compare<-rbind.data.frame(temp_compare,res_compare)
    
    print(i)
    
  }
  write.csv(temp_cox,file = paste0(save_dir,"Association_death_FI.csv"))
  write.csv(temp_compare,file = paste0(save_dir,"Association_death_FI_CIndex_Compare.csv"))
}

### Associations between FI/PFS and mortality ###
save_dir<-"./diseases/"
Cox_model(df_used,Dis_name,"FI_sum_0_pre_z","FIsum_0_z",ConvFactorName_new,save_dir)


############# Visualization ######################
data <- read.csv(file="diseases/C-Association_Disease_FI.csv",row.names = 1)
data<-data[which(data$model=="PredFI"|data$model=="FI"|data$model=="PredFI+Conv"|data$model=="Conv"),]
data$model<-factor(data$model,levels = c("FI","PredFI","Conv","PredFI+Conv"))

### merge mortality results
death<-read.csv(file = "diseases/Association_death_FI.csv",row.names = 1)
death <- death[which(death$model %in% c("FI","PredFI","Conv","PredFI+Conv")),]
death$Disease<-"Mortality"
death$Root<-"Mortality"

data_all<-rbind(data,death)
data_all$Root<-factor(data_all$Root, 
                      levels = unique(data_all$Root)[c(1,5,7,9,2,10,8,4,3,6,12,11,13,14)])

### boxplot
data<-data_all
p <- ggplot(data, aes(x = Root, y = C_index, fill = model, color = model)) +
  stat_boxplot(aes(color = model),geom = "errorbar",width=0.5,size=0.2,alpha=1,position = position_dodge(0.8))+
  geom_boxplot(width=0.7,alpha=0.7,size=0.2,position = position_dodge(0.8),outlier.shape = NA)+
  
  geom_jitter(aes(color = model), size=0.15,position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.7) +  # 添加抖动点，并根据Category填充颜色
  
  scale_fill_manual(values = c("#1c80b9","#404080","#E69F00","#c44438")) +
  scale_color_manual(values = c("#1c80b9","#404080","#E69F00","#c44438")) +
  
  stat_summary(aes(group = model,color=model),fun = "mean", geom = "point", fill="white",shape = 21, alpha = 2.5, size = 1.1, position = position_dodge(0.8)) + # 添加空心圆表示均值
  geom_hline(yintercept=0.8,linetype=2,colour="grey",size=1) +
  
  theme_minimal() +
  labs(title = "Prediction AUC by Chapter", x = "Chapter", y = "Prediction AUC") +
  coord_cartesian(ylim = c(0.45, 1.0)) +  # 设置Y轴范围+
  scale_y_continuous(breaks = seq(0.5, 1.0, by = 0.1), 
                     labels = seq(0.5, 1.0, by = 0.1)) +
  theme(
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.line = element_line(color = "black"),
    axis.ticks.length = unit(0.05, "cm"),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    legend.position = "top"
  )

print(p)
