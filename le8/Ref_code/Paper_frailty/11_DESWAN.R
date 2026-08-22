# ==================================================================================
# DE-SWAN Analysis
# Xueqing Jia, 2025
# ==================================================================================

### Load R packages
library(DEswan)

### Load significant proteins (Model 2)
select <- read.csv(file = "lm_proteins&FI_Model2.csv",row.names = 1)
select <- select[which(select$P_bon<=0.05),"protein"]

### Load data
load(file = "Analysis sample_part1.Rdata")
colnames(temp_data)<-gsub("-","_",colnames(temp_data))

# distribution of Age_XC
temp_data$Age_XC<-round(temp_data$Age_XC,0)
summary(temp_data$Age_XC)

hist(temp_data$Age_XC,main="qt distribution",xlab="Age (years)")

### DEswan
Age.DEswan <- DEswan(data.df = temp_data[,names(temp_data) %in% select], qt = temp_data$Age_XC,
                     window.center = seq(42, 68, 1), covariates = temp_data[,c(3:11,2927:2946)], buckets.size = 3)

Age.DEswan.wide.p <- reshape.DEswan(Age.DEswan, parameter = 1, factor = "qt")
Age.DEswan.wide.q <- q.DEswan(Age.DEswan.wide.p,method="BH")
Age.DEswan.wide.q.signif <- nsignif.DEswan(Age.DEswan.wide.q) 

# Plot
plot_data <- reshape2::melt(Age.DEswan.wide.q.signif)
plot_data$Var2 <- as.numeric(gsub("X","",plot_data$Var2))

p <- ggplot(plot_data[which(plot_data$Var1==0.05),], aes(x = Var2, y = value, group = Var1)) +
  geom_line(size=0.7) +
  geom_point(size=2) +
  theme_bw() +
  theme(
    axis.title = element_text(size=15),
    axis.text=element_text(size = 15)
  )+
  labs(x = "Age (years)", y = "No. of significant proteins (q < 0.05)")+
  scale_x_continuous(limits = c(43, 68), breaks = seq(45,65, by = 5))+
  scale_y_continuous(limits = c(200, 800), breaks = seq(200, 800, by = 200))+
  geom_vline(xintercept = c(50,63),linetype=2,colour=c("#88c4e8","#fdc58f"),size=0.8)
p