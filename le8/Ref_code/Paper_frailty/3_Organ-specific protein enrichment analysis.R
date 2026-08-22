# ============================================================================
# Organ-Specific Protein Enrichment Analysis
# Xueqing Jia, 2025
# ============================================================================

### Load R packages
library(dplyr)
library(fdrci)
library(ggplot2)

### Load data
data<-read.csv(file = "FP_matched_organ_age.csv")

### Perform hypergeometric test
for (i in 1:nrow(data)) {
  inter <- data$N_FI[i]
  m <- data$N_all[i]
  n <- sum(data$N_all)-m
  k <- sum(data$N_FI)
  data$pvalue[i]<-phyper(inter-1,m,n,k,lower.tail = F)
  data$ratio1[i]<-paste0(inter,"/",k)
  data$ratio1_num[i]<-inter/k
  data$ratio2[i]<-paste0(m,"/",sum(data$N_all))
}

# Correct for multiple testing using FDR
data$fdr<-p.adjust(data$pvalue,method = "fdr")

write.csv(data,file = "results_enrichment_organproteins.csv")

### Visualization
data$fdr_log <- -log10(data$fdr)
data<-subset(data,data$fdr<=0.05)
data$fdr_log[which(data$fdr_log>15)]<-"15"
data$fdr_log<-as.numeric(data$fdr_log)
scale_factor <- max(data$fdr_log) / max(data$ratio1_num)

data <- data %>% arrange(fdr_log,data)

p <- ggplot(data, aes(y = reorder(Organ, fdr_log))) +
  geom_bar(aes(x = fdr_log), stat = "identity", fill = "#c3dfed") +
  geom_path(aes(x = ratio1_num * scale_factor, y = reorder(Organ, fdr_log), group = 1),
            color = "black", size = 0.75) +
  geom_point(aes(x = ratio1_num * scale_factor),
             shape = 21,
             fill = "#478ecc",
             color = "black",
             size = 3) +
  scale_x_continuous(
    name = "-log10(pvalue)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Ratio")
  ) +
  labs(y = "Organ") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12),
    axis.ticks = element_line(colour = "black", linewidth = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
p
