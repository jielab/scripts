####################################################################
####            Correlation matrix of the features              ####
####################################################################

# This script generates:
## SFig 6: Correlation matrix of OBSCORE features

rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 1)
setwd("<path>")

## --> packages needed <-- ##
require(data.table)
require(tidyverse)
require(dplyr)
require(magrittr)
require(corrplot)

## import data
ukb.phe <-  fread("<path>")
## import labs
lab.phe <- fread("<path>")
## import LOFO importances
LOFO_res <- fread("<path>")

## proper feature names
LOFO_res %<>% left_join(., lab.phe[, c(6, 9:10)], by = c("ExcludedFeature" =
                                                           "short_name_new"))
## get top 20 features
features <- LOFO_res %>% filter(ExcludedFeature != "None") %>%
  group_by(ExcludedFeature) %>%
  summarise(mean_Cindex_diff = mean(Cindex.diff, na.rm = TRUE)) %>%
  arrange(desc(-mean_Cindex_diff)) %>% slice(1:20) %>% pull(ExcludedFeature)

## subset dataset to top 20
cor_data <- ukb.phe %>% select(all_of(features))

## compute correlation matrix
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

## save plot
pdf("graphics/CorrMatrix.pdf", width = 7, height = 7)
corrplot(
  cor_matrix,
  method = "square",
  addCoef.col = "grey50",
  tl.col = "black",
  tl.cex = .75,
  number.cex = 0.5
)
dev.off()
