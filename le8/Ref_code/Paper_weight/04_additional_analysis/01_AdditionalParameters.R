##################################
####          N death         ####
##################################
rm(list=ls())

options(stringsAsFactors = F)
options(scipen = 1)
setwd("<path")

## packages needed
require(data.table)

## load data
ukb.phe <-  fread("<path>")

## count death
nrow(ukb.phe[death==1])
round(nrow(ukb.phe[death==1])/nrow(ukb.phe),3)