########################################
####          Data with splits      ####
########################################

rm(list = ls())
options(stringsAsFactors = F)
setwd("<path>")

## --> packages needed <-- ##
require(data.table)
require(dplyr)

##--  load, bind data  --##
imp1 <- fread("<path>")
imp2 <- fread("<path>")
imp3 <- fread("<path>")
imp4 <- fread("<path>")
imp5 <- fread("<path>")
imp6 <- fread("<path>")

ukb.phe <- rbind(imp1, imp2, imp3, imp4, imp5, imp6)

##--  create data splits  --##
# seed
set.seed(123)
# training
t.idx <- sample(ukb.phe$f.eid, ceiling(nrow(ukb.phe) * .5))
# optimization
o.idx <- sample(setdiff(ukb.phe$f.eid, t.idx), ceiling(nrow(ukb.phe) * .25))
# remaining samples for validation
v.idx <- setdiff(ukb.phe$f.eid, c(t.idx, o.idx))

## create new column for splits
ukb.phe <- ukb.phe %>%
  mutate(split = case_when(f.eid %in% t.idx ~ "t.idx",
                           f.eid %in% o.idx ~ "o.idx",
                           TRUE ~ "v.idx"))

## save dataset
fwrite(ukb.phe,
       "<path>",
       sep = "\t",
       row.names = F,
       na = NA)