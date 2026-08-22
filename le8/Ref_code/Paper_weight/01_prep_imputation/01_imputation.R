######################################################
####           Impute one-hot encoded data        ####
######################################################
#!/usr/bin/env Rscript

# This script runs imputation for the one-hot encoded dataset
# multiple imputation to impute missing data using the miceRanger package
# variables with missingness above 25% excluded to ensure stable imputation
# generates a single dataset with five iterations
# with submission file, the dataset is split in 6 chunks

rm(list = ls())

## command-line args & options
args <- commandArgs(trailingOnly = T)
## get array task ID
array_task_id <-
  as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
options(stringsAsFactors = F, scipen = 1)

## set working directory
setwd("<path>")

## --> packages needed <-- ##
require(data.table)
require(doMC)
require(miceRanger)

## set seed
set.seed(12345)

## import data
ukb.comb <- fread("<path>")
lab.set  <- fread("<path>")

## missing percentage
lab.set[, miss.per := sapply(short_name_new, function(x)
  nrow(ukb.comb[is.na(eval(as.name(x)))])) / nrow(ukb.comb) * 100]

## convert some date variables
ukb.comb <- as.data.table(ukb.comb)

## take only hour of the day
ukb.comb[, Time_of_blow_measurement := hour(Time_of_blow_measurement)]

## get the data type for each feature
lab.set[, type := sapply(short_name_new, function(x)
  class(unlist(ukb.comb[, eval(x), with = F])))]

## convert to data frame to ease coding
ukb.comb <- as.data.frame(ukb.comb)

## convert integer to numeric
for (j in lab.set[type == "integer"]$short_name_new) {
  ukb.comb[, j] <- as.numeric(ukb.comb[, j])
}

## divide into chunks for feasible computation
## number of chunks
n_chunks <- 6

## total number of rows in the dataset
n_rows <- nrow(ukb.comb)

## calculate chunk sizes
base_chunk_size <- floor(n_rows / n_chunks)
remainder <- n_rows %% n_chunks

## create a vector of chunk sizes
chunk_sizes <- rep(base_chunk_size, n_chunks)
if (remainder > 0) {
  chunk_sizes[1:remainder] <- chunk_sizes[1:remainder] + 1
}

## determine start and end rows for each chunk
chunk_starts <- cumsum(c(1, chunk_sizes[-n_chunks]))
chunk_ends <- cumsum(chunk_sizes)

## set start and end
start_row <- chunk_starts[array_task_id]
end_row <- chunk_ends[array_task_id]

## extract the chunk for this array
ukb_chunk <- ukb.comb[start_row:end_row,]

## print number needed to impute
cat("Need to impute", nrow(lab.set[miss.per > 0]), "variables\n")

## imputation
ukb_imp <-
  miceRanger(ukb_chunk,
             m = 1,
             returnModels = F,
             verbose = T)
gc(reset = T)

## save imputed data chunk
output_file <-
  paste0("<path>",
         array_task_id,
         "<path>")
fwrite(
  completeData(ukb_imp)[[1]],
  output_file,
  sep = "\t",
  row.names = F,
  na = NA
)
