dir0 <- "D:"
gwas_file <- paste0(dir0, "/data/gwas/main/clean/x.bmi.gz")
glist <- "D:/files/glist-hg38"
ld_dir = paste0(dir0, "/data/ldref/cs/ldblk_1kg_eur")
hapref = "1kg"
grch = 38
race = "EUR"
p_min = 3; p_max = 50

source("D:/scripts/f/blockzoom.f.R")

blockzoom(gwas_file = gwas_file, block_id = 105)
