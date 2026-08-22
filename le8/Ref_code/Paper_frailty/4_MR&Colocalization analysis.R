# ============================================================================
# Two-sample Mendelian randomization (MR) and Colocalization analysis
# Weijing Gao, 2025
# ============================================================================

# Load R packages
rm(list = ls());gc()
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(tools)
library(VariantAnnotation)
library(gwasglue)
library(MRPRESSO)
library(coloc)
library(tidyverse)
library(locuscomparer)

#### Clumping ####
# function for picking the first column name that exists in the data from a set of candidates
pick_col <- function(nms, cand) { 
  hit <- intersect(cand, nms)
  if (length(hit) == 0) stop("None of these column names were found: ", paste(cand, collapse=", "))
  hit[1]
}

# function for pQTL data cleaning and clumping
process_one <- function(f, pop = pop) {
  message("Processing: ", f)
  
  # GWAS summary statistics of cis-pQTL data (UKB-ppp)
  df <- fread(f)
  nms <- names(df)
  if("LOG10P" %in% nms){
    df$P <- 10^(-df$LOG10P)
    nms <- names(df)
  }
  # Select the required columns
  col_snp  <- pick_col(nms, c("SNP", "rsids"))
  col_beta <- pick_col(nms, c("BETA","beta", "Beta"))
  col_se   <- pick_col(nms, c("SE","se"))
  col_eaf  <- pick_col(nms, c("A1FREQ","EAF","eaf"))
  col_p    <- pick_col(nms, c("P","pval", "Pval"))
  col_ea   <- pick_col(nms, c("ALLELE1","A1","effect_allele", "effectAllele", "EA"))
  col_oa   <- pick_col(nms, c("ALLELE0","A2","other_allele",  "otherAllele",  "OA"))
  col_chr  <- pick_col(nms, c("CHR", "Chrom", "chr", "chrom"))
  col_pos  <- pick_col(nms, c("BP", "Pos", "position"))
  col_n    <- pick_col(nms, c("N"))
  
  # Standardize to a unified set of column names
  std <- data.frame(
    SNP            = df[[col_snp]],
    beta           = as.numeric(df[[col_beta]]),
    se             = as.numeric(df[[col_se]]),
    eaf            = as.numeric(df[[col_eaf]]),
    pval           = as.numeric(df[[col_p]]),
    effect_allele  = df[[col_ea]],
    other_allele   = df[[col_oa]],
    CHR            = df[[col_chr]],
    BP             = as.numeric(df[[col_pos]]),
    N              = as.numeric(df[[col_n]]),
    stringsAsFactors = FALSE
  )
  
  # Handle chromosome labels that may look like "chr6"
  std$CHR <- gsub("^chr", "", as.character(std$CHR), ignore.case = TRUE)
  suppressWarnings(std$CHR <- as.integer(std$CHR))
  
  # Compute MAF and filter by thresholds
  std$MAF <- pmin(std$eaf, 1 - std$eaf)
  std <- std[is.finite(std$pval) & std$pval < 5e-08 &
               is.finite(std$MAF) & std$MAF > 0.01, ]
  
  if (nrow(std) == 0) {
    message("No variants passed the thresholds: ", basename(f))
    return(invisible(NULL))
  }
  
  # For the same SNP, keep only the row with the smallest p value
  std <- std[order(std$SNP, std$pval), ]
  std <- std[!duplicated(std$SNP), ]
  
  # Format as TwoSampleMR exposure data
  fd_args <- list(
    std,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col   = "se",
    eaf_col  = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    pval_col = "pval",
    samplesize_col = "N",
    chr_col = "CHR",
    pos_col = "BP"
  )
  # If the TwoSampleMR version supports phenotype_col, attach the file name
  if ("phenotype_col" %in% names(formals(TwoSampleMR::format_data))) {
    std$file_name <- basename(f)
    fd_args$phenotype_col <- basename(f)
  } 
  
  exp_dat <- do.call(TwoSampleMR::format_data, fd_args)
  
  # Fallback: ensure exposure/id.exposure and file_name are present
  fname <- basename(f)
  if ("exposure"     %in% names(exp_dat)) exp_dat$exposure     <- fname
  if ("id.exposure"  %in% names(exp_dat)) exp_dat$id.exposure  <- fname
  if (!"file_name"    %in% names(exp_dat)) exp_dat$file_name    <- fname
  
  # Clumping
  clumped <- suppressWarnings(
    clump_data(exp_dat, clump_kb = 10000, clump_r2 = 0.1, pop = "EUR")
  )
  
  if (!is.null(clumped) && nrow(clumped) > 0) {
    # Attach the file_name column (exposure can map directly)
    if (!"file_name" %in% names(clumped)) clumped$file_name <- fname
    out <- paste0(file_path_sans_ext(fname), ".clumped_r2_0.1.csv")
    write.csv(clumped, out, row.names = FALSE)
    message("Completed: ", fname, " -> ", out, "; SNPs retained = ", nrow(clumped))
  } else {
    message("No SNPs remain after clumping: ", fname)
  }
  
  invisible(TRUE)
}

# Batch-process all files in the specified directory
files <- list.files("./pqtl/ukb/", full.names = TRUE)

# path for saving the clumped results
setwd("./IVclumped")

if (length(files) == 0) {
  message("No files found in the current directory.")
} else {
  # You can change pop to "EAS", "AFR", etc.
  for (f in files) {
    try(process_one(f, pop = "EUR"), silent = TRUE)
  }
}

#### Two-sample MR analysis ####

# removed instrumental variants associated with more than five proteins
DATA <- data.frame()
file <- list.files("./IVclumped/")
for(i in file){
  data <- fread(sprintf("./IVclumped/%s",i))
  DATA <- rbind(DATA,data)
}
snpcount <- table(DATA$SNP) %>% as.data.frame()
DATA1 <- DATA[!DATA$SNP %in% snpcount[snpcount$Freq>5,1],]
write.csv(DATA1,"./exposure.csv",row.names = F)

#Exposure
exposureFile = "./exposure.csv" 
expo_rt=read_exposure_data(filename=exposureFile,
                           sep = ",",
                           snp_col = "SNP",
                           beta_col = "beta.exposure",
                           se_col = "se.exposure",
                           pval_col = "pval.exposure",
                           effect_allele_col="effect_allele.exposure",
                           other_allele_col = "other_allele.exposure",
                           eaf_col = "eaf.exposure",
                           phenotype_col = "id.exposure",
                           id_col = "id.exposure",
                           samplesize_col = "samplesize.exposure",
                           chr_col="chr.exposure", pos_col = "pos.exposure",
                           clump=FALSE)
#Outcome
outc_rt <- read_outcome_data(
  snps = expo_rt$SNP,
  filename = "./FI.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA_TG",
  se_col = "SE_TG",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  eaf_col = "EAF_TG",
  pval_col = "P_TG",
  samplesize_col = "N")

#harmonise data
harm_rt <- harmonise_data(
  exposure_dat =  expo_rt, 
  outcome_dat = outc_rt, action=2)

#calculating instrumental variables strength (R2, F-stat), filter F > 10
harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
  (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
     2 * harm_rt$samplesize.exposure * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
harm_rt$F <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
harm_rt <- harm_rt[harm_rt$F > 10, ]
#Steiger filtering
harm_rt <- steiger_filtering(harm_rt)
harm_rt <- harm_rt %>% dplyr::filter(steiger_dir=="TRUE" | steiger_pval>0.05)

# 1) Basic checks
stopifnot(exists("harm_rt"))
setDT(harm_rt)
harm_rt$id.exposure <- harm_rt$exposure
need_cols <- c("id.exposure","id.outcome","SNP")
miss <- setdiff(need_cols, names(harm_rt))
if (length(miss) > 0) stop("harm_rt missing columns: ", paste(miss, collapse = ", "))

# 2) Keep rows marked usable for MR by TwoSampleMR
if ("mr_keep" %in% names(harm_rt)) {
  harm_rt <- harm_rt[is.na(harm_rt$mr_keep) | harm_rt$mr_keep == TRUE,]
}

# 3) Create list of exposure–outcome pairs
pairs_dt <- unique(harm_rt[, c("id.exposure", "id.outcome")])

# 4)function for create safe file name (strip invalid characters)
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

#perform MR
for (k in seq_len(nrow(pairs_dt))) {
  
  ex_id <- pairs_dt$id.exposure[k]
  out_id <- pairs_dt$id.outcome[k]
  
  dat <- copy(harm_rt[id.exposure == ex_id & id.outcome == out_id])
  if (nrow(dat) == 0) next
  
  folder <- paste0(safe_name(ex_id), "_", safe_name(out_id))
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  
  # --- Main MR analysis ---
  mr_result <- NULL
  result_or <- NULL
  try({
    mr_result <- mr(dat)
    result_or <- generate_odds_ratios(mr_result)
  }, silent = TRUE)
  
  # Write harmonised subset
  fwrite(dat, file = file.path(folder, "harmonise.txt"), sep = "\t")
  
  # Write MR results (raw + OR)
  if (!is.null(mr_result)) {
    write.table(mr_result, file = file.path(folder, "mr_result.txt"),
                row.names = FALSE, sep = "\t", quote = FALSE)
  }
  if (!is.null(result_or)) {
    # Your old habit: write columns 5 to end; if fewer columns, write all
    if (ncol(result_or) >= 5) {
      write.table(result_or[, 5:ncol(result_or), drop = FALSE],
                  file = file.path(folder, "OR.txt"),
                  row.names = FALSE, sep = "\t", quote = FALSE)
    } else {
      write.table(result_or, file = file.path(folder, "OR.txt"),
                  row.names = FALSE, sep = "\t", quote = FALSE)
    }
  }
  
  # --- Additional analyses & plots (best-effort; do not interrupt on errors) ---
  # 1) Egger intercept (horizontal pleiotropy)
  try({
    pleiotropy <- mr_pleiotropy_test(dat)
    write.table(pleiotropy, file = file.path(folder, "pleiotropy.txt"),
                row.names = FALSE, sep = "\t", quote = FALSE)
  })
  
  # 2) MR_PRESSO
    set.seed(123)
    try({
    presso <- MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome   = "se.outcome",
      SdExposure  = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = as.data.frame(dat),
      NbDistribution = 1000,
      SignifThreshold = 0.05
    )
    capture.output(presso, file = file.path(folder, "presso.txt"))
    }, silent = TRUE)
    
  # 3) Heterogeneity tests
  try({
    heterogeneity <- mr_heterogeneity(dat)
    write.table(heterogeneity, file = file.path(folder, "heterogeneity.txt"),
                row.names = FALSE, sep = "\t", quote = FALSE)
  }, silent = TRUE)
  
  # 4) Scatter plot
  try({
    if (!is.null(mr_result)) {
      p1 <- mr_scatter_plot(mr_result, dat)
      if (length(p1) > 0) {
        ggsave(filename = file.path(folder, "scatter.pdf"),
               plot = p1[[1]], width = 8, height = 8)
      }
    }
  }, silent = TRUE)
  
  # 5) Single-SNP analysis + forest plot
  singlesnp_res <- NULL
  try({
    singlesnp_res <- mr_singlesnp(dat)
    singlesnpOR <- generate_odds_ratios(singlesnp_res)
    write.table(singlesnpOR, file = file.path(folder, "singlesnpOR.txt"),
                row.names = FALSE, sep = "\t", quote = FALSE)
    p2 <- mr_forest_plot(singlesnp_res)
    if (length(p2) > 0) {
      ggsave(filename = file.path(folder, "forest.pdf"),
             plot = p2[[1]], width = 8, height = 8)
    }
  }, silent = TRUE)
  
  # 6) Leave-one-out + sensitivity plot (>= 2 instruments is more stable)
  try({
    if (uniqueN(dat$SNP) >= 2) {
      sen_res <- mr_leaveoneout(dat)
      p3 <- mr_leaveoneout_plot(sen_res)
      if (length(p3) > 0) {
        ggsave(filename = file.path(folder, "sensitivity-analysis.pdf"),
               plot = p3[[1]], width = 8, height = 8)
      }
    }
  }, silent = TRUE)
  
  # 7) Funnel plot
  try({
    if (!is.null(singlesnp_res)) {
      p4 <- mr_funnel_plot(singlesnp_res)
      if (length(p4) > 0) {
        ggsave(filename = file.path(folder, "funnelplot.pdf"),
               plot = p4[[1]], width = 8, height = 8)
      }
    }
  }, silent = TRUE)
  
  
  message("Done: ", ex_id, " vs ", out_id, " -> ", folder)
}


#### Colocalization analysis ####

# GWAS summary statistics of FI
data1 <- fread("./FI.txt")
data1 <- data1 %>% dplyr::select("SNP","CHR","BP",'EA','OA',
                                 "EAF_TG","BETA_TG","SE_TG","P_TG","N")
colnames(data1) <- c('SNP','chrom',"pos",'effect_allele','other_allele', "eaf","beta","se","P","samplesize")
# computing MAF
data1$MAF <- ifelse(data1$eaf < 0.5, data1$eaf, 1 - data1$eaf)
# de-duplicating SNPs
data1 <- subset(data1, !duplicated(SNP))
# removing NAs and invalid rows
GWASdata0 <- data1 %>% na.omit() %>% filter(P > 0 & MAF > 0)

# load gene position file
pm_sel <- fread('./olink_protein_map_3k_v1.tsv') 
# proteins associated with FI from MR analysis
QTL_list <- read.table('./Protein_OlinkID_list.txt')$V1
# pQTL data
dir <- list.files('./pqtl/ukb/')  
dir <- dir[which(dir %in% QTL_list)]          

# perform colocalization analysis
coloc_result <- NULL
for (i in 1:length(dir)) {

  data1 <- fread(paste0('./pqtl/ukb/', dir[i]))
  data1$end <- data1$BP
  data1$P <- 10^(-data1$LOG10P)
  data1 <- data1 %>% dplyr::select("SNP","CHR","BP","end","ALLELE1","ALLELE0",
                                   "A1FREQ","BETA","SE","P","N")
  colnames(data1) <- c('SNP','chrom',"start","end",'effect_allele','other_allele',
                       "eaf","beta","se","P","samplesize")
  data1$MAF <- ifelse(data1$eaf < 0.5, data1$eaf, 1 - data1$eaf)
  data1$chrom <- str_replace(data1$chrom, 'chr', '')  # remove "chr" prefix
  
  
  # extract cis-pQTLs
  coords <- pm_sel %>% dplyr::filter(OlinkID == OID) %>% dplyr::slice(1)
  if (nrow(coords) == 0 || any(is.na(coords$gene_start), is.na(coords$gene_end))) {
    rm(data1); rm(lead); gc()
    next
  }
  leadstart <- as.numeric(coords$gene_start)
  leadend   <- as.numeric(coords$gene_end)
  
  # gene coordinates ±250 kb
  QTLdata <- data1[data1$start > leadstart - 250000 & data1$end < leadend + 250000, ]  
  QTLdata <- subset(QTLdata, !duplicated(SNP))
  QTLdata <- na.omit(QTLdata) %>% filter(P > 0 & MAF > 0)
  
  sameSNP <- intersect(QTLdata$SNP, GWASdata0$SNP)
  if (length(sameSNP) > 0) {
    QTLdata  <- QTLdata[QTLdata$SNP %in% sameSNP, ]    %>% dplyr::arrange(SNP) %>% na.omit()
    GWASdata <- GWASdata0[GWASdata0$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
    
    # colocalization analysis
    result <- coloc.abf(
      dataset1 = list(pvalues = GWASdata$P, snp = GWASdata$SNP, type = "quant",
                      N = GWASdata$samplesize[1], MAF = GWASdata$MAF),
      dataset2 = list(pvalues = QTLdata$P,  snp = QTLdata$SNP,  type = "quant",
                      N = QTLdata$samplesize[1]),
      MAF = QTLdata$MAF
    )
    
    result <- t(result$summary) %>% data.frame() %>% dplyr::mutate(protein = dir[i])
    coloc_result <- rbind(coloc_result, result)
    
    ## Visualization
    gwas_fn <- GWASdata[, c('SNP', 'P')] %>% dplyr::rename(rsid = SNP, pval = P)
    pqtl_fn <- QTLdata[, c('SNP', 'P')] %>% dplyr::rename(rsid = SNP, pval = P)
    pdf(paste0(dir[i], "_locuscompare.pdf"), width = 5, height = 4)
    print(locuscompare(in_fn1 = gwas_fn,
                       in_fn2 = pqtl_fn,
                       title1 = 'GWAS',
                       title2 = 'pQTL'))
    dev.off()
  }
}
write.csv(coloc_result, "coloc_result_meta.csv", row.names = FALSE)
