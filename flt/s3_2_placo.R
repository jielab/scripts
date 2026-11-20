# G732和mdd_2025的LDSC(P < 0.05)

library(data.table)
library(R.utils)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 PLACO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir0 <- if (Sys.info()[["sysname"]] == "Windows") "D:" else "/mnt/d"
project_dir <- file.path(dir0, "projects/flt")
outdir <- file.path(project_dir, "analysis/03placo")

trait1 <- "G732"
trait2 <- "mdd_2025"
files <- c(file.path(dir0, "data/gwas/flt/common", trait1, "gwas", paste0(trait1, ".gz")),
           file.path(dir0, "data/gwas/brain/common", trait2, "gwas", paste0(trait2, ".gz")))
build <- c("GRCh38", "GRCh38") # 填实际版本；这里不会做liftOver。
method <- "PLACO"             # 可选"PLACO+"；不自动切换方法。
z_sq_max <- 80
p_null <- 1e-4                # 估计参数用的阈值
p_cut <- 1e-5                 # 筛选候选SNP的阈值
maf_min <- 0.01
batch_size <- 10000L
resume <- TRUE               # 相同输入和参数下，续跑已完成的批次。
# BETA对应EA；疾病的BETA必须是log(OR)。

# 主流程 ------------------------------------------------------------
run_placo <- function() {
  stopifnot(trait1 != trait2, all(file.exists(files)), length(build) == 2L, !anyNA(build), build[1] == build[2], method %in% c("PLACO", "PLACO+"), batch_size >= 1L, maf_min >= 0, maf_min < 0.5)
  impl <- load_placo()
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(outdir)) stop("Cannot create results directory: ", outdir)
  prefix <- file.path(outdir, paste0(trait1, "__", trait2, "_", gsub("+", "plus", method, fixed = TRUE)))
  qc <- list()
  qc_add <- function(trait, step, n) {
    qc[[length(qc) + 1L]] <<- data.table(trait, step, n)
    message(sprintf("QC | %-12s | %-36s | n = %s", trait, step, format(n, big.mark = ",")))
    fwrite(rbindlist(qc), paste0(prefix, "_qc.tsv"), sep = "\t")
  }
  # 2. 读取完整GWAS并协调等位基因 -------------------------------------
  a <- read_placo(files[1], trait1, qc_add)
  b <- read_placo(files[2], trait2, qc_add)
  d <- harmonise_placo(a, b, qc_add)
  rm(a, b)
  extreme <- d$Z1^2 > z_sq_max | d$Z2^2 > z_sq_max
  fwrite(d[extreme], paste0(prefix, "_excluded_extremeZ.tsv.gz"), sep = "\t")
  d <- d[!extreme]
  qc_add("pair", "Z2_filtered", nrow(d))
  fwrite(rbindlist(qc), paste0(prefix, "_qc.tsv"), sep = "\t")
  if (nrow(d) < 1000L) {
    print(rbindlist(qc))
    stop("Only ", nrow(d), " SNPs remain after QC. This does NOT identify the cause. ",
         "See the first sharp decrease in the QC table above and ", paste0(prefix, "_qc.tsv"),
         "; check missing frequencies, shared IDs, coordinates, alleles, and Z values. ",
         "PLACO was not run; the small-SNP safeguard has not been relaxed.")
  }
  if (nrow(d) < 100000L) warning("Low SNP coverage; check that the full GWAS, not an HM3 subset, was supplied.")

  # 3. 估计Z方差及误差相关 ---------------------------------------------
  zmat <- as.matrix(d[, .(Z1, Z2)])
  pmat <- 2 * pnorm(abs(zmat), lower.tail = FALSE)
  varz <- impl$env$var.placo(zmat, pmat, p.threshold = p_null)
  null_both <- pmat[, 1] >= p_null & pmat[, 2] >= p_null
  if (sum(null_both) < 1000L) stop("Too few marginally null variants to estimate Z-error correlation.")
  rho <- cor(zmat[null_both, 1], zmat[null_both, 2])
  stopifnot(all(is.finite(varz) & varz > 0), is.finite(rho), abs(rho) < 1)
  message("Testing ", nrow(d), " SNPs; VarZ = ", paste(signif(varz, 5), collapse = ", "),
          "; null-SNP Z correlation = ", signif(rho, 5))
  if (method == "PLACO" && abs(rho) > 0.1) warning("Non-negligible null-SNP Z correlation. 0.1 is a diagnostic flag, not an overlap test. ", "Check cohort overlap and consider method <- 'PLACO+'. Do not substitute LDSC rg for rho.")

  # 4. 分批PLACO检验；保留断点缓存及输入校验 ----------------------------
  tested <- batch_pvalues(d, zmat, varz, rho, impl, prefix)
  pval <- tested$pval
  sig <- tested$signature

  # 5. 导出全部SNP、候选SNP、FUMA及QC ----------------------------------
  # hits是逐SNP候选；独立位点需继续做FUMA。
  d[, `:=`(P_PLACO_RAW = pval, P_PLACO = pmin(pval, 1), T_PLACO = Z1 * Z2)]
  d[, direction := fifelse(BETA1 * BETA2 > 0, "concordant", fifelse(BETA1 * BETA2 < 0, "discordant", "zero"))]
  d[, `:=`(trait1 = trait1, trait2 = trait2, method = method)]
  fwrite(d, paste0(prefix, "_all.tsv.gz"), sep = "\t")
  hits <- d[P_PLACO < p_cut][order(P_PLACO)]
  fwrite(hits, paste0(prefix, "_hits.tsv"), sep = "\t")
  fuma <- d[grepl("^rs[0-9]+$", SNP), .(SNP, A1 = EA, A2 = NEA, P = pmax(P_PLACO, 1e-300))]
  if (nrow(fuma)) {fwrite(fuma, paste0(prefix, "_FUMA.txt.gz"), sep = "\t")} else warning("No rsIDs: no FUMA file exported. Map IDs or supply GRCh37 coordinates first.")
  n_concordant <- sum(hits$direction == "concordant")
  summary <- data.table(trait1, trait2, method, build = build[1], n_tested = nrow(d),
    varZ1 = varz[1], varZ2 = varz[2], rho_null_Z = rho, n_null_both = sum(null_both),
    p_cut, n_hits = nrow(hits), n_concordant,
    pct_concordant = if (nrow(hits)) 100 * n_concordant / nrow(hits) else NA_real_,
    n_P_lt_5e8 = sum(d$P_PLACO < 5e-8), n_p_capped_at_1 = sum(pval > 1),
    n_FUMA = nrow(fuma), n_FUMA_P_floor = sum(d$P_PLACO < 1e-300 & grepl("^rs[0-9]+$", d$SNP)))
  fwrite(summary, paste0(prefix, "_summary.tsv"), sep = "\t")
  writeLines(c(paste("Official source:", impl$url),
               paste("Source MD5:", impl$source_md5),
               paste("Input:", files, "MD5:", sig$inputs_md5),
               "FUMA: SNP2GENE; SNP/P columns; 1000G Phase3 EUR; MAF 0.01; discovery P < 1e-5; independent SNP r2=0.6; lead SNP r2=0.1; merge distance=250 kb.",
               "FUMA input contains PLACO P, not marginal trait P. Do not map a joint PLACO BETA/OR: it does not exist.",
               "FUMA lead/independent SNP and risk-locus counts are not the raw PLACO hit count.",
               capture.output(sessionInfo())), paste0(prefix, "_run_info.txt"))
  print(summary)
  message("Saved: ", prefix)
  invisible(summary)
}

# 数据读取与质控 ----------------------------------------------------
read_placo <- function(file, label, qc_add) {
  hdr <- names(fread(file, nrows = 0L))
  required <- c("SNP", "EA", "NEA", "BETA", "SE", "P")
  if (!all(required %in% hdr)) stop(label, ": missing columns: ", paste(setdiff(required, hdr), collapse = ", "))
  pick <- function(aliases) {
    h <- hdr[toupper(hdr) %in% aliases]
    if (length(h)) h[1L] else NA_character_
  }
  chr <- pick(c("CHR", "CHROM", "CHROMOSOME"))
  pos <- pick(c("POS", "BP", "POSITION"))
  af_cols <- hdr[toupper(hdr) %in% c("EAF", "MAF")]
  d <- fread(file, select = unique(c(required, na.omit(c(chr, pos)), af_cols)))
  qc_add(label, "read", nrow(d))
  message(label, ": input = ", file)
  d[, `:=`(SNP = trimws(as.character(SNP)), EA = toupper(trimws(EA)), NEA = toupper(trimws(NEA)))]
  for (v in c("BETA", "SE", "P")) set(d, j = v, value = suppressWarnings(as.numeric(d[[v]])))
  d <- d[!is.na(SNP) & nzchar(SNP) & SNP != "."]
  qc_add(label, "nonmissing_id", nrow(d))
  d <- d[!duplicated(SNP) & !duplicated(SNP, fromLast = TRUE)]
  qc_add(label, "unique_nonmissing_id", nrow(d))
  d <- d[EA %chin% c("A", "C", "G", "T") & NEA %chin% c("A", "C", "G", "T") & EA != NEA & is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P >= 0 & P <= 1]
  qc_add(label, "valid_snv_and_statistics", nrow(d))
  d <- d[!paste0(EA, NEA) %chin% c("AT", "TA", "CG", "GC")]
  qc_add(label, "valid_nonpalindromic_snv", nrow(d))
  af <- NA_character_
  best_count <- 0L
  for (column in af_cols) {
    raw <- trimws(as.character(d[[column]]))
    supplied <- !is.na(raw) & !toupper(raw) %chin% c("", ".", "NA", "NAN", "NULL")
    f <- suppressWarnings(as.numeric(raw))
    upper <- if (toupper(column) == "MAF") 0.5 else 1
    valid <- is.finite(f) & f >= 0 & f <= upper
    message(label, ": ", column, " valid frequency = ", sum(valid), "/", nrow(d),
            "; missing = ", sum(!supplied))
    if (any(supplied & !valid)) {
      stop(label, ": invalid ", column, " values (expected proportions, not percentages): ",
           paste(head(unique(raw[supplied & !valid]), 5L), collapse = ", "),
           ". See QC above; fix the frequency column before testing.")
    }
    if (sum(valid) > best_count) {
      af <- column
      best_count <- sum(valid)
    }
  }
  if (!is.na(af)) {
    f <- suppressWarnings(as.numeric(d[[af]]))
    upper <- if (toupper(af) == "MAF") 0.5 else 1 - maf_min
    d <- d[is.finite(f) & f >= maf_min & f <= upper]
    qc_add(label, paste0("maf_filtered_using_", af), nrow(d))
  } else {
    warning(label, ": EAF/MAF absent or entirely missing; frequency filtering was NOT applied. ",
            "SNPs were retained. Confirm upstream common-variant QC; missing frequencies do not establish MAF >= 0.01.")
    qc_add(label, "maf_unavailable_not_filtered", nrow(d))
  }
  d[, Z := BETA / SE]
  d <- d[is.finite(Z)]
  chr_values <- if (!is.na(chr)) toupper(sub("^chr", "", trimws(as.character(d[[chr]])), ignore.case = TRUE)) else rep(NA_character_, nrow(d))
  chr_numeric <- suppressWarnings(as.numeric(chr_values))
  is_integer <- is.finite(chr_numeric) & chr_numeric == trunc(chr_numeric)
  chr_values[is_integer] <- as.character(as.integer(chr_numeric[is_integer]))
  pos_values <- if (!is.na(pos)) suppressWarnings(as.numeric(d[[pos]])) else rep(NA_real_, nrow(d))
  set(d, j = "CHR_OUT", value = chr_values)
  set(d, j = "POS_OUT", value = pos_values)
  if (!is.na(chr)) d <- d[CHR_OUT %chin% as.character(1:22)]
  qc_add(label, "usable", nrow(d))
  d[, .(SNP, CHR = CHR_OUT, POS = POS_OUT, EA, NEA, BETA, SE, P, Z)]
}

# 等位基因协调 ------------------------------------------------------
harmonise_placo <- function(a, b, qc_add) {
  d <- merge(a, b, by = "SNP", suffixes = c("1", "2"), sort = TRUE)
  qc_add("pair", "shared_id", nrow(d))
  bad <- (!is.na(d$CHR1) & !is.na(d$CHR2) & d$CHR1 != d$CHR2) | (is.finite(d$POS1) & is.finite(d$POS2) & d$POS1 != d$POS2)
  if (any(bad)) {
    message("Coordinate mismatch examples (verify the actual genome builds):")
    print(head(d[bad, .(SNP, CHR1, POS1, CHR2, POS2)], 5L))
  }
  d <- d[!bad]
  qc_add("pair", "matching_coordinates", nrow(d))
  cEA2 <- chartr("ACGT", "TGCA", d$EA2)
  cNEA2 <- chartr("ACGT", "TGCA", d$NEA2)
  same <- (d$EA1 == d$EA2 & d$NEA1 == d$NEA2) | (d$EA1 == cEA2 & d$NEA1 == cNEA2)
  flip <- (d$EA1 == d$NEA2 & d$NEA1 == d$EA2) | (d$EA1 == cNEA2 & d$NEA1 == cEA2)
  d[, flip2 := flip]
  d <- d[same | flip]
  d[flip2 == TRUE, `:=`(BETA2 = -BETA2, Z2 = -Z2)]
  d[, `:=`(EA = EA1, NEA = NEA1)]
  d[, c("EA1", "NEA1", "EA2", "NEA2") := NULL]
  qc_add("pair", "harmonised", nrow(d))
  d
}

# 官方PLACO / PLACO+函数（离线内嵌，统计公式未改）--------------------
load_placo <- function() {
  env <- new.env(parent = globalenv())
  evalq({
    message("====================================================================")
    message("                  PLACO v0.2.0 is loaded")
    message("====================================================================")
    message("If you use this software, please cite:")
    message("")
    message("Original PLACO paper -->")
    message("Ray and Chatterjee (2020) A powerful method for pleiotropic analysis")
    message("    under composite null hypothesis identifies novel shared")
    message("    loci between type 2 diabetes and prostate cancer.")
    message("    PLoS Genetics 16(12): e1009218")
    message("")
    message("PLACO+ paper -->")
    message("Park and Ray (2025+) A robust pleiotropy method with")
    message("    applications to lipid traits and to inflammatory bowel")
    message("    disease subtypes with sample overlap. Submitted.")
    message("")
    message("**See https://github.com/RayDebashree/PLACO for updated citations")
    message("--------------------------------------------------------------------")
    message("")
    .pdfx<-function(x) besselK(x=abs(x),nu=0)/pi
    .p.bessel<-function(z, varz, AbsTol=1e-13){
    	p1 <- 2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[1])), Inf, abs.tol=AbsTol)$value)
    	p2 <- 2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]/sqrt(varz[2])), Inf, abs.tol=AbsTol)$value)
    	p0 <- 2*as.double(integrate(Vectorize(.pdfx), abs(z[1]*z[2]), Inf, abs.tol=AbsTol)$value)
    	pval.compnull <- p1 + p2 - p0
    	return(pval.compnull)
    }
    .pdfx.cor<-function(x, corz){ 
      f1 <- besselK(x=(abs(x)/(1-corz^2)), nu=0)
      f2 <- exp(corz*x/(1-corz^2))/(pi*sqrt(1-corz^2))
      y<-f1*f2
      for(i in 1:length(y)){
        if(f1[i]==0){ y[i]<-0 }else{ y[i] <- f1[i]*f2[i] }
      }
      return(y)}
    .p.bessel.cor<-function(z, varz, corz, AbsTol=1e-13){
      p1 <- as.double(integrate(Vectorize(.pdfx.cor), abs(z[1]*z[2]/sqrt(varz[1])), Inf, corz=corz, abs.tol=AbsTol)$value + integrate(Vectorize(.pdfx.cor), -Inf, -abs(z[1]*z[2]/sqrt(varz[1])), corz=corz, abs.tol=AbsTol)$value)
      p2 <- as.double(integrate(Vectorize(.pdfx.cor), abs(z[1]*z[2]/sqrt(varz[2])), Inf, corz=corz, abs.tol=AbsTol)$value + integrate(Vectorize(.pdfx.cor), -Inf, -abs(z[1]*z[2]/sqrt(varz[2])), corz=corz, abs.tol=AbsTol)$value)
      p0 <- as.double(integrate(Vectorize(.pdfx.cor), abs(z[1]*z[2]), Inf, corz=corz, abs.tol=AbsTol)$value + integrate(Vectorize(.pdfx.cor), -Inf, -abs(z[1]*z[2]), corz=corz, abs.tol=AbsTol)$value)
      pval.compnull <- p1 + p2 - p0
      return(pval.compnull)
    }
    var.placo<-function(Z.matrix, P.matrix, p.threshold=1e-4){
      k<-ncol(Z.matrix)
      if(k!=2) stop("This method is meant for 2 traits only. Columns correspond to traits.")
      ZP<-cbind(Z.matrix,P.matrix)
      ZP<-na.omit(ZP)
      rows.alt<-which(ZP[,3]<p.threshold & ZP[,4]<p.threshold)
      if(length(rows.alt)>0){
        ZP<-ZP[-rows.alt,]
        if(nrow(ZP)==0) stop(paste("No 'null' variant left at p-value threshold",p.threshold))
        if(nrow(ZP)<30) warning(paste("Too few 'null' variants at p-value threshold",p.threshold))
      }
      varz<-diag(var(ZP[,c(1,2)]))
      return(varz)
    }
    cor.pearson<-function(Z.matrix, P.matrix, p.threshold=1e-4, returnMatrix=TRUE){
      k <- ncol(Z.matrix)
      if(k!=2) stop("This method is meant for 2 traits only.")
        row.exclude <- which( apply(P.matrix, MARGIN = 1, function(x) any(x < p.threshold)) == TRUE )
        if(length(row.exclude)>0) Z.matrix<-Z.matrix[-row.exclude,]
      R <- cor(Z.matrix)
      if(isTRUE(returnMatrix)){ return(R) }else{ return(R[1,2]) }
    }
    placo <- function(Z, VarZ, AbsTol=.Machine$double.eps^0.8){
       k <- length(Z)		
       if(k!=2) stop("This method is meant for 2 traits only.")
       if(length(VarZ)!=k) stop("Provide variance estimates for 2 traits as obtained using var.placo() function.")
       pvalue.b = .p.bessel(z=Z, varz=VarZ, AbsTol=AbsTol)
       return(list(T.placo=prod(Z), p.placo=pvalue.b))
    }
    placo.plus <- function(Z, VarZ, CorZ, AbsTol=.Machine$double.eps^0.8){
      k <- length(Z)		
      if(k!=2) stop("This method is meant for 2 traits only.")
      if(length(VarZ)!=k) stop("Provide variance estimates for 2 traits as obtained using var.placo() function.")
      if(is.matrix(CorZ)) stop("CorZ is estimated correlation parameter (scalar), not a matrix. Use cor.pearson() function with returnMatrix=FALSE to obtain it.")
      if(CorZ < -1 | CorZ > 1) stop("Estimated correlation parameter (CorZ) must be between -1 and 1. Use cor.pearson() function to obtain it.")
      pvalue.b = .p.bessel.cor(z=Z, varz=VarZ, corz=CorZ, AbsTol=AbsTol)
      return(list(T.placo.plus=prod(Z), p.placo.plus=pvalue.b))
    }
  }, envir = env)
  stopifnot(is.function(env$placo), is.function(env$placo.plus), is.function(env$var.placo))
  message("PLACO functions ready (embedded offline version; no download).")
  list(env = env, source_md5 = "ca5b88d24031f02ba9c8f217f52621aa",
       url = "https://github.com/RayDebashree/PLACO/blob/3ba3cae1d323ad117fb4540e620bcefa79f70663/PLACO_v0.2.0.R")
}

# 按原规则校验输入并续跑数值积分；已完成批次可复用。
batch_pvalues <- function(d, zmat, varz, rho, impl, prefix) {
  cache_dir <- paste0(prefix, "_cache")
  dir.create(cache_dir, showWarnings = FALSE)
  sig <- list(pipeline = "s3_3_placo_20261103_offline_v3", inputs_md5 = unname(tools::md5sum(files)),
              implementation_md5 = impl$source_md5, build = build,
              method = method, z_sq_max = z_sq_max, p_null = p_null, maf_min = maf_min,
              batch_size = batch_size, varz = varz, rho = rho, n = nrow(d))
  sig_tmp <- tempfile(tmpdir = cache_dir)
  saveRDS(d[, .(SNP, EA, NEA, Z1, Z2)], sig_tmp, compress = FALSE)
  sig$harmonised_md5 <- unname(tools::md5sum(sig_tmp))
  unlink(sig_tmp)
  manifest <- file.path(cache_dir, "signature.rds")
  if (!resume || !file.exists(manifest) || !identical(readRDS(manifest), sig)) {
    unlink(list.files(cache_dir, pattern = "^batch_[0-9]+\\.rds$", full.names = TRUE))
    saveRDS(sig, manifest)
  }
  test_one <- function(z) {
    for (tol in c(1e-13, 1e-15)) {
      ans <- tryCatch({
        if (method == "PLACO") impl$env$placo(z, VarZ = varz, AbsTol = tol)$p.placo
        else impl$env$placo.plus(z, VarZ = varz, CorZ = rho, AbsTol = tol)$p.placo.plus
      }, error = function(e) NA_real_)
      if (length(ans) == 1L && is.finite(ans) && ans >= 0) return(ans)
    }
    NA_real_
  }
  pval <- rep(NA_real_, nrow(d))
  for (start in seq.int(1L, nrow(d), by = batch_size)) {
    idx <- seq.int(start, min(start + batch_size - 1L, nrow(d)))
    cache_file <- file.path(cache_dir, sprintf("batch_%09d.rds", start))
    ans <- if (resume && file.exists(cache_file)) tryCatch(readRDS(cache_file), error = function(e) NULL) else NULL
    if (is.null(ans) || length(ans) != length(idx) || any(!is.finite(ans) | ans < 0)) {
      ans <- vapply(idx, function(i) test_one(zmat[i, ]), numeric(1))
      if (anyNA(ans)) {
        fwrite(d[idx[is.na(ans)]], paste0(prefix, "_integration_failures.tsv"), sep = "\t")
        stop("PLACO integration failed; SNPs saved. Earlier batches remain cached; no final result was written.")
      }
      tmp <- paste0(cache_file, ".tmp")
      saveRDS(ans, tmp)
      if (!file.rename(tmp, cache_file)) stop("Could not save batch: ", cache_file)
    }
    pval[idx] <- ans
    message("PLACO: ", max(idx), "/", nrow(d), " SNPs completed")
  }
  list(pval = pval, signature = sig)
}

# 执行：结果保存在outdir，控制台打印汇总。
result <- run_placo()
## 结果解读
# n_tested	实际接受 PLACO 检验的 SNP 数
# n_hits	P_PLACO < 1×10⁻⁵的候选共享 SNP
# n_P_lt_5e8	达到更严格阈值的 SNP，仍可能彼此存在 LD
# n_concordant	n_hits个候选中，两性状效应方向一致的 SNP
# pct_concordant	 n_concordant/n_hits；其余n_hits-n_concordant 个方向相反
# varZ1、varZ2	PLACO 用于校准检验的方差参数，不是遗传力
# rho_null_Z	背景 SNP 的 Z 值相关，不是 LDSC/HDL 的遗传相关rg


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 下一步 FUMA ("D:\projects\flt\analysis\03placo\FUMA_job787285")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 将 XX_XX_PLACO_FUMA.txt.gz 文件上传到fuma（ https://fuma.ctglab.nl/snp2gene ）中
## FUMA 参数 
# FUMA 参数 | 设置 | 依据或说明 |
# GWAS summary statistics | 完整 PLACO 文件 | 论文方法 |
# Pre-defined lead SNPs | 不上传 | 建议由当前 PLACO 结果识别 |
# Pre-defined genomic regions | 不上传 | 建议不预先限制区域 |
# Identify additional lead SNPs | 保持自动识别 | 未提供预定义 lead 时适用 |
# Maximum lead SNP P-value | **`1e-5`** | 对应论文共享位点筛选阈值 |
# Maximum GWAS P-value | **`0.05`** | 默认/建议；限制 LD 候选 SNP |
# Independent significant SNPs 的 r² | **`0.6`** | 论文报告 |
# Lead SNPs 的第二个 r² | **`0.1`** | 论文报告 |
# Reference panel | **1000 Genomes Phase 3 EUR** | 与论文使用的欧洲 LD 参考一致 |
# Include variants from reference panel | **Yes** | 默认/建议；作者具体开关待配置文件确认 |
# Minimum MAF | **`0`** | 官方文档默认；作者具体值未报告 |
# Maximum distance of LD blocks to merge | **`250 kb`** | 默认/建议；论文也描述了 250 kb 距离设置 |
# Gene type | **Protein coding** | 默认/建议，与论文基因分析描述一致 |
# Exclude MHC | **勾选** | 论文报告排除注释中的 MHC |
# MHC exclusion scope | **Only from annotations** | 附表支持的设置，见下文 |
# Extended MHC region | 留空 | 使用默认 MOG–COL11A2 区域；未见作者提供自定义范围 |

## 基因映射参数：用于得到 Data 17
# 参数 | 设置 |
# Positional mapping | **开启** |
# Mapping criterion | **按距离映射** |
# Maximum distance to genes | **10 kb** |
# eQTL mapping | **开启** |
# eQTL 数据来源① | **eQTLGen cis-eQTLs** |
# eQTL 数据来源② | **GTEx v8 → Whole Blood** |
# Use only significant SNP–gene pairs | **勾选，对应数据源 FDR≤0.05** |
# Chromatin interaction mapping | **关闭**，本文这部分没有使用 |
# CADD filtering | **不启用**，论文未报告此筛选 |
# RegulomeDB filtering | **不启用**，论文未报告此筛选 |
# Chromatin-state filtering | **不启用**，论文未报告此筛选 |
# 其他功能注释筛选 | 不额外增加 |

## MAGMA 参数：用于得到 Data 18、23 和组织表达结果
# 参数 | 设置 | 说明 |
# **Perform MAGMA** | **勾选** | 当前界面需要明确开启 |
# LD reference | **1000G Phase 3 EUR** | 与前面一致 |
# Gene model | **SNP-wise mean** | 论文使用的方法 |
# **MAGMA gene annotation window** | **0 kb，作为默认起点** | 作者未明确报告，需配置确认 |
# MAGMA gene expression dataset | **GTEx v8** | 论文组织表达分析使用该版本 |
# Gene-set analysis | 保留/运行 | 使用完整 MAGMA 基因结果 |
# Gene-set collection | Curated gene sets＋GO terms | 复现时还需匹配数据库版本 |
# MHC in MAGMA | **保留** | 与 Data 18 的 HLA 结果一致 |

## 下载数据 FUMA 输出的是 GRCh37/hg19

