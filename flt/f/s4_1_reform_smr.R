source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
library(data.table)

## 1、统一处理 Chicken liking（G732）和抑郁症（mdd_2025）。
gwas_dirs <- c(
  G732 = file.path(Sys.getenv("FLT_GWAS"), "flt/common/G732/gwas"),
  mdd_2025 = file.path(Sys.getenv("FLT_GWAS"), "brain/common/mdd_2025/gwas")
)

for (id in names(gwas_dirs)) {
  outdir <- gwas_dirs[[id]]
  infile <- file.path(outdir, paste0(id, ".gz"))
  outfile <- file.path(Sys.getenv("FLT_SMR_INPUT"), paste0(id, ".ma"))
  if (!file.exists(infile)) stop("找不到清洗后的 GWAS：", infile)
  message("正在整理：", id)
  # 将清洗后的 GWAS 整理为 SMR 格式： 左边是 SMR 列名，右边是清洗后的 GWAS 列名。
  cols <- c(SNP = "SNP", A1 = "EA", A2 = "NEA", freq = "EAF", b = "BETA", se = "SE", p = "P", n = "N")
  header <- names(fread(file = infile, nrows = 0L))
  missing_cols <- setdiff(unname(cols), header)
  if (length(missing_cols)) {stop(id, " 输入文件缺少列：", paste(missing_cols, collapse = ", "))}
  d <- fread(file = infile, select = unname(cols))
  setnames(d, old = unname(cols), new = names(cols))
  input_rows <- nrow(d)
  d[, `:=`(SNP = trimws(SNP), A1 = toupper(trimws(A1)),  A2  = toupper(trimws(A2)))]
  for (col in c("freq", "b", "se", "p", "n")) {set(d, j = col, value = suppressWarnings(as.numeric(d[[col]])))}
  # 缺失或超出有效范围的 EAF 和样本量写为 NA
  d[!is.finite(freq) | freq <= 0 | freq >= 1, freq := NA_real_]
  d[!is.finite(n) | n <= 0, n := NA_real_]
  # 检查 rsID、等位基因和效应统计量
  d <- d[grepl("^rs[0-9]+$", SNP) & grepl("^[ACGT]+$", A1) & grepl("^[ACGT]+$", A2) & A1 != A2 & is.finite(b) & is.finite(se) & se > 0 & is.finite(p) & p >= 0 & p <= 1]
  # 完全相同的记录只保留一条。仍有冲突的重复 rsID 整组剔除
  d <- unique(d)
  conflicts <- unique(d$SNP[duplicated(d$SNP)])
  d <- d[!SNP %in% conflicts]
  if (!nrow(d)) stop(id, "：质控后没有可导出的 SNP。")
  setcolorder(d, names(cols))
  tmp <- tempfile(tmpdir=dirname(outfile))
  fwrite(d, file = tmp, sep = "\t", quote = FALSE, na = "NA")
  if (!file.rename(tmp,outfile)) stop("Could not publish ",outfile)
  valid_eaf <- sum(!is.na(d$freq))
  print(c(
    input_rows           = input_rows,
    output_rows          = nrow(d),
    conflicting_rsids    = length(conflicts),
    valid_EAF            = valid_eaf,
    missing_EAF          = nrow(d) - valid_eaf,
    missing_sample_size  = sum(is.na(d$n))
  ))
  cat("输入文件：", infile, "\n输出文件：", outfile, "\n", sep = "")
  if (valid_eaf == 0L) {warning(id, "：导出的 freq 全部为 NA，请核对源文件的 EAF。")}
}
