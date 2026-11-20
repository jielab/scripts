# 任意两组性状：LDSC + HDL；全基因组汇总数据，不使用 MR clump SNP
# trait1 / trait2 可选择整组或单个性状

library(data.table)
library(R.utils)
library(GenomicSEM)
library(HDL)

dir0 <- if (Sys.info()[["sysname"]] == "Windows") "D:" else "/mnt/d"
data_dirs <- c(
  flt = file.path(dir0, "data/gwas/flt"),
  brain = file.path(dir0, "data/gwas/brain"),
  covs = file.path(dir0, "data/gwas/covs")
)
outdir <- file.path(dir0, "projects/flt/analysis/02ldsc_hdl")

# ===== 分析选择：填组名或单个性状 ID =====
# 例如 "covs" × "flt"、"education" × "brain"、"G755" × "ad"。
# 重名性状可写 "covs/education"；每侧只填写一个组或一个性状。
trait1 <- "G732"
trait2 <- "mdd_2025"
disease_ids <- c("dementia_r13", "anystroke_eur", "sleep_disorder_r13", "mdd_2025", "anxiety_2026")

reuse_inputs <- TRUE # 共享预处理文件跨运行复用；若原始 GWAS 或筛选规则改变，请设置 FALSE 重建。
ld_dir <- file.path(dir0, "data/ldref/ldsc/hm3")
hm3 <- file.path(ld_dir, "w_hm3.snplist")
LD.path <- file.path(dir0, "data/ldref/hdl/UKB_imputed_SVD_eigen99_extraction")

resolve_dir <- function(base) {
  if (dir.exists(file.path(base, "common"))) base <- file.path(base, "common")
  normalizePath(base, winslash = "/", mustWork = TRUE)
}
discover_ids <- function(base) {
  ids <- list.dirs(base, recursive = FALSE, full.names = FALSE)
  sort(ids[file.exists(file.path(base, ids, "gwas", paste0(ids, ".gz")))])
}
data_dirs <- vapply(data_dirs, resolve_dir, character(1))
group_ids <- lapply(data_dirs, discover_ids)
missing_diseases <- setdiff(disease_ids, group_ids$brain)
if (length(missing_diseases)) stop("疾病清单中的文件不存在：", paste(missing_diseases, collapse = ", "))
group_ids$brain <- sort(unique(disease_ids))
select_traits <- function(selection) {
  if (!is.character(selection) || length(selection) != 1L || is.na(selection) || !nzchar(selection)) stop("trait1 / trait2 必须是一个组名或性状 ID")
  if (selection %in% names(group_ids)) {
    group <- selection
    ids <- group_ids[[group]]
  } else if (grepl("/", selection, fixed = TRUE)) {
    parts <- strsplit(selection, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2L || !parts[1] %in% names(group_ids)) stop("无效选择：", selection)
    group <- parts[1]
    ids <- parts[2]
    if (!ids %in% group_ids[[group]]) stop("性状不存在或未列入疾病清单：", selection)
  } else {
    groups <- names(group_ids)[vapply(group_ids, function(ids) selection %in% ids, logical(1))]
    if (!length(groups)) stop("性状不存在或未列入疾病清单：", selection)
    if (length(groups) > 1L) stop("性状重名，请使用 组名/性状ID：", selection)
    group <- groups
    ids <- selection
  }
  if (!length(ids)) stop("所选组没有可分析性状：", selection)
  list(group = group, ids = ids, base = unname(data_dirs[[group]]))
}
selection1 <- select_traits(trait1)
selection2 <- select_traits(trait2)
trait1_dir <- selection1$base
trait2_dir <- selection2$base
trait1_ids <- selection1$ids
trait2_ids <- selection2$ids
pairs <- CJ(trait1 = trait1_ids, trait2 = trait2_ids)
# 按输入组合归档，每次运行单独保存，避免不同路径或子集覆盖历史结果。
outdir <- tempfile(paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_"),
  tmpdir = file.path(outdir, paste0(selection1$group, "__", selection2$group)))
required <- c(hm3, file.path(ld_dir, paste0(1:22, ".l2.ldscore.gz")), file.path(ld_dir, paste0(1:22, ".l2.M_5_50")))
if (any(!file.exists(required))) stop("缺少 LDSC 参考：", paste(required[!file.exists(required)], collapse = "; "))
# 输入文件持久保存；日志、逐配对结果和临时解压文件按方法归档。
method_dirs <- setNames(lapply(c("ldsc", "hdl"), function(m) {
  c(inputs = file.path(outdir, m, "inputs"),
    setNames(file.path(outdir, m, c("logs", "results", "temp")),
             c("logs", "results", "temp")))
}), c("ldsc", "hdl"))
for (path in c(outdir, unlist(method_dirs, use.names = FALSE))) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  if (!dir.exists(path)) stop("无法创建输出目录：", path)
}
# 新缓存以组名命名，避免不同组中相同性状 ID 相互覆盖。
# 兼容已有的 <ID>.sumstats.gz / <ID>_hdl.rds（仅用于组间不重名的 ID）。
input_path <- function(method, group, id, reuse = reuse_inputs) {
  suffix <- if (method == "ldsc") ".sumstats.gz" else "_hdl.rds"
  root <- method_dirs[[method]][["inputs"]]
  target <- file.path(root, paste0(group, "__", id, suffix))
  if (reuse) {
    candidates <- target
    if (sum(vapply(group_ids, function(ids) id %in% ids, logical(1))) == 1L)
      candidates <- c(candidates, file.path(root, paste0(id, suffix)))
    valid <- file.exists(candidates) & !is.na(file.info(candidates)$size) & file.info(candidates)$size > 0
    if (any(valid)) return(list(path = candidates[which(valid)[1]], cached = TRUE))
  }
  list(path = target, cached = FALSE)
}
# 先写同目录临时文件，成功后发布；中断留下的临时文件不会被当成输入。
write_input <- function(path, writer) {
  staging <- tempfile(".partial_", tmpdir = dirname(path))
  on.exit(unlink(staging), add = TRUE)
  writer(staging)
  if (file.exists(path) && !file.remove(path)) stop("无法替换缓存：", path)
  if (!file.rename(staging, path)) stop("无法保存缓存：", path)
  path
}
fwrite(pairs, file.path(outdir, "planned_pairs.csv"))
fwrite(rbind(data.table(side = "trait1", id = trait1_ids, base = trait1_dir),
             data.table(side = "trait2", id = trait2_ids, base = trait2_dir)),
       file.path(outdir, "input_manifest.csv"))
# 均报告观察尺度 h2；不将病例占比当作人群患病率。
sample.prev <- population.prev <- c(NA_real_, NA_real_)
ldsc_ref <- fread(hm3, select = c("SNP", "A1", "A2"))
setnames(ldsc_ref, c("A1", "A2"), c("R1", "R2"))
ldsc_ref[, `:=`(SNP = as.character(SNP), R1 = toupper(R1), R2 = toupper(R2))]
stopifnot(!anyDuplicated(ldsc_ref$SNP))
munge_data <- function(d, trait, path) {
  d <- d[SNP %in% ldsc_ref$SNP]  # 先缩小到 HM3，再转换和合并，减少内存开销
  for (v in c("BETA", "SE", "P", "N")) set(d, j = v, value = as.numeric(d[[v]]))
  d[, `:=`(EA = toupper(trimws(EA)), NEA = toupper(trimws(NEA)))]
  d <- d[is.finite(BETA) & is.finite(P) & P >= 0 & P <= 1 & is.finite(N) & N > 0]
  d <- merge(d, ldsc_ref, by = "SNP")
  d <- d[(EA == R1 & NEA == R2) | (EA == R2 & NEA == R1)]
  d[, Z := sign(BETA) * sqrt(qchisq(P, df = 1, lower.tail = FALSE))]
  d[P == 0 & is.finite(SE) & SE > 0, Z := BETA / SE]  # 保留原来的 P 值下溢处理
  d[EA != R1, Z := -Z]
  d <- d[is.finite(Z) & !duplicated(SNP) & !duplicated(SNP, fromLast = TRUE)]
  if (!nrow(d)) stop("没有可用 SNP，请检查列名、逐 SNP N 和 HM3 匹配：", trait)
  write_input(path, function(dest) fwrite(d[, .(SNP, N, Z, A1 = R1, A2 = R2)], dest, sep = "\t", compress = "gzip"))
  message(trait, "：保留 ", nrow(d), " 个 SNP")
  path
}
panel_files <- list.files(LD.path)
bims <- list.files(LD.path, pattern = "\\.bim$", full.names = TRUE)
if (!length(bims) || !any(grepl("_snp_counter", panel_files)) || !any(grepl("_snp_list", panel_files))) stop("请下载并解压完整的 HDL 面板，修改 LD.path。")
hdl_ref <- rbindlist(lapply(bims, function(f) {
  d <- fread(f, header = FALSE, select = c(2, 5, 6))
  setnames(d, c("SNP", "R1", "R2"))
  d
}))
hdl_ref[, `:=`(SNP = as.character(SNP), R1 = toupper(R1), R2 = toupper(R2))]
stopifnot(!anyDuplicated(hdl_ref$SNP))
# 3. 输入 SNP / A1 / A2 / N / Z；N 使用文件中的逐 SNP 样本量 --------------
# EA 是 BETA 的效应等位基因；BETA 为回归系数/log(OR)。Z 显式用 BETA/SE。
read_hdl <- function(d, file) {
  d <- d[SNP %in% hdl_ref$SNP]  # 匹配 HDL 面板，不事先限制为 LDSC 的 HM3 子集
  for (v in c("BETA", "SE", "N")) set(d, j = v, value = as.numeric(d[[v]]))
  d[, `:=`(EA = toupper(trimws(EA)), NEA = toupper(trimws(NEA)))]
  d <- d[is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(N) & N > 0]
  d <- merge(d, hdl_ref, by = "SNP")
  d <- d[(EA == R1 & NEA == R2) | (EA == R2 & NEA == R1)]
  d[, Z := BETA / SE]
  d[EA != R1, Z := -Z]
  d <- d[is.finite(Z) & !duplicated(SNP) & !duplicated(SNP, fromLast = TRUE)]
  if (!nrow(d)) stop("没有可用 SNP：", file)
  as.data.frame(d[, .(SNP, A1 = R1, A2 = R2, N, Z)])
}

run_ldsc <- function(file1, file2, trait1, trait2, prefix) {
old_wd <- setwd(method_dirs$ldsc[["logs"]])
on.exit(setwd(old_wd), add = TRUE)
LDSCoutput <- GenomicSEM::ldsc(traits = c(file1, file2), sample.prev = sample.prev,
  population.prev = population.prev, ld = ld_dir, wld = ld_dir,
  trait.names = c(trait1, trait2), stand = FALSE, ldsc.log = prefix)
S <- LDSCoutput$S; V <- LDSCoutput$V; I <- LDSCoutput$I
# V 的顺序：h2_exp、rgcov、h2_out；下面的 rg_SE 与原 V_Stand 缩放等价。
scale <- if (all(is.finite(diag(S)) & diag(S) > 0)) sqrt(S[1, 1] * S[2, 2]) else NA_real_
rg <- S[1, 2] / scale
rg_se <- sqrt(V[2, 2]) / scale
rg_p <- if (is.finite(rg_se) && rg_se > 0) 2 * pnorm(abs(rg / rg_se), lower.tail = FALSE) else NA_real_
line <- grep("Cross trait Intercept:", readLines(paste0(prefix, "_ldsc.log"), warn = FALSE), value = TRUE)
int_se <- if (length(line) == 1L) suppressWarnings(as.numeric(sub(".*\\(([^)]+)\\).*", "\\1", line))) else NA_real_
result <- data.table(trait1 = trait1, trait2 = trait2,
  h2_exp = S[1, 1], h2_exp_se = sqrt(V[1, 1]), h2_out = S[2, 2], h2_out_se = sqrt(V[3, 3]),
  int_exp = I[1, 1], int_out = I[2, 2], gcov_int = I[1, 2], gcov_int_se = int_se,
  rgcov = S[1, 2], rgcov_se = sqrt(V[2, 2]), rg = rg, rg_se = rg_se, rg_p = rg_p)

result
}
# pairs 按 trait1 排序；仅缓存当前 trait1，避免每个 trait2 配对重复读盘。
# 缓存限定为一个数据框，不将所有 GWAS 常驻内存。
hdl_cache <- new.env(parent = emptyenv())
hdl_cache$file <- NULL
run_hdl <- function(file1, file2, trait1, trait2, prefix) {
old_wd <- setwd(method_dirs$hdl[["logs"]])
on.exit(setwd(old_wd), add = TRUE)
if (!identical(hdl_cache$file, file1)) {
  hdl_cache$data <- readRDS(file1)
  hdl_cache$file <- file1
}
log_file <- paste0(prefix, "_hdl.log")
if (file.exists(log_file)) unlink(log_file)  # 避免包内部 rm 命令影响 Windows 路径
fit <- HDL::HDL.rg(gwas1.df = hdl_cache$data, gwas2.df = readRDS(file2), LD.path = LD.path, output.file = log_file, intercept.output = TRUE)
e <- fit$estimates.df
result <- data.table(trait1 = trait1, trait2 = trait2,
  h2_exp = e["Heritability_1", "Estimate"], h2_exp_se = e["Heritability_1", "se"],
  h2_out = e["Heritability_2", "Estimate"], h2_out_se = e["Heritability_2", "se"],
  int_exp = e["Heritability_1_intercept", "Estimate"], int_out = e["Heritability_2_intercept", "Estimate"],
  gcov_int = e["Genetic_Covariance_intercept", "Estimate"], gcov_int_se = e["Genetic_Covariance_intercept", "se"],
  rgcov = e["Genetic_Covariance", "Estimate"], rgcov_se = e["Genetic_Covariance", "se"],
  rg = fit$rg, rg_se = fit$rg.se, rg_p = fit$P)

result
}

# 先生成各方法 inputs 文件，再从这些路径读取进行配对分析。
# 优先复用共享 inputs；两种方法都命中时，不读取或解压原始 GWAS。
main <- function() {
  on.exit({
    hdl_cache$file <- NULL
    hdl_cache$data <- NULL
  }, add = TRUE)
  # GenomicSEM/readr 使用 vroom 读取压缩文件；其解压缓存也放到 D 盘。
  tmp <- tempfile("run_", tmpdir = method_dirs$ldsc[["temp"]])
  if (!dir.create(tmp)) stop("无法创建解压缓存目录：", tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  old_vroom_temp <- Sys.getenv("VROOM_TEMP_PATH", unset = NA_character_)
  Sys.setenv(VROOM_TEMP_PATH = normalizePath(tmp, winslash = "/", mustWork = TRUE))
  on.exit({
    if (is.na(old_vroom_temp)) Sys.unsetenv("VROOM_TEMP_PATH") else
      Sys.setenv(VROOM_TEMP_PATH = old_vroom_temp)
  }, add = TRUE)
  inputs <- rbind(data.table(id = unique(pairs$trait1), base = trait1_dir, side = "trait1", group = selection1$group),
                  data.table(id = unique(pairs$trait2), base = trait2_dir, side = "trait2", group = selection2$group))
  prepared <- list(ldsc = list(), hdl = list())
  failures <- data.table(method = character(), trait1 = character(), trait2 = character(), error = character())
  results <- list()
  for (i in seq_len(nrow(inputs))) {
    id <- inputs$id[i]
    key <- paste0(inputs$side[i], "__", id)
    file <- file.path(inputs$base[i], id, "gwas", paste0(id, ".gz"))
    paths <- setNames(lapply(c("ldsc", "hdl"), function(m)
      input_path(m, inputs$group[i], id)), c("ldsc", "hdl"))
    raw <- NULL
    if (any(!vapply(paths, function(p) p$cached, logical(1)))) {
      # 仅当至少一种格式缺失时读取一次原始文件。
      raw <- tryCatch({
        d <- fread(file, select = c("SNP", "EA", "NEA", "BETA", "SE", "P", "N"))
        d[, SNP := trimws(as.character(SNP))]
        d
      }, error = function(e) list(error = conditionMessage(e)))
    }
    for (method in c("ldsc", "hdl")) {
      prepared[[method]][[key]] <- tryCatch({
        path <- paths[[method]]$path
        if (paths[[method]]$cached) {
          message("复用 ", method, " 输入：", path)
          path
        } else {
          if (!is.data.table(raw)) stop(raw$error)
          if (method == "ldsc") munge_data(raw, id, path) else
            write_input(path, function(dest) saveRDS(read_hdl(raw, file), dest, compress = FALSE))
        }
      }, error = function(e) list(error = conditionMessage(e)))
    }
  }
  rm(raw)
  for (i in seq_len(nrow(pairs))) {
    x <- pairs$trait1[i]; y <- pairs$trait2[i]
    for (method in c("ldsc", "hdl")) {
      message(i, "/", nrow(pairs), " ", method, ": ", x, " × ", y)
      prefix <- file.path(method_dirs[[method]][["logs"]], paste0(x, "__", y))
      result_file <- file.path(method_dirs[[method]][["results"]], paste0(x, "__", y, "_results.csv"))
      result <- tryCatch({
        a <- prepared[[method]][[paste0("trait1__", x)]]; b <- prepared[[method]][[paste0("trait2__", y)]]
        if (is.list(a)) stop(x, ": ", a$error)
        if (is.list(b)) stop(y, ": ", b$error)
        f <- if (method == "ldsc") run_ldsc else run_hdl
        z <- f(a, b, x, y, prefix)
        z[, method := method]
        fwrite(z, result_file)
        z
      }, error = function(e) {
        failures <<- rbind(failures, data.table(method = method, trait1 = x, trait2 = y, error = conditionMessage(e)))
        message(conditionMessage(e))
        NULL
      })
      results[[length(results) + 1L]] <- list(result)
      fwrite(failures, file.path(outdir, "failures.csv"))
    }
  }
  result <- rbindlist(unlist(results, recursive = FALSE), fill = TRUE)
  if (!nrow(result)) stop("全部配对失败，详见 failures.csv")
  result[, rg_p := fifelse(is.finite(rg_p) & rg_p >= 0 & rg_p <= 1, rg_p, NA_real_)]
  # 每种方法对本次全部 trait1 × trait2 配对统一校正；失败配对仍计入检验数。
  result[, FDR := p.adjust(rg_p, "BH", n = nrow(pairs)), by = method]
  fwrite(result, file.path(outdir, "ldsc_hdl_results.csv"))
  for (m in c("ldsc", "hdl")) fwrite(result[method == m], file.path(outdir, m, paste0(m, "_results.csv")))
  fwrite(result[is.finite(FDR) & FDR < 0.05], file.path(outdir, "ldsc_hdl_hits.csv"))
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
  message("结果已保存：", outdir)
  invisible(result)
}
result <- main()
