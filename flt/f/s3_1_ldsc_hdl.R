source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
# 任意两组性状：LDSC + HDL；全基因组汇总数据，不使用 MR clump SNP
# trait1 / trait2 可选择整组或单个性状

library(data.table)
library(R.utils)
library(GenomicSEM)
Sys.setenv(VROOM_THREADS = "1")
run_hdl_enabled <- Sys.getenv("FLT_RUN_HDL", "TRUE") != "FALSE"
methods <- if (run_hdl_enabled) c("ldsc", "hdl") else "ldsc"
if (run_hdl_enabled) library(HDL) else message("[SKIP] HDL reference unavailable; LDSC remains enabled.")

phase <- Sys.getenv("FLT_PHASE")
focus <- c("G732", "mdd_2025")
foods <- as.character(readxl::read_excel(Sys.getenv("FLT_DICTIONARY"), sheet = "item182")$gwas_id)
stopifnot(length(foods) == 182L, !anyDuplicated(foods), !anyNA(foods))
outcomes <- flt_panel()$id
root <- flt_run("ldsc")
cache <- file.path(Sys.getenv("FLT_TMP"), "ldsc")
logs <- file.path(Sys.getenv("FLT_PROJECT"), "logs", "ldsc")
qc <- file.path(Sys.getenv("FLT_QC"), "ldsc")
for (path in c(root, cache, logs, qc)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
ld_dir <- Sys.getenv("FLT_LDSC_REF", file.path(Sys.getenv("FLT_LDREF"), "ldsc/hm3"))
hm3 <- file.path(ld_dir, "w_hm3.snplist")
LD.path <- Sys.getenv("FLT_HDL_REF", file.path(Sys.getenv("FLT_LDREF"), "hdl/YGB_imputed_SVD_eigen99_extraction"))
sample.prev <- population.prev <- c(NA_real_, NA_real_)
required <- c(hm3, file.path(ld_dir, paste0(1:22, ".l2.ldscore.gz")), file.path(ld_dir, paste0(1:22, ".l2.M_5_50")))
if (any(!file.exists(required))) stop("Missing LDSC references: ", paste(required[!file.exists(required)], collapse = "; "))
write_input <- function(path, writer) {
  staging <- tempfile("partial-", tmpdir = dirname(path))
  on.exit(unlink(staging), add = TRUE)
  writer(staging)
  if (!file.rename(staging, path)) stop("Cannot publish ", path)
  path
}
munge_data <- function(d, trait, path, ldsc_ref) {
  d <- d[SNP %in% ldsc_ref$SNP]  # 先缩小到 HM3，再转换和合并，减少内存开销
  for (v in c("BETA", "SE", "P", "N")) set(d, j = v, value = as.numeric(d[[v]]))
  d[, `:=`(EA = toupper(trimws(EA)), NEA = toupper(trimws(NEA)))]
  d <- d[is.finite(BETA) & is.finite(P) & P >= 0 & P <= 1 & is.finite(N) & N > 0]
  if (!nrow(d)) stop("No HM3 SNPs with usable N, beta and P for ", trait)
  d <- d[N >= (as.numeric(stats::quantile(N, .9)) / 1.5)]
  d <- merge(d, ldsc_ref, by = "SNP")
  d <- d[EA %in% c("A", "C", "G", "T") & NEA %in% c("A", "C", "G", "T") & EA != NEA &
         !paste0(EA, NEA) %in% c("AT", "TA", "CG", "GC")]
  complement <- c(A = "T", T = "A", C = "G", G = "C")
  flip <- which(!((d$EA == d$R1 & d$NEA == d$R2) | (d$EA == d$R2 & d$NEA == d$R1)))
  d[flip, `:=`(EA = unname(complement[EA]), NEA = unname(complement[NEA]))]
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

# 3. 输入 SNP / A1 / A2 / N / Z；N 使用文件中的逐 SNP 样本量 --------------
# EA 是 BETA 的效应等位基因；BETA 为回归系数/log(OR)。Z 显式用 BETA/SE。
read_hdl <- function(d, file, hdl_ref) {
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

run_ldsc <- function(file1, file2, trait1, trait2, prefix, ldsc_fit) {
old_wd <- setwd(logs)
on.exit(setwd(old_wd), add = TRUE)
LDSCoutput <- ldsc_fit(traits = c(file1, file2), sample.prev = sample.prev,
  population.prev = population.prev, ld = ld_dir, wld = ld_dir,
  trait.names = c(trait1, trait2), stand = FALSE, ldsc.log = prefix)
S <- LDSCoutput$S; V <- LDSCoutput$V; I <- LDSCoutput$I
# V 的顺序：h2_exp、rgcov、h2_out；下面的 rg_SE 与原 V_Stand 缩放等价。
scale <- if (all(is.finite(diag(S)) & diag(S) > 0)) sqrt(S[1, 1] * S[2, 2]) else NA_real_
rg <- S[1, 2] / scale
gradient <- c(-rg / (2 * S[1, 1]), 1 / scale, -rg / (2 * S[2, 2]))
rg_se <- if (is.finite(rg)) sqrt(as.numeric(crossprod(gradient, V %*% gradient))) else NA_real_
rg_p <- if (is.finite(rg_se) && rg_se > 0) 2 * pnorm(abs(rg / rg_se), lower.tail = FALSE) else NA_real_
line <- grep("Cross trait Intercept:", readLines(paste0(prefix, "_ldsc.log"), warn = FALSE), value = TRUE)
int_se <- if (length(line) == 1L) suppressWarnings(as.numeric(sub(".*\\(([^)]+)\\).*", "\\1", line))) else NA_real_
result <- data.table(trait1 = trait1, trait2 = trait2,
  h2_exp = S[1, 1], h2_exp_se = sqrt(V[1, 1]), h2_out = S[2, 2], h2_out_se = sqrt(V[3, 3]),
  int_exp = I[1, 1], int_out = I[2, 2], gcov_int = I[1, 2], gcov_int_se = int_se,
  rgcov = S[1, 2], rgcov_se = sqrt(V[2, 2]), rg = rg, rg_se = rg_se, rg_p = rg_p,
  h2_z_flt = S[1, 1] / sqrt(V[1, 1]), h2_z_outcome = S[2, 2] / sqrt(V[3, 3]),
  rg_lower = rg - 1.96 * rg_se, rg_upper = rg + 1.96 * rg_se)

result
}
# pairs 按 trait1 排序；仅缓存当前 trait1，避免每个 trait2 配对重复读盘。
# 缓存限定为一个数据框，不将所有 GWAS 常驻内存。
hdl_cache <- new.env(parent = emptyenv())
hdl_cache$file <- NULL
run_hdl <- function(file1, file2, trait1, trait2, prefix) {
old_wd <- setwd(logs)
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

file_for <- function(group, id) file.path(Sys.getenv("FLT_GWAS"), group, "common", id, "gwas", paste0(id, ".gz"))
cached <- function(group, id, method = "ldsc") {
  ext <- if (method == "ldsc") ".sumstats.gz" else ".rds"
  file.path(cache, paste0(paste(group, id, method, sep = "-"), ext))
}

prepare_inputs <- function() {
  ldsc_ref <- fread(hm3, select = c("SNP", "A1", "A2"))
  setnames(ldsc_ref, c("A1", "A2"), c("R1", "R2"))
  ldsc_ref[, `:=`(SNP = as.character(SNP), R1 = toupper(R1), R2 = toupper(R2))]
  stopifnot(!anyDuplicated(ldsc_ref$SNP))


  registry <- rbind(data.table(group = "flt", id = foods), data.table(group = "brain", id = outcomes))
  registry <- registry[flt_shard(nrow(registry))]
  if (run_hdl_enabled && any(registry$id %in% focus)) {
    bims <- list.files(LD.path, pattern = "[.]bim$", full.names = TRUE)
    hdl_ref <- rbindlist(lapply(bims, function(f) {
      d <- fread(f, header = FALSE, select = c(2, 5, 6))
      setnames(d, c("SNP", "R1", "R2"))
      d
    }))
    hdl_ref[, `:=`(SNP = as.character(SNP), R1 = toupper(R1), R2 = toupper(R2))]
    stopifnot(nrow(hdl_ref) > 0L, !anyDuplicated(hdl_ref$SNP))
  }
  audit <- lapply(seq_len(nrow(registry)), function(i) {
    g <- registry$group[i]; id <- registry$id[i]
    f <- file_for(g, id)
    tryCatch({
      raw <- fread(f, select = c("SNP", "EA", "NEA", "BETA", "SE", "P", "N"))
      raw[, SNP := trimws(as.character(SNP))]
      path <- cached(g, id)
      munge_data(raw, id, path, ldsc_ref)
      hdl_detail <- ""
      if (run_hdl_enabled && id %in% focus) {
        hdl_file <- cached(g, id, "hdl")
        hdl_detail <- tryCatch({
          write_input(hdl_file, function(dest) saveRDS(read_hdl(raw, id, hdl_ref), dest, compress = FALSE))
          ""
        }, error = function(e) {
          unlink(hdl_file)
          conditionMessage(e)
        })
      }
      data.table(group = g, id = id, status = "ok", detail = hdl_detail, input = f)
    }, error = function(e) {
      unlink(c(cached(g, id), cached(g, id, "hdl")))
      data.table(group = g, id = id, status = "unavailable", detail = conditionMessage(e), input = f)
    })
  })
  fwrite(rbindlist(audit), file.path(qc, paste0("inputs-", Sys.getenv("FLT_SHARD"), ".csv")))
}
fit_pairs <- function() {
  y <- Sys.getenv("FLT_OUTCOME")
  # Reuse immutable LD-score and target-GWAS tables across the pair fits in
  # this job. Keep the upstream GenomicSEM namespace and regression unchanged.
  ldsc_fit <- GenomicSEM::ldsc
  readers <- new.env(parent = environment(ldsc_fit))
  cache_inputs <- function(reader, paths, ...) {
    saved <- setNames(lapply(paths, function(path) suppressMessages(reader(path, ...))), paths)
    function(file, ...) {
      if (length(file) == 1L && is.character(file) && !is.null(saved[[file]])) saved[[file]]
      else reader(file, ...)
    }
  }
  keys <- file.path(ld_dir, paste0(1:22, ".l2.ldscore.gz"))
  target <- cached("brain", y)
  if (file.exists(target)) keys <- c(keys, target)
  readers$read_delim <- cache_inputs(readr::read_delim, keys,
      delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE, num_threads = 1L)
  readers$read_csv <- cache_inputs(readr::read_csv,
      file.path(ld_dir, paste0(1:22, ".l2.M_5_50")), col_names = FALSE, num_threads = 1L)
  environment(ldsc_fit) <- readers
  one <- function(x) {
    result <- tryCatch(run_ldsc(cached("flt", x), cached("brain", y), x, y,
                      paste0(x, "-", y), ldsc_fit), error = function(e) data.table(
                        trait1 = x, trait2 = y, rg = NA_real_, rg_p = NA_real_,
                        detail = conditionMessage(e)))
    result[, method := "ldsc"]
    result[, status := ifelse(is.finite(rg_p), "ok", "unavailable")]
    if (run_hdl_enabled && identical(c(x, y), focus)) {
      h <- tryCatch(run_hdl(cached("flt", x, "hdl"), cached("brain", y, "hdl"), x, y,
                    paste0(x, "-", y)), error = function(e) data.table(
                      trait1 = x, trait2 = y, rg = NA_real_, rg_p = NA_real_, detail = conditionMessage(e)))
      h[, `:=`(method = "hdl", status = ifelse(is.finite(rg_p), "ok", "unavailable"))]
      result <- rbindlist(list(result, h), fill = TRUE)
    }
    result
  }
  data.table::setDTthreads(1L)
  cpus <- if (.Platform$OS.type == "unix") flt_cores() else 1L
  fits <- parallel::mclapply(foods, one, mc.cores = cpus, mc.preschedule = TRUE)
  if (any(vapply(fits, inherits, logical(1), "try-error"))) stop("An LDSC worker failed")
  fwrite(rbindlist(fits, fill = TRUE), file.path(root, paste0(y, "-pairs.csv")))
}
collect_pairs <- function() {
  results <- rbindlist(lapply(outcomes, function(y) fread(file.path(root, paste0(y, "-pairs.csv")))), fill = TRUE)
  stopifnot(!anyDuplicated(results[, .(trait1, trait2, method)]))
  results <- merge(results, as.data.table(flt_panel()), by.x = "trait2", by.y = "id", sort = FALSE)
  results[, FDR_outcome := p.adjust(rg_p, "BH", n = .N), by = .(trait2, method)]
  results[, FDR := p.adjust(rg_p, "BH", n = .N), by = .(panel, method)]
  for (v in c("h2_z_flt", "h2_z_outcome")) if (!v %in% names(results)) results[, (v) := NA_real_]
  results[, h2_precision := fifelse(method != "ldsc", NA_character_,
                         fifelse(h2_z_flt >= 4 & h2_z_outcome >= 4, "both_h2_Z_ge_4", "low_h2_precision"))]
  fwrite(results, file.path(root, "ldsc_hdl_results.csv"))
  if (!any(results$status == "ok")) stop("All LDSC/HDL pairs unavailable; see result details and logs.")
  fwrite(results[is.finite(FDR) & FDR < .05], file.path(root, "ldsc_hdl_hits.csv"))
}
switch(phase, "ldsc-input" = prepare_inputs(), "ldsc-pairs" = fit_pairs(), "ldsc" = collect_pairs(),
       stop("Run flt.sh run ldsc"))
