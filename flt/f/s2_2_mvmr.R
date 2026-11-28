source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
# 多变量 MR（MVMR）：FLT -> brain，调整 income、education、physical_activity
# 结局沿用 s2_1_uvmr.R 的二分类设定，BETA 应为 log(OR)。
# Windows Rscript 可能继承无效的 locale；启用 UTF-8 以支持中文项目路径。
if (.Platform$OS.type == "windows") {
  if (!nzchar(Sys.setlocale("LC_CTYPE", ".UTF-8"))) stop("无法启用 UTF-8 locale")
}

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(MVMR)
  library(openxlsx)
  library(robustbase)
})

# 运行设置 ------------------------------------
flt_dir <- file.path(Sys.getenv("FLT_GWAS"), "flt/common")
brain_dir <- file.path(Sys.getenv("FLT_GWAS"), "brain/common")
covs_dir <- file.path(Sys.getenv("FLT_GWAS"), "covs/common")
rootdir <- Sys.getenv("FLT_OUT")
instrument_cache_dir <- file.path(Sys.getenv("FLT_TMP"), "mvmr-instruments")
reuse_clump <- TRUE
plink <- Sys.getenv("PLINK_BIN", Sys.getenv("FLT_PLINK"))
ld_prefix <- Sys.getenv("MR_LD_PREFIX", file.path(Sys.getenv("FLT_LDREF"), "opengwas/1kg.v3/EUR"))
p_cut <- 5e-8
clump_kb <- 10000
clump_r2 <- 0.001

outdir <- flt_run("mvmr")
mvdir <- outdir
clump_dir <- file.path(Sys.getenv("FLT_TMP"), "mvmr-clump", Sys.getenv("FLT_PANEL"))

covariates <- c("income", "edu", "physical_activity")

# 仅分析以下；每个模型同时纳入三个协变量。

# 条件F/Q默认gencov=0；共享样本时需谨慎解释。误差相关矩阵可选。
error_cor_dir <- Sys.getenv("FLT_MVMR_ERROR_COR_DIR")

# 单个配对的主流程 ------------------------------------
twosample_mvmr_function <- function(x, y) {
  # 1. 主流程已读取 GWAS，assoc 保存第一轮工具在各性状中的效应。
  ids <- c(x, covariates)
  key <- paste(x, y, sep = "_")

  # 2. 合并本配对四个暴露各自 clump 后的工具，再联合 clump。
  pair_tools <- instruments[ids]
  first <- rbindlist(lapply(pair_tools, function(tool) tool$data))
  # 联合 clump 按最小 P 去除重复 rsID。候选 SNP 已在第一轮通过 p_cut。
  keep <- clump(first, paste0("mv_", key), p1 = 0.99, log_dir = clump_dir)
  counts <- data.table(
    food = x, outcome = y, exposure = ids,
    n_significant = vapply(pair_tools, function(tool) tool$n_significant, integer(1)),
    n_first_clump = vapply(pair_tools, function(tool) length(tool$snps), integer(1)),
    n_union = length(unique(first$SNP)),
    n_joint_clump = length(keep), n_matched = NA_integer_,
    n_allele_aligned = NA_integer_, n_harmonised = NA_integer_
  )
  mv_instrument_counts[[key]] <<- counts

  # 3. 匹配所有暴露和结局的效应，统一效应等位基因并构造矩阵。
  aligned <- harmonise_mvmr(x, y, ids, keep, counts)
  bx <- aligned$bx
  sx <- aligned$sx
  mvdat <- aligned$mvdat

  # 4. 拟合 MV-IVW 和敏感性模型。食物为第一列，其余列为协变量。
  input <- MendelianRandomization::mr_mvinput(
    bx = bx, bxse = sx, by = mvdat$outcome_beta, byse = mvdat$outcome_se,
    exposure = ids, outcome = y, snps = rownames(bx)
  )
  models <- fit_mvmr(input, x, y, ids, bx, mvdat)

  # 5. 检查条件工具强度和异质性，保存每个配对的输入与模型。
  conditional <- conditional_diagnostics(x, y, ids, bx, sx, mvdat)
  if (!is.null(conditional$F)) mv_conditional_f[[key]] <<- conditional$F
  if (!is.null(conditional$Q)) mv_pleiotropy[[key]] <<- conditional$Q
  saveRDS(
    list(input = input, fit = models$fits[["MV-IVW"]], fits = models$fits,
      harmonised = aligned$harmonised, mvdat = mvdat,
      exposure_order = ids, conditional = conditional),
    file.path(mvdir, "harmonised", paste0("mvmr_", key, ".rds"))
  )
  rbindlist(models$results, fill = TRUE)
}

# 1. GWAS 读取与质控 ------------------------------------
read_gwas <- function(file) {
  cols <- c("SNP", "CHR", "EA", "NEA", "BETA", "SE", "P", "EAF")
  # 仅解压首行获取列名，避免 fread(nrows=0) 再次解压整个 gzip 文件。
  con <- gzfile(file, open = "rt")
  header <- tryCatch(readLines(con, n = 1L, warn = FALSE), finally = close(con))
  if (!length(header)) stop("GWAS 文件为空：", file)
  header_cols <- names(fread(text = paste0(header, "\n"), nrows = 0, showProgress = FALSE))
  missing_cols <- setdiff(cols, header_cols)
  if (length(missing_cols)) stop("GWAS 缺列：", paste(missing_cols, collapse = ", "), "；", file)
  cols <- c(cols, intersect("N", header_cols))
  d <- fread(file, select = cols, showProgress = FALSE)
  if (!"N" %in% names(d)) d[, N := NA_real_]
  for (v in c("BETA", "SE", "P", "EAF", "N")) set(d, j = v, value = as.numeric(d[[v]]))
  d[, `:=`(SNP = as.character(SNP), CHR = as.integer(CHR), EA = toupper(EA), NEA = toupper(NEA))]
  d <- d[grepl("^rs[0-9]+$", SNP) & CHR %in% 1:22 & EA %in% c("A", "C", "G", "T") &
    NEA %in% c("A", "C", "G", "T") & EA != NEA & is.finite(BETA) & is.finite(SE) & SE > 0 &
    is.finite(P) & P >= 0 & P <= 1]
  d <- d[!duplicated(SNP) & !duplicated(SNP, fromLast = TRUE)]
  d[!is.finite(EAF) | EAF <= 0 | EAF >= 1, EAF := NA_real_]
  d[, P := pmax(P, 1e-300)]
  d
}

# 2. 各暴露及联合 LD clumping ------------------------------------
# 两份脚本共用第一轮工具缓存。按输入文件状态和 clumping 参数判定是否可复用。
instrument_stamp <- function(file) {
  inputs <- normalizePath(c(file, plink, ref_files), winslash = "/", mustWork = TRUE)
  stamp <- list(version = 2L, files = file.info(inputs)[, c("size", "mtime", "ctime")],
    p_cut = p_cut, clump_kb = clump_kb, clump_r2 = clump_r2)
  paste0("# inputs: ", paste(deparse(stamp, width.cutoff = 500L, control = "all"), collapse = " "))
}
get_instruments <- function(file, id) {
  cache_file <- file.path(instrument_cache_dir, paste0(id, "_clump.txt"))
  stamp <- instrument_stamp(file)
  if (reuse_clump && file.exists(cache_file)) {
    cached <- tryCatch({
      meta <- readLines(cache_file, n = 3L)
      if (!identical(meta[1L], stamp)) stop("缓存已过期")
      n_significant <- as.integer(sub("^# n_significant=", "", meta[2L]))
      n_clumped <- as.integer(sub("^# n_clumped=", "", meta[3L]))
      if (anyNA(c(n_significant, n_clumped)) || n_clumped < 0L || n_significant < n_clumped) {stop("缓存 SNP 计数无效")}
      d <- fread(cache_file, skip = 3L, sep = "\t", header = TRUE, na.strings = "NA",
        colClasses = list(character = c("SNP", "EA", "NEA"), integer = "CHR",
          numeric = c("BETA", "SE", "P", "EAF", "N")), showProgress = FALSE)
      if (!setequal(names(d), c("SNP", "CHR", "EA", "NEA", "BETA", "SE", "P", "EAF", "N")) ||
          nrow(d) != n_clumped || anyNA(d$SNP) || anyDuplicated(d$SNP)) stop("缓存 SNP 表无效")
      list(n_significant = n_significant, data = d)
    }, error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(cached)) {message("读取工具 SNP TXT 缓存：", id); cached$snps <- cached$data$SNP; return(cached)}
  }
  message("筛选并 clump 工具 SNP：", id)
  significant <- read_gwas(file)[P < p_cut]
  snps <- clump(significant, id)
  result <- list(n_significant = nrow(significant), data = significant[SNP %in% snps])
  if (!identical(stamp, instrument_stamp(file))) stop("筛选工具时 GWAS、PLINK 或 LD 输入发生变化，请重试：", id)
  # 同一个 TXT 保存三行注释及 clump 后的 SNP 表。完整写出后再发布。
  tmp_file <- tempfile(".instruments-", tmpdir = instrument_cache_dir, fileext = ".tmp")
  on.exit(unlink(tmp_file), add = TRUE)
  writeLines(c(stamp, paste0("# n_significant=", result$n_significant), paste0("# n_clumped=", nrow(result$data))), tmp_file)
  fwrite(result$data, tmp_file, sep = "\t", na = "NA", append = TRUE, col.names = TRUE)
  if (!file.rename(tmp_file, cache_file)) stop("无法写入 clump SNP TXT：", cache_file)
  result$snps <- result$data$SNP
  result
}

clump <- function(d, tag, p1 = p_cut, log_dir = clump_dir) {
  if (!nrow(d)) return(character())
  d <- d[order(P)][!duplicated(SNP)]
  prefix <- file.path(log_dir, paste0(tag, "_all"))
  input <- paste0(prefix, ".tsv")
  output <- paste0(prefix, ".clumped")
  fwrite(d[, .(SNP, P)], input, sep = "\t")
  unlink(output)
  status <- system2(plink, c(
    "--bfile", shQuote(ld_prefix), "--clump", shQuote(input),
    "--clump-p1", p1, "--clump-p2", 1, "--clump-r2", clump_r2,
    "--clump-kb", clump_kb, "--threads", "1", "--out", shQuote(prefix)
  ), stdout = paste0(prefix, ".stdout.txt"), stderr = paste0(prefix, ".stderr.txt"))
  if (status != 0) stop("PLINK 失败，请检查：", prefix, ".log")
  if (!file.exists(output)) return(character())

  hits <- fread(output, showProgress = FALSE)
  if (!nrow(hits)) return(character())
  # 保持染色体顺序，确保下游随机方法的输入顺序一致。
  unique(hits[order(CHR), SNP])
}

# 3. 匹配和协调等位基因，构建效应矩阵 ------------------------------------
format_mr_data <- function(id, type, keep) {
  d <- as.data.frame(assoc[[id]][SNP %in% keep])
  if (!nrow(d)) stop("没有可匹配 SNP：", id)
  d$trait <- d$id <- id
  sample_n <- trait_info$N[match(id, trait_info$id)]
  d$N[!is.finite(d$N) | d$N <= 2] <- sample_n
  TwoSampleMR::format_data(d, type = type, phenotype_col = "trait", id_col = "id", snp_col = "SNP",
    beta_col = "BETA", se_col = "SE", effect_allele_col = "EA", other_allele_col = "NEA",
    eaf_col = "EAF", pval_col = "P", samplesize_col = "N", min_pval = 1e-300)
}

harmonise_mvmr <- function(x, y, ids, keep, counts) {
  # 在所有暴露和结局中取共同 SNP；只在任一暴露中显著即可，不要求四者均显著。
  keep <- Reduce(intersect, c(list(keep), lapply(assoc[c(ids, y)], function(d) d$SNP)))
  counts[, n_matched := length(keep)]
  if (!length(keep)) stop("MVMR 联合工具在所有暴露和结局中无交集")
  anchor <- format_mr_data(x, "exposure", keep)
  hh <- lapply(c(covariates, y), function(id) {
    h <- TwoSampleMR::harmonise_data(anchor, format_mr_data(id, "outcome", keep), action = 2)
    h[h$mr_keep %in% TRUE, ]
  })
  shared <- Reduce(intersect, lapply(hh, function(h) h$SNP))
  counts[, n_allele_aligned := length(shared)]
  if (length(shared) <= length(ids)) stop("Insufficient complete instruments: require more than ", length(ids))
  hh <- lapply(hh, function(h) h[match(shared, h$SNP), ])
  # 把已对齐的 covariate 数据转为暴露长表，交回参考代码的 mv_harmonise_data。
  dat_exp_mv <- list(as.data.table(anchor[match(shared, anchor$SNP), ]))
  cols <- c("SNP", paste0(
    c("beta", "se", "effect_allele", "other_allele", "eaf", "pval", "samplesize"), ".outcome"
  ))
  for (j in seq_along(covariates)) {
    d <- as.data.table(hh[[j]])[, ..cols]
    setnames(d, sub("\\.outcome$", ".exposure", names(d)))
    d[, `:=`(id.exposure = covariates[j], exposure = covariates[j])]
    dat_exp_mv[[j + 1L]] <- d
  }
  mvdat <- TwoSampleMR::mv_harmonise_data(
    exposure_dat = as.data.frame(rbindlist(dat_exp_mv, fill = TRUE)),
    outcome_dat = format_mr_data(y, "outcome", shared), harmonise_strictness = 2)
  # 列可按 id 排序；明确食物在第一列，保证 orientate=1 的含义不变。
  order_x <- match(ids, colnames(mvdat$exposure_beta))
  if (anyNA(order_x)) stop("MVMR 暴露列缺失")
  bx <- mvdat$exposure_beta[, order_x, drop = FALSE]
  sx <- mvdat$exposure_se[, order_x, drop = FALSE]
  counts[, n_harmonised := nrow(bx)]
  if (!all(is.finite(c(bx, sx, mvdat$outcome_beta, mvdat$outcome_se))) ||
      any(sx <= 0) || any(mvdat$outcome_se <= 0)) stop("MVMR 效应/SE 缺失或无效")
  if (nrow(bx) <= length(ids) || qr(bx)$rank < length(ids)) stop("MVMR 工具不足或矩阵秩不足")
  list(bx = bx, sx = sx, mvdat = mvdat, harmonised = hh)
}

# 4. MV-IVW 和敏感性模型 ------------------------------------
fit_mvmr <- function(input, x, y, ids, bx, mvdat) {
  # Egger/Lasso 均以第一个暴露（食物喜好）定向；Lasso 使用默认异质性停止规则。
  fitters <- list(
    "MV-IVW" = function() MendelianRandomization::mr_mvivw(
      input, model = "default", correl = FALSE, distribution = "normal", alpha = 0.05
    ),
    "MVMR-median" = function() MendelianRandomization::mr_mvmedian(input,
      distribution = "normal", alpha = 0.05, iterations = 10000, seed = 314159265),
    "MVMR-Egger" = function() {
      if (nrow(bx) <= length(ids) + 1 ||
          qr(cbind(1, bx * ifelse(bx[, 1] < 0, -1, 1)))$rank < length(ids) + 1)
        stop("MVMR-Egger 需 >5 个 SNP，且含截距的暴露矩阵满秩")
      MendelianRandomization::mr_mvegger(
        input, orientate = 1, correl = FALSE, distribution = "normal", alpha = 0.05
      )
    },
    "MVMR-Lasso" = function() MendelianRandomization::mr_mvlasso(input,
      orientate = 1, distribution = "normal", alpha = 0.05, lambda = numeric(0))
  )
  fits <- results <- list()
  for (method in names(fitters)) {
    a <- attempt(paste0("mvmr_", method), x, y, function() {
      fit <- fitters[[method]]()
      egger <- method == "MVMR-Egger"
      se <- if (egger) fit@StdError.Est else fit@StdError
      p <- if (egger) fit@Pvalue.Est else fit@Pvalue
      lo <- if (egger) fit@CILower.Est else fit@CILower
      hi <- if (egger) fit@CIUpper.Est else fit@CIUpper
      if (!all(is.finite(c(fit@Estimate, se, p, lo, hi))) || any(se <= 0))
        stop(method, " 未返回有效估计")
      fits[[method]] <<- fit
      if (method %in% c("MV-IVW", "MVMR-Egger")) {
        mv_diagnostics[[paste(x, y, method, sep = "_")]] <<- data.table(
          food = x, outcome = y, method = method, nsnp = fit@SNPs,
          Q = fit@Heter.Stat[1], Q_df = fit@SNPs - length(ids) - as.integer(egger),
          Q_pval = fit@Heter.Stat[2], intercept = if (egger) fit@Intercept else NA_real_,
          intercept_se = if (egger) fit@StdError.Int else NA_real_,
          intercept_pval = if (egger) fit@Pvalue.Int else NA_real_)
      }
      data.table(food = x, outcome = y, exposure = ids, method = method,
        nsnp = if (method == "MVMR-Lasso") fit@Valid else fit@SNPs,
        nsnp_input = nrow(bx), b = as.numeric(fit@Estimate), se = as.numeric(se),
        pval = as.numeric(p), OR = exp(fit@Estimate), OR_low = exp(lo), OR_high = exp(hi),
        significant = p < 0.05)
    })
    results[[method]] <- a
  }
  # 附件 Robust 为额外敏感性结果，不替代四种核心方法或改变既有稳健性筛选规则。
  results[["MVMR-Robust"]] <- attempt("mvmr_Robust", x, y, function() {
    set.seed(314159265)
    fit <- mvmr_robust(bx, mvdat$outcome_beta, mvdat$outcome_se)
    b <- as.numeric(fit$coefficients)
    se <- as.numeric(fit$se)
    if (length(b) != length(ids) || !all(is.finite(c(b, se))) || any(se <= 0))
      stop("MVMR-Robust 未返回有效估计")
    fits[["MVMR-Robust"]] <<- fit
    p <- 2 * pnorm(abs(b / se), lower.tail = FALSE)
    data.table(food = x, outcome = y, exposure = ids, method = "MVMR-Robust",
      nsnp = nrow(bx), nsnp_input = nrow(bx), b = b, se = se, pval = p,
      OR = exp(b), OR_low = exp(b - 1.96 * se), OR_high = exp(b + 1.96 * se), significant = p < 0.05)
  })
  list(fits = fits, results = results)
}

mvmr_robust <- function(bx, by, seby, k.max = 1000, maxit.scale = 1000) {
  fit <- robustbase::lmrob(by ~ bx - 1, weights = seby^-2,
    k.max = k.max, maxit.scale = maxit.scale)
  s <- summary(fit)
  list(coefficients = s$coef[, 1], se = s$coef[, 2] / min(s$sigma, 1), fit = fit)
}

# 5. 条件 F 与条件 Q ------------------------------------
conditional_diagnostics <- function(x, y, ids, bx, sx, mvdat) {
  # 四个暴露各输出一个条件 F；任一诊断失败都不会丢失已完成的模型结果。
  conditional <- attempt("mvmr_conditional_input", x, y, function() {
    gencov <- 0
    assumption <- "gencov=0_assumed"
    error_cor_file <- file.path(error_cor_dir, paste0(x, ".csv"))
    if (nzchar(error_cor_dir) && file.exists(error_cor_file)) {
      r <- as.matrix(read.csv(error_cor_file, row.names = 1L, check.names = FALSE))
      r <- r[ids, ids, drop = FALSE]
      if (!all(is.finite(r)) || !isSymmetric(r) || any(abs(diag(r) - 1) > 1e-8) ||
          min(eigen(r, symmetric = TRUE)$values) < -1e-8) stop("无效的 MVMR 误差相关矩阵：", x)
      gencov <- MVMR::phenocov_mvmr(r, sx)
      assumption <- "provided_error_correlation"
    }
    fdata <- MVMR::format_mvmr(BXGs = bx, BYG = mvdat$outcome_beta,
      seBXGs = sx, seBYG = mvdat$outcome_se, RSID = rownames(bx))
    list(input = fdata, gencov = gencov, assumption = assumption)
  })
  if (!is.null(conditional)) {
    f_table <- attempt("mvmr_conditional_F", x, y, function() {
      f <- MVMR::strength_mvmr(r_input = conditional$input, gencov = conditional$gencov)
      values <- as.numeric(as.matrix(f))
      if (length(values) != length(ids) || any(!is.finite(values))) stop("条件 F 返回值无效")
      data.table(food = x, outcome = y, exposure = ids, nsnp = nrow(bx),
        conditional_F = values, F_ge_10 = values >= 10,
        covariance_assumption = conditional$assumption, MVMR_version = as.character(utils::packageVersion("MVMR")))
    })
    q_table <- attempt("mvmr_conditional_Q", x, y, function() {
      q <- MVMR::pleiotropy_mvmr(r_input = conditional$input, gencov = conditional$gencov)
      data.table(food = x, outcome = y, nsnp = nrow(bx), Q = as.numeric(q$Qstat),
        Q_pval = as.numeric(q$Qpval), covariance_assumption = conditional$assumption,
        all_conditional_F_ge_10 = if (is.null(f_table)) NA else all(f_table$F_ge_10))
    })
    conditional$F <- f_table
    conditional$Q <- q_table  # 条件工具较弱时，Q 检验不宜解释为多效性证据。
  }
  conditional
}

# 稳健性判定、错误记录与输出 ------------------------------------
flag_mr_robustness <- function(tab, keys, main_methods, sensitivity_methods,
                               main_p_col, main_cutoff = 0.05, sensitivity_cutoff = 0.05) {
  d <- as.data.frame(tab)
  required <- c(keys, "method", "b", "se", "pval", main_p_col)
  if (!all(required %in% names(d)) || anyNA(d[keys])) stop("稳健性筛选字段缺失或配对标识为空")
  if (length(intersect(main_methods, sensitivity_methods))) stop("主方法不能作为敏感性支持")
  defaults <- list(primary_significant = NA, sensitivity_n_available = 0L,
    sensitivity_n_supporting = NA_integer_, supporting_methods = "",
    sensitivity_support = NA, robust = NA, robust_status = "primary_unavailable")
  d[names(defaults)] <- lapply(defaults, rep, times = nrow(d))
  valid <- function(z) is.finite(z$b) & is.finite(z$se) & z$se > 0 &
    is.finite(z$pval) & z$pval >= 0 & z$pval <= 1
  groups <- do.call(paste, c(d[keys], list(sep = "\r")))
  for (ii in split(seq_len(nrow(d)), groups)) {
    z <- d[ii, , drop = FALSE]
    main <- z[z$method %in% main_methods, , drop = FALSE]
    sens <- z[z$method %in% sensitivity_methods, , drop = FALSE]
    sens <- sens[valid(sens), , drop = FALSE]
    d$sensitivity_n_available[ii] <- length(unique(sens$method))
    if (nrow(main) > 1) stop("同一配对存在多个主方法结果")
    if (nrow(main) != 1 || !valid(main)) next
    p <- main[[main_p_col]][1]
    if (!is.finite(p) || p < 0 || p > 1) next
    d$primary_significant[ii] <- p < main_cutoff
    if (nrow(sens)) {
      support <- sort(unique(sens$method[sign(sens$b) == sign(main$b) &
        main$b != 0 & sens$pval < sensitivity_cutoff]))
      d$sensitivity_n_supporting[ii] <- length(support)
      d$supporting_methods[ii] <- paste(support, collapse = "; ")
      d$sensitivity_support[ii] <- length(support) > 0
    }
    if (p >= main_cutoff) {
      d$robust[ii] <- FALSE
      d$robust_status[ii] <- "primary_not_significant"
    } else if (!nrow(sens)) {
      d$robust_status[ii] <- if (main$method == "Wald ratio") "single_snp_not_assessable" else "sensitivity_unavailable"
    } else {
      d$robust[ii] <- d$sensitivity_support[ii]
      d$robust_status[ii] <- if (length(support)) "supported" else "no_concordant_significant_sensitivity"
    }
  }
  d
}

attempt <- function(stage, x, y, fun) {
  tryCatch(fun(), error = function(e) {
    errors <<- rbind(errors, data.table(stage = stage, exposure = x, outcome = y, message = conditionMessage(e)))
    fwrite(errors, file.path(outdir, "failed_models.csv"))
    message(stage, " ", x, " -> ", y, "：", conditionMessage(e))
    NULL
  })
}

label_results <- function(x) {
  x <- copy(as.data.table(x))
  if (!ncol(x)) return(data.table(status = character()))
  for (v in intersect(c("exposure", "outcome", "food"), names(x))) {
    set(x, j = paste0(v, "_name"), value = unname(trait_names[x[[v]]]))
  }
  x
}

save_result <- function(x, name, folder = outdir) {
  x <- label_results(x)
  if (!ncol(x)) x <- data.table(status = character())
  fwrite(x, file.path(folder, paste0(name, ".csv")))
  saveRDS(x, file.path(folder, paste0(name, ".rds")))
}

write_xlsx <- function(tables, filename, folder) {
  tables <- lapply(tables, function(x) as.data.frame(label_results(x)))
  openxlsx::write.xlsx(tables, file.path(folder, filename), overwrite = TRUE)
}

# 执行与输出 ------------------------------------
ref_files <- paste0(ld_prefix, c(".bed", ".bim", ".fam"))
if (!file.exists(plink)) stop("找不到 PLINK：", plink)
if (any(!file.exists(ref_files))) stop("缺少 LD 参考：", paste(ref_files[!file.exists(ref_files)], collapse = "; "))
for (folder in c(outdir, clump_dir, file.path(outdir, "harmonised")))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
dir.create(instrument_cache_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(12345)

# 1. 读取配对与 GWAS 清单。
u <- flt_uv("forward")
hits <- u[u$Method == "IVW" & is.finite(u[["P-FDR"]]) & u[["P-FDR"]] < 0.05, ]
pairs <- data.table(exposure=hits$Exposure, outcome=names(flt_labels)[match(hits$Outcome, flt_labels)])
if (!all(c("exposure", "outcome") %in% names(pairs))) stop("配对表需要exposure和outcome两列")
pairs <- unique(rbind(
  pairs[, .(exposure = as.character(exposure), outcome = as.character(outcome))],
  data.table(exposure = "G732", outcome = "mdd_2025")))
pairs[, selection := ifelse(exposure == "G732" & outcome == "mdd_2025", "focal", "forward_ivw_fdr")]
if (anyNA(pairs)) stop("Unrecognised MVMR outcome")
if (!nrow(pairs)) {
  fwrite(data.table(exposure=character(), food=character(), outcome=character(), method=character(),
    OR=numeric(), OR_low=numeric(), OR_high=numeric(), pval=numeric()), file.path(outdir,"mvmr_all_results.csv"))
  quit(status=0)
}

# 只读取指定配对涉及的食物、疾病及三个协变量。
food <- data.table(id = unique(pairs$exposure))
food[, `:=`(name = id, file = file.path(flt_dir, id, "gwas", paste0(id, ".gz")))]
other <- data.table(id = unique(c(covariates, pairs$outcome)))
other[, `:=`(name = id, file = file.path(
  ifelse(id %in% setdiff(covariates, "edu"), covs_dir, brain_dir), id, "gwas", paste0(id, ".gz")))]
gwas <- rbindlist(list(food, other))
fwrite(pairs, file.path(outdir, "mvmr_pairs.csv"))
trait_names <- setNames(gwas$name, gwas$id)
fwrite(gwas, file.path(outdir, "gwas_manifest.tsv"), sep = "\t")

# N仅用于格式化的可选补充；没有此表也可运行MVMR。
info_file <- file.path(flt_run("uvmr"), "steiger_trait_info_used.csv")
trait_info <- if (file.exists(info_file)) fread(info_file, select = c("id", "N")) else
  data.table(id = gwas$id, N = NA_real_)

if (anyDuplicated(gwas$id)) stop("GWAS 清单的 id 重复，请检查配对和协变量")
if (any(!file.exists(gwas$file))) stop("缺少 GWAS：", paste(gwas$file[!file.exists(gwas$file)], collapse = "；"))

# 2. 食物和三个协变量各自读取第一轮缓存，缺失时才筛选并 clump。
exposure_ids <- unique(c(pairs$exposure, covariates))
instruments <- setNames(lapply(exposure_ids, function(id) {
  get_instruments(gwas$file[match(id, gwas$id)], id)
}), exposure_ids)
instrument_counts <- data.table(
  id = exposure_ids,
  n_significant = vapply(instruments, function(tool) tool$n_significant, integer(1)),
  n_clumped = vapply(instruments, function(tool) length(tool$snps), integer(1))
)
fwrite(instrument_counts, file.path(outdir, "instrument_counts.csv"))

# 3. 从各暴露和结局的 GWAS 提取第一轮工具的交叉效应，供每个配对匹配。
snps <- unique(unlist(lapply(instruments, function(tool) tool$snps), use.names = FALSE))
assoc <- list()
for (i in seq_len(nrow(gwas))) {
  message("提取交叉效应：", gwas$id[i])
  assoc[[gwas$id[i]]] <- read_gwas(gwas$file[i])[SNP %in% snps]
}

# 逐配对执行联合 clumping、匹配、拟合和诊断，再汇总输出。
errors <- data.table(stage = character(), exposure = character(), outcome = character(), message = character())
fwrite(errors, file.path(outdir, "failed_models.csv"))
mv_diagnostics <- mv_conditional_f <- mv_pleiotropy <- mv_instrument_counts <- list()

mv <- lapply(seq_len(nrow(pairs)), function(i) {
  x <- pairs$exposure[i]
  y <- pairs$outcome[i]
  attempt("mvmr", x, y, function() twosample_mvmr_function(x, y))
})
mv_result <- rbindlist(mv, fill = TRUE)
if (!nrow(mv_result)) mv_result <- data.table(food = character(), outcome = character(), exposure = character(),
  method = character(), nsnp = integer(), nsnp_input = integer(), b = numeric(), se = numeric(), pval = numeric(),
  OR = numeric(), OR_low = numeric(), OR_high = numeric(), significant = logical())
# 逐个 exposure 系数判断稳健性，协变量不能为食物效应提供支持。
mv_result <- as.data.table(flag_mr_robustness(mv_result, c("food", "outcome", "exposure"),
  "MV-IVW", c("MVMR-median", "MVMR-Egger", "MVMR-Lasso"), "pval"))

mv_main <- mv_result[method == "MV-IVW" & exposure == food]

diagnostic_table <- rbindlist(mv_diagnostics, fill = TRUE)
conditional_f_table <- rbindlist(mv_conditional_f, fill = TRUE)
conditional_q_table <- rbindlist(mv_pleiotropy, fill = TRUE)
instrument_count_table <- rbindlist(mv_instrument_counts, fill = TRUE)

save_result(mv_result, "mvmr_all_results", mvdir)
save_result(mv_result[method == "MV-IVW"], "mvmr_results", mvdir)
save_result(mv_result[method != "MV-IVW"], "mvmr_sensitivity_results", mvdir)
save_result(mv_result[method == "MVMR-Robust"], "mvmr_robust_regression", mvdir)
save_result(mv_main, "mvmr_robustness", mvdir)
save_result(mv_main[robust %in% TRUE], "mvmr_robust_hits", mvdir)
save_result(diagnostic_table, "mvmr_diagnostics", mvdir)
save_result(conditional_f_table, "mvmr_conditional_F", mvdir)
save_result(conditional_q_table, "mvmr_conditional_Q", mvdir)
save_result(instrument_count_table, "mvmr_instrument_counts", mvdir)

notes <- data.frame(Definition = c(
  "Only the explicitly supplied pairs are analysed; UVMR and reverse MR are not executed.",
  "Instruments: clump each exposure, union retained SNPs, then jointly clump; three covariates retained.",
  "MV-IVW, median, Egger and Lasso match the original fitting settings. Robust is an additional method.",
  "Robust flag: MV-IVW P<0.05 and >=1 concordant median/Egger/Lasso P<0.05 for the same coefficient.",
  "Conditional F/Q use gencov=0 unless an error-correlation matrix is supplied; F_mean is not conditional F.",
  "MVMR P-values are unadjusted. Shared-sample F/Q under zero covariance are assumption-dependent.",
  "UVMR candidate pairs and MVMR robustness are separate outputs."
))

write_xlsx(list(Food_effects = mv_result[exposure == food], All_effects = mv_result,
  Diagnostics = diagnostic_table,
  Robust_regression = mv_result[method == "MVMR-Robust"],
  Conditional_F = conditional_f_table,
  Conditional_Q = conditional_q_table,
  Instrument_counts = instrument_count_table, Notes = notes), "mvmr_summary.xlsx", mvdir)

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

message("MVMR完成：", outdir)
