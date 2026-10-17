# Aggregate publication for PCA/CSx/Disco and the final -> Shiny data contract.
# Only aggregate tables and figures enter the dashboard; no individual IDs/models.
suppressPackageStartupMessages({library(data.table); library(openxlsx); library(ggplot2)})
setDTthreads(4)

result_stamp <- function(paths) {
  paths <- unique(paths[nzchar(paths) & !is.na(paths)])
  z <- file.info(paths)
  data.table(path = paths, size = z$size, modified = as.numeric(z$mtime))
}
same_stamp <- function(a,b) isTRUE(all.equal(as.data.frame(a),as.data.frame(b),tolerance=0,check.attributes=FALSE))
stable_read <- function(path, fun) {
  before <- result_stamp(path)
  value <- fun(path)
  if (!same_stamp(before, result_stamp(path))) stop("Result changed while reading: ", path, "; retry after publication")
  value
}
read_sheets <- function(path) {
  if (!file.exists(path)) return(list())
  stable_read(path, function(p) setNames(lapply(getSheetNames(p), function(s)
    as.data.table(read.xlsx(p, s, check.names = FALSE))), getSheetNames(p)))
}
prefix_tables <- function(z, prefix) {
  if(!length(z)) return(z)
  setNames(z,paste0(prefix,names(z)))
}
atomic_save <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".publish-", tmpdir = dirname(path)); on.exit(unlink(tmp))
  saveRDS(value, tmp, compress = "gzip")
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
}
copy_result <- function(from, to) {
  if (!file.exists(from)) return(FALSE)
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  before <- result_stamp(from)
  tmp <- tempfile(".publish-", tmpdir = dirname(to)); on.exit(unlink(tmp))
  if (!file.copy(from, tmp, copy.date = TRUE) || !same_stamp(before, result_stamp(from)))
    stop("Result changed while copying: ", from)
  if (!file.rename(tmp, to)) stop("Cannot publish ", to)
  TRUE
}
write_book <- function(tables, path) {
  tables <- tables[vapply(tables, ncol, integer(1)) > 0L]
  if (!length(tables)) tables <- list(status = data.frame(message = "No available aggregate results"))
  names(tables) <- make.unique(substr(names(tables), 1, 29))
  tmp <- tempfile(fileext = ".xlsx", tmpdir = dirname(path)); on.exit(unlink(tmp))
  write.xlsx(lapply(tables, as.data.frame), tmp, overwrite = TRUE, asTable = TRUE)
  if (!file.rename(tmp, path)) stop("Cannot publish ", path)
}
score_table <- function(path) {
  d <- as.data.table(stable_read(path, readRDS))
  id <- intersect(c("eid", "IID", "#IID"), names(d))[1]
  if (is.na(id)) stop("Missing sample ID: ", path)
  setnames(d, id, "eid"); d[, eid := as.character(eid)]
  if (anyNA(d$eid) || anyDuplicated(d$eid)) stop("Missing/duplicate sample IDs: ", path)
  d
}
publish_module <- function(module, analysis_root, score_dir, pca_dir, ancestry_file,
                           traits, subset = "") {
  if (module == "pca") {
    dest <- file.path(analysis_root, "pca")
    names_out <- c("projection.xlsx", "ancestry.xlsx", "projection.png", "ancestry.png")
    legacy <- file.path(pca_dir, c("pca_qc.xlsx", "ancestry_auto_qc.xlsx", "Fig1.PCA_QC.png", "Fig2.Ancestry_QC.png"))
    sources <- ifelse(file.exists(file.path(dest,names_out)),file.path(dest,names_out),legacy)
    if (!any(file.exists(sources[1:2]))) return(invisible(NULL))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    projection <- read_sheets(sources[1]); ancestry <- read_sheets(sources[2])
    tables <- c(prefix_tables(projection,"projection_"),prefix_tables(ancestry,"ancestry_"))
    for (i in seq_along(sources)) if(sources[i]!=file.path(dest,names_out[i])) copy_result(sources[i], file.path(dest, names_out[i]))
    atomic_save(list(module = module, trait = "all", tables = tables, inputs = result_stamp(sources),
      figures = names_out[3:4], generated_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)), file.path(dest, "summary.rds"))
    return(invisible(NULL))
  }
  for (trait in traits) {
    srcdir <- file.path(score_dir, trait, subset)
    source <- file.path(srcdir, if (module == "csx") "1csx.scores.rds" else "2disco.scores.rds")
    if (!file.exists(source)) next
    dest <- file.path(analysis_root, if (module == "csx") Sys.getenv("GRID_CSX_MODULE", "csx") else module, trait, subset)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    dependencies <- c(source, ancestry_file, if (module == "disco") file.path(srcdir, c("1csx.scores.rds","2disco.coefficients.xlsx"))
                      else file.path(srcdir, "1csx.scores.provenance.json"))
    stamp <- result_stamp(dependencies)
    old <- file.path(dest, "summary.rds")
    if (file.exists(old) && all(file.exists(file.path(dest, c("scores.xlsx", "scores.png")))) &&
        same_stamp(readRDS(old)$inputs, stamp) &&
        (module!="disco" || !file.exists(tail(dependencies,1)) || file.exists(file.path(dest,"coefficients.xlsx")))) next
    d <- score_table(source)
    columns <- if (module == "csx") grep("^csx[.]", names(d), value = TRUE) else intersect("disco", names(d))
    if (!length(columns)) stop("No score columns: ", source)
    d[, target := "ALL"]
    if (file.exists(ancestry_file)) {
      a <- fread(ancestry_file, select = c("eid", "genetic_ancestry"), colClasses = c(eid = "character"))
      if (anyDuplicated(a$eid)) stop("Duplicate ancestry IDs")
      d[a, on = "eid", target := i.genetic_ancestry]
      d[is.na(target) | target == "ALL", target := "UNASSIGNED"]
    }
    long <- melt(d[, c("eid", "target", columns), with = FALSE], id.vars = c("eid", "target"), variable.name = "score")
    stats <- long[, {
      x <- value[is.finite(value)]; q <- if (length(x)) quantile(x, c(.01,.25,.5,.75,.99), names = FALSE) else rep(NA_real_,5)
      .(N = .N, observed = length(x), missing = .N-length(x), mean = if(length(x)) mean(x) else NA_real_,
        sd = sd(x), q01 = q[1], q25 = q[2], median = q[3], q75 = q[4], q99 = q[5])
    }, by = .(target, score)]
    tables <- list(score_distribution = stats)
    if (module == "disco" && file.exists(file.path(srcdir, "1csx.scores.rds"))) {
      csx <- score_table(file.path(srcdir, "1csx.scores.rds"))
      cs <- grep("^csx[.]", names(csx), value = TRUE)
      joined <- merge(d[, .(eid, target, disco)], csx[, c("eid", cs), with = FALSE], by = "eid")
      tables$csx_correlation <- rbindlist(lapply(cs, function(s) joined[, {
        ok <- is.finite(disco) & is.finite(get(s)); x <- disco[ok]; y <- get(s)[ok]
        .(score = s, N = sum(ok), correlation = if(length(x)>2 && sd(x)>0 && sd(y)>0) cor(x,y) else NA_real_)
      }, by = target]))
      copy_result(file.path(srcdir, "2disco.coefficients.xlsx"), file.path(dest, "coefficients.xlsx"))
    }
    if (module == "csx" && file.exists(tail(dependencies,1))) {
      copy_result(tail(dependencies,1), file.path(dest, "provenance.json"))
    }
    tables$source <- data.table(field = c("module", "trait", "score_file", "scope", "interpretation"),
      value = c(module, trait, source, if(nzchar(subset)) subset else "autosomes",
                "Score distributions describe model output; predictive accuracy is evaluated separately by final on a common test cohort."))
    write_book(tables, file.path(dest, "scores.xlsx"))
    plot <- ggplot(stats, aes(target, median, colour = score)) +
      geom_pointrange(aes(ymin=q25,ymax=q75), position=position_dodge(width=.6)) +
      facet_wrap(~score, scales="free_y") + theme_minimal(base_size=12) +
      theme(legend.position="none") + labs(title=paste(toupper(module), trait, "score distributions"),
        subtitle="Median and interquartile range; score scales differ between models", x=NULL, y="Score")
    ggsave(file.path(dest, "scores.png"), plot, width=11, height=6, dpi=150)
    if (!same_stamp(stamp, result_stamp(dependencies))) stop("Score inputs changed while summarizing; retry after publication")
    atomic_save(list(module=module, trait=trait, tables=tables, inputs=stamp, figures="scores.png",
      generated_at=format(Sys.time(), tz="Asia/Shanghai", usetz=TRUE)), old)
    cat("REPORT", module, trait, dest, "\n")
  }
}

build_dashboard <- function(analysis_root, final_root, score_dir, pca_dir, ancestry_file, traits,
                            grid_root=file.path(analysis_root,"abm"), prs_root=analysis_root) {
  for (module in c("pca","csx","disco")) publish_module(module, analysis_root, score_dir, pca_dir, ancestry_file, traits)
  dir.create(final_root, recursive=TRUE, showWarnings=FALSE)
  version <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "-", Sys.getpid())
  relative <- file.path("dashboard",version)
  snapshot <- file.path(final_root, relative)
  dir.create(snapshot, recursive=TRUE)
  tables <- list(); inventory <- list(); artifacts <- list(); results <- list(); findings <- list()
  add_tables <- function(module, trait, z) {
    for (n in names(z)) tables[[paste(module,trait,n,sep="/")]] <<- as.data.frame(z[[n]])
  }
  add_assets <- function(module, trait, root, files) {
    for (f in files) if (file.exists(file.path(root,f))) {
      dest <- file.path(relative,module,trait,basename(f))
      copy_result(file.path(root,f),file.path(final_root,dest))
      artifacts[[length(artifacts)+1L]] <<- data.table(module=module,trait=trait,title=tools::file_path_sans_ext(basename(f)),path=dest)
    }
  }
  add_status <- function(module,trait,status,note) inventory[[length(inventory)+1L]] <<-
    data.table(module=module,trait=trait,status=status,note=note)
  add_finding <- function(trait,category,text) findings[[length(findings)+1L]] <<-
    data.table(trait=trait,category=category,text=text)
  for (module in c("pca","csx","disco")) for (trait in if(module=="pca") "all" else traits) {
    root <- if(module=="pca") file.path(analysis_root,module) else file.path(analysis_root,module,trait)
    path <- file.path(root,"summary.rds")
    if (!file.exists(path)) { add_status(module,trait,"missing","尚无可用汇总"); next }
    value <- stable_read(path,readRDS)
    add_tables(module,trait,value$tables)
    add_assets(module,trait,root,value$figures)
    add_status(module,trait,"available",if(module=="pca") "参考空间、投影和祖源 QC" else "已有评分分布；准确性以 final 同队列比较为准")
    if(module=="pca") {
      counts <- value$tables$ancestry_Assigned_counts
      if(!is.null(counts) && all(c("N","genetic_ancestry") %in% names(counts)) && sum(counts$N)>0) {
        largest <- which.max(counts$N)
        add_finding("all","pca",sprintf("PCA 祖源分配共 %s 人；%s 占 %.1f%%。各祖源样本量不同，预测表现需分组比较并结合区间宽度解读。",
          format(sum(counts$N),big.mark=","),counts$genetic_ancestry[largest],100*counts$N[largest]/sum(counts$N)))
      }
    }
  }
  prs <- file.path(prs_root,"prsformer")
  if (identical(prs_root,analysis_root) && !file.exists(file.path(prs,"split.rds")) &&
      !dir.exists(file.path(prs,"scores")) && file.exists(file.path(grid_root,"benchmark/prsformer/3.prsformer.split.rds"))) {
    prs_root <- file.path(grid_root,"benchmark"); prs <- file.path(prs_root,"prsformer")
  }
  new_layout <- file.exists(file.path(prs,"split.rds")) || dir.exists(file.path(prs,"scores"))
  pf_prefix <- if(new_layout) "" else "3.prsformer."
  pf_split <- file.path(prs,paste0(pf_prefix,"split.rds"))
  pf_scores <- setNames(vapply(traits,function(trait) {
    if(new_layout) file.path(prs,"scores",trait,"predictions.rds") else if(identical(prs_root,analysis_root))
      file.path(score_dir,trait,"3.prsformer.scores.rds") else file.path(prs_root,"scores",trait,"3.prsformer.scores.rds")
  },character(1)),traits)
  pf_available <- file.exists(pf_split) & file.exists(pf_scores)
  names(pf_available) <- traits
  # Cached training artifacts alone do not make this optional method available.
  if(any(pf_available)) for (kind in c("performance","training")) {
    path <- file.path(prs,paste0(pf_prefix,kind,".xlsx"))
    z <- read_sheets(path)
    add_tables("prsformer","all",prefix_tables(z,paste0(kind,"_")))
    add_assets("prsformer","all",prs,paste0(pf_prefix,kind,".png"))
  }
  for (trait in traits) {
    pf_score <- pf_scores[[trait]]
    if(pf_available[[trait]]) add_status("prsformer",trait,"available",if(new_layout)
        "独立训练指标仅描述本模型；跨方法比较以 final 为准" else "旧版输出：测试名单兼容性须由 final 验证")
    gr <- file.path(grid_root,trait)
    available <- file.exists(file.path(gr,"grid.test_individuals.rds"))
    add_status("abm",trait,if(available) "available" else "missing","ABM：测试集预测、匹配支持和内部置信度审计")
    for(kind in c("training","confidence")) {
      z <- read_sheets(file.path(gr,paste0("grid.",kind,".xlsx")))
      # Workbook exports contain aggregate audit tables only.
      add_tables("abm",trait,prefix_tables(z,paste0(kind,"_")))
      if(kind=="confidence" && !is.null(z$summary) && nrow(z$summary) &&
          all(c("test_high_confidence_n","test_high_confidence_coverage","test_delta_lower95","test_delta_upper95") %in% names(z$summary))) {
        s <- z$summary[1,]
        add_finding(trait,"abm_confidence",sprintf("ABM 内部测试审计：高置信组 %s 人，占全部测试样本 %.1f%%；相对全部测试样本的损失差区间为 [%.4f, %.4f]。这是组层面的内部验证，不代表个体正确性保证或外部验证。",
          format(s$test_high_confidence_n,big.mark=","),100*s$test_high_confidence_coverage,s$test_delta_lower95,s$test_delta_upper95))
      }
    }
    add_assets("abm",trait,gr,"grid.confidence.png")
    report <- file.path(final_root,trait)
    payload <- file.path(report,"final.summary.rds")
    legacy <- file.path(report,"final.comparison.xlsx")
    if (!file.exists(payload) && !file.exists(legacy)) next
    if(file.exists(payload)) {
      value <- stable_read(payload,readRDS)
      data <- value$tables; sources <- value$inputs
      stale <- !same_stamp(sources,result_stamp(sources$path))
    } else {
      model <- read_sheets(file.path(report,"final.models.xlsx"))
      data <- c(list(performance=read_sheets(legacy)$performance),
        list(paired_comparison=read_sheets(file.path(report,"final.paired_improvement.xlsx"))$paired_comparison),
        list(distance=read_sheets(file.path(report,"final.distance_bins.xlsx"))$distance_estimates),model)
      methods <- data$performance$method
      source_paths <- c(file.path(score_dir,trait,c("1csx.scores.rds","2disco.scores.rds")),
        if("ABM" %in% methods) file.path(gr,"grid.test_individuals.rds"),
        if("PRSformer" %in% methods) c(pf_score,pf_split))
      times <- file.info(source_paths)$mtime
      stale <- any(is.na(times)) || any(times > file.info(legacy)$mtime, na.rm=TRUE)
    }
    perf <- as.data.table(data$performance)
    if(!nrow(perf)) next
    perf[, `:=`(trait=trait, report_status=if(stale) "stale" else "current")]
    results[[trait]] <- perf
    add_tables("final",trait,data)
    add_assets("final",trait,report,list.files(report,pattern="[.]png$"))
    found <- unique(perf$method)
    findings[[length(findings)+1L]] <- data.table(trait=trait,category="comparison",text=paste0(
      if(stale) "比较报告早于部分输入或输入已改变，请重跑 final。" else "已有保存的比较报告。",
      "纳入方法：",paste(found,collapse="、"),"。"))
    # A point-estimate winner is descriptive, never a significance or generalization claim.
    if(!stale) for (group in unique(perf$target)) {
      z <- perf[target==group & is.finite(estimate)]
      if(nrow(z)>1 && length(unique(z$metric))==1) {
        z <- z[order(-estimate)]
        findings[[length(findings)+1L]] <- data.table(trait=trait,category="point_estimate",
          text=sprintf("%s：%s 的 %s 点估计最高（%.4f）；显著性需查看配对差值及其区间。",group,z$method[1],z$metric[1],z$estimate[1]))
      }
    }
    paired <- as.data.frame(data$paired_comparison)
    if(!stale && nrow(paired) && all(c("target","method","reference","difference","lower95","upper95") %in% names(paired)))
      for(i in seq_len(nrow(paired))) {
        z <- paired[i,]
        evidence <- if(!is.finite(z$lower95)||!is.finite(z$upper95)) "未提供区间" else if(z$lower95>0) "区间完全大于零" else if(z$upper95<0) "区间完全小于零" else "区间包含零"
        add_finding(trait,"paired",sprintf("%s：%s 相对 %s 的配对差值 %.4f，95%% 区间 [%.4f, %.4f]（%s；区间条件于已拟合模型）。",
          z$target,z$method,z$reference,z$difference,z$lower95,z$upper95,evidence))
      }
  }
  inv <- rbindlist(inventory,fill=TRUE)
  perf <- rbindlist(results,fill=TRUE)
  findings <- rbindlist(findings,fill=TRUE)
  if(!nrow(findings)) findings <- data.table(trait="all",category="pending",text="尚无最终比较；先完成各方法，再运行 grid.sh final。")
  assets <- if(length(artifacts)) rbindlist(artifacts) else data.table(module=character(),trait=character(),title=character(),path=character())
  value <- list(schema_version=1L, generated_at=format(Sys.time(),tz="Asia/Shanghai",usetz=TRUE),
    traits=traits, inventory=inv, findings=findings, performance=perf, tables=tables, artifacts=assets,
    workbook=file.path(relative,"summary.xlsx"))
  write_book(list(modules=inv,findings=findings,performance=perf),file.path(snapshot,"summary.xlsx"))
  copy_result(file.path(snapshot,"summary.xlsx"),file.path(final_root,"summary.xlsx"))
  esc <- function(x) { x<-gsub("&","&amp;",as.character(x),fixed=TRUE); x<-gsub("<","&lt;",x,fixed=TRUE); gsub('"',"&quot;",x,fixed=TRUE) }
  module_labels <- c(pca="PCA",csx="CSx",disco="Disco",prsformer="PRSformer",abm="ABM")
  module_labels <- module_labels[names(module_labels) %in% inv$module]
  html <- c('<!doctype html><html lang="zh"><meta charset="utf-8"><title>GRID 综合结果</title>',
    '<style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:0 24px;color:#183139;background:#f6f8f7}h1{font-size:36px}section{background:white;padding:24px;margin:20px 0;border-radius:12px}img{max-width:100%}li{margin:10px 0}a{color:#14786d}</style>',
    '<h1>GRID · 综合结果</h1>',paste0('<p>',esc(value$generated_at),' · ',paste(module_labels,collapse=" / "),'</p>'),
    '<p><a href="summary.xlsx">下载综合汇总</a> · 交互浏览：运行 <code>./grid.sh shiny</code></p>',
    '<section><h2>结果解读</h2><ul>',paste0('<li>',esc(findings$trait),'：',esc(findings$text),'</li>'),'</ul></section>',
    '<section><h2>模块状态</h2><ul>',paste0('<li>',esc(inv$module),' / ',esc(inv$trait),' — ',esc(inv$status),'：',esc(inv$note),'</li>'),'</ul></section>',
    '<section><h2>结果图</h2>',vapply(seq_len(nrow(assets)),function(i) paste0('<h3>',esc(assets$module[i]),' / ',esc(assets$trait[i]),' / ',esc(assets$title[i]),'</h3><img loading="lazy" src="',esc(assets$path[i]),'">'),character(1)),'</section></html>')
  html_path <- tempfile(tmpdir=final_root); writeLines(html,html_path); if(!file.rename(html_path,file.path(final_root,"index.html"))) stop("Cannot publish index")
  atomic_save(value,file.path(final_root,"dashboard.rds"))
  cat("FINAL SUMMARY:",file.path(final_root,"index.html"),"\n")
  invisible(value)
}

results_main <- function(args) {
  mode <- args[1]; args <- args[-1]; opt <- list()
  if(length(args) %% 2L) stop("Options require values")
  if(length(args)) for(i in seq(1,length(args),2)) opt[[sub("^--","",args[i])]] <- args[i+1]
  get <- function(n,default) if(is.null(opt[[n]])) default else opt[[n]]
  root <- get("analysis-root","/mnt/d/analysis/grid")
  traits <- strsplit(get("traits","height,ldl,t2dm"),",",fixed=TRUE)[[1]]
  common <- list(analysis_root=root,score_dir=get("score-dir","/mnt/f/gen/ukb/pgs"),
    pca_dir=get("pca-dir","/mnt/f/gen/ukb/pca_proj"),ancestry_file=get("ancestry-file","/mnt/f/gen/ukb/pca_proj/ukb.ancestry.auto.tsv.gz"),traits=traits)
  if(mode=="summary") do.call(build_dashboard,c(common,list(final_root=get("final-root",file.path(root,"final")),
    grid_root=get("grid-root",file.path(root,"abm")),prs_root=get("prsformer-root",root))))
  else if(mode %in% c("pca","csx","disco")) do.call(publish_module,c(list(module=mode),common,list(subset=get("subset",""))))
  else stop("Use pca, csx, disco or summary")
}
if (sys.nframe()==0L) results_main(commandArgs(TRUE))
