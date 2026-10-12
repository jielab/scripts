grid_read_snapshot <- function(root) {
  path <- file.path(root,"dashboard.rds")
  if(!file.exists(path)) stop("请先运行 ./grid.sh final，生成网页汇总数据。")
  d <- readRDS(path)
  if(!identical(d$schema_version,1L) || !all(c("inventory","findings","tables","artifacts","performance","traits") %in% names(d)))
    stop("不支持的汇总数据格式；请重新运行 final。")
  if(!all(c("module","trait","status","note") %in% names(d$inventory))) stop("汇总模块清单不完整")
  # Keep all downloadable/rendered files inside final, including through symlinks.
  for(path in c(d$artifacts$path,d$workbook)) grid_snapshot_path(root,path)
  d
}
grid_snapshot_path <- function(root,relative) {
  if(length(relative)!=1L || is.na(relative) || grepl("^(/|[A-Za-z]:)|(^|[/\\\\])[.][.]([/\\\\]|$)",relative)) stop("Invalid snapshot artifact path")
  base <- normalizePath(root,mustWork=TRUE)
  path <- normalizePath(file.path(base,relative),mustWork=TRUE)
  if(!startsWith(path,paste0(base,"/"))) stop("Snapshot artifact escapes final directory")
  path
}
grid_performance <- function(snapshot,trait,targets=NULL,methods=NULL) {
  d <- as.data.frame(snapshot$performance)
  if(!nrow(d)) return(d)
  d <- d[d$trait==trait,,drop=FALSE]
  if(!is.null(targets)) d <- d[d$target %in% targets,,drop=FALSE]
  if(!is.null(methods)) d <- d[d$method %in% methods,,drop=FALSE]
  d
}
