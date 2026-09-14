# Shared data access
source(file.path(app_dir,"..","f","phyml_panel_b.R"),local=FALSE)
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(plotly); library(DT); library(DBI); library(RSQLite); library(data.table); library(scales); library(ape)
})

`%||%` <- function(x,y) if (is.null(x) || length(x)==0 || is.na(x) || !nzchar(x)) y else x

default_final <- Sys.getenv("GU_FINAL_DIR", file.path(Sys.getenv("GU_ANALYSIS_ROOT", "/mnt/d/analysis/gu"), "final"))
default_db <- normalizePath(file.path(default_final, "gu.sqlite"), mustWork = FALSE)
db_path <- Sys.getenv("GU_SQLITE", unset = default_db)
if (!file.exists(db_path)) {
  stop("GU SQLite database not found: ", db_path, "\nRun: ./gu.sh final")
}
db_path <- normalizePath(db_path, mustWork = TRUE)
normalize_root <- dirname(db_path)
.gu_resolve_artifact <- function(path) {
  value <- as.character(path[[1]])
  if (is.na(value) || !nzchar(value)) return(value)
  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\])", value)
  normalizePath(if (is_absolute) value else file.path(normalize_root, value), mustWork = FALSE)
}

chr_order <- c(as.character(1:22), "X", "Y", "MT")
fmt_bp <- function(x) ifelse(x >= 1e6, paste0(round(x/1e6,2)," Mb"), ifelse(x>=1e3,paste0(round(x/1e3,1)," kb"),paste0(x," bp")))
q <- function(con, sql, params = NULL) {
  if (is.null(params) || length(params) == 0L) DBI::dbGetQuery(con, sql)
  else DBI::dbGetQuery(con, sql, params = params)
}
scalar_q <- function(con, sql, params=NULL) { z <- q(con,sql,params); if(nrow(z)) z[[1]][1] else NA }

reduce_intervals <- function(d) {
  if (!nrow(d)) return(d[, .(start,end)])
  x <- as.data.table(d)[order(start,end), .(start=as.numeric(start),end=as.numeric(end))]
  out <- vector("list", nrow(x)); k <- 0L; s <- x$start[1]; e <- x$end[1]
  if (nrow(x)>1) for(i in 2:nrow(x)) {
    if (x$start[i] <= e) e <- max(e,x$end[i]) else { k<-k+1L; out[[k]]<-c(s,e); s<-x$start[i]; e<-x$end[i] }
  }
  k<-k+1L; out[[k]]<-c(s,e)
  m <- do.call(rbind,out[seq_len(k)])
  data.table(start=as.numeric(m[,1]), end=as.numeric(m[,2]))
}

interval_bp <- function(d) { r <- reduce_intervals(as.data.table(d)); if(!nrow(r)) 0 else sum(r$end-r$start) }
intersection_bp <- function(a,b) {
  a <- reduce_intervals(as.data.table(a)); b <- reduce_intervals(as.data.table(b)); i<-1L;j<-1L;z<-0
  while(i<=nrow(a) && j<=nrow(b)) { z<-z+max(0,min(a$end[i],b$end[j])-max(a$start[i],b$start[j])); if(a$end[i]<b$end[j]) i<-i+1L else j<-j+1L }
  z
}

gu_card <- function(...) bslib::card(...,fill=FALSE)

# GU evidence review add-on (read-only generated HTML; no database mutation).
.gu_review_root <- Sys.getenv("GU_PHYML_REPORT_DIR", file.path(dirname(db_path), "review"))
source(file.path(app_dir,"phyml.R"),local=TRUE)
if (file.exists(file.path(.gu_review_root, "review.html"))) {
  shiny::addResourcePath("gu_evidence_review", .gu_review_root)
}

source(file.path(app_dir,"dual_lead.R"),local=TRUE)

source(file.path(app_dir,"methods.R"),local=TRUE)
source(file.path(app_dir,"density.R"),local=TRUE)
source(file.path(app_dir,"summary.R"),local=TRUE)
.gu_summary_land <- as.data.frame(data.table::fread(file.path(app_dir,"www","maps","ne_110m_land.tsv")))
.gu_summary_locations <- as.data.frame(data.table::fread(file.path(app_dir,"www","maps","1kg_populations.tsv")))
shiny::addResourcePath("gu_assets",file.path(app_dir,"www"))
# Use the existing indexed reference without copying multi-GB FASTA files.
.gu_reference_dir <- file.path(Sys.getenv("GU_REF_ROOT","/mnt/e/refGen"),"fasta")
source(file.path(app_dir,"ref.R"),local=TRUE)
.gu_reference_resources <- gu_local_igv_references(.gu_reference_dir)
.gu_local_references <- .gu_reference_resources$references
IGV_region_max <- as.numeric(Sys.getenv("GU_IGV_REGION_MAX", "2000000"))
stopifnot(is.finite(IGV_region_max), IGV_region_max > 0)

