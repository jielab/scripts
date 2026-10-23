# Content-addressed, atomic checkpoints. Preserve the original serial RNG stream.
med_hash <- function(x) {
  p <- tempfile(); on.exit(unlink(p)); saveRDS(x,p,version=3)
  unname(tools::md5sum(p))
}
med_atomic <- function(x,p) {
  dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile('.pending-',tmpdir=dirname(p)); on.exit(unlink(tmp))
  saveRDS(x,tmp,version=3)
  if(!file.rename(tmp,p)) stop('Cannot publish checkpoint: ',p)
}
med_cached <- function(directory,key,compute,replace=FALSE) {
  p <- file.path(directory,paste0(med_hash(key),'.rds'))
  before <- .Random.seed
  if(!replace && file.exists(p)) {
    saved <- tryCatch(readRDS(p),error=function(e) NULL)
    if(is.list(saved) && identical(saved$key,key) &&
       identical(saved$before,before) && is.integer(saved$after) &&
       identical(saved$checksum,med_hash(list(saved$value,saved$after)))) {
      assign('.Random.seed',saved$after,envir=.GlobalEnv)
      return(list(value=saved$value,reused=TRUE))
    }
    message('[mediation] Invalid/incompatible checkpoint; recomputing: ',p)
  }
  value <- compute()
  after <- .Random.seed
  med_atomic(list(key=key,before=before,after=after,value=value,
                  checksum=med_hash(list(value,after))),p)
  list(value=value,reused=FALSE)
}
