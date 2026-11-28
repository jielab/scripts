#!/usr/bin/env Rscript
# Actual glmnet on shared consensus-block scores; no held-out phenotype input.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop("Usage: block_glmnet.R INPUT_DIR OUTPUT_DIR")
for (p in c("glmnet","jsonlite")) if (!requireNamespace(p,quietly=TRUE)) stop("Missing R package: ",p)
input <- args[[1]]; output <- args[[2]]
cfg <- jsonlite::read_json(file.path(input,"config.json"),simplifyVector=TRUE)
set.seed(cfg$seed)
con <- file(file.path(input,"X.bin"),"rb")
raw <- readBin(con,"double",n=as.double(cfg$n)*cfg$p,size=8L,endian="little");close(con)
if(length(raw)!=as.double(cfg$n)*cfg$p) stop("Truncated block feature matrix")
X <- matrix(raw,nrow=cfg$n,ncol=cfg$p,byrow=TRUE);rm(raw)
labels <- read.delim(file.path(input,"labels.tsv"))
stopifnot(nrow(labels)==nrow(X),all(labels$role %in% c("train","tune")),
          all(is.finite(X)),all(is.finite(labels$y)),all(is.finite(labels$base)))
tr <- which(labels$role=="train");tu <- which(labels$role=="tune")
clip <- function(x) pmin(pmax(x,1e-7),1-1e-7)
loss <- function(y,p) if(cfg$binary) mean(-y*log(clip(p))-(1-y)*log(1-clip(p))) else mean((y-p)^2)
best_loss <- loss(labels$y[tu],labels$base[tu])
best_coef <- rep(0,ncol(X)+1L);best_alpha <- NA_real_;best_lambda <- NA_real_;best_blend <- 0
records <- list();nr <- 0L;best_model <- NULL
fit_alpha <- function(a) {
    best_loss <- loss(labels$y[tu],labels$base[tu])
    best_coef <- rep(0,ncol(X)+1L);best_alpha <- NA_real_;best_lambda <- NA_real_;best_blend <- 0
    records <- list();nr <- 0L;best_model <- NULL

    if(cfg$binary) {
        if(length(unique(labels$y[tr]))!=2L || length(unique(labels$y[tu]))!=2L) stop("Both classes required")
        fit <- glmnet::glmnet(X[tr,,drop=FALSE],labels$y[tr],family="binomial",alpha=a,
            offset=qlogis(clip(labels$base[tr])),standardize=FALSE,nlambda=cfg$nlambda,maxit=cfg$maxit)
        correction <- as.matrix(predict(fit,newx=X[tu,,drop=FALSE],type="link",newoffset=rep(0,length(tu))))
    } else {
        fit <- glmnet::glmnet(X[tr,,drop=FALSE],labels$y[tr]-labels$base[tr],
            family="gaussian",alpha=a,standardize=FALSE,nlambda=cfg$nlambda,maxit=cfg$maxit)
        correction <- as.matrix(predict(fit,newx=X[tu,,drop=FALSE],type="response"))
    }
    if(!is.null(fit$jerr) && fit$jerr!=0L) stop("glmnet path did not converge; jerr=",fit$jerr)
    for(b in cfg$blend_grid[cfg$blend_grid>0]) for(j in seq_along(fit$lambda)) {
        predicted <- if(cfg$binary) plogis(qlogis(clip(labels$base[tu]))+b*correction[,j]) else labels$base[tu]+b*correction[,j]
        value <- loss(labels$y[tu],predicted);nr <- nr+1L
        records[[nr]] <- data.frame(alpha=a,lambda=fit$lambda[j],blend=b,tune_loss=value)
        if(is.finite(value) && value<best_loss) {
            best_loss <- value;best_alpha <- a;best_lambda <- fit$lambda[j];best_blend <- b
            best_coef <- as.numeric(stats::coef(fit,s=best_lambda));best_model <- fit
        }
    }

    list(loss=best_loss,coef=best_coef,alpha=best_alpha,lambda=best_lambda,
         blend=best_blend,model=best_model,records=records)
}
workers <- if(is.null(cfg$workers)) 1L else max(1L,min(as.integer(cfg$workers),length(cfg$alpha_grid)))
cat(sprintf("[glmnet] independent alpha processes=%d; threads/process=1\n",workers))
fits <- if(workers>1L && .Platform$OS.type=="unix") {
    parallel::mclapply(cfg$alpha_grid,fit_alpha,mc.cores=workers,mc.set.seed=FALSE)
} else lapply(cfg$alpha_grid,fit_alpha)
if(any(vapply(fits,inherits,logical(1),"try-error"))) stop("A parallel glmnet alpha fit failed")
# Strict improvement and input order preserve the serial tie-breaking rule.
for(result in fits) {
    records <- c(records,result$records)
    if(is.finite(result$loss) && result$loss<best_loss) {
        best_loss <- result$loss;best_coef <- result$coef;best_alpha <- result$alpha
        best_lambda <- result$lambda;best_blend <- result$blend;best_model <- result$model
    }
}
write.table(data.frame(feature=0:ncol(X),coefficient=best_coef),file.path(output,"coefficients.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,records),file.path(output,"tuning.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
saveRDS(list(model=best_model,lambda=best_lambda,alpha=best_alpha,blend=best_blend),file.path(output,"model.rds"))
jsonlite::write_json(list(implementation="actual R glmnet",version=as.character(utils::packageVersion("glmnet")),
    family=if(cfg$binary) "binomial with OOF baseline offset" else "gaussian OOF residual",
    alpha=best_alpha,lambda=best_lambda,blend=best_blend,tune_loss=best_loss,
    selection_metric=if(cfg$binary) "logloss" else "MSE",
    train_n=length(tr),tune_n=length(tu),test_outcomes_used=FALSE),
    file.path(output,"fit.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
