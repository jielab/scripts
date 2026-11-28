#!/usr/bin/env Rscript
# Fixed-horizon IPCW logistic LASSO in actual glmnet, using the shared C1 split.
# The input file contains no test or calibration_audit outcome information.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: c1_unified_glmnet.R INPUT_DIR OUTPUT_DIR")
for (p in c("glmnet", "jsonlite")) if (!requireNamespace(p, quietly=TRUE)) stop("Missing R package: ", p)
input <- args[[1]]; output <- args[[2]]
cfg <- jsonlite::read_json(file.path(input, "config.json"), simplifyVector=TRUE)
set.seed(cfg$seed)
con <- file(file.path(input,"X.bin"), "rb")
raw_x <- readBin(con, "double", n=as.double(cfg$n)*cfg$p, size=8L, endian="little")
close(con)
if (length(raw_x) != as.double(cfg$n)*cfg$p) stop("Truncated numeric feature bridge")
X <- matrix(raw_x, nrow=cfg$n, ncol=cfg$p, byrow=TRUE)
rm(raw_x)
if (any(!is.finite(X))) stop("Nonfinite prepared features")
roles <- read.delim(file.path(input,"roles.tsv"), stringsAsFactors=FALSE)
stopifnot(nrow(roles)==nrow(X), all(roles$role %in% c("build","tune_model","tune_gate","calibration_fit","calibration_audit","test")))
if (any(roles$weight[roles$role %in% c("test","calibration_audit")] != 0)) stop("Held-out labels passed to R fit")
b <- which(roles$role=="build" & roles$weight>0)
t <- which(roles$role=="tune_model" & roles$weight>0)
if (length(unique(roles$y[b]))!=2L || length(unique(roles$y[t]))!=2L) stop("Two classes required in build/tune")
# One warm-started path on the original grid. maxit is the pass budget for
# the entire path. Passing symbols keeps the build matrix out of model$call.
budget <- as.integer(cfg$maxit)
stopifnot(length(budget)==1L, is.finite(budget), budget>0L)
fit_path <- function(...) glmnet::glmnet(x=X[b,,drop=FALSE], y=roles$y[b],
    family="binomial", alpha=1, weights=roles$weight[b], standardize=FALSE,
    intercept=TRUE, nlambda=100L, ...)
message(sprintf("[LE8] glmnet build N=%d cases=%d p=%d; tune N=%d; 100 lambdas, total pass budget=%d",
    length(b), sum(roles$y[b]), ncol(X), length(t), budget))
started <- proc.time()[["elapsed"]]
fit <- if ("control" %in% names(formals(glmnet::glmnet))) {
    fit_path(control=list(maxit=budget, trace.it=1L))
} else {
    fit_path(maxit=budget, trace.it=1L)
}
solver_seconds <- proc.time()[["elapsed"]]-started
jerr <- as.integer(fit$jerr)
if (length(jerr)!=1L || is.na(jerr) || jerr > 0L || jerr <= -10000L)
    stop("glmnet solver failure; jerr=", jerr)
partial <- jerr < 0L
grid <- fit$lambda
audit <- list(solver_jerr=jerr, failed_lambda_index=if (partial) -jerr else NULL,
    path_status=if (partial) "converged_prefix" else "solver_completed",
    lambda_requested=100L, lambda_count=length(grid), maxit_used=budget, fit_attempts=1L,
    coordinate_passes=fit$npasses, solver_seconds=solver_seconds, package_version=as.character(utils::packageVersion("glmnet")))
unavailable <- function(reason) {
    jsonlite::write_json(c(audit,list(status="unavailable_lambda_selection",reason=reason)),
        file.path(output,"glmnet_fit.json"), pretty=TRUE,auto_unbox=TRUE,digits=NA,null="null")
    message("[LE8] glmnet unavailable: ", reason)
    quit(save="no", status=0)
}
if (partial && !length(grid)) unavailable("No converged lambda before the iteration limit")
if (any(!is.finite(grid)) || any(grid<=0) || any(diff(grid)>=0) ||
    !length(grid) || (partial && length(grid)>=-jerr)) stop("Invalid glmnet lambda path")
# Without s, glmnet projects the returned coefficients directly (no interpolation).
p <- as.matrix(predict(fit, newx=X[t,,drop=FALSE], type="response"))
if (!identical(dim(p),c(length(t),length(grid))) || any(!is.finite(p))) stop("Invalid tuning predictions")
p <- pmin(pmax(p, 1e-7), 1-1e-7)
y <- roles$y[t]; w <- roles$weight[t]
loss <- colSums((-y*log(p)-(1-y)*log1p(-p))*w)/sum(w)
if (any(!is.finite(loss))) stop("Nonfinite glmnet tuning loss")
chosen <- which.min(loss)
lambda <- grid[[chosen]]
audit$selected_index <- chosen
audit$weakest_converged_loss <- tail(loss,1L)
audit$selected_tune_loss <- loss[chosen]
write.table(data.frame(index=seq_along(grid),lambda=grid,nonzero=fit$df,
    tune_logloss=loss,selected=seq_along(grid)==chosen),
    file.path(output,"glmnet_tuning.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# A truncated path can have an interior minimum while loss has already
# begun falling again in its weak tail. Do not publish that local minimum.
if (partial) {
    weak_tail <- tail(loss, 5L)
    if (length(grid)-chosen < 5L || any(weak_tail <= loss[chosen]+1e-12) ||
        any(diff(weak_tail) < -1e-12))
        unavailable("Pass budget exhausted without five weaker, non-decreasing tuning losses above the minimum")
}
message(sprintf("[LE8] glmnet %s: jerr=%d, passes=%d in %.1f s, selected %d/%d on tune_model; lambda=%.6g",
    audit$path_status, jerr, fit$npasses, solver_seconds, chosen, length(grid), lambda))
raw <- as.numeric(predict(fit,newx=X,type="response",s=lambda))
if (length(raw)!=cfg$n || any(!is.finite(raw))) stop("Invalid prediction output")
con <- file(file.path(output,"raw_probability_LOCAL_ONLY.bin"),"wb")
writeBin(raw, con, size=8L, endian="little"); close(con)
saveRDS(list(model=fit,lambda=lambda,features=cfg$feature_names,seed=cfg$seed,
    fit_role="build",selection_role="tune_model",family="binomial",alpha=1),
    file.path(output,"glmnet_model_LOCAL_ONLY.rds"))
jsonlite::write_json(c(audit,list(status="fitted",family="binomial",alpha=1,lambda=lambda,
    tune_logloss=loss[[chosen]],build_N=length(b),tune_N=length(t),n_features=ncol(X),
    test_outcomes_used=FALSE,calibration="Python calibration_fit only")),
    file.path(output,"glmnet_fit.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA,null="null")
