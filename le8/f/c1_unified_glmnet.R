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
fit <- glmnet::glmnet(X[b,,drop=FALSE], roles$y[b], family="binomial", alpha=1,
    weights=roles$weight[b], standardize=FALSE, intercept=TRUE, nlambda=100L, maxit=cfg$maxit)
if (!is.null(fit$jerr) && fit$jerr != 0L) stop("glmnet did not complete its lambda path; jerr=",fit$jerr)
p <- pmin(pmax(as.matrix(predict(fit, newx=X[t,,drop=FALSE], type="response")), 1e-7), 1-1e-7)
y <- roles$y[t]; w <- roles$weight[t]
loss <- colSums((-y*log(p)-(1-y)*log(1-p))*w)/sum(w)
chosen <- which.min(loss) # lambda path starts most regularized; deterministic ties
lambda <- fit$lambda[[chosen]]
raw <- as.numeric(predict(fit,newx=X,type="response",s=lambda))
if (length(raw)!=cfg$n || any(!is.finite(raw))) stop("Invalid prediction output")
con <- file(file.path(output,"raw_probability_LOCAL_ONLY.bin"),"wb")
writeBin(raw, con, size=8L, endian="little"); close(con)
saveRDS(list(model=fit,lambda=lambda,features=cfg$feature_names,seed=cfg$seed,
    fit_role="build",selection_role="tune_model",family="binomial",alpha=1),
    file.path(output,"glmnet_model_LOCAL_ONLY.rds"))
jsonlite::write_json(list(package_version=as.character(utils::packageVersion("glmnet")),
    family="binomial",alpha=1,lambda=lambda,tune_logloss=loss[[chosen]],
    build_N=length(b),tune_N=length(t),n_features=ncol(X),lambda_count=length(fit$lambda),
    test_outcomes_used=FALSE,calibration="Python calibration_fit only"),
    file.path(output,"glmnet_fit.json"),pretty=TRUE,auto_unbox=TRUE)
