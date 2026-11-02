#!/usr/bin/env Rscript


# 🚩 gwas.inputs
# A common numeric design matrix gives all methods identical covariates.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
if(!length(a) %in% c(8, 10)) stop('Expected phe cov trait type event covs cats output [continuous_file categorical_file]')
d <- fread(a[1], colClasses = c(FID = 'character', IID = 'character'))
continuous_file <- if(length(a) == 10 && nzchar(a[9])) a[9] else a[2]
categorical_file <- if(length(a) == 10 && nzchar(a[10])) a[10] else a[2]
c <- fread(continuous_file, colClasses = c(FID = 'character', IID = 'character'))
k <- if(a[7] == '' || categorical_file == continuous_file) c else
	fread(categorical_file, colClasses = c(FID = 'character', IID = 'character'))
trait <- a[3]; type <- a[4]; event <- a[5]
csv <- function(x, columns) {
	if(x == '') character() else if(x == '*') setdiff(columns, c('FID', 'IID')) else strsplit(x, ',', fixed = TRUE)[[1]]
}
cats <- unique(csv(a[7], names(k)))
covs <- unique(c(csv(a[6], names(c)), cats))
continuous <- setdiff(covs, cats)
need <- c('FID', 'IID', trait, if(type == 't2e') event)
if(length(setdiff(need, names(d)))) stop('Missing phenotype columns: ', paste(setdiff(need, names(d)), collapse = ','))
check_covariates <- function(data, columns, label) {
	missing <- setdiff(c('FID', 'IID', columns), names(data))
	if(length(missing)) stop('Missing ', label, ' covariate columns: ', paste(missing, collapse = ','))
	if(anyDuplicated(data[, .(FID, IID)])) stop('Duplicate participant keys in ', label, ' covariate file')
}
check_covariates(c, continuous, 'continuous')
if(length(cats)) check_covariates(k, cats, 'categorical')
if(anyDuplicated(d[, .(FID, IID)])) stop('Duplicate participant keys')
if(anyDuplicated(d$IID)) stop('IID must be unique for SAIGE')
d <- merge(d[, ..need], c[, c('FID', 'IID', continuous), with = FALSE], by = c('FID', 'IID'))
if(length(cats)) d <- merge(d, k[, c('FID', 'IID', cats), with = FALSE], by = c('FID', 'IID'))
for(v in c(trait, if(type == 't2e') event, setdiff(covs, cats))) {
	old <- d[[v]]; num <- suppressWarnings(as.numeric(old))
	if(any(!is.na(old) & is.na(num))) stop('Nonnumeric phenotype/covariate: ', v)
	set(d, j = v, value = num)
}
d <- d[complete.cases(d)]
for(v in c(trait, if(type == 't2e') event, setdiff(covs, cats))) d <- d[is.finite(get(v))]
if(type == 'bt' && !all(d[[trait]] %in% 0:1)) stop('Binary phenotype must use 0/1/NA')
if(type == 't2e' && (!all(d[[event]] %in% 0:1) || any(d[[trait]] <= 0))) stop('T2E requires positive duration and 0/1 event')
if(nrow(d) < 10 || uniqueN(d[[trait]]) < 2) stop('Insufficient phenotype variation/sample size')
for(v in cats) set(d, j = v, value = factor(d[[v]]))
variable <- covs[vapply(covs, function(v) uniqueN(d[[v]]) > 1, logical(1))]
if(length(variable)) {
	z <- as.data.frame(d[, ..variable])
	# Keep readable names in analysis.phe and the actual GWAS command/log.
	names(z) <- make.names(variable, unique = TRUE)
	m <- model.matrix( ~ ., z)[, - 1, drop = FALSE]
	# Drop dependent columns after complete-case filtering (e.g. male-only bald).
	full <- cbind(Intercept = 1, m); qr0 <- qr(full); keep <- sort(qr0$pivot[seq_len(qr0$rank)])
	m <- full[, setdiff(keep, 1L), drop = FALSE]
	colnames(m) <- make.names(colnames(m), unique = TRUE)
} else m <- matrix(nrow = nrow(d), ncol = 0)
out <- cbind(d[, ..need], as.data.table(m))
fwrite(out, a[8], sep = '\t', na = 'NA', quote = FALSE)
writeLines(paste(colnames(m), collapse = ','), paste0(a[8], '.covars'))
writeLines(c(paste('trait', trait), paste('type', type), paste('complete_case_n', nrow(d)), paste('design_columns', ncol(m))), paste0(a[8], '.qc.txt'))
