# Paired prediction evaluation. SNP weights are fixed; all outcome fitting is OOF.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(survival);library(pROC);library(patchwork)})
setDTthreads(4)
args <- commandArgs(TRUE)
allowed <- c('trait', 'type', 'method', 'score-dir', 'pgs-file', 'disco-file', 'pt-file', 'pheno-file',
 'ancestry-file', 'group-col', 'covar-name', 'phenotype-col', 'event-col', 'time-col', 'prevalence',
 'pca-file', 'med-file', 'distance-pcs', 'distance-bins', 'min-bin-events', 'folds', 'seed', 'bootstrap', 'min-n',
 'remove', 'out-root', 'dir-gwas', 'dir-gen', 'pt-effect', 'threads', 'check', 'allow-missing-scores', 'disco-tune', 'disco-a', 'min-anchor', 'write-predictions', 'run-dir', 'grid-file', 'posterior-file', 'posterior-mode', 'genetic-variance-file', 'training-centers', 'pca-space', 'allow-chromosome-subset', 'individual-metric', 'individual-max-points', 'distance-source')
opt <- list(); i <- 1L
while(i <= length(args)) {
 key <- sub('^--', '', args[i]); if(!key %in% allowed)stop('Unknown option: ', args[i])
 if(key %in% c('check', 'allow-missing-scores')) {opt[[key]] <- TRUE;i <- i + 1L} else {
	if(i == length(args))stop('Missing value: ', args[i]);opt[[key]] <- args[i + 1L];i <- i + 2L
 }
}
arg <- function(k, default = NULL)if(is.null(opt[[k]]))default else opt[[k]]
Y <- arg('trait'); type <- arg('type')
if(is.null(Y) || !type %in% c('ct', 'dt', 't2e'))stop('--trait and --type ct|dt|t2e required')
if(!arg('method', 'all') %in% c('all', 'csx', 'disco'))stop('Invalid --method')
out <- arg('run-dir', file.path(arg('out-root', '/mnt/d/analysis/grid/Yeval'), Y))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, name)fwrite(x, file.path(out, name), sep = '\t', na = 'NA')
intarg <- function(k, v, minimum) {z <- suppressWarnings(as.integer(arg(k, v)));if(is.na(z) || z < minimum)stop('Invalid --', k);z}
nfold <- intarg('folds', 5, 2); seed <- intarg('seed', 20260904, 0)
nboot <- intarg('bootstrap', 200, 0); minn <- intarg('min-n', 100, 20)
npc <- intarg('distance-pcs', 10, 2); nbins <- intarg('distance-bins', 10, 2)
min_bin_events <- intarg('min-bin-events', 20, 5)
tune <- toupper(arg('disco-tune', 'TRUE'))
if(!tune %in% c('TRUE', 'FALSE'))stop('--disco-tune must be TRUE/FALSE')
tune <- tune == 'TRUE'; min_anchor <- intarg('min-anchor', 100, 10)
quality <- as.numeric(strsplit(arg('disco-a', '1,1,1,1'), ',', fixed = TRUE)[[1]])
if(length(quality) != 4 || any(!is.finite(quality) | quality < 0) || !any(quality > 0))stop('Invalid --disco-a (AFR,EAS,EUR,SAS order)')
partial <- isTRUE(arg('allow-missing-scores', FALSE))
pops <- c('EUR', 'AFR', 'EAS', 'SAS'); base_scores <- paste0('csx.', c('AFR', 'EAS', 'EUR', 'SAS'))
quality <- setNames(quality, c('AFR', 'EAS', 'EUR', 'SAS'))[pops]
main_methods <- c('COJO', 'PRS-CSx-auto-meta', 'PRS-CSx', if(tune)'DiscoDivas-tuned' else 'DiscoDivas-untuned')
score_dir <- file.path(arg('score-dir', '/mnt/d/data/ukb/pgs'), Y)
files <- c(phenotype = arg('pheno-file', '/mnt/d/data/ukb/phe/Rdata/all.rds'),
 csx = arg('pgs-file', file.path(score_dir, 'csx.pgs.gz')),
 disco = arg('disco-file', file.path(score_dir, 'disco.pgs.gz')),
 pt = arg('pt-file', file.path(score_dir, 'pt.pgs.gz')),
 ancestry = arg('ancestry-file', '/mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz'),
 pca = arg('pca-file', '/mnt/d/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz'),
 centers = arg('med-file', '/mnt/d/files/DiscoDivas/med.g1000.4pop.tsv'))
read_table <- function(f) {
 if(!file.exists(f))stop('Missing input: ', f)
 if(grepl('\\.rds$', f, ignore.case = TRUE))return(as.data.table(readRDS(f)))
 header <- names(fread(f, nrows = 0))
 fread(f, colClasses = list(character = intersect(header, c('eid', 'IID', '#IID', 'FID', '#FID'))), na.strings = c('NA', 'NaN', ''))
}
ids <- function(d) {
 nm <- intersect(c('eid', 'IID', '#IID', 'ID_2'), names(d));if(!length(nm))stop('Missing sample ID')
 if(nm[1] != 'eid')setnames(d, nm[1], 'eid')
 d[, eid := as.character(eid)];if(anyNA(d$eid) || anyDuplicated(d$eid))stop('Missing/duplicate sample IDs')
 d
}
num <- function(x, name) {y <- suppressWarnings(as.numeric(as.character(x)));if(any(!is.na(x) & is.na(y)))stop('Non-numeric: ', name);y}


# 🚩 posterior
# Joint posterior propagation and distance definitions; loaded by Yeval.R.
posterior_pops <- c('AFR', 'EAS', 'EUR', 'SAS')
cov_pairs <- t(combn(1:4, 2))
cov_columns <- unlist(lapply(1:4, function(j)paste0('cov.', posterior_pops[j], '.', posterior_pops[j:4])))

load_posterior <- function() {
 mode <- arg('posterior-mode', 'required')
 if(!mode %in% c('required', 'off'))stop('--posterior-mode must be required or off')
 if(mode == 'off')return(NULL)
 f <- arg('posterior-file', file.path(score_dir, 'csx.posterior.tsv.gz'))
 if(!file.exists(f))stop('Missing individual posterior: ', f, '. Run 1csx.sh --posterior TRUE; --posterior-mode off is only for a benchmark without individual posterior panels.')
 meta <- read_table(paste0(f, '.metadata.tsv'))
 md <- setNames(meta$value, meta$field)
 if(md[['schema']] != 'grid_csx_moments_v1' || md[['centering']] != 'discovery_EAF' ||
		!identical(unname(md[['md5']]), unname(tools::md5sum(f))))stop('Invalid/stale posterior metadata or centering')
 if(arg('allow-chromosome-subset', 'FALSE') != 'TRUE' && !setequal(strsplit(md[['chromosomes']], ',', fixed = TRUE)[[1]], as.character(1:22)))stop('Posterior lacks all 22 autosomes; use --allow-chromosome-subset TRUE only for deliberate subset analyses')
 z <- ids(read_table(f));req <- c(paste0('csx.', posterior_pops), cov_columns)
 if(!all(req %in% names(z)) || any(!is.finite(as.matrix(z[, ..req]))))stop('Invalid posterior means/covariance')
 # Every individual covariance must be positive semidefinite. Cholesky-like
 # eigenvalue checks on 4x4 matrices are inexpensive relative to scoring.
 mat <- as.matrix(z[, ..cov_columns]);diag_ix <- c(1L, 5L, 8L, 10L)
 if(any(mat[, diag_ix] < 0))stop('Negative posterior variance')
 for(j in 1:4)for(k in j:4){
	v <- z[[paste0('cov.', posterior_pops[j], '.', posterior_pops[k])]]
	bound <- sqrt(z[[paste0('cov.', posterior_pops[j], '.', posterior_pops[j])]] * z[[paste0('cov.', posterior_pops[k], '.', posterior_pops[k])]])
	if(any(abs(v) > bound + 1e-9 * pmax(1, bound)))stop('Posterior covariance violates Cauchy-Schwarz')
 }
 # The producer constructs PSD covariance from synchronized scores. Full PSD
 # checks on a deterministic spread of rows detect malformed imported files.
 for(i in unique(round(seq(1, nrow(z), length.out = min(1000, nrow(z)))))){
	s <- matrix(0, 4, 4);at <- 0L
	for(j in 1:4)for(k in j:4){at <- at + 1L;s[j, k] <- s[k, j] <- mat[i, at]}
	if(min(eigen(s, symmetric = TRUE, only.values = TRUE)$values) <  - 1e-8 * max(1, max(diag(s))))stop('Non-PSD individual covariance')
 }
 z
}

training_geometry <- function(gd, pm, pc, reference_centers) {
 path <- arg('training-centers', file.path(score_dir, 'training_centers.tsv'))
 distance_source <- arg('distance-source', 'reference')
 if(!distance_source %in% c('discovery', 'reference'))stop('--distance-source must be discovery or reference')
 if(distance_source == 'discovery' && !file.exists(path))stop('Discovery centres missing: ', path, '. Build them with f/training_geometry.py, or use the default --distance-source reference for an exploratory plot.')
 if(distance_source == 'discovery'){
	ct <- read_table(path);need <- c('POP', pc, 'N_GWAS', 'pca_space', 'source', 'kind')
	if(!all(need %in% names(ct)) || anyDuplicated(ct$POP) || !all(posterior_pops %in% ct$POP))stop('Training centres need unique POP, PC1.., N_GWAS, pca_space, source, kind')
	ct <- ct[match(posterior_pops, POP)]
	if(any(ct$kind != 'discovery') || any(!nzchar(ct$source)) || anyNA(ct$source))stop('Training centres must explicitly document discovery provenance')
	space <- arg('pca-space', '')
	if(!nzchar(space) || any(ct$pca_space != space))stop('--pca-space must match discovery centres and the actual target projection basis')
	if(any(!is.finite(ct$N_GWAS) | ct$N_GWAS <= 0) || any(!is.finite(as.matrix(ct[, ..pc]))))stop('Invalid discovery centres/N_GWAS')
	weights <- ct$N_GWAS / sum(ct$N_GWAS)
	status <- 'discovery';label <- 'RMS distance to GWAS training groups'
	description <- 'Sample-size-weighted RMS distance to discovery-group means in the supplied common PC space; a multi-training extension, not the original single-training distance.'
 }else{
	ct <- copy(reference_centers)[match(posterior_pops, POP)]
	weights <- rep(.25, 4);status <- 'reference_proxy';label <- 'RMS distance to four 1KG reference groups'
	description <- 'Equal-weight four-reference distance. Reference centres were selected; this axis is a reference proxy and must not be called distance to GWAS training data.'
 }
 distmat <- vapply(seq_len(4), function(j)sqrt(rowSums(sweep(pm, 2, as.numeric(ct[j, ..pc]), '-') ^ 2)), numeric(nrow(pm)))
 for(j in 1:4)gd[, (paste0('distance.training.', posterior_pops[j])) := distmat[, j]]
 gd[, distance.analysis := sqrt(as.numeric(distmat ^ 2 %*% weights))]
 ct[, mixture_weight := weights]
 list(gd = gd, centers = ct, status = status, label = label, description = description)
}

read_variance <- function() {
 f <- arg('genetic-variance-file', file.path(score_dir, 'genetic_variance.tsv'))
 if(!file.exists(f))return(NULL)
 v <- read_table(f)
 if(!all(c('trait', 'target', 'scale', 'source') %in% names(v)))stop('Genetic variance table needs trait,target,scale,source, and h2 or genetic_variance')
 v <- v[trait == Y]
 if(anyDuplicated(v$target) || anyNA(v$source) || any(!nzchar(v$source)))stop('Duplicate variance target or missing provenance')
 if(!'h2' %in% names(v))v[, h2 := NA_real_]
 if(!'genetic_variance' %in% names(v))v[, genetic_variance := NA_real_]
 if(any(is.finite(v$h2) & (!is.finite(v$h2) | v$h2 <= 0 | v$h2 >= 1)))stop('Require 0 < SNP h2 < 1')
 if(any(is.finite(v$genetic_variance) & v$genetic_variance <= 0))stop('Genetic variance must be positive')
 if(any(!is.finite(v$h2) & !is.finite(v$genetic_variance)))stop('Every supplied variance row needs a finite h2 or genetic_variance; template NA values must be filled')
 if(any(is.finite(v$h2) & is.finite(v$genetic_variance)))stop('Supply h2 OR genetic_variance per row, not both')
 expected <- switch(type, ct = 'residual_phenotype', dt = 'log_odds', t2e = 'log_hazard')
 if(any(v$scale != expected))stop('Variance scale must be ', expected, '; liability h2 cannot be substituted for a log-hazard/log-odds variance')
 if(type != 'ct' && any(is.finite(v$h2)))stop('For dt/t2e supply absolute genetic_variance on log-odds/log-hazard scale, not h2')
 v
}

posterior_variance <- function(z, w) {
 ans <- numeric(nrow(z));diagonal <- numeric(nrow(z))
 for(j in 1:4)for(k in j:4){
	term <- w[j] * w[k] * z[[paste0('cov.', posterior_pops[j], '.', posterior_pops[k])]]
	ans <- ans + if(j == k)term else 2 * term
	if(j == k)diagonal <- diagonal + term
 }
 if(any(ans <  - 1e-8 * pmax(1, diagonal)))stop('Negative combined variance: non-PSD input covariance')
 list(total = pmax(ans, 0), diagonal = diagonal)
}

individual_posterior <- function(te, tr, fit, baseline_fit, ss, g, k) {
 w <- as.numeric(coef(fit)[paste0('z', 1:4)]) / ss
 vv <- posterior_variance(te, w)
 prior <- NA_real_;source <- 'not supplied'
 if(!is.null(variance_spec)){
	row <- variance_spec[target == g]
	if(nrow(row)){
	 prior <- row$genetic_variance;source <- row$source
	 if(type == 'ct' && is.finite(row$h2))prior <- row$h2 * var(residuals(baseline_fit))
	}
 }
 r2 <- if(is.finite(prior) && prior > 0)1 - vv$total / prior else rep(NA_real_, nrow(te))
 # Negative reliability is diagnostic of scale/prior/calibration problems;
 # preserve it. Never clip or force a decreasing curve.
 out <- te[, c('eid', 'target', 'fold', 'proj_PC1', 'proj_PC2', 'distance.analysis', paste0('distance.training.', posterior_pops)), with = FALSE]
 out[, `:=`(posterior_variance = vv$total, posterior_sd = sqrt(vv$total),
					 diagonal_only_variance = vv$diagonal, cross_population_covariance = vv$total - vv$diagonal,
					 prior_genetic_variance = prior, individual_R2 = r2, variance_source = source)]
 for(j in 1:4)out[, (paste0('weight.', posterior_pops[j])) := w[j]]
 out
}


variance_spec <- read_variance()
requested_individual <- arg('individual-metric', 'sd')
if(!requested_individual %in% c('auto', 'reliability', 'sd'))stop('Invalid --individual-metric')
if(requested_individual == 'reliability' && (is.null(variance_spec) || !nrow(variance_spec)))stop('Individual model-based R2 requires --genetic-variance-file; use explicit --individual-metric sd to plot posterior uncertainty without claiming accuracy')
posterior_table <- load_posterior()
individual_results <- list()
fold_coefficients <- list()
cat('Reading full phenotype cohort: ', files['phenotype'], '\n', sep = '')
phe <- ids(read_table(files['phenotype']))
gc <- arg('group-col', 'genetic_ancestry')
if(!gc %in% names(phe)) {
 anc <- ids(read_table(files['ancestry']));if(!gc %in% names(anc))stop('Missing ancestry column: ', gc)
 phe <- merge(phe, anc[, c('eid', gc), with = FALSE], by = 'eid', all.x = TRUE)
}
phe[, target := as.character(get(gc))]
phe[is.na(target) | !nzchar(target), target := 'UNASSIGNED']
rem <- arg('remove', '/mnt/d/files/ukb.exclude.id'); excluded <- character()
if(nzchar(rem)) {
 if(!file.exists(rem))stop('Missing withdrawal file: ', rem)
 if(file.info(rem)$size > 0) {r <- fread(rem, header = FALSE, colClasses = 'character');excluded <- r[[min(2L, ncol(r))]]}
}
phe <- phe[!startsWith(eid, '-') & !eid %in% excluded]
covars <- trimws(strsplit(arg('covar-name', 'age,sex,PC1,PC2'), ',', fixed = TRUE)[[1]])
covars[covars == 'PC'] <- 'PC1';covars <- unique(covars[nzchar(covars) & covars != 'none'])
yc <- arg('phenotype-col', Y); ec <- arg('event-col', paste0(Y, '.Yt2e')); tc <- arg('time-col', paste0(Y, '.t2e'))
outcome_definition <- if(type == 't2e')paste(ec, tc) else yc
# Yt2e is incident-only. Baseline disease includes dated prevalent cases, and
# treats subsequently diagnosed people as baseline non-cases, never as missing.
if(type == 'dt' && Y == 't2dm' && is.null(arg('phenotype-col'))) {
 needed <- c('t2dm.Yr2e', 't2dm.Yt2e');if(!all(needed %in% names(phe)))stop('Need explicit 0/1 --phenotype-col or t2dm.Yr2e/Yt2e')
 phe[, t2dm := fifelse(get('t2dm.Yr2e') == 1, 1,
					fifelse(get('t2dm.Yt2e') %in% c(0, 1) | get('t2dm.Yr2e') == 0, 0, NA_real_), na = NA_real_)]
 # fifelse NA condition above must not erase incident cases (Yr2e=NA).
 phe[is.na(get('t2dm.Yr2e')) & get('t2dm.Yt2e') %in% c(0, 1), t2dm := 0]
 outcome_definition <- 'Baseline ICD10 T2D: Yr2e=1 case; valid nonprevalent Yt2e=0/1 control; undated/invalid NA'
}
required <- unique(c(if(type == 't2e')c(ec, tc) else yc, covars))
if(!all(required %in% names(phe)))stop('Missing phenotype/covariates: ', paste(setdiff(required, names(phe)), collapse = ','))
if(length(intersect(covars, c(yc, ec, tc, base_scores, 'csx.auto', 'csx.meta', 'disco'))))stop('Outcome/score cannot be a covariate')
phe[, outcome := num(get(if(type == 't2e')ec else yc), 'outcome')]
if(type != 'ct' && any(!is.na(phe$outcome) & !phe$outcome %in% c(0, 1)))stop('Binary/event outcome must be 0/1')
if(type == 't2e')phe[, time := num(get(tc), tc)] else phe[, time := NA_real_]
# K is estimated before score matching and covariate complete-case selection.
prevalence <- phe[!is.na(target) & is.finite(outcome), .(population_N = .N, population_cases = sum(outcome == 1), K = mean(outcome)), by = target]
prevalence[, source := 'Full phenotype cohort, before PGS/covariate filtering (UKB cohort assumption)']
if(type == 'dt') {
 spec <- arg('prevalence', 'cohort')
 if(spec != 'cohort') {
	if(grepl('=', spec, fixed = TRUE)) {
	 kv <- strsplit(strsplit(spec, ',', fixed = TRUE)[[1]], '=', fixed = TRUE)
	 if(any(lengths(kv) != 2))stop('Invalid prevalence mapping')
	 ks <- setNames(vapply(kv, function(z)as.numeric(z[2]), numeric(1)), vapply(kv, `[`, character(1), 1))
	 if(anyDuplicated(names(ks)))stop('Duplicate prevalence target')
	 prevalence[, K := ks[target]]
	} else prevalence[, K := as.numeric(spec)]
	prevalence[, source := 'User-specified population prevalence']
 }
 if(any(!is.finite(prevalence[target %in% pops]$K) | prevalence[target %in% pops]$K <= 0 | prevalence[target %in% pops]$K >= 1))stop('Each target needs prevalence 0<K<1')
}
keep <- unique(c('eid', 'target', 'outcome', 'time', covars))
d <- phe[, ..keep]; rm(phe); invisible(gc(verbose = FALSE))
for(v in covars)if(is.character(d[[v]]) || is.factor(d[[v]]) || v == 'sex')set(d, j = v, value = factor(d[[v]]))
audit <- list(data.table(stage = 'phenotype_after_withdrawal', target = d$target)[, .N, by = .(stage, target)])
pgs <- ids(read_table(files['csx']))
if(!is.null(posterior_table)) {
 # Means and variances must describe the SAME scores. These are discovery-
 # centred scores with mean-imputed missing genotypes, not the old uncentred sums.
 oldcols <- intersect(base_scores, names(pgs));pgs[, (oldcols) := NULL]
 pgs <- merge(pgs, posterior_table, by = 'eid', all = FALSE)
 if(!nrow(pgs))stop('No IDs shared by scores and posterior moments')
}
missing <- setdiff(c(base_scores, 'csx.auto', 'csx.meta'), names(pgs))
if(length(missing) && !partial)stop('Missing CSx columns: ', paste(missing, collapse = ', '), '. Complete upstream scores or explicitly use --allow-missing-scores.')
available <- intersect(c(base_scores, 'csx.auto', 'csx.meta'), names(pgs))
d <- merge(d, pgs[, c('eid', available, if(!is.null(posterior_table))cov_columns), with = FALSE], by = 'eid');rm(pgs, posterior_table)
for(kind in c('disco', 'pt')) {
 sc <- if(kind == 'disco')'disco' else paste0('pt.', pops)
 if(file.exists(files[kind])) {
	z <- ids(read_table(files[kind]));bad <- setdiff(sc, names(z))
	if(length(bad) && !partial)stop('Missing ', kind, ' columns: ', paste(bad, collapse = ','))
	use <- intersect(sc, names(z));available <- c(available, use);missing <- c(missing, bad)
	d <- merge(d, z[, c('eid', use), with = FALSE], by = 'eid', all.x = TRUE)
 } else {
	missing <- c(missing, sc)
	if(!partial && !isTRUE(arg('check', FALSE)) && !(kind == 'disco' && tune))stop('Missing ', files[kind])
 }
}
grid_scores <- character()
if(!is.null(arg('grid-file'))){
 z <- ids(read_table(arg('grid-file')));grid_scores <- intersect(c(paste0('GRID_', c('AFR', 'EAS', 'EUR', 'SAS')), 'GRID_shared', 'GRID_posterior', 'GRID_matched'), names(z))
 if(!length(grid_scores))stop('No recognized GRID columns')
 d <- merge(d, z[, c('eid', grid_scores), with = FALSE], by = 'eid', all.x = TRUE);available <- c(available, grid_scores)
}
for(s in available)set(d, j = s, value = num(d[[s]], s))
# Use reference-projected PCs only for genetic distance; phenotype PCs remain covariates.
pc <- paste0('PC', seq_len(npc)); pca <- ids(read_table(files['pca'])); centers <- read_table(files['centers'])
if(!all(pc %in% names(pca)) || !all(c('POP', pc) %in% names(centers)))stop('Missing reference-projected PCs/centers')
centers <- centers[match(pops, POP)]
if(anyNA(centers$POP) || anyDuplicated(read_table(files['centers'])$POP))stop('Need one center for each EUR/AFR/EAS/SAS')
pm <- as.matrix(pca[, ..pc]); cm <- as.matrix(centers[, ..pc])
if(any(!is.finite(pm)) || any(!is.finite(cm)))stop('Nonfinite projected PCs or centers')
pc_names <- paste0('ancPC', seq_len(npc))
gd <- data.table(eid = pca$eid, proj_PC1 = pm[, 1], proj_PC2 = pm[, 2])
for(j in seq_len(npc))gd[, (pc_names[j]) := pm[, j]]
for(j in seq_along(pops))gd[, (paste0('distance.', pops[j])) := sqrt(rowSums(sweep(pm, 2, cm[j, ], '-') ^ 2))]
gd[, nearest_distance := do.call(pmin, .SD), .SDcols = paste0('distance.', pops)]
geometry <- training_geometry(gd, pm, pc, centers)
gd <- geometry$gd
d <- merge(d, gd, by = 'eid');rm(gd, pca, pm);invisible(gc(verbose = FALSE))
models <- list(COJO = 'pt.TARGET', `PRS-CSx-auto-meta` = 'csx.auto', `PRS-CSx` = base_scores, `DiscoDivas-untuned` = 'disco', `PRS-CSx-fixed-meta` = 'csx.meta')
if(tune){models[['DiscoDivas-tuned']] <- 'disco.cv';available <- c(available, 'disco.cv');d[, disco.cv := 0]}
for(s in base_scores)models[[s]] <- s
if(length(grid_scores)){
 g4 <- paste0('GRID_', c('AFR', 'EAS', 'EUR', 'SAS'))
 if(all(g4 %in% grid_scores)){models[['GRID-tuned']] <- g4;main_methods <- c(main_methods, 'GRID-tuned')}
 for(g in intersect(c('GRID_shared', 'GRID_posterior', 'GRID_matched'), grid_scores))models[[g]] <- g
}
model_map <- rbindlist(lapply(names(models), function(m)data.table(method = m, scores = paste(models[[m]], collapse = ','))))
manifest <- data.table(field = c('trait', 'type', names(files), 'covariates', 'outcome', 'folds', 'seed', 'bootstrap', 'distance_PCs', 'distance_bins_max', 'min_bin_events', 'group_column', 'missing_scores', 'disco_tuned', 'grid_file', 'uncertainty'),
 value = c(Y, type, files, paste(covars, collapse = ','), outcome_definition, nfold, seed, nboot, npc, nbins, min_bin_events, gc, paste(missing, collapse = ','), tune, arg('grid-file', 'not supplied'), 'Paired subject bootstrap conditional on fixed OOF fits; no discovery/fit uncertainty'))
if(length(missing))cat('Unavailable scores: ', paste(missing, collapse = ', '), '\n', sep = '')
cat('Ancestry-matched candidates:\n');print(d[, .N, by = target])
if(isTRUE(arg('check', FALSE))) {
 print(manifest);print(model_map);print(rbindlist(audit))
 cat('CHECK complete. Missing PT can be generated by a normal run.\n');quit(status = 0)
}
quote_name <- function(x)paste0('`', x, '`')
formula_for <- function(cv, ns = 0L, surv = FALSE) {
 terms <- c(quote_name(cv), if(ns)paste0('z', seq_len(ns)))
 as.formula(paste(if(surv)'Surv(time,outcome)' else 'outcome', '~', if(length(terms))paste(terms, collapse = '+') else '1'))
}
fit_model <- function(f, x, kind) {
 m <- switch(kind, ct = lm(f, data = x), dt = glm(f, data = x, family = binomial()), t2e = coxph(f, data = x, ties = 'efron'))
 if(kind == 'dt' && !m$converged)stop('Logistic model did not converge')
 if(any(!is.finite(coef(m))))stop('Singular/nonfinite fit; check covariates and score collinearity')
 m
}
predict_model <- function(m, x, kind)as.numeric(if(kind == 'ct')predict(m, x) else predict(m, x, type = if(kind == 'dt')'response' else 'lp'))
cindex <- function(y, time, p, fold) {
 num <- den <- 0
 for(k in unique(fold)) {
	take <- which(fold == k);z <- data.frame(y = y[take], time = time[take], p = p[take])
	if(nrow(z) < 2 || !any(z$y == 1 & z$time < max(z$time)))next
	cc <- survival::concordancefit(Surv(z$time, z$y), z$p, reverse = TRUE, std.err = FALSE);np <- sum(cc$count[1:3])
	if(np > 0 && is.finite(cc$concordance)){num <- num + np * cc$concordance;den <- den + np}
 }
 if(den > 0)num / den else NA_real_
}
# All methods share bootstrap indices. Calibration/SSE metrics are saved separately.
summarize_predictions <- function(x, pred, base, linear, linear_base, Kpop = NA_real_, B = nboot) {
 nm <- colnames(pred);n <- nrow(x);P <- mean(x$outcome)
 if(type == 'ct') {
	yres <- x$outcome - base;pres <- sweep(linear, 1, base, '-')
	stat <- function(ix) {
	 counts <- tabulate(ix, nbins = n);w <- counts / sum(counts)
	 ym <- sum(w * yres);pm <- as.numeric(crossprod(w, pres))
	 vy <- sum(w * yres ^ 2) - ym ^ 2
	 vp <- as.numeric(crossprod(w, pres ^ 2)) - pm ^ 2
	 cp <- as.numeric(crossprod(w * yres, pres)) - ym * pm
	 value <- cp ^ 2 / (vy * vp);value[vy <= 0 | vp <= 0] <- NA_real_
	 value <- pmin(1, pmax(0, value));names(value) <- nm;value
	}
 } else if(type == 'dt') {
	auc <- function(y, p)if(length(unique(y)) < 2)NA_real_ else as.numeric(pROC::auc(pROC::roc(y, p, levels = 0:1, direction = '<', quiet = TRUE)))
	stat <- function(ix)c(vapply(nm, function(m)auc(x$outcome[ix], pred[ix, m]), numeric(1)), .baseline = auc(x$outcome[ix], base[ix]))
 } else stat <- function(ix)c(vapply(nm, function(m)cindex(x$outcome[ix], x$time[ix], pred[ix, m], x$fold[ix]), numeric(1)),
				.baseline = cindex(x$outcome[ix], x$time[ix], base[ix], x$fold[ix]))
 point <- stat(seq_len(n));boot_names <- names(point)
 boots <- matrix(NA_real_, B, length(boot_names), dimnames = list(NULL, boot_names))
 for(b in seq_len(B)) {
	ix <- if(type == 'ct')sample.int(n, n, replace = TRUE) else unlist(lapply(split(seq_len(n), x$outcome), function(z)sample(z, length(z), replace = TRUE)), use.names = FALSE)
	boots[b, ] <- stat(ix)
	if(B >= 100L && b %% max(1L, B %/% 4L) == 0L){cat('  bootstrap ', b, '/', B, '\n', sep = '');flush.console()}
 }
 ci <- function(v)if(sum(is.finite(v)) >= max(10, .8 * B))quantile(v, c(.025, .975), na.rm = TRUE, names = FALSE) else c(NA_real_, NA_real_)
 intervals <- if(B > 0)t(apply(boots, 2, ci)) else matrix(NA_real_, length(boot_names), 2, dimnames = list(boot_names, NULL))
 res <- data.table(method = nm, estimate = as.numeric(point[nm]), lower95 = intervals[nm, 1], upper95 = intervals[nm, 2], N = n, events = if(type == 'ct')NA_integer_ else sum(x$outcome), K = Kpop, P = if(type == 'ct')NA_real_ else P)
 if(type == 't2e') {
	res[, `:=`(baseline_C = unname(point['.baseline']), baseline_C_lower95 = intervals['.baseline', 1],
						baseline_C_upper95 = intervals['.baseline', 2], delta_C = estimate - unname(point['.baseline']))]
	dc <- if(B > 0)t(vapply(nm, function(m)ci(boots[, m] - boots[, '.baseline']), numeric(2))) else matrix(NA_real_, length(nm), 2)
	res[, `:=`(delta_C_lower95 = dc[, 1], delta_C_upper95 = dc[, 2])]
 }
 if(type == 'dt') {
	res[, `:=`(baseline_AUC = unname(point['.baseline']), delta_AUC = estimate - unname(point['.baseline']))]
	dc <- if(B > 0)t(vapply(nm, function(m)ci(boots[, m] - boots[, '.baseline']), numeric(2))) else matrix(NA_real_, length(nm), 2)
	res[, `:=`(delta_AUC_lower95 = dc[, 1], delta_AUC_upper95 = dc[, 2])]
 }
 res[, metric := switch(type, ct = 'OOF_prediction_R2', dt = 'OOF_AUC', t2e = 'OOF_Harrell_C')]
 list(performance = res, bootstrap = boots[, nm, drop = FALSE])
}
set.seed(seed)


# 🚩 disco
# Phenotype-tuned DiscoDivas within each OUTER training fold.
# The four input anchor models each combine all four raw PRS-CSx scores.
disco_weights <- function(pc, centers, quality = rep(1, 4)) {
 G <- as.matrix(dist(centers))
 if(!is.finite(kappa(G)) || kappa(G) > 1e12)stop('Singular Disco anchor geometry')
 correction <- as.numeric(solve(G, rep(1, nrow(G)))) * quality
 dist <- vapply(seq_len(nrow(centers)), function(j)sqrt(rowSums(sweep(pc, 2, centers[j, ], '-') ^ 2)), numeric(nrow(pc)))
 eps <- max(G, 1) * 1e-12
 w <- sweep(1 / pmax(dist, eps), 2, correction, '*');den <- rowSums(w)
 if(any(abs(den) <= 1e-12 * pmax(rowSums(abs(w)), 1e-300)))stop('Undefined Disco interpolation denominator')
 w <- w / den
 exact <- which(rowSums(dist <= eps) > 0)
 for(i in exact){hit <- which(dist[i, ] <= eps & correction != 0);if(length(hit) == 1){w[i, ] <- 0;w[i, hit] <- 1}}
 w
}
build_disco_folds <- function(d, eligible, pc_names, quality, min_anchor = 100L) {
 answer <- vector('list', nfold)
 for(k in seq_len(nfold)) {
	training <- which(eligible & d$fold != k)
	models <- matrix(NA_real_, nrow(d), 4, dimnames = list(NULL, pops))
	centers <- matrix(NA_real_, 4, length(pc_names), dimnames = list(pops, pc_names))
	anchors <- list()
	for(j in seq_along(pops)) {
	 ix <- training[d$target[training] == pops[j]]
	 if(length(ix) < min_anchor)stop('Too few Disco anchor training samples: ', pops[j], '; use --disco-tune FALSE for the saved-score diagnostic')
	 tr <- copy(d[ix]);sc <- base_scores
	 mu <- vapply(tr[, ..sc], mean, numeric(1));ss <- vapply(tr[, ..sc], sd, numeric(1))
	 if(any(!is.finite(ss) | ss <= 0))stop('Constant anchor PRS: ', pops[j])
	 for(q in seq_along(sc))tr[, (paste0('z', q)) := (get(sc[q]) - mu[q]) / ss[q]]
	 cv <- covars[vapply(tr[, ..covars], function(z)uniqueN(z) > 1, logical(1))]
	 fit <- fit_model(formula_for(cv, 4L, type == 't2e'), tr, type)
	 beta <- coef(fit)[paste0('z', 1:4)]
	 models[, j] <- as.numeric(sweep(sweep(as.matrix(d[, ..sc]), 2, mu, '-'), 2, ss, '/') %*% beta)
	 centers[j, ] <- vapply(tr[, ..pc_names], median, numeric(1))
	 # Balanced TRAINING reference for ancestry residualization and scaling.
	 anchors[[j]] <- ix
	}
	balanced_n <- min(10000L, lengths(anchors))
	harmonize <- unlist(lapply(anchors, function(ix)sample(ix, balanced_n)), use.names = FALSE)
	pc <- as.matrix(d[, ..pc_names]);hx <- cbind(1, pc[harmonize, , drop = FALSE]);allx <- cbind(1, pc)
	hf <- lm.fit(hx, models[harmonize, , drop = FALSE])
	if(any(!is.finite(hf$coefficients)))stop('Singular training PC harmonization')
	residual <- models - allx %*% hf$coefficients
	mu <- colMeans(residual[harmonize, , drop = FALSE]);ss <- apply(residual[harmonize, , drop = FALSE], 2, sd)
	if(any(!is.finite(ss) | ss <= 0))stop('Zero residual PRS variance')
	scaled <- sweep(sweep(residual, 2, mu, '-'), 2, ss, '/')
	weights <- disco_weights(pc, centers, quality)
	answer[[k]] <- rowSums(scaled * weights)
	cat('Disco tuned: fold ', k, '/', nfold, ' DONE; training=', length(training), '\n', sep = '');flush.console()
 }
 answer
}


groups <- c(pops, setdiff(sort(unique(na.omit(d$target))), pops))
d[, `:=`(fold = 0L, eligible = FALSE, row_index = .I)]
for(g in groups) {
 ix <- which(d$target == g);mm <- lapply(models, function(z)sub('TARGET', g, z, fixed = TRUE))
 mm <- mm[vapply(mm, function(z)all(z %in% available), logical(1))]
 req <- unique(c('outcome', if(type == 't2e')'time', covars, unlist(mm)))
 valid <- complete.cases(d[ix, ..req]);for(v in req)if(is.numeric(d[[v]]))valid <- valid & is.finite(d[[v]][ix])
 if(type == 't2e')valid <- valid & d$time[ix] > 0
 good <- ix[valid];d$eligible[good] <- TRUE
 strata <- if(type == 'ct')list(good) else split(good, d$outcome[good])
 for(j in strata)if(length(j))d$fold[j] <- sample(rep(seq_len(nfold), length.out = length(j)))
}
disco_folds <- if(tune)build_disco_folds(d, d$eligible, pc_names, quality, min_anchor) else NULL
all_predictions <- list()
perf <- differences <- distance_results <- skips <- list()
groups <- c(pops, setdiff(sort(unique(na.omit(d$target))), pops))
for(g in groups) {
 x <- copy(d[target == g]);n_start <- nrow(x)
 mm <- lapply(models, function(s)sub('TARGET', g, s, fixed = TRUE))
 ok_models <- vapply(mm, function(s)all(s %in% available), logical(1));mm <- mm[ok_models]
 for(m in names(models)[!ok_models])skips[[length(skips) + 1L]] <- data.table(target = g, method = m, reason = 'Missing score input or no ancestry-matched PT')
 used <- unique(unlist(mm));required <- c('outcome', if(type == 't2e')'time', covars, used)
 valid <- complete.cases(x[, ..required]);for(v in required)if(is.numeric(x[[v]]))valid <- valid & is.finite(x[[v]])
 if(type == 't2e')valid <- valid & x$time > 0
 x <- droplevels(x[valid]);n <- nrow(x)
 audit[[length(audit) + 1L]] <- data.table(stage = c('matched_before_complete_cases', 'paired_complete_cases'), target = g, N = c(n_start, n))
 if(n < minn || !length(mm) || (type != 'ct' && min(table(factor(x$outcome, levels = 0:1))) < nfold * 2)) {
	skips[[length(skips) + 1L]] <- data.table(target = g, method = 'ALL', reason = paste('Insufficient samples or cases/controls:', n));next
 }
 cv <- covars[vapply(x[, ..covars], function(v)uniqueN(v) > 1, logical(1))]
 if(any(x$fold == 0L))stop('Internal fold assignment mismatch')
 cat('START ', g, ': N=', n, '; folds=', nfold, '\n', sep = '');flush.console()
 pred <- lin <- matrix(NA_real_, n, length(mm), dimnames = list(NULL, names(mm)))
 base <- lb <- rep(NA_real_, n)
 for(k in seq_len(nfold)) {
	tr <- copy(x[fold != k]);te <- copy(x[fold == k]);it <- which(x$fold == k)
	if(tune){tr[, disco.cv := disco_folds[[k]][row_index]];te[, disco.cv := disco_folds[[k]][row_index]]}
	f0 <- formula_for(cv, surv = type == 't2e');bm <- fit_model(f0, tr, type);base[it] <- predict_model(bm, te, type)
	lb[it] <- if(type == 'dt')predict_model(fit_model(formula_for(cv), tr, 'ct'), te, 'ct') else base[it]
	for(m in names(mm)) {
	 sc <- mm[[m]];mu <- vapply(sc, function(s)mean(tr[[s]]), numeric(1));ss <- vapply(sc, function(s)sd(tr[[s]]), numeric(1))
	 if(any(!is.finite(ss) | ss <= 0))stop('Constant training score: ', g, ' ', m)
	 for(j in seq_along(sc)){tr[, (paste0('z', j)) := (get(sc[j]) - mu[j]) / ss[j]];te[, (paste0('z', j)) := (get(sc[j]) - mu[j]) / ss[j]]}
	 f <- formula_for(cv, length(sc), type == 't2e');fm <- fit_model(f, tr, type);pred[it, m] <- predict_model(fm, te, type)
	 lmfit <- if(type == 'dt')fit_model(formula_for(cv, length(sc)), tr, 'ct') else fm
	 lin[it, m] <- if(type == 'dt')predict_model(lmfit, te, 'ct') else pred[it, m]
	 if(m == 'PRS-CSx') {
		fold_coefficients[[length(fold_coefficients) + 1L]] <- data.table(target = g, fold = k, score = sc,
			coefficient = as.numeric(coef(fm)[paste0('z', seq_along(sc))]), training_mean = mu, training_sd = ss,
			raw_weight = as.numeric(coef(fm)[paste0('z', seq_along(sc))]) / ss)
		if(all(cov_columns %in% names(te)))individual_results[[length(individual_results) + 1L]] <- individual_posterior(te, tr, fm, bm, ss, g, k)
	 }
	}
	cat('  fold ', k, '/', nfold, ' DONE\n', sep = '');flush.console()
 }
 if(toupper(arg('write-predictions', 'FALSE')) == 'TRUE')all_predictions[[g]] <- cbind(x[, .(eid, target, fold, outcome, time)], data.table(baseline = base), as.data.table(pred))
 if(any(!is.finite(pred)) || any(!is.finite(lin)))stop('Nonfinite held-out prediction: ', g)
 kp <- if(type == 'dt')prevalence[target == g]$K else NA_real_
 if(type == 'dt' && (length(kp) != 1 || !is.finite(kp) || kp <= 0 || kp >= 1)){skips[[length(skips) + 1L]] <- data.table(target = g, method = 'ALL', reason = 'Missing valid K');next}
 sm <- summarize_predictions(x, pred, base, lin, lb, kp);pp <- sm$performance;pp[, target := g]
 # Additional metrics retain the same cohort and fold assignments.
 for(m in names(mm)) {
	if(type == 'ct') {
	 den <- sum((x$outcome - mean(x$outcome)) ^ 2)
	 sse0 <- sum((x$outcome - base) ^ 2);sse1 <- sum((x$outcome - pred[, m]) ^ 2)
	 pp[method == m, `:=`(baseline_R2 = 1 - sse0 / den, full_R2 = 1 - sse1 / den,
							delta_R2 = (sse0 - sse1) / den, baseline_SSE = sse0, full_SSE = sse1, RMSE = sqrt(sse1 / n),
							SSE_partial_R2 = 1 - sse1 / sse0, prediction_r = cor(x$outcome - base, pred[, m] - base))]
	} else if(type == 'dt') {
	 auc <- function(p)as.numeric(pROC::auc(pROC::roc(x$outcome, p, levels = 0:1, direction = '<', quiet = TRUE)))
	 pp[method == m, `:=`(AUC = auc(pred[, m]), baseline_AUC = auc(base), delta_AUC = auc(pred[, m]) - auc(base), Brier = mean((x$outcome - pred[, m]) ^ 2), baseline_Brier = mean((x$outcome - base) ^ 2))]
	}
 }
 perf[[length(perf) + 1L]] <- pp
 disco_method <- if(tune)'DiscoDivas-tuned' else 'DiscoDivas-untuned'
 if(nboot > 0 && all(c('PRS-CSx', disco_method) %in% names(mm))) {
	bs <- sm$bootstrap;dd <- bs[, disco_method] - bs[, 'PRS-CSx'];ci <- quantile(dd, c(.025, .975), na.rm = TRUE)
	reference <- pp[method == 'PRS-CSx', estimate];val <- pp[method == disco_method, estimate] - reference
	rel <- if(reference > 0 && mean(bs[, 'PRS-CSx'] > 0, na.rm = TRUE) >= .975)100 * dd / bs[, 'PRS-CSx'] else rep(NA_real_, nboot)
	ri <- if(any(is.finite(rel)))quantile(rel, c(.025, .975), na.rm = TRUE) else c(NA_real_, NA_real_)
	differences[[length(differences) + 1L]] <- data.table(target = g, metric = pp$metric[1], difference = val, lower95 = ci[1], upper95 = ci[2], relative_percent = if(reference > 0)100 * val / reference else NA_real_, relative_lower95 = ri[1], relative_upper95 = ri[2], N = n)
 }
 # PRS-CSx only. Quantile boundaries depend on distance; for binary/survival
 # outcomes reduce the bin count until every bin meets the event-count guard.
 # Never select bins for high performance or refit a model within a bin.
 # Distance resampling cannot change the next ancestry's main bootstrap stream.
 rng_before_distance <- .Random.seed
 if('PRS-CSx' %in% names(mm)) {
	axis <- 'distance.analysis'; selected <- 'PRS-CSx'
	bins <- NULL; max_bins <- min(nbins, n %/% minn)
	if(type != 'ct')max_bins <- min(max_bins, sum(x$outcome == 1) %/% min_bin_events, sum(x$outcome == 0) %/% min_bin_events)
	if(max_bins >= 2L)for(q in seq.int(max_bins, 2L)) {
	 breaks <- unique(quantile(x[[axis]], seq(0, 1, length.out = q + 1L), na.rm = TRUE))
	 if(length(breaks) < 3L)next
	 candidate <- cut(x[[axis]], breaks, include.lowest = TRUE, labels = FALSE)
	 counts <- x[, .(N = .N, events = sum(outcome == 1), non_events = sum(outcome == 0)), by = .(bin = candidate)]
	 if(all(counts$N >= minn) &&
			(type == 'ct' || all(counts$events >= min_bin_events & counts$non_events >= min_bin_events))) {bins <- candidate;break}
	}
	if(!is.null(bins))for(bin in sort(unique(bins))) {
	 ix <- which(bins == bin)
	 zz <- summarize_predictions(x[ix], pred[ix, selected, drop = FALSE], base[ix], lin[ix, selected, drop = FALSE], lb[ix], kp, B = min(nboot, 100L))$performance
	 zz[, `:=`(target = g, distance_axis = axis, bin = bin, bins_in_group = uniqueN(bins), distance = median(x[[axis]][ix]),
						distance_min = min(x[[axis]][ix]), distance_max = max(x[[axis]][ix]))]
	 distance_results[[length(distance_results) + 1L]] <- zz
	}
 }
 .Random.seed <- rng_before_distance
 cat('Evaluated ', g, ': N=', n, '; methods=', paste(names(mm), collapse = ', '), '\n', sep = '');flush.console()
}
if(!length(perf))stop('No evaluable ancestry groups')
if(length(all_predictions))write_tsv(rbindlist(all_predictions, fill = TRUE), 'predictions.tsv.gz')
individual <- rbindlist(individual_results, fill = TRUE)
if(nrow(individual))write_tsv(individual, 'individual_posterior.tsv.gz')
write_tsv(rbindlist(fold_coefficients), 'fold_coefficients.tsv')
write_tsv(geometry$centers, 'distance_centers.tsv')
manifest <- rbind(manifest, data.table(field = c('posterior_mode', 'distance_status', 'distance_definition', 'genetic_variance_file'),
 value = c(arg('posterior-mode', 'required'), geometry$status, geometry$description, arg('genetic-variance-file', file.path(score_dir, 'genetic_variance.tsv')))))
performance <- rbindlist(perf, fill = TRUE);write_tsv(performance, 'performance.tsv')
write_tsv(rbindlist(audit, fill = TRUE), 'cohort.tsv')
comparison <- rbindlist(differences)
distperf <- rbindlist(distance_results)
if(!ncol(distperf))distperf <- data.table(method = character(), target = character(), estimate = numeric(), lower95 = numeric(), upper95 = numeric(), N = integer(), distance = numeric())
write_tsv(distperf, 'distance_performance.tsv')


# 🚩 plots
# Report layer. Scores and performance estimates are computed by Yeval.R.
colours <- c('GRID-tuned' = '#278C78', GRID_shared = '#76A88C', GRID_posterior = '#5D9B91', GRID_matched = '#86A9AA',
 COJO = '#7C8798', 'PRS-CSx-auto-meta' = '#E1A63B', 'PRS-CSx' = '#3275B4', 'DiscoDivas-untuned' = '#CD8B95',
 'DiscoDivas-tuned' = '#D45260', 'PRS-CSx-fixed-meta' = '#7755A2', 'csx.AFR' = '#C6813B', 'csx.EAS' = '#60A67A', 'csx.EUR' = '#6593CC', 'csx.SAS' = '#9E83BD')
ancestry_colours <- c(EUR = '#437CB3', AFR = '#D59421', EAS = '#269B78', SAS = '#AE6BA5', OTH = '#8B9098', UNASSIGNED = '#BABEC5')
extra <- setdiff(groups, names(ancestry_colours))
if(length(extra))ancestry_colours <- c(ancestry_colours, setNames(scales::hue_pal()(length(extra)), extra))
metric_label <- switch(type, ct = 'Prediction R²', dt = 'AUC (covariates + PRS)', t2e = 'Harrell C-index (covariates + PRS)')
fmt_n <- function(x)format(x, big.mark = ',', scientific = FALSE, trim = TRUE)
fmt <- function(x, digits = 4)ifelse(is.finite(x), formatC(x, digits = digits, format = 'f'), 'NA')
ci_text <- paste0(nfold, '-fold out-of-fold predictions; ', if(nboot > 0)paste0(nboot, ' paired subject bootstrap resamples') else 'intervals disabled')
counts <- d[eligible == TRUE, .(N = .N, events = if(type == 'ct')NA_integer_ else sum(outcome)), by = target]
target_label <- setNames(vapply(pops, function(g){
 z <- performance[target == g][1]
 if(!nrow(z) || is.na(z$N))return(paste0(g, '\nUnavailable'))
 paste0(g, '\nN = ', fmt_n(z$N), if(type != 'ct')paste0('\n', if(type == 't2e')'Events' else 'Cases', ' = ', fmt_n(z$events)))
}, character(1)), pops)
axis_methods <- function(z){
 labels <- c(COJO = 'COJO', 'PRS-CSx-auto-meta' = 'Auto\nmeta', 'PRS-CSx-fixed-meta' = 'Fixed\nmeta', 'PRS-CSx' = 'Four-score\nfit',
				'DiscoDivas-tuned' = 'Disco\ntuned', 'DiscoDivas-untuned' = 'Disco\nuntuned', 'GRID-tuned' = 'GRID\ntuned')
 ans <- unname(labels[z]);ans[is.na(ans)] <- z[is.na(ans)];ans
}
plot_theme <- theme_classic(base_size = 12) + theme(legend.position = 'bottom', legend.title = element_text(face = 'bold'),
 strip.background = element_rect(fill = '#F2F4F7', colour = NA), strip.text = element_text(face = 'bold', margin = margin(8, 5, 8, 5)),
 plot.title = element_text(face = 'bold', size = 17), plot.subtitle = element_text(colour = '#526071', size = 11, lineheight = 1.15),
 plot.caption = element_text(hjust = 0, colour = '#626A76', size = 10), panel.spacing = grid::unit(1.2, 'lines'), plot.margin = margin(12, 14, 10, 12))
primary_caption <- switch(type,
 ct = 'Prediction R² = squared correlation of covariate-residualized phenotype and PRS prediction, evaluated on held-out people.',
 dt = 'AUC evaluates held-out case-control discrimination; covariates are included. ΔAUC and Brier are saved separately.',
 t2e = 'C-index measures ranking of event times, including covariates. Dashed lines show covariates-only C; 0.5 is chance ranking.')
make_comparison <- function(methods, title, subtitle, caption) {
 z <- merge(CJ(target = pops, method = methods, unique = TRUE), performance[target %in% pops], by = c('target', 'method'), all.x = TRUE)
 z[, method := factor(method, levels = methods)];z[, target := factor(target, levels = pops)]
 p <- ggplot(z, aes(method, estimate, fill = method)) +
	geom_hline(yintercept = if(type != 'ct').5 else 0, colour = '#AEB7C2', linewidth = .35) +
	geom_col(width = .66, na.rm = TRUE, alpha = .95) +
	geom_errorbar(aes(ymin = lower95, ymax = upper95), width = .16, linewidth = .5, na.rm = TRUE) +
	geom_text(data = z[is.na(estimate)], aes(x = method, y = if(type == 't2e').51 else 0, label = 'Unavailable'),
				inherit.aes = FALSE, angle = 90, hjust = 0, size = 2.8, colour = '#777777') +
	facet_grid(. ~ target, labeller = labeller(target = target_label), drop = FALSE) +
	scale_fill_manual(name = 'Predictor', values = colours, limits = methods, breaks = methods, drop = FALSE) +
	scale_x_discrete(labels = axis_methods, drop = FALSE) + scale_y_continuous(expand = expansion(mult = c(.04, .08))) +
	guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + plot_theme +
	theme(axis.text.x = element_text(size = 10, lineheight = .95), legend.text = element_text(size = 11)) +
	labs(title = title, subtitle = subtitle, x = NULL, y = metric_label, caption = caption)
 if(type == 't2e') {
	b <- unique(z[is.finite(baseline_C), .(target, baseline_C)])
	p <- p + geom_hline(data = b, aes(yintercept = baseline_C), linetype = 2, colour = '#303842', linewidth = .6) +
	 coord_cartesian(ylim = c(min(.5, z$lower95, z$estimate, z$baseline_C, na.rm = TRUE),
						min(1, max(.82, z$upper95, z$estimate, z$baseline_C, na.rm = TRUE) + .025)))
 }
 p
}
p1 <- make_comparison(main_methods, paste(Y, '| Prediction by target ancestry'),
 'Overall benchmark: COJO, auto-meta, fitted PRS-CSx and DiscoDivas within each ancestry.\nThe same participants are used for every available predictor in a panel.', primary_caption)
p2 <- make_comparison(c('PRS-CSx-auto-meta', 'PRS-CSx-fixed-meta', 'PRS-CSx', if(tune)'DiscoDivas-tuned' else 'DiscoDivas-untuned'),
 paste(Y, '| Combined-score comparison'),
 'How CSx scores are combined: posterior meta-analysis (auto/fixed phi), four-score regression (PRS-CSx),\nor genetic-distance interpolation (DiscoDivas). COJO is omitted; shared bars repeat the first figure.',
 if(type == 't2e')primary_caption else 'auto/fixed meta combine SNP effects. Four-score regression learns target-ancestry weights in training folds.')

# Four distinct questions: ancestry, empirical prediction, geometry, posterior.
set.seed(seed + 1L)
landscape <- d[eligible == TRUE, .SD[sample.int(.N, min(.N, 2000L))], by = target]
shown_groups <- groups[groups %in% unique(landscape$target)]
category <- merge(data.table(target = shown_groups), performance[method == 'PRS-CSx'], by = 'target', all.x = TRUE)
category[, target := factor(target, levels = shown_groups)]
small_theme <- plot_theme + theme(plot.title = element_text(size = 13, face = 'bold'), plot.subtitle = element_text(size = 10),
 legend.text = element_text(size = 10), axis.title = element_text(size = 11), plot.tag = element_text(face = 'bold', size = 18))
ancestry_scale <- function()scale_colour_manual(name = 'Ancestry', values = ancestry_colours, limits = shown_groups, drop = FALSE)
pa <- ggplot(landscape, aes(proj_PC1, proj_PC2, colour = target)) + geom_point(size = .65, alpha = .55) +
 ancestry_scale() + small_theme + coord_equal() +
 guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 1)) +
 labs(title = 'Ancestry groups', subtitle = paste0('Original labels: ', gc), x = 'Projected PC1', y = 'Projected PC2')

category[, `:=`(value = estimate, lo = lower95, hi = upper95)]
if(type == 't2e')category[, `:=`(value = delta_C, lo = delta_C_lower95, hi = delta_C_upper95)]
if(type == 'dt')category[, `:=`(value = delta_AUC, lo = delta_AUC_lower95, hi = delta_AUC_upper95)]
category_labels <- setNames(vapply(shown_groups, function(g){
 z <- counts[target == g]
 paste0(g, '\nN=', fmt_n(z$N), if(type != 'ct')paste0('\nEvents=', fmt_n(z$events)))
}, character(1)), shown_groups)
gain_label <- switch(type, ct = 'Prediction R²', dt = 'ΔAUC from PRS', t2e = 'ΔC from PRS')
pb <- ggplot(category, aes(target, value, colour = target)) +
 geom_hline(yintercept = 0, colour = '#CDD3DA', linewidth = .4) +
 geom_errorbar(aes(ymin = lo, ymax = hi), width = .12, linewidth = .6, na.rm = TRUE) + geom_point(size = 3.2, na.rm = TRUE) +
 geom_text(aes(label = ifelse(is.finite(value), fmt(value, 3), 'Unavailable')), vjust =  - 1.3, size = 3, na.rm = TRUE) +
 ancestry_scale() + scale_x_discrete(labels = category_labels, drop = FALSE) +
 scale_y_continuous(expand = expansion(mult = c(.12, .24))) + small_theme + theme(legend.position = 'none') +
 labs(title = 'PRS-CSx prediction by group', subtitle = if(type == 'ct')'Covariate-adjusted prediction in held-out participants' else 'Improvement beyond covariates in held-out participants', x = NULL, y = gain_label)

# c uses a small number of real people to explain geometry, not another cloud
# recoloured by distance. All four source centres are shown, with no EUR default.
examples <- landscape[, .SD[which.min(abs(distance.analysis - median(distance.analysis)))], by = target]
examples[, point_id := paste0('i', seq_len(.N))]
gc_plot <- copy(geometry$centers)
gc_plot[, label := paste0(POP, if(geometry$status == 'discovery')' training' else ' reference')]
segments <- rbindlist(lapply(seq_len(nrow(examples)), function(i){
 data.table(x = gc_plot$PC1, y = gc_plot$PC2, xend = examples$proj_PC1[i], yend = examples$proj_PC2[i], weight = gc_plot$mixture_weight)
}))
pcp <- ggplot() +
 geom_segment(data = segments, aes(x = x, y = y, xend = xend, yend = yend, linewidth = weight), colour = '#AAB2BD', alpha = .65,
			arrow = grid::arrow(length = grid::unit(.07, 'inches'), type = 'closed')) +
 scale_linewidth_continuous(range = c(.25, 1.1), guide = 'none') +
 geom_point(data = gc_plot, aes(PC1, PC2), shape = 17, size = 3.3, colour = '#273646') +
 geom_text(data = gc_plot, aes(PC1, PC2, label = label), nudge_y = .04 * diff(range(landscape$proj_PC2)), size = 3, check_overlap = TRUE) +
 geom_point(data = examples, aes(proj_PC1, proj_PC2, colour = target), size = 3.2) +
 geom_text(data = examples, aes(proj_PC1, proj_PC2, label = point_id, colour = target), nudge_y =  - .04 * diff(range(landscape$proj_PC2)), size = 3.4, show.legend = FALSE) +
 ancestry_scale() + scale_x_continuous(expand = expansion(mult = .18)) + scale_y_continuous(expand = expansion(mult = .18)) + small_theme + coord_equal() + theme(legend.position = 'none') +
 labs(title = if(geometry$status == 'discovery')'Distance to GWAS training groups' else 'Reference geometry (exploratory)',
			subtitle = paste0('Triangles: centres; dots: example people; distances use ', npc, ' PCs'), x = 'Projected PC1', y = 'Projected PC2')

individual_metric <- arg('individual-metric', 'sd')
if(!individual_metric %in% c('auto', 'reliability', 'sd'))stop('--individual-metric must be auto, reliability or sd')
if(individual_metric == 'auto')individual_metric <- if(nrow(individual) && any(is.finite(individual$individual_R2)))'reliability' else 'sd'
if(individual_metric == 'reliability' && (!nrow(individual) || !any(is.finite(individual$individual_R2))))stop('Individual reliability requires a valid --genetic-variance-file on the correct outcome scale')
if(nrow(individual)){
 individual[, plot_value := if(individual_metric == 'reliability')individual_R2 else posterior_sd]
 plotted <- individual[is.finite(plot_value)]
 plot_limit <- intarg('individual-max-points', 5000, 100)
 set.seed(seed + 2L)
 dots <- plotted[, .SD[sample.int(.N, min(.N, plot_limit))], by = target]
 ylab <- if(individual_metric == 'reliability')'Individual model-based R²' else switch(type,
		 ct = 'Posterior SD of PRS prediction', dt = 'Posterior SD of PRS log-odds', t2e = 'Posterior SD of PRS log-hazard')
 pd <- ggplot(dots, aes(distance.analysis, plot_value, colour = target)) + geom_point(size = .55, alpha = .3) +
	ancestry_scale() + scale_x_continuous(trans = 'sqrt', labels = scales::label_number()) +
	small_theme + theme(legend.position = 'none') +
	labs(title = if(individual_metric == 'reliability')'Individual PRS-CSx reliability' else 'Individual PRS-CSx uncertainty',
			 subtitle = 'Each dot is one person; no bins or imposed decay curve', x = geometry$label, y = ylab)
 if(individual_metric == 'reliability' && min(dots$plot_value) < 0)pd <- pd + geom_hline(yintercept = 0, colour = '#9DA8B4', linewidth = .4)
 missing_prior <- nrow(individual) - nrow(plotted)
 posterior_caption <- paste0('d: ', fmt_n(nrow(plotted)), ' individual estimates; up to ', fmt_n(plot_limit), ' dots per group displayed. ',
	if(individual_metric == 'reliability')paste0('Model-based reliability is distinct from empirical prediction R² in b; ', missing_prior, ' people lack a variance scale.') else 'Posterior SD is shown; lower values mean less uncertainty. This is not a prediction R².')
}else{
 pd <- ggplot() + theme_void() + labs(title = 'Individual posterior analysis disabled', subtitle = 'Run 1csx.sh --posterior TRUE, then Yeval with --posterior-mode required') + small_theme
 posterior_caption <- 'Individual posterior estimates disabled explicitly.'
}
distance_caption <- paste0('a: original ancestry labels. c: lines illustrate PC1/PC2 only; numerical distances use all ', npc, ' PCs.\n',
 if(geometry$status == 'discovery')'Multi-training distance: sqrt(sum(N_group / N_total × distance_to_group²)).' else 'Reference proxy: equal weights over AFR/EAS/EUR/SAS 1KG centres; these are not the GWAS training centres.',
 '\n', posterior_caption)
p_distance <- (pa + pb) / (pcp + pd) + plot_annotation(title = paste(Y, '| Prediction along genetic distance'),
 subtitle = 'PRS-CSx: observed prediction by ancestry and posterior information for each person.', tag_levels = 'a', caption = distance_caption,
 theme = theme(plot.title = element_text(face = 'bold', size = 18), plot.subtitle = element_text(size = 12, colour = '#526071'), plot.caption = element_text(hjust = 0, size = 9.5), plot.margin = margin(12, 14, 12, 12)))

# The empirical distance bins remain a separate validation plot, not fake
# individual observations. Survival panels emphasize incremental discrimination.
binplot <- copy(distperf)
if(nrow(binplot)){
 binplot[, `:=`(value = estimate, lo = lower95, hi = upper95)]
 if(type == 't2e')binplot[, `:=`(value = delta_C, lo = delta_C_lower95, hi = delta_C_upper95)]
 if(type == 'dt')binplot[, `:=`(value = delta_AUC, lo = delta_AUC_lower95, hi = delta_AUC_upper95)]
 p_bins <- ggplot(binplot, aes(distance, value, colour = target)) + geom_hline(yintercept = 0, colour = '#ADB5C0') +
	geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, alpha = .6) + geom_point(size = 2) +
	ancestry_scale() + facet_wrap( ~ target, scales = 'free_x', nrow = 1) + scale_x_continuous(trans = 'sqrt') + plot_theme +
	labs(title = paste(Y, '| Empirical prediction across distance bins'), subtitle = 'Each point summarizes held-out participants in a distance bin.', x = geometry$label, y = gain_label,
			 caption = 'These group estimates validate prediction on observed outcomes; they are separate from individual posterior reliability or uncertainty.')
}else p_bins <- ggplot() + theme_void() + labs(title = 'Insufficient samples/events for empirical distance bins')
if(nrow(comparison)){
 comparison[, target := factor(target, levels = groups)]
 p_paired <- ggplot(comparison, aes(target, difference)) + geom_hline(yintercept = 0, linetype = 2, colour = '#777777') +
	geom_errorbar(aes(ymin = lower95, ymax = upper95), width = .12, colour = colours[disco_method], linewidth = .7) +
	geom_point(size = 3.5, colour = colours[disco_method]) + plot_theme +
	labs(title = paste(Y, '|', disco_method, 'versus PRS-CSx'), subtitle = 'Positive differences favour DiscoDivas; negative differences favour PRS-CSx.',
			 x = 'Target ancestry', y = paste0('Difference in ', switch(type, ct = 'prediction R²', dt = 'AUC', t2e = 'C-index')))
}else p_paired <- ggplot() + theme_void() + labs(title = 'DiscoDivas versus PRS-CSx', subtitle = 'Comparison unavailable or intervals disabled')


# 🚩 Named PNG figures
figures <- list(comparison = p1, combined_scores = p2, distance_performance = p_distance, paired_improvement = p_paired, distance_bins = p_bins)
for(nm in names(figures)){
	ggsave(file.path(out, paste0(nm, '.png')), figures[[nm]], width = 14,
		height = if(nm == 'distance_performance')11 else 6.8, dpi = 170, bg = 'white')
}

methods_text <- c(paste0('# ', Y, ' PRS evaluation'), '', paste('Outcome:', outcome_definition), paste('Covariates:', paste(covars, collapse = ', ')),
 paste('Evaluation:', ci_text), '',
 '## Prediction metrics',
 '- Continuous: Prediction R² = cor(y - baseline_OOF, full_OOF - baseline_OOF)^2. This is covariate-adjusted squared prediction correlation. All residualization, score standardization and combination fitting use training folds.',
 '- OLS full_OOF - baseline_OOF equals the score-weighted, training-covariate-residualized score. No regression is fitted within a held-out distance bin.',
 '- Squared correlation does not assess calibration or direction. prediction_r, RMSE, full_R2, baseline_R2, delta_R2 and the former SSE_partial_R2 remain in performance.tsv.',
 '- Binary: main metric is logistic AUC; baseline_AUC, delta_AUC, paired delta intervals and Brier are also saved. No unvalidated liability conversion is applied.',
 '- Survival: main benchmark is covariates + PRS Harrell C. Comparisons use within-fold comparable pairs. The ancestry and empirical distance panels show delta_C to isolate PRS increment.',
 '- These metrics follow common PRS reporting conventions but are not a numerical reproduction of any publication, because training data, covariates, splits and outcomes differ.', '',
 '## Four panels',
 '- a: original ancestry labels in PC space. b: empirical group prediction (R² or delta discrimination). c: separate geometric illustration with four source centres and representative real people. d: one individual per posterior point.',
 '- The first bar chart is the overall method benchmark including COJO. The second isolates combined-score strategies, omits COJO and adds fixed-meta; repeated bars are identical.',
 '- The DiscoDivas difference chart follows the four-panel figure. Empirical bins are a separate validation figure after that.', '',
 '## Posterior model and interpretation',
 '- All four population beta draws share a retained iteration within a chromosome. Individual score covariance includes SNP LD and all cross-population covariance terms.',
 '- Scores are centred at harmonized discovery EAF. Missing genotypes are mean-imputed at the same EAF. The posterior means replace the four old uncentred scores for fitting, so means, scales and covariance refer to one predictor.',
 '- Chromosome means and covariance matrices are summed under PRS-CSx chromosome independence. Independent chromosomes are not artificially coupled by matching iteration numbers.',
 '- With training-fold weights w = regression_coefficient / training_score_SD, individual prediction variance is w^T Sigma_i w. Both diagonal-only and cross-population contributions are saved.',
 '- If an external genetic variance Vg on the correct scale is supplied, model-based individual R² = 1 - w^T Sigma_i w / Vg. For a continuous phenotype, supplied residual SNP h² is multiplied by training-fold covariate-residual phenotype variance.',
 '- This stacked-CSx reliability is an exploratory extension conditional on fitted combination/covariate weights and the supplied Vg, not an established equivalence to LDpred2 reliability. It requires model calibration and compatible priors/scales. It is not empirical phenotype prediction R².',
 '- Negative model-based reliability is retained as a diagnostic of variance scaling, uncertainty or model mismatch. No clipping, artificial dots or forced decay is applied.',
 '- If genetic variance is absent, d reports posterior SD explicitly; it is uncertainty, not accuracy. For survival this is SD of the PRS log-hazard contribution; for binary outcomes SD of log-odds. Liability h² cannot be used in their place.',
 '- Bootstrap intervals condition on fixed OOF fits. Posterior covariance conditions on GWAS summary statistics, LD reference and fitted combination weights; neither accounts for all sources of model misspecification.', '',
 '## Genetic distance', geometry$description,
 '- Discovery centres must be means in the same PCA coordinate system, with source and pca_space identifiers. Multiple source groups use N_GWAS-weighted RMS distance, retaining distances to every group in the individual table.',
 '- A mixture centroid alone can conceal large distance from every actual training group. RMS distance avoids that cancellation, but is a descriptive extension, not a proven predictor of multi-ancestry accuracy.',
 '- Without discovery centres, a clearly labelled equal-weight four-reference proxy is used. It is not the former EUR-only distance and is never called a GWAS training distance.',
 paste0('- At most ', nbins, ' empirical quantile bins per ancestry; minimum ', minn, ' people, and for dt/t2e ', min_bin_events, ' events/cases plus non-events/controls. Sparse bins are not used to infer precise trends.'), '',
 '## References',
 '- PRS-CS: https://doi.org/10.1038/s41467-019-09718-5',
 '- PRS-CSx: https://doi.org/10.1038/s41588-022-01054-7',
 '- Individual accuracy: https://doi.org/10.1038/s41586-023-06079-4',
 '- Source method: https://github.com/yidingdd/individual-pgs-accuracy/tree/main/pgs-accuracy')
writeLines(methods_text, file.path(out, 'methods.md'))
escape <- function(x){x <- gsub('&', '&amp;', as.character(x), fixed = TRUE);x <- gsub('<', '&lt;', x, fixed = TRUE);gsub('>', '&gt;', x, fixed = TRUE)}
html_table <- function(z)paste0('<div class="table-wrap"><table><tr>', paste0('<th>', escape(names(z)), '</th>', collapse = ''), '</tr>',
 paste(apply(as.data.frame(z), 1, function(r)paste0('<tr>', paste0('<td>', escape(r), '</td>', collapse = ''), '</tr>')), collapse = ''), '</table></div>')
csx_summary <- performance[method == 'PRS-CSx']
if(type == 'ct'){
 overview <- csx_summary[, .(Ancestry = target, N = fmt_n(N), `Prediction R²` = fmt(estimate), `Prediction r` = fmt(prediction_r), `SSE-based partial R²` = fmt(SSE_partial_R2), RMSE = fmt(RMSE, 3))]
 metric_explanation <- '<p><strong>Prediction R²</strong> is the squared correlation between covariate-residualized observed and predicted phenotypes in held-out participants. It differs from the previous SSE-based partial R², which remains available for comparison.</p>'
}else if(type == 't2e'){
 overview <- csx_summary[, .(Ancestry = target, N = fmt_n(N), Events = fmt_n(events), `Full C` = fmt(estimate), `Covariates C` = fmt(baseline_C), `ΔC from PRS` = fmt(delta_C), `ΔC 95% CI` = paste0(fmt(delta_C_lower95), '–', fmt(delta_C_upper95)))]
 metric_explanation <- '<p><strong>C-index</strong> measures discrimination of event times. The main comparison bars show covariates + PRS; panel b and empirical distance bins show <strong>ΔC</strong> from adding PRS. Posterior SD in d measures uncertainty of the individual PRS log-hazard contribution; there is no individual C-index.</p>'
}else{
 overview <- csx_summary[, .(Ancestry = target, N = fmt_n(N), Cases = fmt_n(events), AUC = fmt(estimate), `Covariates AUC` = fmt(baseline_AUC), `ΔAUC` = fmt(delta_AUC), Brier = fmt(Brier))]
 metric_explanation <- '<p><strong>AUC</strong> measures held-out case-control discrimination. Panel b and empirical distance bins show ΔAUC from adding PRS. Posterior SD measures uncertainty of individual PRS log-odds.</p>'
}
shown <- copy(performance);for(nm in names(shown))if(is.numeric(shown[[nm]]))set(shown, j = nm, value = signif(shown[[nm]], 5))
writeLines(paste0('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>', escape(Y), ' PRS prediction</title>',
 '<style>body{font:16px system-ui;max-width:1440px;margin:32px auto;padding:0 24px;color:#263446;line-height:1.55}img{width:100%;height:auto;margin:18px 0}.table-wrap{overflow-x:auto}table{border-collapse:collapse;font-size:13px;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:right;white-space:nowrap}th{background:#f2f4f7}td:first-child,th:first-child{text-align:left}.metric{background:#f2f6fa;padding:16px 22px;border-radius:8px}pre{white-space:pre-wrap;font:14px system-ui}a{color:#245f9b}details{margin:18px 0}</style>',
 '<h1>', escape(Y), ' | PRS prediction comparison</h1><p><strong>Outcome:</strong> ', escape(outcome_definition), ' · <strong>Covariates:</strong> ', escape(paste(covars, collapse = ', ')), ' · ', nfold, ' held-out folds.</p>',
 '<div class="metric">', metric_explanation, '<p>', escape(geometry$description), '</p></div>',
 '<p><a href="performance.tsv">Performance</a> · <a href="distance_performance.tsv">Empirical bins</a> · <a href="methods.md">Methods</a>',
 if(nrow(individual))' · <a href="individual_posterior.tsv.gz">Individual posterior estimates</a>' else '', '</p>',
 '<h2>PRS-CSx overview</h2>', html_table(overview),
 paste0('<img src="', names(figures), '.png" alt="', names(figures), '">', collapse = ''),
 '<h2>All predictors</h2>', html_table(shown),
 if(length(skips))paste0('<details><summary>Skipped evaluations</summary>', html_table(rbindlist(skips)), '</details>') else '',
 '<details><summary>Evaluation inputs</summary>', html_table(manifest), '</details>',
 '<details><summary>Score definitions</summary>', html_table(model_map), '</details>',
 '<details><summary>Methods and references</summary><pre>', escape(paste(methods_text, collapse = '\n')), '</pre></details></html>'), file.path(out, 'report.html'))


cat('DONE: ', file.path(out, 'report.html'), '\n', sep = '')
