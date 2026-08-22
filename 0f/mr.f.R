#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Packages, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IV.filter <- function(dat, CHR_col="CHR", POS_col="POS", AF_col="EAF", BETA_col="BETA", N_col="N") {
	if(!all(c(CHR_col, POS_col, AF_col, BETA_col, N_col) %in% names(dat))) return(dat)
	filter1 <- dat[[CHR_col]] %in% c(6, "6") & dat[[POS_col]] >= 28510120 & dat[[POS_col]] <= 33480577
	dat$R2 <- 2 * dat[[AF_col]] * (1 - dat[[AF_col]]) * dat[[BETA_col]]^2
	dat$F_stat <- dat$R2 * (dat[[N_col]] - 2) / (1 - dat$R2)
	dat[!(filter1 | is.na(dat$F_stat) | dat$F_stat < 10), ]
}

rb <- function(x) signif(as.numeric(x), 4)
rp <- function(x) signif(as.numeric(x), 3)

read_gwas_mr <- function(file, N_default = 100000) {
	dat <- data.table::fread(file, header = TRUE) %>% as.data.frame()
	names(dat) <- stringi::stri_replace_all_regex(toupper(names(dat)), pattern = toupper(pattern), replacement = replacement, vectorize_all = FALSE)
	if("CHR" %in% names(dat)) dat$CHR <- ifelse(dat$CHR %in% c(23, "23"), "X", as.character(dat$CHR))
	if(!"N" %in% names(dat)) dat$N <- N_default
	dat
}

read_iv_file <- function(file) {
	dat <- data.table::fread(file, header = TRUE) %>% as.data.frame()
	names(dat) <- toupper(names(dat)); if(!"SNP" %in% names(dat)) names(dat)[1] <- "SNP"
	dat %>% dplyr::select(SNP) %>% dplyr::distinct()
}

get_iv_snps <- function(X, dir.X, dat.X.raw, IV.file = "NO", X.use.top = FALSE) {
	if(identical(IV.file, "TOP") || X.use.top) return(dat.X.raw %>% dplyr::select(SNP) %>% dplyr::distinct())
	if(!identical(IV.file, "NO")) return(read_iv_file(IV.file))
	f <- file.path(dir.X, X, paste0(X, ".clump.snp"))
	if(file.exists(f)) return(read_iv_file(f))
	dat.X.raw %>% dplyr::filter(P <= 5e-08) %>% dplyr::select(SNP) %>% dplyr::distinct()
}

extract_snps_from_gz <- function(gz, snps, outfile, N_default = 100000) {
	if(file.exists(outfile) && file.size(outfile) > 0) return(read_gwas_mr(outfile, N_default))
	tmp <- tempfile(); writeLines(unique(snps), tmp); on.exit(unlink(tmp), add = TRUE)
	head <- system(paste("zcat", shQuote(gz), "| head -n 1"), intern = TRUE)
	nm <- toupper(strsplit(head, "\\s+")[[1]]); scol <- which(nm %in% c("SNP", "RSID", "ID"))[1]
	if(is.na(scol)) stop("Cannot find SNP column in ", gz)
	cmd <- paste0("zcat ", shQuote(gz), " | awk 'NR==FNR{a[$1];next} FNR==1||($", scol, " in a)' ", shQuote(tmp), " -")
	dat <- data.table::fread(cmd = cmd, header = TRUE) %>% as.data.frame()
	names(dat) <- stringi::stri_replace_all_regex(toupper(names(dat)), pattern = toupper(pattern), replacement = replacement, vectorize_all = FALSE)
	if("CHR" %in% names(dat)) dat$CHR <- ifelse(dat$CHR %in% c(23, "23"), "X", as.character(dat$CHR))
	if(!"N" %in% names(dat)) dat$N <- N_default
	data.table::fwrite(dat, outfile, sep = "\t")
	dat
}

main_mr_p <- function(fit) {
	m0 <- c("Inverse variance weighted", "Wald ratio", "Inverse variance weighted (multiplicative random effects)")
	z <- fit %>% dplyr::filter(method %in% m0) %>% dplyr::mutate(ii = match(method, m0)) %>% dplyr::arrange(ii)
	if(nrow(z) == 0) NA_real_ else z$pval[1]
}

plot_radial_fallback <- function(dat, fit) {
	b.ivw <- fit %>% dplyr::filter(method %in% c("Inverse variance weighted", "Wald ratio", "Inverse variance weighted (multiplicative random effects)")) %>% dplyr::slice(1) %>% dplyr::pull(b)
	b.egger <- fit %>% dplyr::filter(grepl("Egger", method)) %>% dplyr::slice(1) %>% dplyr::pull(b)
	d <- dat %>% dplyr::mutate(ratio = beta.outcome / beta.exposure, Wj = pmax((beta.exposure / se.outcome)^2, 1e-8), x = sqrt(Wj), y = ratio * sqrt(Wj), q = (y - b.ivw[1] * x)^2) %>% dplyr::arrange(dplyr::desc(q)) %>% dplyr::mutate(lab = ifelse(dplyr::row_number() <= pmin(10, n()), SNP, NA))
	p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) + ggplot2::geom_point(size = 1.3, alpha = .75) + ggplot2::geom_abline(slope = b.ivw[1], intercept = 0, color = "steelblue", linewidth = .7) + ggrepel::geom_text_repel(ggplot2::aes(label = lab), size = 2.5, max.overlaps = Inf, na.rm = TRUE) + ggplot2::labs(x = expression(sqrt(W[j])), y = expression(hat(beta)[j]*sqrt(W[j])), title = "Radial estimates") + ggplot2::theme_bw(11) + ggplot2::theme(aspect.ratio = 1)
	if(length(b.egger) > 0 && !is.na(b.egger[1])) p <- p + ggplot2::geom_abline(slope = b.egger[1], intercept = 0, color = "tomato", linewidth = .7)
	p
}

plot_mr_forest_radial <- function(dat, fit, outfile, title = NULL) {
	ord <- c("Inverse variance weighted", "Simple median", "Weighted median", "MR Egger", "Maximum likelihood", "Contamination mixture", "Wald ratio")
	df <- fit %>% dplyr::filter(method %in% ord, !is.na(b), !is.na(se)) %>% dplyr::mutate(lci = b - 1.96 * se, uci = b + 1.96 * se, method = factor(method, levels = rev(ord[ord %in% method])))
	p1 <- if(nrow(df) > 0) ggplot2::ggplot(df, ggplot2::aes(b, method)) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lci, xmax = uci), height = .08, linewidth = .45) + ggplot2::geom_point(size = 1.9) + ggplot2::labs(x = "Causal estimate (95% CI)", y = NULL, title = title) + ggplot2::theme_bw(11) + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank()) else ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = title)
	p2 <- if(nrow(dat) > 2 && requireNamespace("RadialMR", quietly = TRUE)) tryCatch({
		rd <- RadialMR::format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
		ivw <- RadialMR::ivw_radial(rd, alpha = 0.05, weights = 1, tol = 0.0001, summary = FALSE)
		eg <- tryCatch(RadialMR::egger_radial(rd, alpha = 0.05, weights = 1, summary = FALSE), error = function(e) NULL)
		(RadialMR::plot_radial(if(is.null(eg)) ivw else c(ivw, eg), radial_scale = TRUE, show_outliers = FALSE, scale_match = TRUE) + ggplot2::labs(title = "Radial estimates") + ggplot2::theme_bw(11) + ggplot2::theme(aspect.ratio = 1))
	}, error = function(e) plot_radial_fallback(dat, fit)) else plot_radial_fallback(dat, fit)
	p <- patchwork::wrap_plots(p1, p2, nrow = 1, widths = c(1, 1.25))
	ggplot2::ggsave(outfile, p, width = 10, height = 5.8, dpi = 300, bg = "white")
}

plot_coloc_locus <- function(dat, fit, outfile, X, Y, H4 = NA) {
	pp <- tryCatch(fit$results %>% as.data.frame(), error = function(e) NULL)
	if(!is.null(pp) && all(c("snp", "SNP.PP.H4") %in% names(pp))) dat <- dplyr::left_join(dat, pp %>% dplyr::transmute(SNP = snp, SNP.PP.H4), by = "SNP")
	if(!"SNP.PP.H4" %in% names(dat)) dat$SNP.PP.H4 <- 0
	d <- dat %>% dplyr::mutate(pos_plot = ifelse(is.na(POS), dplyr::row_number(), POS) / ifelse(all(is.na(POS)), 1, 1e6)) %>% dplyr::transmute(SNP, POS = pos_plot, SNP.PP.H4, `Protein pQTL` = -log10(P.X), `Disease GWAS` = -log10(P.Y)) %>% tidyr::pivot_longer(c(`Protein pQTL`, `Disease GWAS`), names_to = "Trait", values_to = "log10P")
	xlab <- if(all(is.na(dat$POS))) "SNP index" else paste0("Position on chr ", unique(stats::na.omit(dat$CHR))[1], " (Mb)")
	p <- ggplot2::ggplot(d, ggplot2::aes(POS, log10P)) + ggplot2::geom_point(ggplot2::aes(size = SNP.PP.H4), alpha = .75) + ggplot2::facet_wrap(~Trait, ncol = 1, scales = "free_y") + ggplot2::scale_size_continuous(range = c(.5, 3), guide = "none") + ggplot2::labs(title = paste0(X, " - ", Y, " coloc; PP.H4 = ", signif(as.numeric(H4), 3)), x = xlab, y = expression(-log[10](P))) + ggplot2::theme_bw(11)
	ggplot2::ggsave(outfile, p, width = 7, height = 5, dpi = 300, bg = "white")
}

run_mr2s <- function(analysis, label, X, dir.X, dat.X.raw, IV.file = "NO", IV_filter = TRUE, Y, dir.Y, dat.Y.raw, X.use.top = FALSE, Pth_MR = 0.05, plot_yes = TRUE) {
	log_file <- paste0(label, ".mr.log")
	dat.X.iv <- get_iv_snps(X, dir.X, dat.X.raw, IV.file, X.use.top)
	if(nrow(dat.X.iv) == 0) { write(paste("SKIP:", analysis, X, Y, "X has no IV"), log_file, append = TRUE); return(NULL) }
	dat.X <- merge(dat.X.raw, dat.X.iv, by = "SNP")
	if(!"N" %in% names(dat.X)) dat.X$N <- 100000
	if(isTRUE(IV_filter)) dat.X <- IV.filter(dat.X, "CHR", "POS", "EAF", "BETA", "N")
	if(nrow(dat.X) == 0) { write(paste("SKIP:", analysis, X, Y, "X is empty after IV filter"), log_file, append = TRUE); return(NULL) }
	dat.X <- dat.X %>% TwoSampleMR::format_data(type = "exposure", snp_col = "SNP", chr_col = "CHR", pos_col = "POS", effect_allele_col = "EA", other_allele_col = "NEA", eaf_col = "EAF", samplesize_col = "N", beta_col = "BETA", se_col = "SE", pval_col = "P") %>% dplyr::mutate(exposure = X, id.exposure = X)
	dat.Y.4x <- dat.Y.raw %>% dplyr::filter(SNP %in% dat.X$SNP)
	if(!"N" %in% names(dat.Y.4x)) dat.Y.4x$N <- 100000
	if(nrow(dat.Y.4x) == 0) { write(paste("SKIP:", analysis, X, Y, "IV SNPs do not exist in Y"), log_file, append = TRUE); return(NULL) }
	dat.Y.4x <- dat.Y.4x %>% TwoSampleMR::format_data(type = "outcome", snp_col = "SNP", chr_col = "CHR", pos_col = "POS", effect_allele_col = "EA", other_allele_col = "NEA", eaf_col = "EAF", samplesize_col = "N", beta_col = "BETA", se_col = "SE", pval_col = "P") %>% dplyr::mutate(outcome = Y, id.outcome = Y)
	dat <- TwoSampleMR::harmonise_data(dat.X, dat.Y.4x, action = 1) %>% dplyr::filter(mr_keep)
	write.table(dat, paste0(label, ".", analysis, ".dat"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
	if(nrow(dat) == 0) { write(paste("SKIP:", analysis, X, Y, "harmonized data is empty"), log_file, append = TRUE); return(NULL) }
	if(nrow(dat) > 1) {
		mths <- c("mr_ivw", "mr_ivw_mre", "mr_simple_median", "mr_weighted_median", "mr_two_sample_ml", "mr_egger_regression")
		if("mr_contamination_mixture" %in% TwoSampleMR::mr_method_list()$obj) mths <- c(mths, "mr_contamination_mixture")
		fit <- TwoSampleMR::mr(dat, method_list = mths)
		p.hetero <- tryCatch(TwoSampleMR::mr_heterogeneity(dat)$Q_pval[1], error = function(e) NA)
		p.pleio <- tryCatch(TwoSampleMR::mr_pleiotropy_test(dat)$pval[1], error = function(e) NA)
	} else { fit <- TwoSampleMR::mr(dat, method_list = c("mr_wald_ratio")); p.hetero <- NA; p.pleio <- NA }
	if(isTRUE(plot_yes) && !is.na(main_mr_p(fit)) && main_mr_p(fit) < Pth_MR) tryCatch(plot_mr_forest_radial(dat, fit, paste0(label, ".", analysis, ".MR.png"), paste0(X, " -> ", Y, " (", analysis, ")")), error = function(e) write(paste("PLOT_ERROR:", analysis, X, Y, conditionMessage(e)), log_file, append = TRUE))
	res <- fit %>% dplyr::mutate(analysis = analysis, exposure = X, outcome = Y, n_IV = nrow(dat), heterogeneity_p = p.hetero, pleiotropy_p = p.pleio, P.Y.min = min(dat$pval.outcome, na.rm = TRUE))
	write(paste(analysis, X, Y, nrow(dat), paste(fit$method, rb(fit$b), rb(fit$se), rp(fit$pval), collapse = " | ")), log_file, append = TRUE)
	res
}

run_coloc <- function(label, X, dir.X, Y, dat.Y.raw, Y.type = "quant", Y.s = NULL, N_default = 100000, H4_thres = 0.8, plot_yes = TRUE) {
	log_file <- paste0(label, ".coloc.log")
	f4 <- file.path(dir.X, X, paste0(X, ".4gcta"))
	if(!file.exists(f4) || file.size(f4) == 0) { write(paste("SKIP:", X, Y, "no .4gcta for coloc"), log_file, append = TRUE); return(NULL) }
	dat.X <- read_gwas_mr(f4, N_default)
	dat.Y <- dat.Y.raw %>% dplyr::filter(SNP %in% dat.X$SNP)
	dat.X <- dat.X %>% dplyr::arrange(P) %>% dplyr::distinct(SNP, .keep_all = TRUE)
	dat.Y <- dat.Y %>% dplyr::arrange(P) %>% dplyr::distinct(SNP, .keep_all = TRUE)
	dat <- merge(dat.X, dat.Y, by = "SNP", suffixes = c(".X", ".Y"))
	dat$CHR <- if("CHR.X" %in% names(dat)) as.character(dat$CHR.X) else if("CHR" %in% names(dat)) as.character(dat$CHR) else NA
	if("CHR.Y" %in% names(dat)) dat$CHR <- ifelse(is.na(dat$CHR), as.character(dat$CHR.Y), dat$CHR)
	dat$POS <- if("POS.X" %in% names(dat)) as.numeric(dat$POS.X) else if("POS" %in% names(dat)) as.numeric(dat$POS) else NA
	if("POS.Y" %in% names(dat)) dat$POS <- ifelse(is.na(dat$POS), as.numeric(dat$POS.Y), dat$POS)
	dat <- dat %>% dplyr::mutate(MAF.X = pmin(EAF.X, 1 - EAF.X), MAF.Y = pmin(EAF.Y, 1 - EAF.Y)) %>% dplyr::filter(!is.na(BETA.X), !is.na(SE.X), !is.na(P.X), !is.na(BETA.Y), !is.na(SE.Y), !is.na(P.Y), MAF.X > 0, MAF.X < 0.5, MAF.Y > 0, MAF.Y < 0.5)
	if(nrow(dat) < 50 || min(dat$P.Y, na.rm = TRUE) > 1e-2) { write(paste("SKIP:", X, Y, "coloc n/minP.Y", nrow(dat), signif(min(dat$P.Y, na.rm = TRUE), 3)), log_file, append = TRUE); return(NULL) }
	d1 <- list(type = "quant", snp = dat$SNP, MAF = dat$MAF.X, N = dat$N.X, beta = dat$BETA.X, varbeta = dat$SE.X^2, pvalues = dat$P.X)
	d2 <- list(type = Y.type, snp = dat$SNP, MAF = dat$MAF.Y, N = dat$N.Y, beta = dat$BETA.Y, varbeta = dat$SE.Y^2, pvalues = dat$P.Y)
	if(!all(is.na(dat$POS))) { d1$position <- dat$POS; d2$position <- dat$POS }
	if(Y.type == "cc" && !is.null(Y.s)) d2$s <- Y.s
	fit <- tryCatch(coloc::coloc.abf(d1, d2), error = function(e) { write(paste("ERROR:", X, Y, conditionMessage(e)), log_file, append = TRUE); NULL })
	if(is.null(fit)) return(NULL)
	ss <- fit$summary
	res <- data.frame(X = X, Y = Y, nSNP = ss["nsnps"], PP.H0.abf = ss["PP.H0.abf"], PP.H1.abf = ss["PP.H1.abf"], PP.H2.abf = ss["PP.H2.abf"], PP.H3.abf = ss["PP.H3.abf"], PP.H4.abf = ss["PP.H4.abf"])
	if(isTRUE(plot_yes) && as.numeric(res$PP.H4.abf) >= H4_thres) tryCatch(plot_coloc_locus(dat, fit, paste0(label, ".coloc.png"), X, Y, res$PP.H4.abf), error = function(e) write(paste("PLOT_ERROR:", X, Y, conditionMessage(e)), log_file, append = TRUE))
	write(paste(X, Y, paste(signif(unlist(res[-c(1:2)]), 3), collapse = " ")), log_file, append = TRUE)
	res
}

# Mediation analysis
run_mrMed <- function(label, X, dir.X, dat.X.raw, file.M, Y, dat.Y.raw, X.use.top = FALSE) {
	log_file <- paste0(label, ".mrMed.log")
	dir.M <- dirname(file.M); M <- sub("\\.(gz|top)$", "", basename(file.M))
	dat.X.iv <- get_iv_snps(X = X, dir.X = dir.X, dat.X.raw = dat.X.raw, IV.file = "NO", X.use.top = X.use.top); IV <- dat.X.iv
	dat.M.raw <- read_gwas_mr(file.M)
	if(file.exists(paste0(dir.M, "/", M, ".top.snp"))) { dat.M.iv <- read.table(paste0(dir.M, "/", M, ".top.snp"), header = TRUE); names(dat.M.iv) <- "SNP" } else if(file.exists(paste0(dir.M, "/", M, ".NEW.top.snp"))) { dat.M.iv <- read.table(paste0(dir.M, "/", M, ".NEW.top.snp"), header = TRUE); if(ncol(dat.M.iv) == 1) names(dat.M.iv) <- "SNP" } else if(grepl("\\.top$", file.M)) { dat.M.iv <- unique(dat.M.raw["SNP"]) } else { dat.M.sig <- dat.M.raw %>% dplyr::filter(P <= 5e-08) %>% dplyr::mutate(mb = ceiling(POS / 1e+05)); dat.M.iv <- dat.M.sig %>% dplyr::group_by(mb) %>% dplyr::slice(which.min(P)) %>% dplyr::ungroup() %>% dplyr::select(SNP); write.table(dat.M.iv, paste0(dir.M, "/", M, ".NEW.top.snp"), append = FALSE, quote = FALSE, row.names = FALSE, col.names = TRUE) }
	dat.XnM.iv <- unique(rbind(IV, dat.M.iv))
	dat.X.mv <- dat.X.raw %>% merge(dat.XnM.iv) %>% format_data(type = "exposure", snp_col = "SNP", effect_allele_col = "EA", other_allele_col = "NEA", beta_col = "BETA", se_col = "SE", pval_col = "P") %>% dplyr::mutate(id.exposure = X)
	dat.M.mv <- dat.M.raw %>% merge(dat.XnM.iv)
	if(nrow(dat.M.mv) == 0) { write(paste("SKIP:", X, Y, M, "no SNP in dat.M.mv"), file = log_file, append = TRUE); return(NULL) }
	dat.M.mv <- dat.M.mv %>% format_data(type = "outcome", snp_col = "SNP", effect_allele_col = "EA", other_allele_col = "NEA", beta_col = "BETA", se_col = "SE", pval_col = "P") %>% dplyr::mutate(id.outcome = M)
	dat.Y.mv <- dat.Y.raw %>% merge(dat.XnM.iv) %>% format_data(type = "outcome", snp_col = "SNP", effect_allele_col = "EA", other_allele_col = "NEA", beta_col = "BETA", se_col = "SE", pval_col = "P") %>% dplyr::mutate(id.outcome = Y)
	dat.mv <- harmonise_data(dat.X.mv, dat.Y.mv, action = 1)
	dat.XM.mv <- harmonise_data(dat.X.mv, dat.M.mv, action = 1); names(dat.XM.mv) <- gsub("outcome", "mediator", names(dat.XM.mv))
	dat <- merge(dat.XM.mv, dat.mv, by = "SNP")
	if(nrow(dat) == 0) { write(paste("SKIP:", X, Y, M, "X-M-Y harmonized data is empty"), file = log_file, append = TRUE); return(NULL) }
	bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
	if(bad_row != 0) { write(paste(X, Y, M, "ERR: X-M-Y have inconsistent alleles"), file = log_file, append = TRUE); return(NULL) }
	names(dat) <- gsub("\\.x$", "", names(dat)); dat <- subset(dat, select = !grepl("\\.y", names(dat)))
	dat1 <- dat; names(dat1) <- stringi::stri_replace_all_regex(names(dat1), pattern = c("exposure", "mediator", "outcome"), replacement = c("X", "M", "Y"), vectorize_all = FALSE)
	dat1 <- dat1 %>% dplyr::mutate(Gx = ifelse(dat$SNP %in% IV$SNP, 1, 0), Gx_plum = Gx, Gm = 1, Gm_plum = Gm, G_mvmr = ifelse(Gx_plum == 0 & Gm_plum == 0, 0, 1))
	if(max(dat1$Gx) * max(dat1$Gm, na.rm = TRUE) == 0) { write(paste(X, Y, M, "mrMed data has 0 rows for Gx or Gm"), file = log_file, append = TRUE); return(NULL) }
	res.xy <- if(nrow(dat.mv) > 1) mr(dat.mv, method_list = c("mr_ivw_mre")) else mr(dat.mv, method_list = c("mr_wald_ratio"))
	beta.X2Y <- res.xy$b[1]; p.X2Y <- res.xy$pval[1]
	res <- try({ mrMed(dat_mrMed = dat1, method_list = "Prod_IVW_0") }, silent = TRUE)
	if(inherits(res, "try-error")) { write(paste("ERROR:", X, Y, M, "mrMed gives ERROR message"), file = log_file, append = TRUE); return(NULL) }
	X_str <- paste0(X, "(", nrow(IV), ")"); M_str <- paste0(M, "(", nrow(dat.M.iv), ")")
	write(paste(X_str, M_str, Y, nrow(dat), rb(beta.X2Y), rp(p.X2Y), paste(rb(res$TE$b), rp(res$TE$p), rb(res$IE$b), rp(res$IE$p), rb(res$DE$b), rp(res$DE$p), rb(res$rho$b), rp(res$rho$p))), file = log_file, append = TRUE)
}

run_mr_all <- function(Y, dir.Y, dat.Y.raw, X, dir.X, gene_pos_file = NULL, N_default = 100000, flank_region = 100000, IV_filter = TRUE, mr2s_yes = TRUE, cisMr_yes = FALSE, coloc_yes = TRUE, mrMed_yes = FALSE, Ms.full = NULL, Y.type = "quant", Y.s = NULL, outdir = "MR_coloc", Pth_MR = 0.05, H4_thres = 0.8, plot_yes = TRUE) {
	dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
	label <- paste0(outdir, "/", X, "-", Y); cat(paste("\n-->> NOW PROCESS", X, "-", Y, "\n"))
	write("Analysis X Y n-IV method BETA SE P", file = paste0(label, ".mr.log"), append = FALSE)
	write("X Y nSNP H0 H1 H2 H3 H4", file = paste0(label, ".coloc.log"), append = FALSE)
	res.mr <- list(); res.coloc <- NULL

	fx <- file.path(dir.X, X, paste0(X, ".clump.assoc"))
	if(mr2s_yes && file.exists(fx) && file.size(fx) > 0) { dat.X.raw <- read_gwas_mr(fx, N_default); res.mr[["gwIV.way1"]] <- run_mr2s("gwIV.way1", label, X, dir.X, dat.X.raw, "TOP", IV_filter, Y, dir.Y, dat.Y.raw, Pth_MR = Pth_MR, plot_yes = plot_yes) } else write(paste("SKIP:", X, Y, "no X.clump.assoc"), file = paste0(label, ".mr.log"), append = TRUE)

	if(mr2s_yes) {
		fy <- c(file.path(dir.Y, Y, paste0(Y, ".clump.assoc")), file.path(dir.Y, paste0(Y, ".clump.assoc"))); fy <- fy[file.exists(fy) & file.size(fy) > 0][1]
		fgz <- file.path(dir.X, X, paste0(X, ".gz"))
		if(!is.na(fy) && file.exists(fgz)) { dat.Y.iv <- read_gwas_mr(fy, N_default); dat.X.4Y <- extract_snps_from_gz(fgz, dat.Y.iv$SNP, paste0(outdir, "/", X, ".for_", Y, ".clump.assoc"), N_default); res.mr[["gwIV.way2"]] <- run_mr2s("gwIV.way2", label, Y, dir.Y, dat.Y.iv, "TOP", IV_filter, X, dir.X, dat.X.4Y, Pth_MR = Pth_MR, plot_yes = FALSE) } else write(paste("SKIP:", X, Y, "no Y.clump.assoc or X.gz for gwIV.way2"), file = paste0(label, ".mr.log"), append = TRUE)
		f4 <- file.path(dir.X, X, paste0(X, ".4gcta")); fj <- file.path(dir.X, X, paste0(X, ".jma.cojo"))
		if(file.exists(f4) && file.size(f4) > 0 && file.exists(fj) && file.size(fj) > 0) { dat.X.cis <- read_gwas_mr(f4, N_default); res.mr[["cis.jma"]] <- run_mr2s("cis.jma", label, X, dir.X, dat.X.cis, fj, IV_filter, Y, dir.Y, dat.Y.raw, Pth_MR = Pth_MR, plot_yes = plot_yes) } else write(paste("SKIP:", X, Y, "no X.4gcta or X.jma.cojo for cis.jma"), file = paste0(label, ".mr.log"), append = TRUE)
	}

	if(isTRUE(cisMr_yes)) write("NOTE: cisMr_yes ignored. cis.jma is implemented by run_mr2s(). Add real cisMr() later.", file = paste0(label, ".mr.log"), append = TRUE)
	if(coloc_yes) res.coloc <- run_coloc(label, X, dir.X, Y, dat.Y.raw, Y.type, Y.s, N_default, H4_thres, plot_yes)
	if(mrMed_yes && !is.null(Ms.full) && exists("dat.X.raw")) for(file.M in Ms.full) run_mrMed(label, X, dir.X, dat.X.raw, file.M, Y, dat.Y.raw)
	list(mr = dplyr::bind_rows(res.mr), cisMr = NULL, coloc = res.coloc)
}
