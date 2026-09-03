pacman::p_load(tidyverse, data.table, stringi, patchwork, qqman, CMplot, topr, TwoSampleMR, plotly)	

dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:/", "/work/sph-huangj")
dirgwa = paste0(dir0, '/data/gwas/stroke/clean')
source(paste0(dir0, '/scripts/0f/0phe.f.R'))
source(paste0(dir0, '/scripts/0f/plot.f.R'))

runs = NA # c("qt", "camsis")
gwass = c("bald.qt", "bald12.bt", "bald13.bt", "bald14.bt")


# 🚩 读取数据
for (run in runs) {
for (gwas in gwass) {
	if (!is.na(run)) fn <- paste0(gwas, ".", run) else fn <- gwas
	dat <- fread(paste0(dirgwa, "/", fn, ".gz"), header = TRUE) %>%
	thinP0(., P = 1e-3, p_col = "P") %>% filter(EAF > 0.005, EAF < 0.995) %>% as.data.frame()			  
	names(dat) <- stri_replace_all_regex(toupper(names(dat)), pattern = toupper(pattern), replacement = replacement, vectorize_all = FALSE)
	dat <- dat %>% dplyr::select(SNP, CHR, POS, EA, NEA, EAF, BETA, SE, P) %>% mutate(
		CHR = ifelse(CHR == "X", 23, CHR), CHR = as.numeric(CHR), 
		across(c(EAF, BETA, SE, P), as.numeric), P = thinP1(P)
	)
	assign(paste0(fn, '.dat.mh'), dat)
	# png(paste0(fn, '.png'), w=1200, h=800); qqman::manhattan(get(paste0(gwas, ".", run, '.dat.mh')), chr = "CHR", bp = "POS", p = "P", snp = "SNP", col = c("blue4", "orange3")); dev.off()
	assign(paste0(fn, ".dat.cm"), dat %>% dplyr::select(SNP, CHR, POS, P) %>% { nm <- paste0("P.", gwas); setNames(., c("SNP","CHR","POS", nm)) })
	assign(paste0(fn, '.dat.tr'), dat %>% rename(REF = NEA, ALT = EA))
	assign(paste0(fn, '.dat.mr'), dat %>% TwoSampleMR::format_data(type = 'exposure', snp_col = 'SNP', chr_col = 'CHR', pos_col = 'POS', effect_allele_col = 'EA', other_allele_col = 'NEA', eaf_col = 'EAF', beta_col = 'BETA', se_col = 'SE', pval_col = 'P'))
	rm(dat); gc()
}
}


# 🚩 曼哈顿图比较
datCM_list <- mget(paste0(gwass, ".dat.cm"))
	datCM <- Reduce(function(x, y) merge(x, y, by = c("SNP", "CHR", "POS"), all = TRUE), datCM_list)
	datCM <- subset(datCM, select = c('SNP', 'CHR', 'POS', grep('^P\\.', names(datCM), value = TRUE)))
	CMplot(datCM, plot.type = "m", multracks = TRUE, cex = 0.2, amplify = FALSE, file.output = TRUE, file = "jpg", file.name = paste(gwas, "cmplot", sep = "."), width = 20, height = 5, dpi = 300) 
	plt.topr <-	topr::manhattan(mget(paste0(gwass, ".dat.tr")), annotate = 5e-20, color = c("darkgray", "blue"), legend_labels=gwass, ntop=1, title=gwas)
	png('topr.png', width = 3508, height = 2480, res = 300); plt.topr; dev.off()


# 🚩 EAF和BETA比较
for (gwas in gwass) {
	X <- paste0(gwas, ".", runs[1])
	Y <- paste0(gwas, ".", runs[2])
	datX <- get(paste0(X, ".dat.mr")); nrow(datX) 
	datY <- get(paste0(Y, ".dat.mr")); nrow(datY); names(datY) <- gsub('exposure', 'outcome', names(datY))
	dat <- harmonise_data(datX, datY, action = 1) %>% filter(pval.exposure < 1e-03)
	dat$logP.exposure <- -log10(dat$pval.exposure); dat$logP.outcome <- -log10(dat$pval.outcome)

	plots <- list()
	for (v in c("beta", "eaf", "logP")) {
		dat$exposure <- dat[[paste0(v, ".exposure")]]
		dat$outcome <- dat[[paste0(v, ".outcome")]]
		cor1 <- cor(dat$exposure, dat$outcome, use = "complete.obs")
		plt <- ggplot(dat, aes(x = exposure, y = outcome)) +
			geom_point(alpha = 0.35, size = 0.6) + geom_abline(slope = 1, intercept = 0) +
			labs(x = paste(X, v), y = paste(Y, v), title = sprintf("r = %.3f; n=%d", cor1, nrow(dat))) +theme_bw()
		plots[[v]] <- plt
	}		
	plt <- (plots$beta | plots$eaf) / (plots$logP | patchwork::plot_spacer()); plt
	ggsave(paste0("comp.", gwas, ".jpeg"), plot = plt, width = 12, height = 12, dpi = 300)
}


# 🚩 查找发表的GWAS结果 http://www.phenoscanner.medschl.cam.ac.uk
pacman::p_load(phenoscanner, ieugwasr)
snps <- read.table("D:/files/stroke_i.snp", header = FALSE)
res <- phenoscanner(snpquery="rs123,rs456,rs789", catalogue = "GWAS", pvalue = 5e-8, proxies = "None", build = 37)
head(res$results)
res_list <- lapply(snps, function(s) ieugwasr::phewas(variant=s, pval=5e-8))
phewas(variant = c("rs123","rs456"), pval = 5e-8)


# 🚩 LDSC热力图
pacman::p_load(ComplexHeatmap)
gwass <- list.files("D:/data/gwas/stroke/clean", pattern = "\\.gz$", full.names = FALSE) %>% sub("\\.gz$", "", .)
dat <- read.table("D:/analysis/ldsc/stroke/all.rg.res", header = TRUE, as.is = TRUE) %>%
#	filter(if_all(c(p1, p2), ~ grepl("cad|self_ihd|self_stroke|self_bmi|st_.*_EUR|move|self_camsis|self_siops", .x) & !grepl("\\.t2e", .x))) %>%
	filter(if_all(c(p1, p2), ~ .x %in% gwass & !grepl("AMR$|AFR$|SAS$|\\.t2e", .x))) %>%
	mutate(rg = pmax(-1, pmin(1, rg)))
plot_rg(dat, "p1", "p2", "rg", "p", alpha = 0.05)
