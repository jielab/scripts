pacman::p_load(tidyverse, data.table, survival)

parse_args <- function(args) {
	out <- list(resume = TRUE)
	i <- 1L
	while (i <= length(args)) {
		key <- args[[i]]
		if (key %in% c("-h", "--help")) {
			cat("Run ./assoc.sh --help for usage.\n"); quit(save = "no")
		}
		if (key == "--no-resume") { out$resume <- FALSE; i <- i + 1L; next }
		if (!startsWith(key, "--") || i == length(args)) stop("Invalid or incomplete option: ", key)
		name <- sub("^--", "", key)
		if (!name %in% c("method", "action", "Xs", "varX", "Ys", "Y", "type", "mode", "scale-X", "scale_X", "out-prefix")) stop("Unknown option: ", key)
		out[[name]] <- args[[i + 1L]]
		i <- i + 2L
	}
	out
}

split_csv <- function(x) {
	if (is.null(x) || !nzchar(trimws(x))) return(character())
	y <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
	unique(y[nzchar(y)])
}
parse_bool <- function(x, name) {
	y <- toupper(trimws(x))
	if (y %in% c("TRUE", "T", "1", "YES", "Y")) return(TRUE)
	if (y %in% c("FALSE", "F", "0", "NO", "N")) return(FALSE)
	stop(name, " must be TRUE or FALSE")
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
method <- opt$method %||% stop("--method is required")
action <- opt$action %||% stop("--action is required")
dir0 <- Sys.getenv("DIR0", unset = ifelse(Sys.info()[["sysname"]] == "Windows", "D:/", "/mnt/d"))
indir <- Sys.getenv("PHEDIR", unset = file.path(dir0, "data/ukb/phe"))
outdir <- Sys.getenv("UKB_OUT", unset = file.path(dir0, "analysis/ukb/assoc"))
script_dir <- Sys.getenv("SCRIPT_DIR", unset = file.path(dir0, "scripts"))
outdir <- file.path(outdir, method)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

invisible(lapply(c("0phe.f.R", "assoc.f.R", "plot.f.R"), function(f) source(file.path(script_dir, "0f", f))))
read_rds <- function(name) readRDS(file.path(indir, "Rdata", paste0(name, ".rds")))
analysis_file <- function(...) file.path(outdir, ...)

gen <- read_rds("gen")

Xs <- split_csv(opt$Xs)
if (!length(Xs)) Xs <- names(gen)[-1]
unknown_x <- setdiff(Xs, names(gen))
if (length(unknown_x)) stop("Unknown --Xs variable(s) in gen.rds: ", paste(unknown_x, collapse = ", "))

mode <- tolower(opt$mode %||% "")
if (nzchar(mode) && !mode %in% c("res", "add", "dom")) stop("--mode must be one of: res, add, dom")
if (nzchar(mode) && !length(split_csv(opt$Xs))) stop("--mode requires an explicit --Xs list of genotype variables")
apply_genetic_mode <- function(dat, Xs, mode) {
	if (!nzchar(mode)) return(dat)
	for (X in Xs) {
		x <- dat[[X]]
		if (!is.numeric(x)) stop("--mode can only be used with numeric genotype data; ", X, " is not numeric")
		finite <- x[is.finite(x)]
		if (!length(finite)) stop("--mode cannot be used because ", X, " has no finite genotype values")
		if (any(finite < 0 | finite > 2)) stop("--mode requires genotype/dosage values in [0, 2]; invalid values found in ", X)
		has_dosage <- any(abs(finite - round(finite)) > 1e-8)
		if (mode == "add") {
			if (has_dosage) message(X, ": additive mode retains imputed dosage values")
		} else {
			if (has_dosage) message(X, ": ", mode, " mode hard-calls imputed dosage by rounding to 0/1/2")
			hard_call <- round(x)
			dat[[X]] <- if (mode == "res") as.numeric(hard_call == 2) else as.numeric(hard_call >= 1)
			dat[[X]][is.na(x)] <- NA_real_
		}
	}
	dat
}
invisible(apply_genetic_mode(gen[, c("eid", Xs), drop = FALSE], Xs, mode))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Phewas
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (method %in% c("phewas", "medwas")) {
	non_numeric_x <- Xs[!vapply(gen[Xs], is.numeric, logical(1))]
	if (length(non_numeric_x)) message("Skipping non-numeric gen variables for ", method, ": ", paste(non_numeric_x, collapse = ", "))
	Xs <- setdiff(Xs, non_numeric_x)
	if (!length(Xs)) stop("No numeric gen variables selected for ", method)
}
varX <- split_csv(opt$varX)
if (!length(varX)) varX <- c(vars.basic, "le8.sco")
out_prefix <- opt[["out-prefix"]] %||% "gen"
if (!grepl("^[A-Za-z0-9._-]+$", out_prefix)) stop("--out-prefix may contain only letters, numbers, dot, underscore, and hyphen")
if (nzchar(mode)) out_prefix <- paste0(out_prefix, ".", mode)
saved_phewas_file <- analysis_file(paste0(out_prefix, ".phewas.rds"))
prior_phewas_res <- if (method == "phewas" && file.exists(saved_phewas_file)) readRDS(saved_phewas_file) else NULL
reuse_phewas <- method == "phewas" && isTRUE(opt$resume) && !is.null(prior_phewas_res) &&
	all(Xs %in% unique(as.character(prior_phewas_res$snp)))
needs_data <- (method == "phewas" && action %in% c("all", "save_rds") && !reuse_phewas) ||
	(method %in% c("forest", "medwas") && action %in% c("all", "save_rds"))
if (needs_data) {
	dat1 <- read_rds("all") %>% filter(ethnic.c == "White")
	missing_x <- setdiff(Xs, names(dat1))
	if (length(missing_x)) stop("Selected gen variables missing from all.rds: ", paste(missing_x, collapse = ", "))
	missing_varX <- setdiff(varX, names(dat1))
	if (length(missing_varX)) stop("Covariates missing from all.rds: ", paste(missing_varX, collapse = ", "))
	dat1 <- apply_genetic_mode(dat1, Xs, mode)
}

cat("METHOD=", method, " ACTION=", action, "\n", sep = "")
cat("Xs=", paste(Xs, collapse = ","), "\n", sep = "")
if (nzchar(mode)) cat("mode=", mode, "\n", sep = "")
if (needs_data) cat("varX=", paste(varX, collapse = ","), "\n", sep = "")

if (method == "forest") {
	result_file <- analysis_file(paste0(out_prefix, ".forest.rds"))
	if (action %in% c("all", "save_rds")) {
		Ys <- split_csv(opt$Ys %||% opt$Y)
		if (!length(Ys)) Ys <- sub("^fod_icd10_", "", grep("^fod_icd10_", names(dat1), value = TRUE))
		if (!length(Ys)) stop("No outcomes selected or discovered")
		type <- opt$type %||% "t2e"
		if (!type %in% c("bt", "qt", "t2e", "t2e.tdc", "ordinal")) stop("Unsupported --type: ", type)
		scale_arg <- opt[["scale-X"]] %||% opt$scale_X
		scale_X <- if (is.null(scale_arg)) TRUE else parse_bool(scale_arg, "--scale-X/--scale_X")
		if (type %in% c("t2e", "t2e.tdc")) {
			for (Y in Ys) {
				dat1[grep(paste0("^", Y, "\\.Y[?]t2e"), names(dat1))] <- NULL
				dat1 <- t2e(dat1, NA, paste0("fod_icd10_", Y), "birth_date", "date_attend", "date_lost", "date_death", date_follow_end, Y, "year")
			}
		}
		cat("Ys=", paste(Ys, collapse = ","), " type=", type, " scale_X=", scale_X, "\n", sep = "")
		res <- assoc_reg(dat1, Xs = Xs, varX = varX, Y = Ys, type = type, scale_X = scale_X) %>%
			mutate(lab.Y = coalesce(unname(dx.lst[Outcome]), as.character(Outcome)), color = "black") %>% as.data.frame()
		saveRDS(res, result_file)
		message("Wrote: ", result_file)
	}
	if (action %in% c("all", "plot")) {
		if (action == "plot") {
			if (!file.exists(result_file)) stop("Saved forest result not found: ", result_file, "\nRun ACTION save_rds first.")
			res <- readRDS(result_file)
		}
		plot_xs <- intersect(Xs, unique(as.character(res$Exposure)))
		if (!length(plot_xs)) stop("None of the selected --Xs are present in: ", result_file)
		forest_plots <- plot_forest(res[res$Exposure %in% plot_xs, , drop = FALSE], lab.X = plot_xs, n_col = "N_event")
		for (X in names(forest_plots)) {
			p <- forest_plots[[X]]
			if (!is.null(p)) {
				png <- analysis_file(paste0(X, if (nzchar(mode)) paste0(".", mode) else "", ".forest.png"))
				ggsave(png, p, width = 9, height = 12, dpi = 300, bg = "white", limitsize = FALSE)
				message("Wrote: ", png)
			}
		}
	}
} else if (method == "phewas") {
	if (!requireNamespace("PheWAS", quietly = TRUE)) stop("R package PheWAS is required")
	suppressPackageStartupMessages(library(PheWAS))
	result_file <- saved_phewas_file
	res <- NULL
	if (action %in% c("all", "save_rds")) {
		if (reuse_phewas) {
			message("Using existing PheWAS result: ", result_file)
			res <- prior_phewas_res
		} else {
			if (!file.exists(file.path(indir, "Rdata", "PheWAS.rds"))) stop("Missing PheWAS.rds; generate it in the phenotype pipeline first")
			checkpoint_outdir <- if (nzchar(mode)) file.path(outdir, paste0("mode_", mode)) else outdir
			phewas_res <- plot_phewas(dat1, phecode = NA, Xs = Xs, varX = varX,
				output_dir = checkpoint_outdir, resume = opt$resume, make_plots = FALSE)
			res <- phewas_res$res
			if (!is.null(prior_phewas_res)) {
				res <- dplyr::bind_rows(prior_phewas_res[!prior_phewas_res$snp %in% Xs, , drop = FALSE], res)
			}
			saveRDS(res, result_file)
			message("Wrote: ", result_file)
		}
	}
	if (action %in% c("all", "plot")) {
		if (is.null(res)) {
			if (!file.exists(result_file)) stop("Saved PheWAS result not found: ", result_file, "\nRun ACTION save_rds first.")
			message("Reading existing PheWAS result: ", result_file)
			res <- if (is.null(prior_phewas_res)) readRDS(result_file) else prior_phewas_res
		}
		plot_xs <- intersect(Xs, unique(as.character(res$snp)))
		if (!length(plot_xs)) stop("None of the selected --Xs are present in: ", result_file)
		for (X in plot_xs) {
			dat2 <- res[res$snp == X & is.finite(res$p), , drop = FALSE]
			if (!nrow(dat2)) { warning("No finite PheWAS results for ", X, call. = FALSE); next }
			p <- suppressMessages(phewasManhattan(dat2,
				annotate.phenotype.description = dat2[, c("phenotype", "description")],
				title = X, OR.direction = FALSE, size.x.labels = 14, size.y.labels = 14)) +
				theme(text = element_text(size = 12))
			png <- analysis_file(paste0(X, if (nzchar(mode)) paste0(".", mode) else "", ".phewas.png"))
			ggsave(png, p, width = 14, height = 8, dpi = 300, bg = "white", limitsize = FALSE)
			message("Wrote: ", png)
		}
	}
} else if (method == "medwas") {
	med_names <- c(img = "Imaging", met = "Metabolites", prot = "Proteome", bbc = "Biochemistry and blood counts")
	result_file <- analysis_file(paste0(out_prefix, ".medwas.rds"))
	res <- NULL
	if (action %in% c("all", "save_rds")) {
		res <- setNames(vector("list", length(med_names)), names(med_names))
		for (M in names(med_names)) {
			med_file <- file.path(indir, "Rdata", paste0(M, ".rds"))
			if (!file.exists(med_file)) stop("Missing MedWAS input: ", med_file)
			dd <- readRDS(med_file)
			if (!"eid" %in% names(dd)) stop(med_file, " does not contain eid")
			Ys <- setdiff(names(dd), "eid")
			Ys <- Ys[vapply(dd[Ys], is.numeric, logical(1))]
			if (!length(Ys)) { warning("No numeric traits in ", med_file, call. = FALSE); next }
			base <- dat1[, unique(c("eid", Xs, varX)), drop = FALSE]
			med_dat <- merge(base, dd[, c("eid", Ys), drop = FALSE], by = "eid")
			med_dat[Ys] <- lapply(med_dat[Ys], function(y) {
				y <- as.numeric(y)
				if (sum(is.finite(y)) < 2L || stats::sd(y, na.rm = TRUE) == 0) return(rep(NA_real_, length(y)))
				as.numeric(scale(y))
			})
			message("Running MedWAS group ", M, " (", length(Ys), " traits)")
			res[[M]] <- assoc_reg(med_dat, Xs = Xs, varX = varX, Y = Ys, type = "qt", scale_X = TRUE)
			res[[M]]$group <- M
			saveRDS(res, result_file)
		}
	}
	if (action %in% c("all", "plot")) {
		if (is.null(res)) {
			if (!file.exists(result_file)) stop("Saved MedWAS result not found: ", result_file, "\nRun ACTION save_rds first.")
			res <- readRDS(result_file)
		}
		if (!requireNamespace("cowplot", quietly = TRUE) || !requireNamespace("ggrepel", quietly = TRUE)) stop("R packages cowplot and ggrepel are required for MedWAS plots")
		for (X in Xs) {
			panels <- lapply(names(med_names), function(M) {
				d <- res[[M]]
				if (is.null(d)) d <- data.frame()
				if (nrow(d)) d <- d[d$Exposure == X & is.finite(d$estimate) & is.finite(d$p.value) & d$p.value > 0, , drop = FALSE]
				if (!nrow(d)) return(ggplot() + theme_void() + labs(title = paste(med_names[[M]], "(no results)")))
				plot_volcano(med_names[[M]], d, "Outcome", "estimate", "p.value", 0.05 / nrow(d), "Standardized beta", expression(-log[10](P)), topN = min(10L, nrow(d)))
			})
			p <- cowplot::plot_grid(plotlist = panels, ncol = 2, labels = names(med_names))
			png <- analysis_file(paste0(X, if (nzchar(mode)) paste0(".", mode) else "", ".medwas.png"))
			ggsave(png, p, width = 18, height = 14, dpi = 300, bg = "white", limitsize = FALSE)
			message("Wrote: ", png)
		}
	}
} else stop("Unknown method: ", method)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 MedWAS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (FALSE) {
Xs <- names(gen)[-1]
vars.base <- setdiff(vars.basic, "sex")
res <- list()
for (M in c("img", "met", "prot", "bbc")) {
	message("Do 👉 ", X, " 👈 ", Sys.time())
	tryCatch({
		dd <- readRDS(paste0(indir, "/Rdata/", X, ".rds"))
		tmp <- merge_check(dat_all = dat_all, dd = dd, by = "eid")
		dat1 <- merge(tmp$dat_all, tmp$dd, by = "eid") %>% filter(sex == 1)
		dat1$bald.o <- droplevels(factor(dat1$bald.c, ordered = TRUE))
		Xs <- setdiff(names(dd), "eid")
		res[[paste0(X, ".bald12")]] <- assoc_reg(dat1, Xs, vars.base, Y = "bald12", type = "bt")
		saveRDS(res, paste0(X, ".bald.rds"))
	}, error = function(e) message("⚠️ skip ", X, ": ", conditionMessage(e)))
}
bind_rows(res, .id = "analysis") %>% filter(p.value < 1e-06) %>% group_by(analysis) %>% arrange(p.value) %>% filter(term == "img_p26536") %>% ungroup() %>% dplyr::select(analysis, term, estimate, p.value)
}
