#!/usr/bin/env Rscript
# Unified R entry point generated from prep_input.R, prep_input_archaic.R, make_hap.R, and make_phy.R.
# Commands: prep_input, add_positive_loci, prep_archaic, make_hap, make_phy.

usage <- function() {
    stop(paste("Usage:",
        "  Rscript locus.R prep_input --dirgwas DIR [--dircojo DIR] --dirout DIR --traits trait ...",
        "  Rscript locus.R prep_input_local --dirgwas DIR --dirout DIR --dirmod DIR [--ref_pop ALL|EUR] --traits trait ...",
        "  Rscript locus.R add_positive_loci --dirgwas DIR --dirmod DIR --dirout DIR --bed BED --trait TRAIT [--ref_pop ALL|EUR]",
        "  Rscript locus.R prep_archaic kg.vcf.gz archaic.raw.vcf.gz out.vcf.gz sample_name",
        "  Rscript locus.R make_hap RES_DIR [trait | trait id | merge]",
        "  Rscript locus.R filter_hap RES_DIR sample.tsv POP MAX_COUNT",
        "  Rscript locus.R make_phy RES_DIR",
        "  Rscript locus.R make_tree RES_DIR",
        "  Rscript locus.R viome --dirgwas DIR --dirout DIR --dirmod DIR --asnp aSNPs.haplotypes.v1.tsv --traits trait ...",
        sep = "\n"), call. = FALSE)
}

run_prep_input <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) script_args else base::commandArgs(FALSE)
    }
    library(data.table)

    args <- commandArgs(TRUE)
    arg <- function(k, multi = FALSE){
    	i <- which(args %in% k)[1]
    	if(is.na(i) || i == length(args)) return(NULL)
    	j <- i + 1
    	if(!multi) return(args[j])
    	z <- args[j:length(args)]
    	z[seq_len(which(c(grepl("^-", z[-1]), TRUE))[1])]
    }
    pick1 <- function(nm, cand){ x <- cand[cand %chin% nm]; if(length(x)) x[1] else NA_character_ }
    norm_chr <- function(x){ x <- toupper(sub("^CHR", "", as.character(x))); x[x == "X"] <- "23"; suppressWarnings(as.integer(x)) }

    dirgwas <- arg(c("--dirgwas", "-dirgwas"))
    dircojo <- arg(c("--dircojo", "-dircojo"))
    dirout  <- arg(c("--dirout", "-dirout"))
    traits  <- arg(c("--traits", "-traits"), TRUE)
    if(is.null(dirgwas) || is.null(dirout) || is.null(traits)) stop("usage: Rscript locus.R prep_input --dirgwas DIR [--dircojo DIR] --dirout DIR --traits trait ...")
    traits <- unlist(strsplit(paste(traits, collapse = " "), "[, ]+")); traits <- traits[nzchar(traits)]

    lead0 <- file.path(dirout, "lead")
    dir.create(lead0, recursive = TRUE, showWarnings = FALSE)

    read_cojo <- function(f){
    	x <- fread(f)
    	cc <- c(chr = pick1(names(x), c("Chr", "CHR", "chr", "#CHROM", "CHROM")),
    	        snp = pick1(names(x), c("SNP", "snp", "lead_snp", "rsid", "ID")),
    	        bp  = pick1(names(x), c("bp", "BP", "pos", "POS", "lead_bp")))
    	if(anyNA(cc)) stop("COJO needs chr/SNP/bp: ", f)
    	unique(x[, .(lead_chr = norm_chr(get(cc["chr"])), lead_snp = as.character(get(cc["snp"])), lead_bp = as.integer(get(cc["bp"])))])[lead_chr %between% c(1L, 23L) & !is.na(lead_bp) & nzchar(lead_snp)]
    }

    read_gwas <- function(f){
    	nm <- names(fread(f, nrows = 0))
    	cc <- c(snp = pick1(nm, c("SNP", "rsid", "RSID", "ID")),
    	        chr = pick1(nm, c("CHR", "Chr", "chr")),
    	        bp  = pick1(nm, c("POS", "BP", "bp", "pos")),
    	        ea  = pick1(nm, c("EA", "A1", "effect_allele", "EffectAllele")),
    	        oa  = pick1(nm, c("NEA", "A2", "other_allele", "OtherAllele")),
    	        b   = pick1(nm, c("BETA", "beta", "b")))
    	if(anyNA(cc)) stop("GWAS needs SNP/CHR/POS/EA/NEA/BETA: ", f)
    	x <- fread(f, select = unname(cc), showProgress = FALSE)
    	setnames(x, unname(cc), c("lead_snp", "gwas_chr", "gwas_bp", "effect_allele", "other_allele", "beta"))
    	x[, `:=`(lead_snp = as.character(lead_snp), gwas_chr = norm_chr(gwas_chr), gwas_bp = as.integer(gwas_bp), effect_allele = toupper(effect_allele), other_allele = toupper(other_allele), beta = as.numeric(beta))]
    	x[nzchar(lead_snp) & !is.na(beta) & nzchar(effect_allele) & nzchar(other_allele)]
    }

    one <- function(tr){
    	# Local layout used by locus.local.sh: DIR/cojo/<trait>/<trait>.jma.cojo and DIR/clean/<trait>/<trait>.gz.
    	# Older COJO layouts are kept as fallbacks when --dircojo is supplied.
    	cj0 <- c(
    		file.path(dirgwas, "cojo", tr, paste0(tr, ".jma.cojo")),
    		file.path(dirgwas, "cojo", tr, "jma.cojo"),
    		if(!is.null(dircojo)) c(
    			file.path(dircojo, tr, "37.jma.cojo"),
    			file.path(dircojo, tr, "jma.cojo"),
    			file.path(dircojo, tr, paste0(tr, ".jma.cojo")),
    			file.path(dircojo, paste0(tr, "37.jma.cojo")),
    			file.path(dircojo, paste0(tr, ".37.jma.cojo")),
    			file.path(dircojo, paste0(tr, ".jma.cojo"))
    		) else character()
    	)
    	cj <- cj0[file.exists(cj0)][1]
    	gw0 <- c(
    		file.path(dirgwas, "clean", tr, paste0(tr, ".gz")),
    		file.path(dirgwas, tr, paste0(tr, ".gz")),
    		file.path(dirgwas, paste0(tr, ".gz"))
    	)
    	gw <- gw0[file.exists(gw0)][1]
    	if(is.na(cj)) stop("missing COJO: ", paste(cj0, collapse = "; "))
    	if(is.na(gw) || !file.exists(gw)) stop("missing GWAS: ", paste(gw0, collapse = "; "))

    	lead <- read_cojo(cj)
    	if(tr == "bald.qt" && !any(lead$lead_snp == "rs35044562")) {
    		lead <- rbind(lead, data.table(lead_chr = 3L, lead_snp = "rs35044562", lead_bp = 45909024L))
    	}
    	setorder(lead, lead_chr, lead_bp, lead_snp)
    	fwrite(lead, file.path(lead0, paste0(tr, ".lead.tsv")), sep = "\t")

    	g <- read_gwas(gw)
    	z <- merge(lead[, row_id := .I], g, by = "lead_snp", all.x = TRUE, sort = FALSE, allow.cartesian = TRUE)
    	z[, match_rank := fifelse(!is.na(gwas_chr) & gwas_chr == lead_chr, 0L, 1L)]
    	setorder(z, row_id, match_rank)
    	z <- z[, .SD[1], by = row_id]
    	z[, match_method := fifelse(!is.na(effect_allele), "SNP", NA_character_)]

    	miss <- is.na(z$effect_allele) | is.na(z$other_allele) | is.na(z$beta)
    	if(any(miss)){
    		fb <- merge(z[miss, .(row_id, lead_chr, lead_bp)], g[, .(lead_chr = gwas_chr, lead_bp = gwas_bp, effect_allele, other_allele, beta)], by = c("lead_chr", "lead_bp"), all.x = TRUE, sort = FALSE, allow.cartesian = TRUE)
    		fb <- fb[, .SD[1], by = row_id]
    		m <- match(z$row_id, fb$row_id)
    		fill <- miss & !is.na(m) & !is.na(fb$effect_allele[m])
    		z[fill, `:=`(gwas_chr = lead_chr, gwas_bp = lead_bp, effect_allele = fb$effect_allele[m[fill]], other_allele = fb$other_allele[m[fill]], beta = fb$beta[m[fill]], match_method = "CHR_BP")]
    	}

    	z[, `:=`(trait = tr, risk_allele = fifelse(beta >= 0, effect_allele, other_allele))]
    	if(tr == "bald.qt" && any(z$lead_snp == "rs35044562")) {
    		z[lead_snp == "rs35044562", `:=`(gwas_chr = 3L, gwas_bp = 45909024L, effect_allele = "ALT", other_allele = "REF", beta = NA_real_, risk_allele = "ALT", match_method = "Nature_rs35044562_ALT")]
    	}
    	bad <- z[is.na(risk_allele) | !nzchar(risk_allele)]
    	if(nrow(bad)) {
    		fwrite(bad, file.path(lead0, paste0(tr, ".lead.fail.tsv")), sep = "\t")
    		stop("missing risk allele for ", tr)
    	}
    	out <- z[, .(trait, lead_chr, lead_snp, lead_bp, effect_allele, other_allele, beta, risk_allele, match_method)]
    	fwrite(out, file.path(lead0, paste0(tr, ".lead.assoc")), sep = "\t")
    	out
    }

    lead_assoc <- rbindlist(lapply(traits, one), fill = TRUE)
    setorder(lead_assoc, trait, lead_chr, lead_bp, lead_snp)
    fwrite(lead_assoc, file.path(lead0, "lead.assoc"), sep = "\t")
    message("Done: ", file.path(lead0, "lead.assoc"))
}
run_prep_input_local <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table)
    fphe <- "/mnt/d/scripts/f/phe.f.R"
    if (file.exists(fphe)) 
        source(fphe)
    args <- commandArgs(TRUE)
    arg <- function(k, multi = FALSE) {
        i <- which(args %in% k)[1]
        if (is.na(i) || i == length(args)) 
            return(NULL)
        j <- i + 1
        if (!multi) 
            return(args[j])
        z <- args[j:length(args)]
        z[seq_len(which(c(grepl("^-", z[-1]), TRUE))[1])]
    }
    pick1 <- function(nm, z) {
        x <- z[z %chin% nm]
        if (length(x)) 
            x[1]
        else NA_character_
    }
    chr_int <- function(x) {
        x <- toupper(sub("^CHR", "", as.character(x)))
        x[x == "X"] <- "23"
        x[x == "Y"] <- "24"
        suppressWarnings(as.integer(x))
    }
    id_mode <- function(x) {
        x <- unique(na.omit(as.character(x)))
        x <- x[nzchar(x) & x != "."]
        x <- head(x, 5000)
        if (!length(x)) 
            return("missing")
        y <- toupper(x)
        rs <- mean(grepl("^RS[0-9]+$", y))
        cp <- mean(grepl("^(CHR)?([0-9]+|X|Y|MT|M):[0-9]+:[ACGTN]+:[ACGTN,]+$", y))
        if (rs > 0.8) 
            "rsid"
        else if (cp > 0.8) 
            "chrpos"
        else "other"
    }
    usable_id <- function(x) {
        !is.na(x) & nzchar(x) & x != "."
    }
    allele_score <- function(ea, nea, ref, alt) {
        ea <- toupper(as.character(ea))
        nea <- toupper(as.character(nea))
        ref <- toupper(as.character(ref))
        alt <- toupper(as.character(alt))
        pair <- ((ea == ref & nea == alt) | (ea == alt & nea == ref))
        one <- (ea == ref | ea == alt | nea == ref | nea == alt)
        fifelse(pair, 2L, fifelse(one, 1L, 0L))
    }
    rd_gwas <- function(f, keep_snps = NULL) {
        h <- names(fread(f, nrows = 0))
        cc <- c(SNP = pick1(h, c("SNP", "rsid", "RSID", "ID")), EA = pick1(h, c("EA", "A1", "effect_allele", 
            "ALT", "alt")), NEA = pick1(h, c("NEA", "A2", "other_allele", "REF", "ref")), BETA = pick1(h, 
            c("BETA", "beta", "b")))
        if (anyNA(cc)) 
            stop("Missing GWAS columns in ", f, ": ", paste(h, collapse = ","))
        if (!is.null(keep_snps)) {
            keep_snps <- unique(as.character(keep_snps[!is.na(keep_snps) & nzchar(keep_snps)]))
            if (!length(keep_snps)) 
                return(data.table(SNP = character(), EA = character(), NEA = character(), BETA = numeric()))
            tmp <- tempfile()
            fwrite(data.table(SNP = keep_snps), tmp, col.names = FALSE, sep = "\t")
            idx <- match(unname(cc), h)
            cmd <- sprintf("gzip -dc %s | awk 'BEGIN{FS=OFS=\"\\t\"} NR==FNR{k[$1]=1; next} FNR==1{next} ($%d in k){print $%d,$%d,$%d,$%d}' %s -", 
                shQuote(f), idx[1], idx[1], idx[2], idx[3], idx[4], shQuote(tmp))
            x <- fread(cmd = cmd, col.names = names(cc), showProgress = FALSE)
            unlink(tmp)
        }
        else {
            x <- fread(f, select = unname(cc), showProgress = FALSE)
            setnames(x, unname(cc), names(cc))
        }
        unique(x[, .(SNP = as.character(SNP), EA = toupper(EA), NEA = toupper(NEA), BETA = as.numeric(BETA))], 
            by = "SNP")
    }
    pvar_files <- function(dirmod, chr) {
        if (chr == 23L) {
            if (ref_pop == "ALL") paste0(dirmod, "/chrX.pvar")
            else unique(c(paste0(dirmod, "/", ref_pop, ".male.chrX.par.pvar"), paste0(dirmod, "/", ref_pop, ".male.chrX.nonPar.pvar"), 
                paste0(dirmod, "/", ref_pop, ".chrX.pvar")))
        } else {
            if (ref_pop == "ALL") paste0(dirmod, "/chr", chr, ".pvar")
            else paste0(dirmod, "/", ref_pop, ".chr", chr, ".pvar")
        }
    }
    pvar_ids <- function(f, n = 2000) {
        if (!file.exists(f)) 
            return(character())
        cmd <- sprintf("awk '/^#/{next} $3!=\".\" && $3!=\"\"{print $3; if(++n==%d) exit}' %s", 
            as.integer(n), shQuote(f))
        system(cmd, intern = TRUE)
    }
    pvar_mode <- function(fs) {
        id_mode(unlist(lapply(fs, pvar_ids), use.names = FALSE))
    }
    rd_pvar <- function(dirmod, lead) {
        out <- list()
        k <- 0L
        for (chr in sort(unique(lead$lead_chr))) {
            fs <- pvar_files(dirmod, chr)
            fs <- fs[file.exists(fs)]
            if (!length(fs)) 
                stop("No 1000G .pvar found for chr ", chr, " under ", dirmod)
            key <- unique(lead[lead_chr == chr, .(lead_bp, lead_snp)])
            tmp <- tempfile()
            fwrite(key, tmp, col.names = FALSE, sep = "\t")
            for (f in fs) {
                awk_script <- paste0(
                  "BEGIN{FS=OFS=\"\\t\"} ",
                  "NR==FNR{p[$1]=1; id[$2]=1; next} ",
                  "/^##/{next} ",
                  "/^#CHROM/{print \"CHR\",$2,$3,$4,$5; next} ",
                  "FNR==1 && ($1==\"CHR\" || $1==\"CHROM\")",
                  "{print \"CHR\",$2,$3,$4,$5; next} ",
                  "(($2 in p) || ($3 in id))",
                  "{gsub(/^chr/,\"\",$1); print $1,$2,$3,$4,$5}"
                )
                cmd <- sprintf("awk %s %s %s", shQuote(awk_script),
                  shQuote(tmp), shQuote(f))
                x <- fread(cmd = cmd, header = TRUE, showProgress = FALSE)
                if (nrow(x) <= 0) 
                  next
                cc <- c(CHR = pick1(names(x), c("CHR", "CHROM")), POS = pick1(names(x), c("POS", 
                  "BP")), ID = pick1(names(x), c("ID", "SNP")), REF = pick1(names(x), c("REF")), 
                  ALT = pick1(names(x), c("ALT")))
                if (anyNA(cc)) 
                  stop("Bad pvar columns after skipping ## lines: ", f, "; names=", paste(names(x), 
                    collapse = ","))
                k <- k + 1L
                out[[k]] <- x[, .(CHR = chr_int(get(cc["CHR"])), POS = as.integer(get(cc["POS"])), 
                  ID = as.character(get(cc["ID"])), REF = toupper(get(cc["REF"])), ALT = toupper(get(cc["ALT"])), 
                  pvar_file = f)]
            }
            unlink(tmp)
        }
        if (!length(out)) 
            return(data.table(CHR = integer(), POS = integer(), ID = character(), REF = character(), 
                ALT = character(), pvar_file = character()))
        unique(rbindlist(out, fill = TRUE)[!is.na(CHR) & !is.na(POS)])
    }
    align_1kg <- function(z, dirmod, lead0, tr) {
        z[, `:=`(rid, .I)]
        if (is.null(dirmod)) {
            z[, `:=`(lead_snp0 = lead_snp, pvar_id = lead_snp, match_type_1kg = "not_checked", 
                pvar_file = NA_character_, pvar_bp = lead_bp)]
            return(z)
        }
        fs <- unique(unlist(lapply(sort(unique(z$lead_chr)), function(chr) pvar_files(dirmod, 
            chr))))
        fs <- fs[file.exists(fs)]
        mode <- pvar_mode(fs)
        cojo_mode <- id_mode(z$lead_snp)
        p <- rd_pvar(dirmod, z)
        fwrite(data.table(source = c("COJO_lead", "1000G_pvar"), id_mode = c(cojo_mode, mode), 
            example = c(paste(head(z$lead_snp, 8), collapse = ","), paste(head(p$ID, 8), collapse = ","))), 
            paste0(lead0, "/", tr, ".id_sanity.tsv"), sep = "\t")
        if (!mode %chin% c("rsid", "chrpos")) {
            message("WARN: 1000G pvar ID mode is ", mode, " for ", tr,
                "; continuing with position/allele matching and downstream --set-missing-var-ids.")
        }
        z[, `:=`(lead_snp0 = lead_snp, lead_bp0 = lead_bp, pvar_id = NA_character_, pvar_bp = NA_integer_, 
            pvar_file = NA_character_, match_type_1kg = NA_character_)]
        p[, `:=`(pvar_rank, fifelse(grepl("male\\.chrX\\.(par|nonPar)\\.pvar$", pvar_file), 1L, 
            2L))]
        m <- data.table(rid = integer(), pvar_id = character(), pvar_bp = integer(), pvar_file = character(), 
            match_type_1kg = character(), priority = integer())
        idm <- merge(z[usable_id(lead_snp), .(rid, lead_chr, lead_snp, lead_bp, EA, NEA)], p[usable_id(ID)], 
            by.x = "lead_snp", by.y = "ID", allow.cartesian = TRUE)
        idm <- idm[lead_chr == CHR]
        if (nrow(idm)) {
            idm[, `:=`(allele_match, allele_score(EA, NEA, REF, ALT))]
            idm <- idm[allele_match > 0]
            idm[, `:=`(pos_same, as.integer(lead_bp == POS))]
            if (nrow(idm)) {
                setorder(idm, rid, -allele_match, -pos_same, pvar_rank)
                idm <- idm[, .SD[1], by = rid]
                m <- rbind(m, idm[, .(rid, pvar_id = lead_snp, pvar_bp = POS, pvar_file, match_type_1kg = fifelse(pos_same == 
                  1L, "ID_allele_POS_same", "ID_allele_POS_from_pvar"), priority = 0L)], fill = TRUE)
            }
        }
        unmatched <- z[!(rid %in% m$rid)]
        if (nrow(unmatched)) {
            pm <- merge(unmatched[, .(rid, lead_chr, lead_snp, lead_bp, EA, NEA)], p, by.x = c("lead_chr", 
                "lead_bp"), by.y = c("CHR", "POS"), allow.cartesian = TRUE)
            if (nrow(pm)) {
                pm[, `:=`(allele_match, allele_score(EA, NEA, REF, ALT))]
                pm <- pm[allele_match > 0]
                if (nrow(pm)) {
                  setorder(pm, rid, -allele_match, pvar_rank)
                  pm <- pm[, .SD[1], by = rid]
                  pm[, `:=`(out_id, fifelse(usable_id(ID), ID, lead_snp))]
                  m <- rbind(m, pm[, .(rid, pvar_id = out_id, pvar_bp = lead_bp, pvar_file, match_type_1kg = fifelse(usable_id(ID), "POS_allele_ID_mismatch", "POS_allele_missing_ID"), 
                    priority = 1L)], fill = TRUE)
                }
            }
        }
        h <- m
        h <- h[!is.na(pvar_id) & nzchar(pvar_id) & pvar_id != "."]
        if (nrow(h)) {
            setorder(h, rid, priority)
            h <- h[, .SD[1], by = rid]
            z[h, `:=`(pvar_id = i.pvar_id, pvar_bp = i.pvar_bp, pvar_file = i.pvar_file, match_type_1kg = i.match_type_1kg), 
                on = .(rid)]
        }
        z[, `:=`(lead_snp, fifelse(!is.na(pvar_id), pvar_id, lead_snp))]
        z[, `:=`(lead_bp, fifelse(!is.na(pvar_bp), pvar_bp, lead_bp))]
        z[]
    }
    dirgwas <- arg(c("--dirgwas", "-dirgwas"))
    dirout <- arg(c("--dirout", "-dirout"))
    dirmod <- arg(c("--dirmod", "-dirmod"))
    ref_pop <- arg(c("--ref_pop", "-ref_pop"))
    if (is.null(ref_pop)) ref_pop <- "ALL"
    traits <- arg(c("--traits", "-traits"), TRUE)
    if (is.null(dirgwas) || is.null(dirout) || is.null(traits)) 
        stop("Usage: Rscript locus.R prep_input --dirgwas DIR --dirout DIR --dirmod 1KG_DIR --traits bald bald12 ...")
    traits <- unlist(strsplit(paste(traits, collapse = " "), "[, ]+"))
    traits <- traits[nzchar(traits)]
    lead0 <- paste0(normalizePath(dirout, winslash = "/", mustWork = FALSE), "/lead")
    dir.create(lead0, recursive = TRUE, showWarnings = FALSE)
    all_match <- list()
    prep_version <- "prep_input_v3_id_first_pos_if_missing_id"
    for (tr in traits) {
        cj_f <- paste0(dirgwas, "/cojo/", tr, "/", tr, ".jma.cojo")
        gw_f <- paste0(dirgwas, "/clean/", tr, "/", tr, ".gz")
        cj <- fread(cj_f)
        cc <- c(chr = pick1(names(cj), c("Chr", "CHR", "chr", "#CHROM", "CHROM")), bp = pick1(names(cj), 
            c("bp", "BP", "POS", "pos")))
        if (anyNA(cc) || !"SNP" %chin% names(cj)) 
            stop("COJO needs Chr/SNP/bp: ", cj_f)
        lead <- unique(cj[, .(lead_chr = chr_int(get(cc["chr"])), lead_snp = as.character(SNP), 
            lead_bp = as.integer(get(cc["bp"])))])[lead_chr %between% c(1L, 23L) & !is.na(lead_bp)]
        g <- rd_gwas(gw_f, lead$lead_snp)
        z <- merge(lead, g, by.x = "lead_snp", by.y = "SNP", all.x = TRUE, sort = FALSE)
        z[, `:=`(fail, is.na(EA) | !nzchar(EA) | is.na(NEA) | !nzchar(NEA) | is.na(BETA))]
        fail <- z[fail == TRUE, .(trait = tr, lead_chr, lead_snp, lead_bp, effect_allele = EA, 
            other_allele = NEA, beta = BETA, fail_reason = "not_in_gwas_or_missing_allele_beta")]
        z <- z[fail == FALSE]
        if (nrow(z)) 
            z <- align_1kg(z, dirmod, lead0, tr)
        else z[, `:=`(lead_snp0 = lead_snp, lead_bp0 = lead_bp, pvar_id = NA_character_, pvar_bp = NA_integer_, 
            pvar_file = NA_character_, match_type_1kg = NA_character_)]
        all_match[[tr]] <- z[, .(trait = tr, lead_chr, lead_snp0, lead_bp0, lead_snp, lead_bp, 
            EA, NEA, BETA, pvar_id, pvar_bp, pvar_file, match_type_1kg)]
        fail <- rbind(fail, z[is.na(pvar_id), .(trait = tr, lead_chr, lead_snp = lead_snp0, lead_bp = lead_bp0, 
            effect_allele = EA, other_allele = NEA, beta = BETA, fail_reason = "not_in_1000G_pvar_or_allele_mismatch")], 
            fill = TRUE)
        z <- z[!is.na(pvar_id)][, .(trait = tr, lead_chr, lead_snp, lead_bp, effect_allele = EA, 
            other_allele = NEA, beta = BETA, risk_allele = fifelse(BETA >= 0, EA, NEA), lead_snp0, 
            match_type_1kg)]
        fwrite(z[, .(lead_chr, lead_snp, lead_bp)], paste0(lead0, "/", tr, ".lead.3col"), sep = "\t")
        if (nrow(fail)) 
            fwrite(fail, paste0(lead0, "/", tr, ".lead.fail.tsv"), sep = "\t")
        else unlink(paste0(lead0, "/", tr, ".lead.fail.tsv"))
        fwrite(z, paste0(lead0, "/", tr, ".lead.assoc"), sep = "\t")
    }
    fwrite(rbindlist(lapply(traits, function(tr) fread(paste0(lead0, "/", tr, ".lead.assoc"))), 
        fill = TRUE), paste0(lead0, "/lead.assoc"), sep = "\t")
    fwrite(rbindlist(all_match, fill = TRUE), paste0(lead0, "/lead.match.tsv"), sep = "\t")
    writeLines(prep_version, paste0(lead0, "/prep.version"))
}


run_prep_archaic <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) script_args else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table)
    args <- commandArgs(TRUE)
    if(length(args) < 4) {
    	stop("Usage: Rscript input_prep_archaic.R kg.vcf.gz archaic.raw.vcf.gz out.vcf.gz sample_name")
    }

    kg_vcf  <- args[1]
    arc_vcf <- args[2]
    out_vcf <- args[3]
    sample  <- args[4]

    if(!file.exists(kg_vcf))  stop("Missing kg_vcf: ", kg_vcf)
    if(!file.exists(arc_vcf)) stop("Missing arc_vcf: ", arc_vcf)

    bgzip <- Sys.which("bgzip")
    tabix <- Sys.which("tabix")
    if(bgzip == "") stop("bgzip not found in PATH")
    if(tabix == "") stop("tabix not found in PATH")

    message("kg_vcf  = ", kg_vcf)
    message("arc_vcf = ", arc_vcf)
    message("out_vcf = ", out_vcf)
    message("sample  = ", sample)

    read_query <- function(vcf, fmt, cols, label){
    	tmp <- tempfile()
    	cmd <- paste("bcftools query -f", shQuote(fmt), shQuote(vcf), ">", shQuote(tmp))
    	rc <- system(cmd)
    	if(rc != 0) stop("bcftools query failed for ", label, ": ", vcf)
    	if(!file.exists(tmp) || file.info(tmp)$size == 0){
    		unlink(tmp)
    		return(as.data.table(setNames(rep(list(character()), length(cols)), cols)))
    	}
    	x <- fread(tmp, header=FALSE, col.names=cols, showProgress=FALSE)
    	unlink(tmp)
    	x
    }

    kg <- read_query(
    	kg_vcf,
    	"%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n",
    	c("CHR", "POS", "ID", "REF_kg", "ALT_kg"),
    	"1000G template"
    )

    arc <- read_query(
    	arc_vcf,
    	"%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n",
    	c("CHR", "POS", "REF_arc", "ALT_arc", "GT_arc"),
    	"archaic raw"
    )

    if(nrow(kg) == 0) stop("kg_vcf has 0 variants: ", kg_vcf)
    if(nrow(arc) == 0) warning("arc_vcf has 0 records: ", arc_vcf)

    kg[, `:=`(
    	idx = .I,
    	CHR = as.character(CHR),
    	POS = as.integer(POS),
    	ID = fifelse(is.na(ID) | ID == "", ".", as.character(ID)),
    	REF_kg = toupper(REF_kg),
    	ALT_kg = toupper(ALT_kg)
    )]

    arc[, `:=`(
    	CHR = as.character(CHR),
    	POS = as.integer(POS),
    	REF_arc = toupper(REF_arc),
    	ALT_arc = toupper(ALT_arc),
    	GT_arc = gsub("\\|", "/", as.character(GT_arc))
    )]

    x <- merge(kg, arc, by = c("CHR", "POS"), all.x = TRUE, sort = FALSE)
    setorder(x, idx)

    map_gt <- function(ref_kg, alt_kg, ref_arc, alt_arc, gt_arc){
    	if(is.na(ref_arc) || is.na(gt_arc) || gt_arc %chin% c(".", "./.", ".|.")) return("./.")
    	if(is.na(alt_arc) || alt_arc == "") alt_arc <- "."
    	alts <- if(alt_arc == ".") character(0) else strsplit(alt_arc, ",", fixed=TRUE)[[1]]
    	alleles <- c(ref_arc, alts)
    	gt_arc <- gsub("\\|", "/", gt_arc)
    	if(gt_arc == "0" || gt_arc == "1") gt_arc <- paste(gt_arc, gt_arc, sep="/")
    	g <- strsplit(gt_arc, "/", fixed=TRUE)[[1]]
    	if(length(g) != 2L || any(g == ".")) return("./.")
    	a <- suppressWarnings(as.integer(g)) + 1L
    	if(any(is.na(a)) || any(a < 1L) || any(a > length(alleles))) return("./.")
    	b <- alleles[a]
    	out <- ifelse(b == ref_kg, "0", ifelse(b == alt_kg, "1", NA_character_))
    	if(any(is.na(out))) return("./.")
    	paste(out, collapse="/")
    }

    x[, GT := mapply(map_gt, REF_kg, ALT_kg, REF_arc, ALT_arc, GT_arc)]

    # Keep only standard A/C/G/T SNP template from 1000G
    x <- x[REF_kg %chin% c("A", "C", "G", "T") & ALT_kg %chin% c("A", "C", "G", "T")]

    vcf <- x[, .(
    	CHR,
    	POS,
    	ID,
    	REF = REF_kg,
    	ALT = ALT_kg,
    	QUAL = ".",
    	FILTER = ".",
    	INFO = ".",
    	FORMAT = "GT",
    	GT
    )]

    dir.create(dirname(out_vcf), recursive = TRUE, showWarnings = FALSE)

    tmp_vcf <- tempfile(fileext = ".vcf")
    cat("##fileformat=VCFv4.2\n", file = tmp_vcf)
    cat("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", file = tmp_vcf, append = TRUE)
    cat(paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample, sep = "\t"), "\n", file = tmp_vcf, append = TRUE)
    fwrite(vcf, tmp_vcf, sep = "\t", append = TRUE, col.names = FALSE)

    if(file.exists(out_vcf)) file.remove(out_vcf)
    if(file.exists(paste0(out_vcf, ".tbi"))) file.remove(paste0(out_vcf, ".tbi"))

    cmd_bgzip <- paste(shQuote(bgzip), "-c", shQuote(tmp_vcf), ">", shQuote(out_vcf))
    rc1 <- system(cmd_bgzip)
    if(rc1 != 0) stop("bgzip failed: ", cmd_bgzip)
    if(!file.exists(out_vcf) || file.info(out_vcf)$size == 0) stop("bgzip did not create output: ", out_vcf)

    rc2 <- system2(tabix, c("-f", "-p", "vcf", out_vcf))
    if(rc2 != 0) stop("tabix failed for: ", out_vcf)
    if(!file.exists(paste0(out_vcf, ".tbi"))) stop("tabix did not create index: ", paste0(out_vcf, ".tbi"))

    message("Wrote: ", out_vcf)
    message("n_template=", nrow(kg))
    message("n_archaic_raw=", nrow(arc))
    message("n_output=", nrow(vcf))
    message("GT_0_0=", sum(vcf$GT == "0/0"))
    message("GT_0_1=", sum(vcf$GT %chin% c("0/1", "1/0")))
    message("GT_1_1=", sum(vcf$GT == "1/1"))
    message("GT_missing=", sum(vcf$GT == "./."))
}
run_make_hap <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table)
    args0 <- commandArgs(TRUE)
    has_root <- length(args0) && grepl("[/\\\\]", args0[1])
    root <- normalizePath(if (has_root) 
        args0[1]
    else Sys.getenv("BALD_RES", "/data/sph-zhaor/analysis/bald/res"), winslash = "/", mustWork = FALSE)
    args <- if (has_root) 
        args0[-1]
    else args0
    mat0 <- file.path(root, "mat")
    hap0 <- file.path(root, "hap")
    ld0 <- file.path(root, "ld")
    corevcf0 <- file.path(root, "coreVcf")
    report0 <- file.path(root, "report")
    riskf <- file.path(root, "lead", "lead.assoc")
    dir.create(hap0, recursive = TRUE, showWarnings = FALSE)
    dir.create(report0, recursive = TRUE, showWarnings = FALSE)
    base_set <- c("A", "C", "G", "T")
    pick_lineage_method <- Sys.getenv("PICK_LINEAGE_METHOD", "simple")
    hap_filter_method <- Sys.getenv("HAP_FILTER_METHOD", "simple")
    yri_freq_th <- suppressWarnings(as.numeric(Sys.getenv("YRI_FREQ_TH", "0.05")))
    if (!is.finite(yri_freq_th) || yri_freq_th < 0) 
        stop("YRI_FREQ_TH must be a non-negative number")
    diagnostic_delta_th <- suppressWarnings(as.numeric(Sys.getenv("DIAGNOSTIC_DELTA_TH", "0.5")))
    if (!is.finite(diagnostic_delta_th)) 
        stop("DIAGNOSTIC_DELTA_TH must be a number")
    if (pick_lineage_method == "hpc") pick_lineage_method <- "simple"
    if (hap_filter_method == "hpc") hap_filter_method <- "simple"
    pick1 <- function(nm, cand) {
        x <- cand[cand %chin% nm]
        if (length(x)) 
            x[1]
        else NA_character_
    }
    arch_label <- function(f) {
        x <- tools::file_path_sans_ext(basename(f))
        paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    }
    arch_files <- function(d0) {
        fs <- list.files(d0, pattern = "\\.tsv$", full.names = TRUE)
        fs <- fs[!basename(fs) %in% c("kg.tsv", "kg.samples.tsv")]
        setNames(fs, arch_label(fs))
    }
    lineage_group <- function(x) {
        y <- tolower(x)
        fifelse(grepl("vindija|altai|chagyr|neand", y), "Neanderthal", fifelse(grepl("denisova|denisovan", 
            y), "Denisovan", "Archaic"))
    }
    lineage_group_fun <- lineage_group
    read_risk <- function(f) {
        x <- fread(f)
        nm <- names(x)
        tr_col <- pick1(nm, c("trait", "Trait"))
        sn_col <- pick1(nm, c("lead_snp", "SNP", "snp", "rsid", "RSID", "rsID", "ID", "id"))
        rk_col <- pick1(nm, c("risk_allele", "RiskAllele", "riskAllele"))
        ea_col <- pick1(nm, c("effect_allele", "EA", "ea", "A1", "a1", "ALT", "alt", "Allele1", 
            "allele1", "tested_allele", "TestedAllele"))
        oa_col <- pick1(nm, c("other_allele", "OA", "oa", "NEA", "nea", "A2", "a2", "REF", "ref", 
            "Allele2", "allele2", "non_effect_allele", "Non_Effect_Allele"))
        bt_col <- pick1(nm, c("beta", "BETA", "Beta", "effect", "Effect", "Estimate", "estimate", 
            "B", "b"))
        or_col <- pick1(nm, c("OR", "or", "OddsRatio", "odds_ratio", "oddsratio"))
        if (is.na(tr_col) || is.na(sn_col)) 
            stop("lead.assoc must contain trait and lead_snp")
        keep <- unique(na.omit(c(tr_col, sn_col, rk_col, ea_col, oa_col, bt_col, or_col)))
        x <- x[, ..keep]
        map <- c(trait = tr_col, lead_snp = sn_col, risk_allele = rk_col, effect_allele = ea_col, 
            other_allele = oa_col, beta = bt_col, OR = or_col)
        for (nm0 in names(map)) if (!is.na(map[nm0])) 
            setnames(x, map[nm0], nm0)
        if (!"risk_allele" %in% names(x)) {
            if (!all(c("effect_allele", "other_allele") %in% names(x))) 
                stop("lead.assoc missing risk_allele and effect_allele/other_allele")
            if ("beta" %in% names(x)) 
                x[, `:=`(risk_allele, fifelse(as.numeric(beta) >= 0, effect_allele, other_allele))]
            else if ("OR" %in% names(x)) 
                x[, `:=`(risk_allele, fifelse(as.numeric(OR) >= 1, effect_allele, other_allele))]
            else stop("lead.assoc missing risk_allele and beta/OR")
        }
        x[, `:=`(trait = as.character(trait), lead_snp = as.character(lead_snp), risk_allele = toupper(as.character(risk_allele)))]
        unique(x[, .(trait, lead_snp, risk_allele)])
    }
    risk <- read_risk(riskf)
    read_panel <- function() {
        env_panel <- unique(c(Sys.getenv("SAMPLE_FILE", ""), Sys.getenv("sample_file", "")))
        env_panel <- env_panel[nzchar(env_panel)]
        cand <- c(env_panel, file.path(root, "integrated_call_samples_v3.20130502.ALL.panel"), 
            file.path(root, "1kg.v3.sample.txt"), "/mnt/d/files/1kg.v3.sample.txt", 
            "/mnt/d/refGen/1kg/ucsc/vcf/samples_v3.ALL.panel",
            "/data/sph-zhaor/refGen/1kg_phase3/integrated_call_samples_v3.20130502.ALL.panel", 
            "/data/sph-zhaor/analysis/bald/data/1kg/integrated_call_samples_v3.20130502.ALL.panel")
        f <- cand[file.exists(cand)][1]
        if (is.na(f)) {
            message("WARN: no 1000G sample panel found for YRI filter; set SAMPLE_FILE to samples_v3.ALL.panel")
            return(data.table(sample = character(), pop = character(), super_pop = character()))
        }
        x <- fread(f, fill = TRUE)
        setnames(x, names(x)[1:3], c("sample", "pop", "super_pop"))
        x <- x[, .(sample = as.character(sample), pop = as.character(pop), super_pop = as.character(super_pop))]
        message("Sample panel for YRI filter: ", f, "; samples=", nrow(x), "; YRI=", sum(x$pop == "YRI", na.rm = TRUE))
        x
    }
    panel <- read_panel()
    afr_pops <- c("YRI", "LWK", "GWD", "MSL", "ESN")
    afr_samples <- setNames(lapply(afr_pops, function(pp) panel[pop == pp, sample]), afr_pops)
    yri_samples <- afr_samples[["YRI"]]
    read_vcf_samples <- function(vcf) {
        if (!file.exists(vcf)) return(character(0))
        con <- gzfile(vcf, "rt")
        on.exit(close(con), add = TRUE)
        repeat {
            x <- readLines(con, n = 1000, warn = FALSE)
            if (!length(x)) return(character(0))
            h <- grep("^#CHROM", x, value = TRUE)
            if (length(h)) return(strsplit(h[1], "\t", fixed = TRUE)[[1]][-(1:9)])
        }
    }
    chr_id <- function(x) {
        x <- toupper(sub("^CHR", "", as.character(x)))
        x[x == "X"] <- "23"
        suppressWarnings(as.integer(x))
    }
    parse_id <- function(id) {
        z <- strsplit(id, "\\.")[[1]]
        list(chr = chr_id(z[1]), snp = paste(z[2:(length(z) - 1)], collapse = "."), bp = suppressWarnings(as.integer(z[length(z)])))
    }
    gt_hap <- function(x) {
        x <- gsub("\\|", "/", x)
        sp <- tstrsplit(x, "/", fixed = TRUE)
        list(a = sp[[1]], b = sp[[2]])
    }
    roman_n <- function(x) vapply(x, function(i) as.character(as.roman(i)), character(1))
    to_base <- function(h, ref, alt) {
        out <- rep("N", length(h))
        out[h == "0"] <- ref[h == "0"]
        out[h == "1"] <- alt[h == "1"]
        out
    }
    arch_base <- function(gt, ref, alt) {
        gt <- gsub("\\|", "/", gt)
        sp <- strsplit(gt, "/", fixed = TRUE)
        vapply(seq_along(sp), function(i) {
            z <- sp[[i]]
            z <- z[z != "."]
            if (length(z) != 2L || z[1] != z[2]) 
                return(NA_character_)
            a <- z[1]
            alts <- if (alt[i] %in% c(".", "")) 
                character(0)
            else strsplit(alt[i], ",", fixed = TRUE)[[1]]
            if (a == "0") 
                return(ref[i])
            ai <- suppressWarnings(as.integer(a))
            if (is.na(ai) || ai < 1L || ai > length(alts)) 
                return(NA_character_)
            alts[ai]
        }, character(1))
    }
    score_seq <- function(x, y) {
        ok <- x != "N" & y != "N" & !is.na(x) & !is.na(y)
        c(match = sum(x[ok] == y[ok]), n = sum(ok))
    }
    ils_p <- function(size_bp, recomb_cM_Mb = 0.53, split_years = 550000, archaic_age_years = 50000, 
        gen_years = 29) {
        r <- recomb_cM_Mb * 1e-08
        L <- 1/(r * ((split_years/gen_years) + ((split_years - archaic_age_years)/gen_years)))
        1 - pgamma(size_bp, shape = 2, rate = 1/L)
    }
    read_mat <- function(f, type) {
        x <- fread(f, header = FALSE)
        if (!nrow(x)) 
            return(NULL)
        if (type == "kg") 
            setnames(x, c("chr", "pos", "ref", "alt", "aa", paste0("s", seq_len(ncol(x) - 5))))
        else setnames(x, c("chr", "pos", "ref", "alt", "gt"))
        x[, `:=`(chr, sub("^chr", "", as.character(chr), ignore.case = TRUE))]
        x[, `:=`(ref = toupper(trimws(ref)), alt = toupper(trimws(alt)))]
        if ("aa" %in% names(x)) {
            x[, `:=`(aa, toupper(sub("\\|.*$", "", trimws(aa))))]
        }
        x
    }
    read_arch <- function(f, lab) {
        x <- read_mat(f, "arch")
        if (is.null(x)) 
            return(NULL)
        x[, `:=`(allele, arch_base(gt, ref, alt))]
        x <- x[, .(chr, pos, allele)]
        setnames(x, "allele", lab)
        x
    }
    merge_arch <- function(kg, fs) {
        x <- kg
        keep <- character(0)
        for (lab in names(fs)) {
            a <- read_arch(fs[[lab]], lab)
            if (is.null(a)) 
                next
            x <- merge(x, a, by = c("chr", "pos"), all.x = TRUE)
            keep <- c(keep, lab)
        }
        setorder(x, pos)
        list(x = x, arch = keep)
    }
    comp_allele <- function(x) {
        x <- toupper(x)
        ifelse(nchar(x) == 1L, chartr("ACGT", "TGCA", x), x)
    }
    harmonize_lead <- function(a, ref, alt) {
        a <- toupper(trimws(a))
        if (!nzchar(a)) 
            return(NA_character_)
        if (a == ref || a == alt) 
            return(a)
        b <- comp_allele(a)
        if (b == ref || b == alt) 
            b
        else NA_character_
    }
    lead_risk_base <- function(a, ref, alt) {
        a <- toupper(trimws(a))
        if (a == "ALT") return(alt)
        if (a == "REF") return(ref)
        harmonize_lead(a, ref, alt)
    }
    pick_lineage_local <- function(arch_stat, p, p_cut = 0.1, min_snp = 2L, min_prop = 0.5) {
        if (!is.finite(p) || p >= p_cut || !nrow(arch_stat)) 
            return(list(best_lineage = NA_character_, keep_arch = character(0)))
        z <- copy(arch_stat)
        z[, `:=`(group, lineage_group(archaic))]
        z[, `:=`(ok, n_compared_risk >= min_snp & n_match_risk >= min_snp & prop_match_risk >= 
            min_prop)]
        z <- z[ok == TRUE]
        if (!nrow(z)) 
            return(list(best_lineage = NA_character_, keep_arch = character(0)))
        grp <- z[, .(prop = max(prop_match_risk, na.rm = TRUE), match = max(n_match_risk, na.rm = TRUE)), 
            by = group][order(-prop, -match)][1, group]
        list(best_lineage = grp, keep_arch = z[group == grp][order(-prop_match_risk, -n_match_risk), 
            archaic])
    }
    pick_lineage_simple <- function(arch_stat, p, p_cut = 0.1) {
        if (!is.finite(p) || p >= p_cut || !nrow(arch_stat)) 
            return(list(best_lineage = NA_character_, keep_arch = character(0)))
        z <- copy(arch_stat)
        if (!"lineage_group" %in% names(z)) z[, `:=`(lineage_group, lineage_group_fun(archaic))]
        g <- z[lineage_group %chin% c("Neanderthal", "Denisovan"), .(n_compared_risk = sum(n_compared_risk, 
            na.rm = TRUE), n_match_risk = sum(n_match_risk, na.rm = TRUE)), by = lineage_group]
        g[, `:=`(prop_match_risk, fifelse(n_compared_risk > 0, n_match_risk/n_compared_risk, 
            0))]
        nean <- g[lineage_group == "Neanderthal"]
        den <- g[lineage_group == "Denisovan"]
        ns <- if (nrow(nean)) nean$prop_match_risk[1] else 0
        ds <- if (nrow(den)) den$prop_match_risk[1] else 0
        nm <- if (nrow(nean)) nean$n_match_risk[1] else 0L
        dm <- if (nrow(den)) den$n_match_risk[1] else 0L
        if (ns <= 0 && ds <= 0) 
            return(list(best_lineage = NA_character_, keep_arch = character(0)))
        if (ns > ds || (ns == ds && nm > dm)) 
            return(list(best_lineage = "Neanderthal", keep_arch = z[lineage_group == "Neanderthal" & 
                n_match_risk > 0, archaic]))
        if (ds > ns || (ds == ns && dm > nm)) 
            return(list(best_lineage = "Denisovan", keep_arch = z[lineage_group == "Denisovan" & 
                n_match_risk > 0, archaic]))
        list(best_lineage = NA_character_, keep_arch = character(0))
    }
    pick_lineage <- function(arch_stat, p) {
        if (pick_lineage_method == "simple") pick_lineage_simple(arch_stat, p) else pick_lineage_local(arch_stat, 
            p)
    }
    write_one <- function(x, f) {
        if (nrow(x)) 
            fwrite(x, f, sep = "\t")
    }
    clean_one <- function(tr0, id0) {
        od <- file.path(hap0, tr0)
        dir.create(od, recursive = TRUE, showWarnings = FALSE)
        unlink(file.path(od, paste0(id0, c(".hap.tsv", ".site.tsv", ".arch.tsv", ".region.tsv", 
            ".core.tsv"))), force = TRUE)
    }
    process_one <- function(tr0, id0) {
        od <- file.path(hap0, tr0)
        dir.create(od, recursive = TRUE, showWarnings = FALSE)
        p0 <- parse_id(id0)
        chr0 <- p0$chr
        snp0 <- p0$snp
        bp0 <- p0$bp
        ldf <- file.path(ld0, tr0, "ld.tsv")
        bdf <- file.path(ld0, tr0, "block.tsv")
        if (!file.exists(ldf) || !file.exists(bdf)) 
            return(invisible(NULL))
        ld <- fread(ldf)
        ld[, `:=`(lead_chr = as.integer(lead_chr), lead_bp = as.integer(lead_bp), pos = as.integer(pos), 
            R2 = as.numeric(R2))]
        ld <- ld[!is.na(lead_chr) & !is.na(lead_bp) & !is.na(pos)]
        blk <- fread(bdf)
        blk[, `:=`(lead_chr = as.integer(lead_chr), lead_bp = as.integer(lead_bp), start = as.integer(start), 
            end = as.integer(end), n = as.integer(n), size_bp = as.integer(size_bp))]
        blk <- blk[lead_chr == chr0 & lead_snp == snp0 & lead_bp == bp0]
        if (!nrow(blk)) 
            return(invisible(NULL))
        blk <- blk[1]
        rk_raw <- risk[trait == tr0 & lead_snp == snp0, risk_allele]
        rk_raw <- if (length(rk_raw)) 
            rk_raw[1]
        else NA_character_
        p <- ils_p(as.numeric(blk$size_bp))
        region <- data.table(trait = tr0, id = id0, lead_chr = chr0, lead_snp = snp0, lead_bp = bp0, 
            core_start = blk$start, core_end = blk$end, n_ld_snp = blk$n, core_size_bp = blk$size_bp, 
            p_ils = p, best_lineage = NA_character_, matched_archaics = "")
        site <- data.table(trait = tr0, id = id0, n_site_raw = 0L, n_site_called = 0L, n_site_keep = 0L, 
            n_hap_raw = 0L, n_hap_keep = 0L, best_lineage = NA_character_, matched_archaics = "")
        d0 <- file.path(mat0, tr0, id0)
        kgf <- file.path(d0, "kg.tsv")
        afs <- arch_files(d0)
        if (!file.exists(kgf) || !length(afs)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        kg <- read_mat(kgf, "kg")
        if (is.null(kg)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        ma <- merge_arch(kg, afs)
        x <- ma$x
        arch_names <- ma$arch
        if (!length(arch_names)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        site[, `:=`(n_site_raw = nrow(x), n_site_called = nrow(x))]
        core_pos <- ld[lead_chr == chr0 & lead_snp == snp0 & lead_bp == bp0, sort(unique(pos))]
        xcore <- x[pos %in% core_pos]
        setorder(xcore, pos)
        if (!nrow(xcore)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        samp_core <- setdiff(names(xcore), c("chr", "pos", "ref", "alt", "aa", arch_names))
        Hcore <- do.call(cbind, c(lapply(samp_core, function(s) to_base(gt_hap(xcore[[s]])$a, 
            xcore$ref, xcore$alt)), lapply(samp_core, function(s) to_base(gt_hap(xcore[[s]])$b, 
            xcore$ref, xcore$alt))))
        colnames(Hcore) <- c(paste0(samp_core, "_1"), paste0(samp_core, "_2"))
        lead_i <- which(xcore$pos == bp0)
        rk <- if (length(lead_i) == 1L) 
            lead_risk_base(rk_raw, xcore$ref[lead_i], xcore$alt[lead_i])
        else NA_character_
        carry <- if (length(lead_i) == 1L && !is.na(rk)) 
            Hcore[lead_i, ] == rk
        else rep(FALSE, ncol(Hcore))
        core <- copy(xcore)[, .(trait = tr0, id = id0, lead_chr = chr0, lead_snp = snp0, lead_bp = bp0, 
            pos, ref, alt)]
        core[, `:=`(lead_risk_raw = rk_raw, lead_risk_harmonized = rk, n_lead_risk_haps = sum(carry, 
            na.rm = TRUE), risk_core_allele = NA_character_, carry_freq = NA_real_, noncarry_freq = NA_real_)]
        if (sum(carry, na.rm = TRUE) > 0L) {
            rka <- vapply(seq_len(nrow(Hcore)), function(i) {
                z <- Hcore[i, carry, drop = TRUE]
                z <- z[!is.na(z) & z != "N"]
                if (!length(z)) 
                  NA_character_
                else names(sort(table(z), decreasing = TRUE))[1]
            }, character(1))
            core[, `:=`(risk_core_allele, rka)]
            core[, `:=`(carry_freq, vapply(seq_len(nrow(Hcore)), function(i) {
                z <- Hcore[i, carry, drop = TRUE]
                z <- z[!is.na(z) & z != "N"]
                if (!length(z) || is.na(rka[i])) 
                  NA_real_
                else mean(z == rka[i])
            }, numeric(1)))]
            core[, `:=`(noncarry_freq, vapply(seq_len(nrow(Hcore)), function(i) {
                z <- Hcore[i, !carry, drop = TRUE]
                z <- z[!is.na(z) & z != "N"]
                if (!length(z) || is.na(rka[i])) 
                  NA_real_
                else mean(z == rka[i])
            }, numeric(1)))]
        }
        core[, `:=`(risk_core_side = fifelse(risk_core_allele == ref, "REF", fifelse(risk_core_allele == 
            alt, "ALT", NA_character_)), risk_freq_delta = carry_freq - noncarry_freq)]
        vcf_samples <- read_vcf_samples(file.path(corevcf0, tr0, id0, "kg.vcf.gz"))
        yri_idx <- match(yri_samples, vcf_samples)
        yri_samp <- paste0("s", yri_idx[!is.na(yri_idx)])
        yri_haps <- intersect(c(paste0(yri_samp, "_1"), paste0(yri_samp, "_2")), colnames(Hcore))
        afr_haps <- setNames(lapply(afr_pops, function(pp) {
            pop_idx <- match(afr_samples[[pp]], vcf_samples)
            pop_samp <- paste0("s", pop_idx[!is.na(pop_idx)])
            intersect(c(paste0(pop_samp, "_1"), paste0(pop_samp, "_2")), colnames(Hcore))
        }), afr_pops)
        if (length(yri_haps)) {
            yri_called <- vapply(seq_len(nrow(Hcore)), function(i) sum(Hcore[i, yri_haps] %chin% 
                base_set, na.rm = TRUE), integer(1))
            yri_count <- vapply(seq_len(nrow(Hcore)), function(i) {
                if (is.na(core$risk_core_allele[i])) return(NA_integer_)
                sum(Hcore[i, yri_haps] == core$risk_core_allele[i], na.rm = TRUE)
            }, integer(1))
        }
        else {
            yri_called <- rep(NA_integer_, nrow(Hcore))
            yri_count <- rep(NA_integer_, nrow(Hcore))
        }
        core[, `:=`(yri_called_alleles = yri_called, yri_risk_count = yri_count, 
            yri_risk_freq = fifelse(yri_called > 0, yri_count/yri_called, NA_real_))]
        for (pp in afr_pops) {
            pop_haps <- afr_haps[[pp]]
            if (length(pop_haps)) {
                pop_called <- vapply(seq_len(nrow(Hcore)), function(i) sum(Hcore[i, pop_haps] %chin%
                    base_set, na.rm = TRUE), integer(1))
                pop_count <- vapply(seq_len(nrow(Hcore)), function(i) {
                    if (is.na(core$risk_core_allele[i])) return(NA_integer_)
                    sum(Hcore[i, pop_haps] == core$risk_core_allele[i], na.rm = TRUE)
                }, integer(1))
            } else {
                pop_called <- rep(NA_integer_, nrow(Hcore))
                pop_count <- rep(NA_integer_, nrow(Hcore))
            }
            core[, (paste0(pp, "_called_alleles")) := pop_called]
            core[, (paste0(pp, "_risk_count")) := pop_count]
            core[, (paste0(pp, "_risk_freq")) := fifelse(pop_called > 0, pop_count/pop_called, NA_real_)]
        }
        arch_stat <- rbindlist(lapply(arch_names, function(an) {
            ca <- xcore[[an]]
            ok1 <- ca %chin% base_set
            ok2 <- ok1 & !is.na(core$risk_core_allele)
            z1 <- core[ok1]
            z2 <- core[ok2]
            a2 <- ca[ok2]
            data.table(trait = tr0, id = id0, lead_chr = chr0, lead_snp = snp0, lead_bp = bp0, 
                archaic = an, lineage_group = lineage_group(an), n_core_ld_snp = length(core_pos), 
                n_risk_defined = sum(!is.na(core$risk_core_allele)), n_callable = nrow(z1), n_compared_risk = nrow(z2), 
                n_match_risk = sum(a2 == z2$risk_core_allele, na.rm = TRUE), n_match_ref = sum(ca[ok1] == 
                  z1$ref, na.rm = TRUE), n_match_alt = sum(ca[ok1] == z1$alt, na.rm = TRUE), 
                prop_match_risk = fifelse(nrow(z2) > 0L, sum(a2 == z2$risk_core_allele, na.rm = TRUE)/nrow(z2), 
                  NA_real_))
        }), fill = TRUE)
        cls <- pick_lineage(arch_stat, p)
        keep_arch <- cls$keep_arch
        best_lineage_val <- cls$best_lineage
        matched_archaics_val <- paste(keep_arch, collapse = ";")
        core[, `:=`(diagnostic_archaic_allele = NA_character_, n_diagnostic_archaics = 0L,
            is_diagnostic_archaic = FALSE, yri_diagnostic_count = NA_integer_,
            yri_diagnostic_freq = NA_real_)]
        for (pp in afr_pops) {
            core[, (paste0(pp, "_diagnostic_count")) := NA_integer_]
            core[, (paste0(pp, "_diagnostic_freq")) := NA_real_]
        }
        if (length(keep_arch)) {
            diag_arch <- if (identical(best_lineage_val, "Neanderthal")) {
                intersect(c("Altai", "Chagyr", "Vindija"), arch_names)
            } else keep_arch
            diag_arch <- intersect(diag_arch, names(xcore))
            if (length(diag_arch)) {
                arch_mat <- as.matrix(xcore[, ..diag_arch])
                diag_allele <- vapply(seq_len(nrow(arch_mat)), function(i) {
                  z <- as.character(arch_mat[i, ])
                  if (all(z %chin% base_set) && length(unique(z)) == 1L) z[1] else NA_character_
                }, character(1))
                diag_ok <- !is.na(diag_allele) & !is.na(core$risk_core_allele) &
                  diag_allele == core$risk_core_allele & is.finite(core$risk_freq_delta) &
                  core$risk_freq_delta > diagnostic_delta_th
                diag_count <- rep(NA_integer_, nrow(core))
                if (length(yri_haps)) {
                  diag_count <- vapply(seq_len(nrow(Hcore)), function(i) {
                    if (!diag_ok[i]) return(NA_integer_)
                    sum(Hcore[i, yri_haps] == diag_allele[i], na.rm = TRUE)
                  }, integer(1))
                }
                core[, `:=`(diagnostic_archaic_allele = diag_allele,
                  n_diagnostic_archaics = fifelse(!is.na(diag_allele), length(diag_arch), 0L),
                  is_diagnostic_archaic = diag_ok,
                  yri_diagnostic_count = diag_count,
                  yri_diagnostic_freq = fifelse(diag_ok & yri_called_alleles > 0,
                    diag_count / yri_called_alleles, NA_real_))]
                for (pp in afr_pops) {
                    pop_haps <- afr_haps[[pp]]
                    if (length(pop_haps)) {
                        pop_diag_count <- vapply(seq_len(nrow(Hcore)), function(i) {
                            if (!diag_ok[i]) return(NA_integer_)
                            sum(Hcore[i, pop_haps] == diag_allele[i], na.rm = TRUE)
                        }, integer(1))
                    } else {
                        pop_diag_count <- rep(NA_integer_, nrow(Hcore))
                    }
                    pop_called_col <- paste0(pp, "_called_alleles")
                    core[, (paste0(pp, "_diagnostic_count")) := pop_diag_count]
                    core[, (paste0(pp, "_diagnostic_freq")) := fifelse(diag_ok & get(pop_called_col) > 0,
                        pop_diag_count / get(pop_called_col), NA_real_)]
                }
            }
        }
        diag_i <- which(core$is_diagnostic_archaic == TRUE & !is.na(core$diagnostic_archaic_allele))
        for (pp in afr_pops) {
            pop_haps <- afr_haps[[pp]]
            hf <- NA_real_
            if (length(diag_i) && length(pop_haps)) {
                subH <- Hcore[diag_i, pop_haps, drop = FALSE]
                called_haps <- apply(subH, 2, function(z) all(z %chin% base_set))
                if (any(called_haps)) {
                    diag_vec <- core$diagnostic_archaic_allele[diag_i]
                    hf <- mean(apply(subH[, called_haps, drop = FALSE], 2, function(z) all(z == diag_vec)))
                }
            }
            region[, (paste0(pp, "_haplo_freq")) := hf]
        }
        region[, `:=`(best_lineage = best_lineage_val, matched_archaics = matched_archaics_val)]
        site[, `:=`(best_lineage = best_lineage_val, matched_archaics = matched_archaics_val)]
        write_one(core, file.path(od, paste0(id0, ".core.tsv")))
        write_one(arch_stat, file.path(od, paste0(id0, ".arch.tsv")))
        if (!length(keep_arch)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        x <- x[nchar(ref) == 1L & nchar(alt) == 1L & ref %chin% base_set & alt %chin% base_set]
        samp <- setdiff(names(x), c("chr", "pos", "ref", "alt", "aa", arch_names))
        H <- do.call(cbind, c(lapply(samp, function(s) to_base(gt_hap(x[[s]])$a, x$ref, x$alt)), 
            lapply(samp, function(s) to_base(gt_hap(x[[s]])$b, x$ref, x$alt))))
        colnames(H) <- c(paste0(samp, "_1"), paste0(samp, "_2"))
        keep_modern <- apply(H, 1, function(z) {
            z <- z[!is.na(z) & z %chin% base_set]
            if (!length(z)) 
                return(FALSE)
            tab <- table(z)
            length(tab) >= 2L && min(tab) >= 2L
        })
        keep_arch_called <- Reduce(`&`, lapply(keep_arch, function(an) x[[an]] %chin% base_set))
        keep <- keep_modern & keep_arch_called
        x2 <- x[keep]
        H2 <- H[keep, , drop = FALSE]
        site[, `:=`(n_site_keep, nrow(x2))]
        if (!nrow(x2)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        hap <- data.table(copy = colnames(H2), seq = apply(H2, 2, paste0, collapse = ""))[, .(n = .N, 
            copies = paste(copy, collapse = ";")), by = seq][order(-n, seq)]
        site[, `:=`(n_hap_raw, nrow(hap))]
        hap <- hap[n >= 1L][order(-n, seq)]
        if (!nrow(hap)) {
            write_one(site, file.path(od, paste0(id0, ".site.tsv")))
            write_one(region, file.path(od, paste0(id0, ".region.tsv")))
            return(invisible(NULL))
        }
        hap[, `:=`(hap_id, roman_n(.I))]
        site[, `:=`(n_hap_keep, nrow(hap))]
        for (an in keep_arch) {
            aseq <- paste0(x2[[an]], collapse = "")
            hap[, `:=`((paste0(an, "_match")), vapply(seq, function(sq) score_seq(strsplit(sq, 
                "", fixed = TRUE)[[1]], strsplit(aseq, "", fixed = TRUE)[[1]])["match"], numeric(1)))]
        }
        cols <- paste0(keep_arch, "_match")
        hap[, `:=`(best_arch, keep_arch[max.col(.SD, ties.method = "first")]), .SDcols = cols]
        hap[, `:=`(best_match, apply(.SD, 1, max)), .SDcols = cols]
        risk_i <- which(x2$pos == bp0)
        if (length(risk_i) == 1L) {
            risk_a <- lead_risk_base(rk_raw, x2$ref[risk_i], x2$alt[risk_i])
            hap[, `:=`(carry_risk, if (!is.na(risk_a)) 
                substring(seq, risk_i, risk_i) == risk_a
            else NA)]
        }
        else hap[, `:=`(carry_risk, NA)]
        hap[, `:=`(trait = tr0, id = id0, best_lineage = best_lineage_val, matched_archaics = matched_archaics_val)]
        setcolorder(hap, c("trait", "id", "hap_id", "n", "best_lineage", "matched_archaics", 
            "best_arch", "best_match", "carry_risk", "seq", "copies", cols))
        write_one(hap, file.path(od, paste0(id0, ".hap.tsv")))
        write_one(site, file.path(od, paste0(id0, ".site.tsv")))
        write_one(region, file.path(od, paste0(id0, ".region.tsv")))
    }
    merge_all <- function() {
        gather <- function(pat) {
            fs <- list.files(hap0, pattern = pat, recursive = TRUE, full.names = TRUE)
            fs <- fs[file.info(fs)$size > 0]
            if (!length(fs)) 
                return(data.table())
            rbindlist(lapply(fs, fread), fill = TRUE)
        }
        hap_dt <- gather("\\.hap\\.tsv$")
        site_dt <- gather("\\.site\\.tsv$")
        arch_dt <- gather("\\.arch\\.tsv$")
        region_dt <- gather("\\.region\\.tsv$")
        core_dt <- gather("\\.core\\.tsv$")
        out_files <- c("hap_match.csv", "site_count.csv", "archaic_match.csv", "region_summary.csv", 
            "risk_core.csv", "inherited_segments.csv", "hap_match.tsv", "site_count.tsv", "archaic_match.tsv", 
            "region_summary.tsv", "risk_core.tsv", "inherited_segments.tsv", "hap_site_count.tsv", 
            "core_archaic_match.tsv", "core_risk.tsv", "selected_region.tsv", "yri_filter.tsv")
        unlink(c(file.path(root, out_files), file.path(report0, out_files)), force = TRUE)
        write_report <- function(dt, stem) {
            if (!nrow(dt)) 
                return(invisible(NULL))
            fwrite(dt, file.path(report0, paste0(stem, ".tsv")), sep = "\t")
        }
        if (hap_filter_method == "simple") {
            write_report(hap_dt, "hap_match")
            write_report(site_dt, "hap_site_count")
            write_report(arch_dt, "core_archaic_match")
            write_report(core_dt, "core_risk")
            afr_pops_out <- c("YRI", "LWK", "GWD", "MSL", "ESN")
            afr_summary_cols <- as.vector(rbind(paste0(afr_pops_out, "_max_freq"),
                paste0(afr_pops_out, "_haplo_freq")))
            summary_cols <- c("trait", "id", "lead_chr", "lead_snp", "lead_bp", "core_start",
                "core_end", "n_ld_snp", "core_size_bp", "p_ils", "best_lineage",
                "matched_archaics", "yri_freq", "yri_freq_all", "n_yri_filter_sites",
                afr_summary_cols)
            empty_summary <- data.table(trait = character(), id = character(),
                lead_chr = integer(), lead_snp = character(), lead_bp = integer(),
                core_start = integer(), core_end = integer(), n_ld_snp = integer(),
                core_size_bp = integer(), p_ils = numeric(), best_lineage = character(),
                matched_archaics = character(), yri_freq = numeric(), yri_freq_all = numeric(),
                n_yri_filter_sites = integer())
            for (nm in afr_summary_cols) empty_summary[, (nm) := numeric()]
            if (nrow(region_dt)) {
                reg_cols <- c("trait", "id", "lead_chr", "lead_snp", "lead_bp", "core_start",
                  "core_end", "n_ld_snp", "core_size_bp", "p_ils", "best_lineage", "matched_archaics",
                  paste0(afr_pops_out, "_haplo_freq"))
                reg <- unique(region_dt[, intersect(reg_cols, names(region_dt)), with = FALSE])
                have_yri <- nrow(core_dt) && all(c("risk_core_allele", "yri_risk_freq", 
                  "yri_risk_count", "yri_called_alleles") %in% names(core_dt))
                have_diag_yri <- have_yri && all(c("is_diagnostic_archaic", "yri_diagnostic_freq",
                  "yri_diagnostic_count") %in% names(core_dt))
                if (have_yri) {
                  if (have_diag_yri) {
                    yri <- core_dt[!is.na(risk_core_allele), {
                      yy_all <- if (any(is.finite(yri_risk_freq))) max(yri_risk_freq, na.rm = TRUE) else NA_real_
                      yy <- if (any(is_diagnostic_archaic == TRUE & is.finite(yri_diagnostic_freq))) 
                        max(yri_diagnostic_freq[is_diagnostic_archaic == TRUE], na.rm = TRUE) else NA_real_
                      list(yri_freq = yy, yri_freq_all = yy_all,
                        n_yri_filter_sites = sum(is_diagnostic_archaic == TRUE, na.rm = TRUE),
                        yri_missing = !any(is_diagnostic_archaic == TRUE) ||
                        any(is.na(yri_diagnostic_count[is_diagnostic_archaic == TRUE]) |
                        yri_called_alleles[is_diagnostic_archaic == TRUE] <= 0))
                    }, by = .(trait, id)]
                  } else {
                    yri <- core_dt[!is.na(risk_core_allele), {
                    yy_all <- if (any(is.finite(yri_risk_freq))) max(yri_risk_freq, na.rm = TRUE) else NA_real_
                    informative <- is.finite(risk_freq_delta) & risk_freq_delta > diagnostic_delta_th
                    yy <- if (any(informative & is.finite(yri_risk_freq))) max(yri_risk_freq[informative], 
                      na.rm = TRUE) else NA_real_
                    list(yri_freq = yy, yri_freq_all = yy_all, n_yri_filter_sites = sum(informative, 
                      na.rm = TRUE), yri_missing = any(is.na(yri_risk_count[informative]) | 
                      yri_called_alleles[informative] <= 0) || !any(informative))
                    }, by = .(trait, id)]
                  }
                  reg <- merge(reg, yri, by = c("trait", "id"), all.x = TRUE)
                }
                else reg[, `:=`(yri_freq = NA_real_, yri_freq_all = NA_real_, 
                  n_yri_filter_sites = NA_integer_, yri_missing = NA)]
                for (pp in afr_pops_out) {
                  diag_col <- paste0(pp, "_diagnostic_freq")
                  risk_col <- paste0(pp, "_risk_freq")
                  max_col <- paste0(pp, "_max_freq")
                  if (nrow(core_dt) && all(c("risk_core_allele", "is_diagnostic_archaic", diag_col) %in% names(core_dt))) {
                    pop_max <- core_dt[!is.na(risk_core_allele), {
                      z <- get(diag_col)[is_diagnostic_archaic == TRUE]
                      list(pop_max_freq = if (any(is.finite(z))) max(z, na.rm = TRUE) else NA_real_)
                    }, by = .(trait, id)]
                  } else if (pp == "YRI" && "yri_freq" %in% names(reg)) {
                    pop_max <- reg[, .(trait, id, pop_max_freq = yri_freq)]
                  } else if (nrow(core_dt) && all(c("risk_core_allele", "risk_freq_delta", risk_col) %in% names(core_dt))) {
                    pop_max <- core_dt[!is.na(risk_core_allele), {
                      informative <- is.finite(risk_freq_delta) & risk_freq_delta > diagnostic_delta_th
                      z <- get(risk_col)[informative]
                      list(pop_max_freq = if (any(is.finite(z))) max(z, na.rm = TRUE) else NA_real_)
                    }, by = .(trait, id)]
                  } else {
                    pop_max <- data.table(trait = character(), id = character(), pop_max_freq = numeric())
                  }
                  if (nrow(pop_max)) {
                    setnames(pop_max, "pop_max_freq", max_col)
                    reg <- merge(reg, unique(pop_max, by = c("trait", "id")), by = c("trait", "id"), all.x = TRUE)
                  } else reg[, (max_col) := NA_real_]
                  hap_col <- paste0(pp, "_haplo_freq")
                  if (!hap_col %in% names(reg)) reg[, (hap_col) := NA_real_]
                }
                cand <- reg[!is.na(matched_archaics) & nzchar(matched_archaics)]
                selected <- if (!have_yri || yri_freq_th >= 1) copy(cand) else cand[yri_missing == FALSE & 
                  !is.na(YRI_max_freq) & YRI_max_freq <= yri_freq_th]
                yri_filter <- data.table(
                  yri_freq_th = yri_freq_th,
                  diagnostic_delta_th = diagnostic_delta_th,
                  yri_filter_enabled = have_yri && yri_freq_th < 1,
                  yri_filter_method = if (have_diag_yri) "diagnostic_archaic_allele" else "risk_core_delta_fallback",
                  filter_variable = "YRI_max_freq",
                  matched_loci = nrow(cand),
                  selected_loci = nrow(selected)
                )
                cand[, `:=`(yri_missing = NULL)]
                selected[, `:=`(yri_missing = NULL)]
                for (nm in setdiff(summary_cols, names(cand))) cand[, (nm) := NA]
                for (nm in setdiff(summary_cols, names(selected))) selected[, (nm) := NA]
                fwrite(cand[, ..summary_cols], file.path(report0, "region_summary.tsv"), sep = "\t")
                fwrite(selected[, ..summary_cols], file.path(report0, "selected_region.tsv"), sep = "\t")
                fwrite(yri_filter, file.path(report0, "yri_filter.tsv"), sep = "\t")
            } else {
                fwrite(empty_summary, file.path(report0, "region_summary.tsv"), sep = "\t")
                fwrite(empty_summary, file.path(report0, "selected_region.tsv"), sep = "\t")
                fwrite(data.table(yri_freq_th = yri_freq_th, yri_filter_enabled = FALSE,
                  diagnostic_delta_th = diagnostic_delta_th,
                  yri_filter_method = "none",
                  matched_loci = 0L, selected_loci = 0L), file.path(report0, "yri_filter.tsv"), sep = "\t")
            }
        }
        else {
            write_report(hap_dt, "hap_match")
            write_report(site_dt, "site_count")
            write_report(arch_dt, "archaic_match")
            write_report(region_dt, "region_summary")
            write_report(core_dt, "risk_core")
            if (nrow(region_dt)) {
                inherited <- region_dt[!is.na(best_lineage) & nzchar(best_lineage) & nzchar(matched_archaics)]
                write_report(inherited, "inherited_segments")
            }
        }
    }
    run_trait <- function(tr0) {
        td <- file.path(mat0, tr0)
        if (!dir.exists(td)) 
            stop("trait not found: ", tr0)
        for (id0 in list.dirs(td, recursive = FALSE, full.names = FALSE)) {
            clean_one(tr0, id0)
            process_one(tr0, id0)
        }
    }
    run_all <- function() {
        for (tr0 in list.dirs(mat0, recursive = FALSE, full.names = FALSE)) run_trait(tr0)
        merge_all()
    }
    if (!length(args)) 
        run_all()
    else if (length(args) == 1L && args[1] == "merge") 
        merge_all()
    else if (length(args) == 1L) 
        run_trait(args[1])
    else if (length(args) == 2L) {
        clean_one(args[1], args[2])
        process_one(args[1], args[2])
    }
    else stop("usage: Rscript locus.R make_hap <res_dir> [trait | trait id | merge]")
}

run_filter_hap <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table, writexl)
    args <- commandArgs(TRUE)
    root <- normalizePath(if (length(args) >= 1L) 
        args[1]
    else "/mnt/d/analysis/gu/locus", winslash = "/", mustWork = FALSE)
    sample_file <- normalizePath(if (length(args) >= 2L) 
        args[2]
    else "/mnt/d/files/1kg.v3.sample.txt", winslash = "/", mustWork = FALSE)
    target_pop <- if (length(args) >= 3L) 
        args[3]
    else Sys.getenv("FILTER_POP", "YRI")
    max_count <- as.integer(if (length(args) >= 4L) args[4] else Sys.getenv("FILTER_MAX_COUNT", 
        "1"))
    if (is.na(max_count)) 
        stop("FILTER_MAX_COUNT must be an integer")
    report0 <- file.path(root, "report")
    mat0 <- file.path(root, "mat")
    lead0 <- file.path(root, "lead")
    dir.create(report0, recursive = TRUE, showWarnings = FALSE)
    hapf <- file.path(report0, "hap_match.tsv")
    segf <- file.path(report0, "inherited_segments.tsv")
    if (!file.exists(hapf)) 
        stop("missing hap_match.tsv: ", hapf)
    if (!file.exists(segf)) 
        stop("missing inherited_segments.tsv: ", segf)
    if (!file.exists(sample_file)) 
        stop("missing sample file: ", sample_file)
    hap <- fread(hapf)
    seg <- fread(segf)
    sample_meta <- fread(sample_file)
    if (!all(c("sample", "pop") %in% names(sample_meta))) 
        stop("sample file must contain columns: sample,pop")
    if (!"super_pop" %in% names(sample_meta)) 
        sample_meta[, `:=`(super_pop, NA_character_)]
    sample_meta[, `:=`(pop = as.character(pop), super_pop = as.character(super_pop))]
    pop_levels <- sort(unique(na.omit(sample_meta$pop)))
    super_levels <- sort(unique(na.omit(sample_meta$super_pop)))
    pop_cols <- paste0("pop_", pop_levels)
    super_cols <- paste0("super_", super_levels)
    main_hap <- hap[n > 10]
    main_hap <- main_hap[trait %chin% seg$trait & id %chin% seg$id]
    filter_scope <- "n_gt_10"
    if ("carry_risk" %in% names(main_hap)) {
        main_hap[, `:=`(carry_risk, as.logical(carry_risk))]
        filter_hap <- main_hap[carry_risk == TRUE]
        filter_scope <- "carry_risk_TRUE_and_n_gt_10"
    }
    else {
        filter_hap <- main_hap
    }
    sample_cache <- new.env(parent = emptyenv())
    get_sample_map <- function(trait, id) {
        key <- paste(trait, id, sep = "\r")
        if (exists(key, envir = sample_cache, inherits = FALSE)) 
            return(get(key, envir = sample_cache))
        f <- file.path(mat0, trait, id, "kg.samples.tsv")
        if (!file.exists(f)) 
            stop("missing kg sample list: ", f)
        x <- fread(f, header = FALSE, col.names = "sample")
        x[, `:=`(sidx, paste0("s", .I))]
        x <- merge(x, sample_meta[, .(sample, pop, super_pop)], by = "sample", all.x = TRUE, 
            sort = FALSE)
        assign(key, x, envir = sample_cache)
        x
    }
    copy_table <- function(trait, id, hap_id, n, copies) {
        cp <- unlist(strsplit(copies, ";", fixed = TRUE), use.names = FALSE)
        cp <- cp[nzchar(cp)]
        if (!length(cp)) {
            return(data.table(trait = trait, id = id, hap_id = hap_id, hap_n = n, copy = character(), 
                sidx = character(), sample = character(), pop = character(), super_pop = character()))
        }
        x <- data.table(copy = cp, sidx = sub("_[12]$", "", cp))
        x <- merge(x, get_sample_map(trait, id)[, .(sidx, sample, pop, super_pop)], by = "sidx", 
            all.x = TRUE, sort = FALSE)
        x[, `:=`(trait = trait, id = id, hap_id = hap_id, hap_n = n)]
        setcolorder(x, c("trait", "id", "hap_id", "hap_n", "copy", "sidx", "sample", "pop", "super_pop"))
        x
    }
    copies_long <- if (nrow(filter_hap)) {
        rbindlist(lapply(seq_len(nrow(filter_hap)), function(i) {
            copy_table(filter_hap$trait[i], filter_hap$id[i], filter_hap$hap_id[i], filter_hap$n[i], 
                filter_hap$copies[i])
        }), fill = TRUE)
    }
    else data.table(trait = character(), id = character(), hap_id = character(), hap_n = integer(), 
        copy = character(), sidx = character(), sample = character(), pop = character(), super_pop = character())
    wide_counts <- function(x, by_cols) {
        base <- unique(seg[, .(trait, id)])
        if (!nrow(x)) {
            for (cc in c(pop_cols, super_cols)) base[, `:=`((cc), 0L)]
            return(base)
        }
        p <- dcast(x[!is.na(pop), .N, by = c(by_cols, "pop")], paste(paste(by_cols, collapse = " + "), 
            "~ pop"), value.var = "N", fill = 0)
        s <- dcast(x[!is.na(super_pop), .N, by = c(by_cols, "super_pop")], paste(paste(by_cols, 
            collapse = " + "), "~ super_pop"), value.var = "N", fill = 0)
        for (z in list(p, s)) {
            if (!nrow(z)) 
                next
        }
        if (nrow(p)) 
            setnames(p, setdiff(names(p), by_cols), paste0("pop_", setdiff(names(p), by_cols)))
        if (nrow(s)) 
            setnames(s, setdiff(names(s), by_cols), paste0("super_", setdiff(names(s), by_cols)))
        out <- Reduce(function(a, b) merge(a, b, by = by_cols, all = TRUE, sort = FALSE), Filter(nrow, 
            list(p, s)))
        if (!length(out)) 
            out <- unique(x[, ..by_cols])
        for (cc in c(pop_cols, super_cols)) if (!cc %in% names(out)) 
            out[, `:=`((cc), 0L)]
        for (cc in c(pop_cols, super_cols)) set(out, which(is.na(out[[cc]])), cc, 0L)
        out
    }
    locus_counts <- wide_counts(copies_long, c("trait", "id"))
    hap_counts <- if (nrow(copies_long)) 
        wide_counts(copies_long, c("trait", "id", "hap_id"))
    else data.table()
    locus_meta <- filter_hap[, .(n_filter_hap = .N, n_filter_copies = sum(n)), by = .(trait, 
        id)]
    locus_summary <- merge(seg, locus_meta, by = c("trait", "id"), all.x = TRUE, sort = FALSE)
    locus_summary <- merge(locus_summary, locus_counts, by = c("trait", "id"), all.x = TRUE, 
        sort = FALSE)
    locus_summary[is.na(n_filter_hap), `:=`(n_filter_hap, 0L)]
    locus_summary[is.na(n_filter_copies), `:=`(n_filter_copies, 0L)]
    for (cc in c(pop_cols, super_cols)) {
        if (!cc %in% names(locus_summary)) 
            locus_summary[, `:=`((cc), 0L)]
        set(locus_summary, which(is.na(locus_summary[[cc]])), cc, 0L)
    }
    target_col <- paste0("pop_", target_pop)
    if (!target_col %in% names(locus_summary)) 
        locus_summary[, `:=`((target_col), 0L)]
    locus_summary[, `:=`(keep, n_filter_hap > 0L & get(target_col) <= max_count)]
    front <- c("trait", "id", "lead_chr", "lead_snp", "lead_bp", "core_start", "core_end", "n_ld_snp", 
        "core_size_bp", "p_ils", "best_lineage", "matched_archaics", "n_filter_hap", "n_filter_copies", 
        target_col, "keep")
    setcolorder(locus_summary, c(front[front %in% names(locus_summary)], setdiff(names(locus_summary), 
        front)))
    if (nrow(hap_counts)) {
        hap_counts <- merge(filter_hap[, setdiff(names(filter_hap), "copies"), with = FALSE], 
            hap_counts, by = c("trait", "id", "hap_id"), all.x = TRUE, sort = FALSE)
        for (cc in c(pop_cols, super_cols)) {
            if (!cc %in% names(hap_counts)) 
                hap_counts[, `:=`((cc), 0L)]
            set(hap_counts, which(is.na(hap_counts[[cc]])), cc, 0L)
        }
    }
    filtered_loci <- locus_summary[keep == TRUE]
    filtered_keys <- filtered_loci[, .(trait, id)]
    filtered_hap <- merge(filter_hap, filtered_keys, by = c("trait", "id"), all = FALSE, sort = FALSE)
    filtered_seg <- merge(seg, filtered_keys, by = c("trait", "id"), all = FALSE, sort = FALSE)
    criteria <- data.table(filter_pop = target_pop, max_count = max_count, filter_scope = filter_scope, 
        input_inherited_loci = nrow(seg), input_filter_hap = nrow(filter_hap), kept_hap = nrow(filtered_hap), 
        kept_loci = nrow(filtered_seg))
    read_optional <- function(stem) {
        f <- file.path(report0, paste0(stem, ".tsv"))
        if (file.exists(f) && file.info(f)$size > 0) 
            fread(f)
        else data.table()
    }
    filter_dt <- function(x) {
        if (!nrow(x) || !all(c("trait", "id") %in% names(x)) || !nrow(filtered_keys)) 
            return(x[0])
        merge(x, filtered_keys, by = c("trait", "id"), all = FALSE, sort = FALSE)
    }
    write_book <- function(sheets, path) {
        sheets <- sheets[vapply(sheets, function(x) is.data.frame(x) && nrow(x) >= 0, logical(1))]
        write_xlsx(sheets, path)
    }
    all_sheets <- list(inherited_loci_counts = locus_summary, inherited_haplotype_counts = hap_counts, 
        inherited_segments = seg, archaic_match = read_optional("archaic_match"), risk_core = read_optional("risk_core"), 
        site_count = read_optional("site_count"), region_summary = read_optional("region_summary"), 
        filter_criteria = criteria)
    filtered_sheets <- list(filtered_loci = filtered_loci, filtered_haplotypes = hap_counts[filtered_keys, 
        on = .(trait, id), nomatch = 0], filtered_inherited_segments = filtered_seg, filtered_archaic_match = filter_dt(read_optional("archaic_match")), 
        filtered_risk_core = filter_dt(read_optional("risk_core")), filtered_site_count = filter_dt(read_optional("site_count")), 
        filter_criteria = criteria)
    selected_sheets <- list(selected_loci = if (file.exists(file.path(lead0, "pick.tsv"))) fread(file.path(lead0, 
        "pick.tsv")) else data.table(), single_snp_loci = if (file.exists(file.path(lead0, "pick.single.tsv"))) fread(file.path(lead0, 
        "pick.single.tsv")) else data.table(), lead_assoc = if (file.exists(file.path(lead0, 
        "lead.assoc"))) fread(file.path(lead0, "lead.assoc")) else data.table())
    fwrite(locus_summary, file.path(report0, "hap_pop_counts.tsv"), sep = "\t")
    fwrite(filtered_hap, file.path(report0, "filtered_hap_match.tsv"), sep = "\t")
    fwrite(filtered_seg, file.path(report0, "filtered_inherited_segments.tsv"), sep = "\t")
    write_book(all_sheets, file.path(report0, "all.xlsx"))
    write_book(filtered_sheets, file.path(report0, "filtered.xlsx"))
    write_book(selected_sheets, file.path(report0, "selected.xlsx"))
    cat(sprintf("filter_pop=%s max_count=%d input_inherited_loci=%d input_filter_hap=%d kept_hap=%d kept_loci=%d\n", 
        target_pop, max_count, nrow(seg), nrow(filter_hap), nrow(filtered_hap), nrow(filtered_seg)))
}

run_make_phy <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    library(data.table)
    args0 <- commandArgs(TRUE)
    has_root <- length(args0) && grepl("[/\\\\]", args0[1])
    root <- normalizePath(if (has_root) 
        args0[1]
    else Sys.getenv("BALD_RES", "/data/sph-zhaor/analysis/bald/res"), winslash = "/", mustWork = FALSE)
    args <- if (has_root) 
        args0[-1]
    else args0
    target_trait <- if (length(args) >= 1L) 
        args[1]
    else NA_character_
    target_id <- if (length(args) >= 2L) 
        args[2]
    else NA_character_
    make_phy_method <- Sys.getenv("MAKE_PHY_METHOD", "strict")
    report0 <- file.path(root, "report")
    if (make_phy_method == "hpc") make_phy_method <- "simple"
    hapf <- if (make_phy_method == "simple") file.path(report0, "hap_match.tsv") else file.path(report0, 
        "filtered_hap_match.tsv")
    if (!file.exists(hapf)) hapf <- file.path(report0, "hap_match.tsv")
    if (!file.exists(hapf)) 
        hapf <- file.path(root, "hap_match.tsv")
    selectedf <- file.path(report0, "selected_region.tsv")
    positivef <- file.path(root, "lead", "positive_pick.tsv")
    mat0 <- file.path(root, "mat")
    phy0 <- file.path(root, "phy")
    dir.create(phy0, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(hapf) || file.info(hapf)$size == 0) 
        quit(save = "no", status = 1)
    base_set <- c("A", "C", "G", "T")
    hap <- fread(hapf)
    if (make_phy_method == "simple") {
        if (!file.exists(selectedf) || file.info(selectedf)$size == 0) 
            quit(save = "no", status = 1)
        selected <- fread(selectedf)
        chr_cols <- c("trait", "id", "lead_snp", "best_lineage", "matched_archaics")
        int_cols <- c("lead_chr", "lead_bp", "core_start", "core_end", "n_ld_snp", "core_size_bp")
        num_cols <- c("p_ils", "yri_freq")
        for (nm in intersect(chr_cols, names(selected))) selected[, (nm) := as.character(get(nm))]
        for (nm in intersect(int_cols, names(selected))) selected[, (nm) := as.integer(get(nm))]
        for (nm in intersect(num_cols, names(selected))) selected[, (nm) := as.numeric(get(nm))]
        phy_regions <- copy(selected)
        if (file.exists(positivef) && file.info(positivef)$size > 0) {
            positive <- fread(positivef)
            if (nrow(positive)) {
                positive[, `:=`(trait = as.character(trait), lead_snp = as.character(lead_snp))]
                positive[, id := paste(lead_chr, lead_snp, lead_bp, sep = ".")]
                pos_region <- positive[, .(trait, id, lead_chr, lead_snp, lead_bp, 
                  core_start = start, core_end = end, n_ld_snp = n, 
                  core_size_bp = size_bp, p_ils = NA_real_, 
                  best_lineage = "forced_control", matched_archaics = "", 
                  yri_freq = NA_real_)]
                phy_regions <- unique(rbindlist(list(phy_regions, pos_region), 
                  fill = TRUE), by = c("trait", "id"))
            }
        }
        fwrite(phy_regions, file.path(report0, "phy_region.tsv"), sep = "\t")
        if (!nrow(phy_regions)) 
            quit(save = "no", status = 0)
        hap <- merge(hap, unique(phy_regions[, .(trait, id)]), by = c("trait", "id"))
    }
    if (!is.na(target_trait)) 
        hap <- hap[trait == target_trait]
    if (!is.na(target_id)) 
        hap <- hap[id == target_id]
    if (!nrow(hap)) 
        quit(save = "no", status = 0)
    old_root <- if (is.na(target_trait)) 
        phy0
    else file.path(phy0, target_trait)
    oldf <- list.files(old_root, pattern = "\\.(phy|meta\\.tsv|phy_phyml_tree\\.txt|phy_phyml_stats\\.txt)$", 
        recursive = TRUE, full.names = TRUE)
    if (is.na(target_id) && length(oldf)) 
        unlink(oldf, force = TRUE)
    arch_label <- function(f) {
        x <- tools::file_path_sans_ext(basename(f))
        paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    }
    arch_files <- function(d0) {
        fs <- list.files(d0, pattern = "\\.tsv$", full.names = TRUE)
        fs <- fs[!basename(fs) %in% c("kg.tsv", "kg.samples.tsv")]
        setNames(fs, arch_label(fs))
    }
    gt_hap <- function(x) {
        x <- gsub("\\|", "/", x)
        sp <- tstrsplit(x, "/", fixed = TRUE)
        list(a = sp[[1]], b = sp[[2]])
    }
    to_base <- function(h, ref, alt) {
        out <- rep("N", length(h))
        out[h == "0"] <- ref[h == "0"]
        out[h == "1"] <- alt[h == "1"]
        out
    }
    arch_base <- function(gt, ref, alt) {
        gt <- gsub("\\|", "/", gt)
        sp <- strsplit(gt, "/", fixed = TRUE)
        vapply(seq_along(sp), function(i) {
            z <- sp[[i]]
            z <- z[z != "."]
            if (length(z) != 2L || z[1] != z[2]) 
                return(NA_character_)
            a <- z[1]
            alts <- if (alt[i] %in% c(".", "")) 
                character(0)
            else strsplit(alt[i], ",", fixed = TRUE)[[1]]
            if (a == "0") 
                return(ref[i])
            ai <- suppressWarnings(as.integer(a))
            if (is.na(ai) || ai < 1L || ai > length(alts)) 
                return(NA_character_)
            alts[ai]
        }, character(1))
    }
    read_mat <- function(f, type) {
        x <- fread(f, header = FALSE)
        if (!nrow(x)) 
            return(NULL)
        if (type == "kg") 
            setnames(x, c("chr", "pos", "ref", "alt", "aa", paste0("s", seq_len(ncol(x) - 5))))
        else setnames(x, c("chr", "pos", "ref", "alt", "gt"))
        x[, `:=`(chr, sub("^chr", "", as.character(chr), ignore.case = TRUE))]
        x[, `:=`(ref = toupper(trimws(ref)), alt = toupper(trimws(alt)))]
        if ("aa" %in% names(x)) {
            x[, `:=`(aa, toupper(sub("\\|.*$", "", trimws(aa))))]
        }
        x
    }
    read_arch <- function(f, lab) {
        x <- read_mat(f, "arch")
        if (is.null(x)) 
            return(NULL)
        x[, `:=`(allele, arch_base(gt, ref, alt))]
        x <- x[, .(chr, pos, allele)]
        setnames(x, "allele", lab)
        x
    }
    merge_arch <- function(kg, fs, keep_arch) {
        x <- kg
        keep <- character(0)
        for (lab in keep_arch) {
            if (!lab %in% names(fs)) 
                next
            a <- read_arch(fs[[lab]], lab)
            if (is.null(a)) 
                next
            x <- merge(x, a, by = c("chr", "pos"), all.x = TRUE)
            keep <- c(keep, lab)
        }
        setorder(x, pos)
        list(x = x, arch = keep)
    }
    clean_phy <- function(tr0, id0) {
        od <- file.path(phy0, tr0)
        dir.create(od, recursive = TRUE, showWarnings = FALSE)
        unlink(file.path(od, paste0(id0, c(".full.phy", ".full.meta.tsv", ".full.phy_phyml_tree.txt", 
            ".full.phy_phyml_stats.txt", ".main.phy", ".main.meta.tsv", ".main.phy_phyml_tree.txt", 
            ".main.phy_phyml_stats.txt"))), force = TRUE)
    }
    keys <- unique(hap[, .(trait, id, matched_archaics)])
    for (i in seq_len(nrow(keys))) {
        tr0 <- keys$trait[i]
        id0 <- keys$id[i]
        clean_phy(tr0, id0)
        keep_arch <- unlist(strsplit(keys$matched_archaics[i], ";", fixed = TRUE))
        keep_arch <- keep_arch[nzchar(keep_arch)]
        if (!length(keep_arch)) 
            next
        h <- hap[trait == tr0 & id == id0]
        if (!nrow(h)) 
            next
        d0 <- file.path(mat0, tr0, id0)
        kgf <- file.path(d0, "kg.tsv")
        fs <- arch_files(d0)
        if (!file.exists(kgf) || !length(fs)) 
            next
        kg <- read_mat(kgf, "kg")
        if (is.null(kg)) 
            next
        ma <- merge_arch(kg, fs, keep_arch)
        x <- ma$x
        keep_arch <- ma$arch
        if (!length(keep_arch)) 
            next
        samp <- setdiff(names(x), c("chr", "pos", "ref", "alt", "aa", names(fs)))
        H <- do.call(cbind, c(lapply(samp, function(s) to_base(gt_hap(x[[s]])$a, x$ref, x$alt)), 
            lapply(samp, function(s) to_base(gt_hap(x[[s]])$b, x$ref, x$alt))))
        keep_modern <- apply(H, 1, function(z) {
            z <- z[!is.na(z) & z %chin% base_set]
            if (!length(z)) 
                return(FALSE)
            tab <- table(z)
            length(tab) >= 2L && min(tab) >= 2L
        })
        keep_arch_called <- Reduce(`&`, lapply(keep_arch, function(an) x[[an]] %chin% base_set))
        x2 <- x[keep_modern & keep_arch_called]
        if (!nrow(x2)) 
            next
        anc <- x2$aa
        anc[!anc %in% base_set] <- "N"
        od <- file.path(phy0, tr0)
        dir.create(od, recursive = TRUE, showWarnings = FALSE)
        make_one <- function(sub, tag) {
            if (!nrow(sub)) 
                return(invisible(NULL))
            seqs <- setNames(sub$seq, sub$hap_id)
            for (an in keep_arch) {
                z <- x2[[an]]
                z[!z %in% base_set] <- "N"
                seqs <- c(seqs, setNames(paste0(z, collapse = ""), an))
            }
            seqs <- c(seqs, Ancestral = paste0(anc, collapse = ""))
            len <- unique(nchar(seqs))
            if (length(len) != 1L) 
                return(invisible(NULL))
            lab0 <- names(seqs)
            lab <- if (make_phy_method == "simple") substr(lab0, 1, 12) else sprintf("S%09d", 
                seq_along(seqs))
            phyf <- file.path(od, paste0(id0, ".", tag, ".phy"))
            con <- file(phyf, "w")
            writeLines(sprintf("%d %d", length(seqs), len), con)
            if (make_phy_method == "simple") {
                for (j in seq_along(seqs)) writeLines(sprintf("%-12s%s", lab[j], seqs[j]), con)
            }
            else {
                for (j in seq_along(seqs)) writeLines(sprintf("%-10s %s", lab[j], seqs[j]), con)
            }
            close(con)
            meta <- data.table(label = lab, original_label = lab0, type = fifelse(lab0 %in% 
                keep_arch, "archaic", fifelse(lab0 == "Ancestral", "ancestral", "modern")))
            meta <- merge(meta, sub[, .(original_label = hap_id, n, best_lineage, best_arch, 
                best_match)], by = "original_label", all.x = TRUE)
            fwrite(meta, file.path(od, paste0(id0, ".", tag, ".meta.tsv")), sep = "\t")
        }
        if (make_phy_method == "simple") make_one(h, "full")
        make_one(if (make_phy_method == "simple") h[n > 10] else h, "main")
    }
}

run_make_tree <- function(script_args) {
    commandArgs <- function(trailingOnly = FALSE, ...) {
        if (isTRUE(trailingOnly)) 
            script_args
        else base::commandArgs(FALSE)
    }
    pacman::p_load(data.table, ape, eoffice)
    args <- commandArgs(TRUE)
    root <- normalizePath(if (length(args)) 
        args[1]
    else Sys.getenv("BALD_RES", "D:/bald/res"), winslash = "/", mustWork = FALSE)
    phy0 <- file.path(root, "phy")
    out0 <- file.path(root, "plot")
    dir.create(out0, recursive = TRUE, showWarnings = FALSE)
    keep_existing <- Sys.getenv("MAKE_TREE_KEEP_EXISTING", "0") == "1"
    write_pptx <- Sys.getenv("MAKE_TREE_PPTX", "1") != "0"
    if (!keep_existing) 
        unlink(list.files(out0, pattern = "^s8_tree_.*\\.png$", full.names = TRUE), force = TRUE)
    ppt_full <- file.path(out0, "s8_tree_full.pptx")
    ppt_main <- file.path(out0, "s8_tree_main.pptx")
    if (write_pptx && !keep_existing && file.exists(ppt_full)) 
        file.remove(ppt_full)
    if (write_pptx && !keep_existing && file.exists(ppt_main)) 
        file.remove(ppt_main)
    min_boot <- 70
    scale_bar_len <- 0.005
    arch_col <- "#1f78b4"
    anc_col <- "#33a02c"
    hi_col <- "#d73027"
    style <- function(tag) {
        if (tag == "main") {
            list(dot_r = 0.0035, cex_tip = 0.36, cex_boot = 0.36, lwd_tree = 0.45, lwd_guide = 0.22, 
                r_lab = 1.1, r_text = 1.18, r_lim = 1.3, ppt_w = 3.6, ppt_h = 3.6)
        }
        else {
            list(dot_r = 0.0028, cex_tip = 0.27, cex_boot = 0.28, lwd_tree = 0.34, lwd_guide = 0.16, 
                r_lab = 1.08, r_text = 1.15, r_lim = 1.32, ppt_w = 4, ppt_h = 3.8)
        }
    }
    topptx_safe <- function(p, filename, width, height, append = FALSE) {
        if ("append" %in% names(formals(topptx))) {
            topptx(p, filename = filename, width = width, height = height, append = append)
        }
        else {
            if (append) 
                warning("topptx() has no append argument; output may be overwritten")
            topptx(p, filename = filename, width = width, height = height)
        }
    }
    read_phyml_tree <- function(f) {
        x <- readLines(f, warn = FALSE)
        x <- x[nzchar(trimws(x))]
        i <- grep("^\\(", x)
        if (!length(i)) 
            stop("no tree found: ", f)
        read.tree(text = x[i[1]])
    }
    desc_tips <- function(tr, node) {
        nt <- Ntip(tr)
        kid <- tr$edge[tr$edge[, 1] == node, 2]
        unlist(lapply(kid, function(k) if (k <= nt) 
            k
        else desc_tips(tr, k)), use.names = FALSE)
    }
    parent_node <- function(tr, node) {
        x <- tr$edge[tr$edge[, 2] == node, 1]
        if (length(x)) 
            x[1]
        else NA_integer_
    }
    sector_range <- function(a) {
        a <- sort((a + 2 * pi)%%(2 * pi))
        if (length(a) <= 1) 
            return(c(a, a))
        g <- c(diff(a), a[1] + 2 * pi - a[length(a)])
        j <- which.max(g)
        s <- a[(j%%length(a)) + 1]
        e <- a[j]
        if (e < s) 
            e <- e + 2 * pi
        c(s, e)
    }
    tip_type <- function(meta, labs) {
        out <- rep(NA_character_, length(labs))
        names(out) <- labs
        if (!is.null(meta) && all(c("label", "type") %in% names(meta))) {
            idx <- match(labs, meta$label)
            out[!is.na(idx)] <- as.character(meta$type[idx[!is.na(idx)]])
        }
        out[grepl("vindija|altai|chagyr|denisova|neand|archaic", labs, ignore.case = TRUE)] <- "archaic"
        out[labs == "Ancestral"] <- "ancestral"
        out[is.na(out)] <- "modern"
        out
    }
    highlight_node <- function(tr, meta) {
        nt <- Ntip(tr)
        labs <- tr$tip.label
        tp <- tip_type(meta, labs)
        arch <- labs[tp == "archaic"]
        modern <- labs[tp == "modern"]
        arch <- arch[!is.na(arch) & nzchar(arch)]
        modern <- modern[!is.na(modern) & nzchar(modern)]
        if (!length(arch) || !length(modern)) 
            return(NA_integer_)
        candidates <- function(all_arch = TRUE) {
            rbindlist(lapply((nt + 1):(nt + tr$Nnode), function(node) {
                tips <- labs[desc_tips(tr, node)]
                tips <- tips[!is.na(tips) & nzchar(tips)]
                ok_arch <- if (all_arch) 
                  isTRUE(all(arch %in% tips))
                else isTRUE(any(tips %in% arch))
                has_modern <- isTRUE(any(tips %in% modern))
                has_ancestral <- isTRUE(any(tips == "Ancestral"))
                support <- suppressWarnings(as.numeric(tr$node.label[node - nt]))
                if (length(support) != 1L) 
                  support <- NA_real_
                if (!ok_arch || !has_modern || has_ancestral || is.na(support) || support < min_boot) 
                  return(NULL)
                data.table(node = node, support = support, n_tips = length(tips), n_mod = sum(tips %in% 
                  modern))
            }), fill = TRUE)
        }
        z <- candidates(TRUE)
        if (!nrow(z)) 
            z <- candidates(FALSE)
        if (!nrow(z)) 
            return(NA_integer_)
        setorder(z, n_tips, -support, -n_mod)
        z$node[1]
    }
    draw_tree <- function(tr, meta = NULL, tag = c("full", "main")) {
        tag <- match.arg(tag)
        st <- style(tag)
        tr$tip.label <- trimws(tr$tip.label)
        if ("Ancestral" %in% tr$tip.label) {
            tr <- tryCatch(root(tr, outgroup = "Ancestral", resolve.root = TRUE), error = function(e) tr)
        }
        tr <- ladderize(tr)
        nt <- Ntip(tr)
        ia <- which(tr$tip.label == "Ancestral")
        dep <- node.depth.edgelength(tr)[1:nt]
        r_tree <- max(dep[setdiff(seq_len(nt), ia)], na.rm = TRUE)
        if (!is.finite(r_tree)) 
            r_tree <- max(dep, na.rm = TRUE)
        r_lab <- r_tree * st$r_lab
        r_text <- r_tree * st$r_text
        r_lim <- r_tree * st$r_lim
        par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA, pty = "s")
        plot.phylo(tr, type = "fan", use.edge.length = TRUE, show.tip.label = FALSE, no.margin = TRUE, 
            edge.width = st$lwd_tree, x.lim = c(-r_lim, r_lim), y.lim = c(-r_lim, r_lim))
        pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        xx <- pp$xx
        yy <- pp$yy
        xt <- xx[1:nt]
        yt <- yy[1:nt]
        ang <- atan2(yt, xt)
        labs <- tr$tip.label
        tp <- tip_type(meta, labs)
        hi <- highlight_node(tr, meta)
        if (!is.na(hi)) {
            idx <- desc_tips(tr, hi)
            aa <- sector_range(atan2(yy[idx], xx[idx]))
            th <- seq(aa[1], aa[2], length.out = 500)
            r0 <- sqrt(xx[hi]^2 + yy[hi]^2)
            polygon(c(r0 * cos(th), rev(r_lab * cos(th))), c(r0 * sin(th), rev(r_lab * sin(th))), 
                col = adjustcolor(hi_col, alpha.f = 0.22), border = NA)
        }
        segments(xt, yt, r_lab * cos(ang), r_lab * sin(ang), lty = 3, col = "grey45", lwd = st$lwd_guide)
        symbols(xt, yt, circles = rep(r_tree * st$dot_r, nt), inches = FALSE, add = TRUE, bg = "black", 
            fg = "black", lwd = 0.15)
        lab_col <- rep("black", nt)
        lab_col[tp == "archaic"] <- arch_col
        lab_col[tp == "ancestral"] <- anc_col
        deg <- ang * 180/pi
        flip <- deg < -90 | deg > 90
        srt <- ifelse(flip, deg + 180, deg)
        for (i in seq_len(nt)) {
            text(r_text * cos(ang[i]), r_text * sin(ang[i]), labs[i], srt = srt[i], adj = if (flip[i]) 
                c(1, 0.5)
            else c(0, 0.5), cex = st$cex_tip, col = lab_col[i])
        }
        if (!is.na(hi)) {
            boot <- suppressWarnings(as.numeric(tr$node.label))
            nd <- nt + seq_len(tr$Nnode)
            show <- unique(c(parent_node(tr, hi), hi))
            idx <- match(show[show > nt], nd)
            idx <- idx[!is.na(idx) & !is.na(boot[idx])]
            if (length(idx)) 
                nodelabels(boot[idx], node = nd[idx], frame = "n", cex = st$cex_boot, col = "black")
        }
        usr <- par("usr")
        add.scale.bar(x = usr[2] - 0.14 * diff(usr[1:2]), y = usr[3] + 0.06 * diff(usr[3:4]), 
            length = scale_bar_len, lwd = 0.8, cex = 0.65)
    }
    assign("draw_tree", draw_tree, envir = .GlobalEnv)
    make_editable_plot <- function(tr, meta, tag) {
        assign(".locus_tree_tr", tr, envir = .GlobalEnv)
        assign(".locus_tree_meta", meta, envir = .GlobalEnv)
        assign(".locus_tree_tag", tag, envir = .GlobalEnv)
        assign(".locus_tree_draw_once", function() {
            draw_tree(get(".locus_tree_tr", envir = .GlobalEnv), get(".locus_tree_meta", envir = .GlobalEnv), 
                get(".locus_tree_tag", envir = .GlobalEnv))
        }, envir = .GlobalEnv)
        convertplot(.locus_tree_draw_once())
    }
    files <- list.files(phy0, pattern = "\\.main\\.phy_phyml_tree\\.txt$", recursive = TRUE, 
        full.names = TRUE)
    if (!length(files)) 
        quit(save = "no", status = 0)
    n_ppt <- c(full = 0L, main = 0L)
    for (treef in files) {
        base <- sub("\\.phy_phyml_tree\\.txt$", "", basename(treef))
        tag <- if (grepl("\\.main$", base)) 
            "main"
        else "full"
        st <- style(tag)
        metaf <- file.path(dirname(treef), paste0(base, ".meta.tsv"))
        meta <- if (file.exists(metaf)) 
            fread(metaf)
        else NULL
        tr <- tryCatch(read_phyml_tree(treef), error = function(e) NULL)
        if (is.null(tr)) 
            next
        trait <- basename(dirname(treef))
        pngf <- file.path(out0, paste0("s8_tree_", tag, "_", trait, "_", base, ".png"))
        if (!keep_existing || !file.exists(pngf) || file.info(pngf)$size == 0) {
            png(pngf, width = st$ppt_w, height = st$ppt_h, units = "in", res = 300)
            draw_tree(tr, meta, tag)
            dev.off()
        }
        p <- make_editable_plot(tr, meta, tag)
        if (write_pptx) {
            p <- make_editable_plot(tr, meta, tag)
            topptx_safe(p, if (tag == "main") 
                ppt_main
            else ppt_full, width = st$ppt_w, height = st$ppt_h, append = keep_existing || n_ppt[tag] > 
                0L)
        }
        n_ppt[tag] <- n_ppt[tag] + 1L
    }
}


run_add_positive_loci <- function(script_args){
    library(data.table)
    args <- script_args
    arg <- function(k){ i <- which(args %in% k)[1]; if(is.na(i) || i == length(args)) return(NULL); args[i+1] }
    dirgwas <- arg(c('--dirgwas','-dirgwas'))
    dirmod  <- arg(c('--dirmod','-dirmod'))
    dirout  <- arg(c('--dirout','-dirout'))
    bed     <- arg(c('--bed','-bed'))
    trait   <- arg(c('--trait','-trait'))
    ref_pop <- arg(c('--ref_pop','-ref_pop')); if(is.null(ref_pop)) ref_pop <- 'ALL'
    if(is.null(dirgwas) || is.null(dirout) || is.null(bed) || is.null(trait)) stop('usage: Rscript locus.R add_positive_loci --dirgwas DIR --dirmod DIR --dirout DIR --bed BED --trait TRAIT [--ref_pop ALL|EUR]')
    if(!file.exists(bed) || !isTRUE(file.info(bed)$size > 0)) {
        message('Positive-control BED not found or empty; skip adding positive loci: ', bed)
        return(invisible(NULL))
    }
    lead0 <- file.path(dirout, 'lead'); dir.create(lead0, recursive=TRUE, showWarnings=FALSE)
    assocf <- file.path(lead0, paste0(trait, '.lead.assoc'))
    threef <- file.path(lead0, paste0(trait, '.lead.3col'))
    if(!file.exists(assocf)) stop('missing lead.assoc: ', assocf)
    x <- fread(assocf)
    if(!'trait' %in% names(x)) x[, trait := trait]
    need <- c('trait','lead_chr','lead_snp','lead_bp','effect_allele','other_allele','beta','risk_allele')
    for(nm in need) if(!nm %in% names(x)) x[, (nm) := NA]
    # BED-like file may mix tabs and spaces; parse as general whitespace.
    bed_lines <- readLines(bed, warn = FALSE)
    bed_lines <- trimws(bed_lines)
    bed_lines <- bed_lines[nzchar(bed_lines) & !grepl('^#', bed_lines)]
    sp <- strsplit(bed_lines, '[[:space:]]+')
    sp <- sp[lengths(sp) >= 4L]
    if(!length(sp)) { message('No valid positive loci with at least 4 columns; skip adding positive loci: ', bed); return(invisible(NULL)) }
    loc <- rbindlist(lapply(sp, function(z) data.table(chr=z[1], start=z[2], end=z[3], snp=z[4])), fill=TRUE)
    loc[, chr := toupper(sub('^CHR','',as.character(chr)))]
    loc[, lead_chr := fifelse(chr == 'X', 23L, as.integer(chr))]
    loc[, `:=`(start=as.integer(start), end=as.integer(end), snp=as.character(snp))]
    loc <- loc[!is.na(lead_chr) & nzchar(snp)]
    if(!nrow(loc)){ message('No positive loci to add.'); return(invisible(NULL)) }
    pick1 <- function(nm, cand){ z <- cand[cand %in% nm]; if(length(z)) z[1] else NA_character_ }
    read_gwas <- function(tr){
        cand <- c(file.path(dirgwas,'clean',tr,paste0(tr,'.gz')), file.path(dirgwas,'clean',tr,paste0(tr,'.tsv.gz')), file.path(dirgwas,paste0(tr,'.gz')), file.path(dirgwas,paste0(tr,'.tsv.gz')))
        f <- cand[file.exists(cand)][1]
        if(is.na(f)) return(data.table())
        nm <- names(fread(f, nrows=0))
        cc <- c(SNP=pick1(nm,c('SNP','snp','rsid','RSID','ID','id')), EA=pick1(nm,c('EA','A1','effect_allele','EffectAllele','ALT','alt')), NEA=pick1(nm,c('NEA','A2','other_allele','OtherAllele','REF','ref')), BETA=pick1(nm,c('BETA','beta','b','effect','Effect')))
        if(anyNA(cc)) return(data.table())
        g <- fread(f, select=unname(cc), showProgress=FALSE)
        setnames(g, unname(cc), names(cc))
        g[, .(lead_snp=as.character(SNP), effect_allele=toupper(as.character(EA)), other_allele=toupper(as.character(NEA)), beta=as.numeric(BETA))]
    }
    pvar_files <- function(chr, bp){
        if(chr == 23L){
            par <- if(!is.na(bp) && ((bp>=60001 && bp<=2699520) || (bp>=154931044 && bp<=155260560))) 'par' else 'nonPar'
            if(ref_pop == 'ALL') {
                file.path(dirmod, 'chrX.pvar')
            } else {
                unique(c(file.path(dirmod, paste0(ref_pop,'.male.chrX.',par,'.pvar')), file.path(dirmod, paste0(ref_pop,'.chrX.pvar'))))
            }
        } else {
            if(ref_pop == 'ALL') {
                file.path(dirmod, paste0('chr',chr,'.pvar'))
            } else {
                file.path(dirmod, paste0(ref_pop,'.chr',chr,'.pvar'))
            }
        }
    }
    query_pvar <- function(chr, snp, bp){
        fs <- pvar_files(chr, bp); fs <- fs[file.exists(fs)]
        for(f in fs){
            cmd <- sprintf("awk 'BEGIN{FS=OFS=\"\\t\"} $1!~/^#/ && $3==%s{print $2,$3,$4,$5; exit}' %s", shQuote(snp), shQuote(f))
            z <- tryCatch(system(cmd, intern=TRUE), error=function(e) character())
            if(length(z) && nzchar(z[1])){
                sp <- strsplit(z[1], '\\t')[[1]]
                return(list(bp=as.integer(sp[1]), id=sp[2], ref=toupper(sp[3]), alt=toupper(sp[4]), pvar=f))
            }
        }
        for(f in fs){
            cmd <- sprintf("awk 'BEGIN{FS=OFS=\"\\t\"} $1!~/^#/ && $2==%s{print $2,$3,$4,$5; exit}' %s", shQuote(as.character(bp)), shQuote(f))
            z <- tryCatch(system(cmd, intern=TRUE), error=function(e) character())
            if(length(z) && nzchar(z[1])){
                sp <- strsplit(z[1], '\\t')[[1]]
                id <- if(!is.na(sp[2]) && nzchar(sp[2]) && sp[2] != '.') sp[2] else snp
                return(list(bp=as.integer(sp[1]), id=id, ref=toupper(sp[3]), alt=toupper(sp[4]), pvar=f))
            }
        }
        list(bp=as.integer(bp), id=snp, ref=NA_character_, alt=NA_character_, pvar=NA_character_)
    }
    g <- read_gwas(trait)
    rows <- list(); k <- 0L
    for(i in seq_len(nrow(loc))){
        snp <- loc$snp[i]
        if(any(as.character(x$lead_snp) == snp)) next
        pv <- query_pvar(loc$lead_chr[i], snp, loc$end[i])
        gg <- if(nrow(g)) g[lead_snp == snp][1] else data.table()
        ea <- if(nrow(gg)) gg$effect_allele[1] else pv$alt
        oa <- if(nrow(gg)) gg$other_allele[1] else pv$ref
        be <- if(nrow(gg)) gg$beta[1] else 1
        rk <- if(!is.na(be) && be >= 0) ea else oa
        k <- k + 1L
        row <- as.list(rep(NA, length(names(x)))); names(row) <- names(x)
        row[['trait']] <- trait; row[['lead_chr']] <- loc$lead_chr[i]; row[['lead_snp']] <- snp; row[['lead_bp']] <- pv$bp
        row[['effect_allele']] <- ea; row[['other_allele']] <- oa; row[['beta']] <- be; row[['risk_allele']] <- rk
        if('lead_snp0' %in% names(x)) row[['lead_snp0']] <- snp
        if('match_type_1kg' %in% names(x)) row[['match_type_1kg']] <- 'POSITIVE_CONTROL'
        rows[[k]] <- as.data.table(row)
    }
    if(length(rows)){
        add <- rbindlist(rows, fill=TRUE)
        x <- rbindlist(list(x, add), fill=TRUE)
        fwrite(x, assocf, sep='\t')
        fwrite(x[, .(lead_chr, lead_snp, lead_bp)], threef, sep='\t')
        message('Added positive control loci to ', assocf, ': ', nrow(add))
    } else {
        message('No new positive control loci added to ', assocf)
    }
    invisible(NULL)
}

run_hap_sample_map <- function(script_args) {
    suppressPackageStartupMessages(library(data.table))
    args <- script_args
    if (length(args) < 3) stop("usage: Rscript locus.R hap_sample_map DIR_OUT SAMPLE_FILE OUT_TSV")
    root <- args[1]
    sample_file <- args[2]
    outf <- args[3]
    hapf <- file.path(root, "report", "hap_match.tsv")
    hap <- fread(hapf, select = c("trait", "id", "hap_id", "n", "copies"))
    panel <- if (file.exists(sample_file)) fread(sample_file, fill = TRUE) else data.table()
    if (nrow(panel)) {
        setnames(panel, names(panel)[seq_len(min(3, ncol(panel)))], c("sample", "pop", "super_pop")[seq_len(min(3, ncol(panel)))])
        panel[, sample := as.character(sample)]
    } else {
        panel <- data.table(sample = character(), pop = character(), super_pop = character())
    }
    rows <- rbindlist(lapply(seq_len(nrow(hap)), function(i) {
        cp <- unlist(strsplit(as.character(hap$copies[i]), ";", fixed = TRUE), use.names = FALSE)
        cp <- cp[nzchar(cp)]
        if (!length(cp)) return(NULL)
        data.table(trait = hap$trait[i], id = hap$id[i], hap_id = hap$hap_id[i], hap_n = hap$n[i], copy = cp)
    }), fill = TRUE)
    if (!nrow(rows)) {
        fwrite(data.table(), outf, sep = "\t")
        return(invisible(NULL))
    }
    rows[, `:=`(
        sample_index = suppressWarnings(as.integer(sub("^s([0-9]+)_[12]$", "\\1", copy))),
        haplotype = suppressWarnings(as.integer(sub("^s[0-9]+_([12])$", "\\1", copy)))
    )]
    sample_cache <- new.env(parent = emptyenv())
    get_samples <- function(tr, id) {
        key <- paste(tr, id, sep = "\t")
        if (exists(key, sample_cache, inherits = FALSE)) return(get(key, sample_cache))
        f <- file.path(root, "mat", tr, id, "kg.samples.tsv")
        x <- if (file.exists(f)) fread(f, header = FALSE)[[1]] else character()
        if (!length(x)) {
            vcf <- file.path(root, "coreVcf", tr, id, "kg.vcf.gz")
            if (file.exists(vcf)) {
                con <- gzfile(vcf, "rt")
                repeat {
                    z <- readLines(con, n = 1000, warn = FALSE)
                    if (!length(z)) break
                    h <- grep("^#CHROM", z, value = TRUE)
                    if (length(h)) {
                        x <- strsplit(h[1], "\t", fixed = TRUE)[[1]][-(1:9)]
                        break
                    }
                }
                close(con)
            }
        }
        assign(key, as.character(x), sample_cache)
        as.character(x)
    }
    rows[, sample := {
        ss <- get_samples(trait[1], id[1])
        out <- rep(NA_character_, .N)
        ok <- !is.na(sample_index) & sample_index >= 1L & sample_index <= length(ss)
        out[ok] <- ss[sample_index[ok]]
        out
    }, by = .(trait, id)]
    rows <- merge(rows, panel[, intersect(c("sample", "pop", "super_pop"), names(panel)), with = FALSE],
        by = "sample", all.x = TRUE, sort = FALSE)
    front <- intersect(c("trait", "id", "hap_id", "hap_n", "copy", "sample", "haplotype", "pop", "super_pop"), names(rows))
    rows <- rows[, c(front, setdiff(names(rows), front)), with = FALSE]
    setorder(rows, trait, id, hap_id, sample, haplotype)
    fwrite(rows, outf, sep = "\t")
    invisible(NULL)
}

run_positive_loci_fate <- function(script_args) {
    suppressPackageStartupMessages(library(data.table))
    args <- script_args
    if (length(args) < 3) stop("usage: Rscript locus.R positive_loci_fate DIR_OUT BED OUT_TSV")
    root <- args[1]
    bedf <- args[2]
    outf <- args[3]
    read_tsv <- function(f) if (file.exists(f) && file.info(f)$size > 0) fread(f) else data.table()
    lines <- trimws(readLines(bedf, warn = FALSE))
    lines <- lines[nzchar(lines) & !grepl("^#", lines)]
    sp <- strsplit(lines, "[[:space:]]+")
    sp <- sp[lengths(sp) >= 4L]
    bed <- rbindlist(lapply(sp, function(z) data.table(chr = toupper(sub("^CHR", "", z[1])),
        start = as.integer(z[2]), end = as.integer(z[3]), lead_snp = z[4])), fill = TRUE)
    bed[, lead_chr := fifelse(chr == "X", 23L, suppressWarnings(as.integer(chr)))]
    bed[, locus_id_from_bed := paste(lead_chr, lead_snp, end, sep = ".")]
    pick <- read_tsv(file.path(root, "lead", "pick.tsv"))
    single <- read_tsv(file.path(root, "lead", "pick.single.tsv"))
    pospick <- read_tsv(file.path(root, "lead", "positive_pick.tsv"))
    region <- read_tsv(file.path(root, "report", "region_summary.tsv"))
    selected <- read_tsv(file.path(root, "report", "selected_region.tsv"))
    phy <- read_tsv(file.path(root, "report", "phy_region.tsv"))
    core <- read_tsv(file.path(root, "report", "core_risk.tsv"))
    arch <- read_tsv(file.path(root, "report", "core_archaic_match.tsv"))
    mkid <- function(x) {
        if (!nrow(x)) return(x)
        if (!"id" %in% names(x) && all(c("lead_chr","lead_snp","lead_bp") %in% names(x))) {
            x[, id := paste(lead_chr, lead_snp, lead_bp, sep = ".")]
        }
        x
    }
    pick <- mkid(pick); single <- mkid(single); pospick <- mkid(pospick)
    region <- mkid(region); selected <- mkid(selected); phy <- mkid(phy); core <- mkid(core); arch <- mkid(arch)
    out <- copy(bed)
    out[, trait := if (nrow(pick) && "trait" %in% names(pick)) pick$trait[1] else NA_character_]
    all_loci <- rbindlist(list(
        if (nrow(pick)) pick[, .(id, lead_chr = as.integer(lead_chr), lead_snp = as.character(lead_snp), lead_bp = as.integer(lead_bp))] else NULL,
        if (nrow(single)) single[, .(id, lead_chr = as.integer(lead_chr), lead_snp = as.character(lead_snp), lead_bp = as.integer(lead_bp))] else NULL,
        if (nrow(region)) region[, .(id, lead_chr = as.integer(lead_chr), lead_snp = as.character(lead_snp), lead_bp = as.integer(lead_bp))] else NULL,
        if (nrow(core)) unique(core[, .(id, lead_chr = as.integer(lead_chr), lead_snp = as.character(lead_snp), lead_bp = as.integer(lead_bp))]) else NULL
    ), fill = TRUE)
    all_loci <- unique(all_loci[!is.na(lead_chr) & !is.na(lead_bp) & nzchar(lead_snp)], by = c("lead_chr", "lead_snp", "lead_bp", "id"))
    out[, id := vapply(seq_len(.N), function(i) {
        hit <- all_loci[lead_chr == out$lead_chr[i] & lead_snp == out$lead_snp[i] &
            lead_bp >= out$start[i] & lead_bp <= out$end[i]]
        if (nrow(hit)) hit$id[1] else out$locus_id_from_bed[i]
    }, character(1))]
    for (nm in c("in_pick", "in_single_snp", "in_positive_pick", "has_hap_region", "in_region_summary", "in_selected_region", "in_phy_region")) out[, (nm) := FALSE]
    flag <- function(ids) out$id %in% ids
    out[, in_pick := flag(unique(pick$id))]
    out[, in_single_snp := flag(unique(single$id))]
    out[, in_positive_pick := flag(unique(pospick$id))]
    hap_region_ids <- sub("\\.region\\.tsv$", "", basename(list.files(file.path(root, "hap"), pattern = "\\.region\\.tsv$", recursive = TRUE, full.names = TRUE)))
    out[, has_hap_region := flag(hap_region_ids)]
    out[, in_region_summary := flag(unique(region$id))]
    out[, in_selected_region := flag(unique(selected$id))]
    out[, in_phy_region := flag(unique(phy$id))]
    yri <- if (nrow(region) && "YRI_max_freq" %in% names(region)) region[, .(id, YRI_max_freq)] else
        if (nrow(region) && "yri_freq" %in% names(region)) region[, .(id, YRI_max_freq = yri_freq)] else data.table(id=character(), YRI_max_freq=numeric())
    best <- if (nrow(region)) region[, intersect(c("id","p_ils","best_lineage","matched_archaics"), names(region)), with=FALSE] else data.table(id=character())
    out <- merge(out, unique(best, by = "id"), by = "id", all.x = TRUE, sort = FALSE)
    out <- merge(out, unique(yri, by = "id"), by = "id", all.x = TRUE, sort = FALSE)
    arch_sum <- if (nrow(arch)) arch[, .(max_prop_match_risk = max(prop_match_risk, na.rm = TRUE),
        max_n_match_risk = max(n_match_risk, na.rm = TRUE)), by = id] else data.table(id=character())
    core_sum <- if (nrow(core) && "yri_risk_freq" %in% names(core)) core[!is.na(risk_core_allele),
        .(max_yri_risk_freq = max(yri_risk_freq, na.rm = TRUE)), by = id] else data.table(id=character())
    out <- merge(out, arch_sum, by = "id", all.x = TRUE, sort = FALSE)
    out <- merge(out, core_sum, by = "id", all.x = TRUE, sort = FALSE)
    yri_freq_th <- suppressWarnings(as.numeric(Sys.getenv("YRI_FREQ_TH", "0.05")))
    yri_msg <- if (is.finite(yri_freq_th) && yri_freq_th >= 1) "s6 simple YRI filter disabled" else
        paste0("s6 simple YRI filter: archaic-like region found but YRI_max_freq > ", 
            if (is.finite(yri_freq_th)) yri_freq_th else 0.05, " or missing")
    out[, filter_step := fifelse(!in_pick & in_single_snp, "s2/s3: single-SNP LD block excluded before haplotype analysis",
        fifelse(!in_pick, "s2/s3: not selected as a multi-SNP LD block",
        fifelse(in_pick & !has_hap_region, "s6: forced into pick.tsv but absent from LD block.tsv, so make_hap wrote no region",
        fifelse(has_hap_region & !in_region_summary, "s6 merge: no matched archaic haplotype/lineage retained",
        fifelse(in_region_summary & !in_selected_region, yri_msg,
        fifelse(in_selected_region & !in_phy_region, "s7: selected but not sent to phylogeny", "kept for plot/phylogeny"))))))]
    setcolorder(out, c("id","chr","start","end","lead_snp","filter_step","in_pick","in_single_snp","in_positive_pick","has_hap_region","in_region_summary","in_selected_region","in_phy_region",
        setdiff(names(out), c("id","chr","start","end","lead_snp","filter_step","in_pick","in_single_snp","in_positive_pick","has_hap_region","in_region_summary","in_selected_region","in_phy_region"))))
    fwrite(out, outf, sep = "\t")
    invisible(NULL)
}


run_viome <- function(script_args) {
    suppressPackageStartupMessages(library(data.table))
    args <- script_args
    arg <- function(k, multi = FALSE) {
        i <- which(args %in% k)[1]
        if (is.na(i) || i == length(args)) return(NULL)
        j <- i + 1
        if (!multi) return(args[j])
        z <- args[j:length(args)]
        z[seq_len(which(c(grepl("^-", z[-1]), TRUE))[1])]
    }
    pick1 <- function(nm, cand) {
        x <- cand[cand %chin% nm]
        if (length(x)) x[1] else NA_character_
    }
    chr_int <- function(x) {
        x <- toupper(sub("^CHR", "", as.character(x)))
        x[x == "X"] <- "23"
        x[x == "Y"] <- "24"
        suppressWarnings(as.integer(x))
    }
    `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a
    clean_chr <- function(x) {
        y <- toupper(sub("^CHR", "", as.character(x)))
        y[y == "23"] <- "X"
        y[y == "24"] <- "Y"
        y
    }
    norm_variant_id <- function(chr, pos, ref, alt, with_chr = TRUE) {
        chr <- clean_chr(chr)
        if (with_chr) chr <- paste0("chr", chr)
        paste(chr, as.integer(pos), toupper(as.character(ref)), toupper(as.character(alt)), sep = ":")
    }
    parse_variant_id <- function(x) {
        x0 <- as.character(x)
        x0 <- sub("^chr", "", x0, ignore.case = TRUE)
        sp <- tstrsplit(x0, ":", fixed = TRUE)
        if (length(sp) < 4) {
            return(data.table(var_chr = NA_integer_, var_pos = NA_integer_, var_ref = NA_character_, var_alt = NA_character_))
        }
        data.table(var_chr = chr_int(sp[[1]]), var_pos = suppressWarnings(as.integer(sp[[2]])),
                   var_ref = toupper(sp[[3]]), var_alt = toupper(sp[[4]]))
    }
    first_existing <- function(x) {
        x[file.exists(x)][1]
    }
    gwas_file <- function(dirgwas, tr) {
        first_existing(c(
            file.path(dirgwas, "clean", tr, paste0(tr, ".gz")),
            file.path(dirgwas, "clean", paste0(tr, ".gz")),
            file.path(dirgwas, tr, paste0(tr, ".gz")),
            file.path(dirgwas, paste0(tr, ".gz")),
            file.path(dirgwas, "clean", tr, paste0(tr, ".tsv.gz")),
            file.path(dirgwas, paste0(tr, ".tsv.gz")),
            file.path(dirgwas, paste0(tr, ".txt.gz"))
        ))
    }
    read_gwas_header <- function(f) names(fread(f, nrows = 0, showProgress = FALSE))
    read_gwas_viome <- function(f) {
        h <- read_gwas_header(f)
        cc <- c(
            SNP = pick1(h, c("SNP", "rsid", "RSID", "ID", "id", "variant", "variant_id", "MarkerName")),
            CHR = pick1(h, c("CHR", "Chr", "chr", "#CHROM", "CHROM")),
            POS = pick1(h, c("POS", "BP", "bp", "pos", "position")),
            REF = pick1(h, c("REF", "ref", "A2", "NEA", "other_allele", "Allele2")),
            ALT = pick1(h, c("ALT", "alt", "A1", "EA", "effect_allele", "Allele1")),
            BETA = pick1(h, c("BETA", "beta", "b", "Effect", "effect", "Estimate", "estimate")),
            P = pick1(h, c("P_BOLT_LMM", "P_BOLT_LMM_INF", "P", "p", "PVAL", "pval", "p_value", "P_VALUE")),
            AF = pick1(h, c("A1FREQ", "A1_FREQ", "AF", "ALT_FREQ", "EAF", "effect_allele_frequency", "freq"))
        )
        need <- c("SNP", "CHR", "POS", "REF", "ALT", "P")
        if (any(is.na(cc[need]))) stop("GWAS missing required columns SNP/CHR/POS/REF/ALT/P in ", f, ". Header: ", paste(h, collapse = ","))
        keep <- unique(na.omit(unname(cc)))
        x <- fread(f, select = keep, showProgress = FALSE)
        for (nm in names(cc)) if (!is.na(cc[nm]) && cc[nm] %chin% names(x)) setnames(x, cc[nm], nm)
        if (!"BETA" %in% names(x)) x[, BETA := NA_real_]
        if (!"AF" %in% names(x)) x[, AF := NA_real_]
        x[, `:=`(SNP = as.character(SNP), CHR = chr_int(CHR), POS = as.integer(POS),
                 REF = toupper(as.character(REF)), ALT = toupper(as.character(ALT)),
                 BETA = suppressWarnings(as.numeric(BETA)), P = suppressWarnings(as.numeric(P)),
                 AF = suppressWarnings(as.numeric(AF)))]
        x[!is.na(CHR) & !is.na(POS) & !is.na(P) & nzchar(REF) & nzchar(ALT),
          `:=`(varid_chr = norm_variant_id(CHR, POS, REF, ALT, TRUE),
               varid_nochr = norm_variant_id(CHR, POS, REF, ALT, FALSE))]
        x
    }
    standardize_asnp <- function(f) {
        a <- fread(f, showProgress = FALSE)
        nm <- names(a)
        id_col <- pick1(nm, c("aSNP", "asnp", "variant", "variant_id", "id", "ID", "SNP", "snp"))
        if (is.na(id_col) && ncol(a) >= 3L) id_col <- nm[3]
        hap_col <- pick1(nm, c("chr_st_end", "haplotype", "hap_id", "hapID", "segment", "region", "introgressed_haplotype"))
        sup_col <- pick1(nm, c("SupPop", "super_pop", "superpop", "population", "Population", "POP", "pop"))
        freq_col <- pick1(nm, c("haplo_freq", "haplotype_freq", "archaic_allele_frequency", "archaic_freq", "freq", "AF", "Frequency"))
        lin_col <- pick1(nm, c("lineage", "Lineage", "archaic", "archaic_genome", "closest", "best_lineage"))
        n_col <- pick1(nm, c("n_asnp", "n_aSNP", "nSNP", "n_snp", "n_marker"))
        if (is.na(id_col)) stop("Cannot identify aSNP variant-id column in ", f)
        a[, asnp_id_raw := as.character(get(id_col))]
        p <- parse_variant_id(a$asnp_id_raw)
        a <- cbind(a, p)
        a[, `:=`(asnp_id_chr = norm_variant_id(var_chr, var_pos, var_ref, var_alt, TRUE),
                 asnp_id_nochr = norm_variant_id(var_chr, var_pos, var_ref, var_alt, FALSE))]
        if (is.na(hap_col)) a[, viome_haplotype := paste0("hap_", asnp_id_chr)] else a[, viome_haplotype := as.character(get(hap_col))]
        if (is.na(sup_col)) a[, SupPop_viome := NA_character_] else a[, SupPop_viome := as.character(get(sup_col))]
        if (is.na(freq_col)) a[, viome_freq := NA_real_] else a[, viome_freq := suppressWarnings(as.numeric(get(freq_col)))]
        if (is.na(lin_col)) a[, viome_lineage := NA_character_] else a[, viome_lineage := as.character(get(lin_col))]
        if (is.na(n_col)) {
            n_by_hap <- a[!is.na(var_pos), .(viome_n_asnp_haplotype = uniqueN(asnp_id_chr)), by = viome_haplotype]
            a <- merge(a, n_by_hap, by = "viome_haplotype", all.x = TRUE, sort = FALSE)
        } else {
            a[, viome_n_asnp_haplotype := suppressWarnings(as.integer(get(n_col)))]
        }
        a[is.na(viome_n_asnp_haplotype), viome_n_asnp_haplotype := 1L]
        a
    }
    read_existing_ld <- function(root, tr) {
        fs <- list.files(file.path(root, "ld", tr), pattern = "ld\\.tsv$", recursive = TRUE, full.names = TRUE)
        if (!length(fs)) return(data.table())
        rbindlist(lapply(fs, function(f) tryCatch(fread(f, showProgress = FALSE), error = function(e) data.table())), fill = TRUE)
    }
    attach_ld <- function(hit, ld) {
        hit[, viome_ld_r2_with_lead := NA_real_]
        hit[, viome_ld_source := "not_available"]
        if (!nrow(hit) || !nrow(ld)) return(hit)
        nm <- names(ld)
        if (!all(c("lead_snp", "snp", "R2") %chin% nm)) return(hit)
        z <- ld[, .(lead_snp = as.character(lead_snp), POS = suppressWarnings(as.integer(pos)), snp = as.character(snp), R2 = suppressWarnings(as.numeric(R2)))]
        hit <- merge(hit, z[, .(lead_snp, POS, viome_ld_r2_with_lead_tmp = R2)], by = c("lead_snp", "POS"), all.x = TRUE, sort = FALSE)
        hit[!is.na(viome_ld_r2_with_lead_tmp), `:=`(viome_ld_r2_with_lead = viome_ld_r2_with_lead_tmp, viome_ld_source = "existing_locus_ld")]
        hit[, viome_ld_r2_with_lead_tmp := NULL]
        hit
    }
    write_empty <- function(report0) {
        empty_hits <- data.table(trait=character(), lead_chr=integer(), lead_snp=character(), lead_bp=integer(), SNP=character(), CHR=integer(), POS=integer(), REF=character(), ALT=character(), BETA=numeric(), P=numeric(), AF=numeric(), asnp_id=character(), viome_haplotype=character(), viome_n_asnp_haplotype=integer(), viome_freq=numeric(), SupPop_viome=character(), viome_lineage=character(), viome_ld_r2_with_lead=numeric(), viome_ld_source=character())
        fwrite(empty_hits, file.path(report0, "viome_aSNP_hits.tsv"), sep = "\t")
        fwrite(empty_hits, file.path(report0, "viome_inherited_segments.tsv"), sep = "\t")
        fwrite(data.table(), file.path(report0, "viome_region_summary.tsv"), sep = "\t")
    }
    make_network_one <- function(root, dirmod, row, a_sub, ld_r2) {
        if (!requireNamespace("pegas", quietly = TRUE) || !requireNamespace("ape", quietly = TRUE)) return("skip: R packages pegas/ape not installed")
        chr_lab <- clean_chr(row$lead_chr[1])
        st <- max(1L, suppressWarnings(as.integer(min(a_sub$var_pos, na.rm = TRUE))))
        en <- suppressWarnings(as.integer(max(a_sub$var_pos, na.rm = TRUE)))
        if (!is.finite(st) || !is.finite(en) || st >= en) return("skip: invalid aSNP region")
        vcf_candidates <- c(file.path(dirmod, "vcf", paste0("chr", chr_lab, ".vcf.gz")), file.path(dirmod, paste0("chr", chr_lab, ".vcf.gz")), file.path(dirname(dirmod), "vcf", paste0("chr", chr_lab, ".vcf.gz")))
        vcf <- vcf_candidates[file.exists(vcf_candidates)][1]
        if (is.na(vcf)) return("skip: 1000G VCF not found")
        if (Sys.which("bcftools") == "") return("skip: bcftools not found")
        dir.create(file.path(root, "viome"), recursive = TRUE, showWarnings = FALSE)
        prefix <- file.path(root, "viome", paste(row$trait[1], row$lead_chr[1], row$lead_snp[1], row$lead_bp[1], sep = "."))
        gt_file <- paste0(prefix, ".network_genotypes.tsv")
        cmd <- sprintf("bcftools view -v snps -r %s:%d-%d %s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n' > %s", shQuote(paste0("chr", chr_lab)), st, en, shQuote(vcf), shQuote(gt_file))
        system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
        if (!file.exists(gt_file) || file.info(gt_file)$size == 0) return("skip: no genotype rows extracted")
        gen <- tryCatch(fread(gt_file, header = FALSE, showProgress = FALSE), error = function(e) data.table())
        if (nrow(gen) < 2 || ncol(gen) < 6) return("skip: too few variants/samples")
        gen <- gen[nchar(V3) == 1 & nchar(V4) == 1 & V3 %chin% c("A","C","G","T") & V4 %chin% c("A","C","G","T")]
        if (nrow(gen) < 2) return("skip: too few biallelic SNVs")
        gt <- as.matrix(gen[, -(1:4), with = FALSE])
        hap_strings <- character(0)
        for (j in seq_len(ncol(gt))) {
            g <- gt[, j]
            sp <- tstrsplit(gsub("/", "|", g), "|", fixed = TRUE)
            if (length(sp) < 2) next
            h1 <- ifelse(sp[[1]] == "1", gen$V4, ifelse(sp[[1]] == "0", gen$V3, "N"))
            h2 <- ifelse(sp[[2]] == "1", gen$V4, ifelse(sp[[2]] == "0", gen$V3, "N"))
            hap_strings <- c(hap_strings, paste(h1, collapse = ""), paste(h2, collapse = ""))
        }
        hap_strings <- hap_strings[nchar(hap_strings) == nrow(gen) & !grepl("N", hap_strings)]
        if (length(unique(hap_strings)) < 2) return("skip: fewer than two haplotypes")
        tab <- sort(table(hap_strings), decreasing = TRUE)
        core <- data.table(core_haplotype = names(tab), n = as.integer(tab))
        fwrite(core, paste0(prefix, ".core_haplotypes.tsv"), sep = "\t")
        mat <- do.call(rbind, strsplit(names(tab), ""))
        rownames(mat) <- paste0("H", seq_len(nrow(mat)), "_n", as.integer(tab))
        pdf_file <- paste0(prefix, ".haplotype_network.pdf")
        tryCatch({
            grDevices::pdf(pdf_file, width = 7, height = 6)
            h <- pegas::haplotype(ape::as.DNAbin(mat))
            net <- pegas::haploNet(h)
            plot(net, size = sqrt(as.integer(tab)) + 1, show.mutation = 1, main = paste0(row$trait[1], " ", row$lead_snp[1]))
            grDevices::dev.off()
        }, error = function(e) {
            try(grDevices::dev.off(), silent = TRUE)
        })
        if (file.exists(pdf_file)) paste0("ok: ", pdf_file) else "skip: pegas plotting failed"
    }

    dirgwas <- arg(c("--dirgwas", "-dirgwas"))
    dircojo <- arg(c("--dircojo", "-dircojo"))
    dirout <- arg(c("--dirout", "-dirout"))
    dirmod <- arg(c("--dirmod", "-dirmod"))
    asnp_file <- arg(c("--asnp", "-asnp"))
    traits <- arg(c("--traits", "-traits"), TRUE)
    p_th <- suppressWarnings(as.numeric(arg(c("--p-th", "-p-th")) %||% "1e-9"))
    freq_th <- suppressWarnings(as.numeric(arg(c("--freq-th", "-freq-th")) %||% "0.01"))
    min_asnp <- suppressWarnings(as.integer(arg(c("--min-asnp", "-min-asnp")) %||% "5"))
    window_kb <- suppressWarnings(as.numeric(arg(c("--window-kb", "-window-kb")) %||% "1000"))
    ld_r2 <- suppressWarnings(as.numeric(arg(c("--ld-r2", "-ld-r2")) %||% "0.9"))
    make_network <- suppressWarnings(as.integer(arg(c("--make-network", "-make-network")) %||% "1"))
    if (is.null(dirgwas) || is.null(dirout) || is.null(asnp_file) || is.null(traits)) stop("usage: Rscript locus.R viome --dirgwas DIR --dirout DIR --dirmod DIR --asnp FILE --traits trait ...")
    traits <- unlist(strsplit(paste(traits, collapse = " "), "[, ]+")); traits <- traits[nzchar(traits)]
    report0 <- file.path(dirout, "report")
    dir.create(report0, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(dirout, "viome"), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(asnp_file)) stop("viome aSNP map not found: ", asnp_file)
    lead_file <- file.path(dirout, "lead", "lead.assoc")
    if (!file.exists(lead_file)) stop("viome requires prepared lead.assoc: ", lead_file)
    lead <- fread(lead_file, showProgress = FALSE)
    if (!nrow(lead)) { write_empty(report0); return(invisible(NULL)) }
    if (!all(c("trait", "lead_chr", "lead_snp", "lead_bp") %chin% names(lead))) stop("lead.assoc must contain trait, lead_chr, lead_snp, lead_bp")
    lead[, `:=`(lead_chr = as.integer(lead_chr), lead_bp = as.integer(lead_bp), lead_snp = as.character(lead_snp))]
    lead <- lead[trait %chin% traits]
    a <- standardize_asnp(asnp_file)
    a_key <- unique(a[!is.na(var_chr) & !is.na(var_pos), .(asnp_id_chr, asnp_id_nochr, viome_haplotype, viome_n_asnp_haplotype, viome_freq, SupPop_viome, viome_lineage, var_chr, var_pos, var_ref, var_alt)])
    all_hits <- list(); region_sum <- list(); network_log <- list()
    for (tr in traits) {
        gf <- gwas_file(dirgwas, tr)
        if (is.na(gf) || !file.exists(gf)) {
            region_sum[[length(region_sum)+1L]] <- data.table(trait = tr, status = "missing_gwas", message = paste("No GWAS file found under", dirgwas))
            next
        }
        message("viome: reading GWAS ", tr, " from ", gf)
        g <- read_gwas_viome(gf)
        ld <- read_existing_ld(dirout, tr)
        leads_tr <- lead[trait == tr]
        for (ii in seq_len(nrow(leads_tr))) {
            l <- leads_tr[ii]
            st <- max(1L, as.integer(l$lead_bp - window_kb * 1000))
            en <- as.integer(l$lead_bp + window_kb * 1000)
            reg <- g[CHR == l$lead_chr & POS >= st & POS <= en]
            if (!nrow(reg)) {
                region_sum[[length(region_sum)+1L]] <- data.table(trait=tr, lead_chr=l$lead_chr, lead_snp=l$lead_snp, lead_bp=l$lead_bp, start=st, end=en, n_gwas=0L, n_asnp_hits=0L, n_significant_asnp=0L, best_p=NA_real_, best_asnp=NA_character_, best_haplotype=NA_character_, status="no_gwas_rows")
                next
            }
            h1 <- merge(reg, a_key, by.x = "varid_chr", by.y = "asnp_id_chr", all = FALSE, sort = FALSE, allow.cartesian = TRUE)
            h2 <- merge(reg, a_key, by.x = "varid_nochr", by.y = "asnp_id_nochr", all = FALSE, sort = FALSE, allow.cartesian = TRUE)
            hit <- unique(rbindlist(list(h1, h2), fill = TRUE), by = c("SNP", "CHR", "POS", "REF", "ALT", "viome_haplotype", "SupPop_viome"))
            if (nrow(hit)) {
                hit[, `:=`(trait = tr, lead_chr = l$lead_chr, lead_snp = l$lead_snp, lead_bp = l$lead_bp,
                           region_start = st, region_end = en,
                           asnp_id = norm_variant_id(CHR, POS, REF, ALT, TRUE),
                           viome_test_freq = fifelse(!is.na(viome_freq), viome_freq, AF),
                           pass_p = P <= p_th,
                           pass_freq = is.na(viome_test_freq) | viome_test_freq >= freq_th,
                           pass_min_asnp = viome_n_asnp_haplotype >= min_asnp)]
                hit <- attach_ld(hit, ld)
                hit[, pass_viome := pass_p & pass_freq & pass_min_asnp]
                viome_cols <- c("trait","lead_chr","lead_snp","lead_bp","region_start","region_end","SNP","CHR","POS","REF","ALT","BETA","P","AF","asnp_id","viome_haplotype","viome_n_asnp_haplotype","viome_freq","viome_test_freq","SupPop_viome","viome_lineage","viome_ld_r2_with_lead","viome_ld_source","pass_p","pass_freq","pass_min_asnp","pass_viome")
                setcolorder(hit, c(intersect(viome_cols, names(hit)), setdiff(names(hit), viome_cols)))
                all_hits[[length(all_hits)+1L]] <- hit
                best <- hit[order(P)][1]
                sig <- hit[pass_viome == TRUE]
                region_sum[[length(region_sum)+1L]] <- data.table(trait=tr, lead_chr=l$lead_chr, lead_snp=l$lead_snp, lead_bp=l$lead_bp, start=st, end=en,
                    n_gwas=nrow(reg), n_asnp_hits=nrow(hit), n_significant_asnp=nrow(sig), best_p=best$P, best_asnp=best$asnp_id, best_haplotype=best$viome_haplotype,
                    best_lineage=best$viome_lineage, best_freq=best$viome_freq, best_n_asnp_haplotype=best$viome_n_asnp_haplotype,
                    max_ld_r2_with_lead=suppressWarnings(max(hit$viome_ld_r2_with_lead, na.rm = TRUE)), status=ifelse(nrow(sig)>0, "viome_candidate", "asnp_overlap_no_pass"))
                if (make_network == 1L && nrow(sig)) {
                    a_sub <- a[viome_haplotype %chin% unique(sig$viome_haplotype)]
                    msg <- tryCatch(make_network_one(dirout, dirmod %||% "", l, a_sub, ld_r2), error = function(e) paste("skip:", conditionMessage(e)))
                    network_log[[length(network_log)+1L]] <- data.table(trait=tr, lead_snp=l$lead_snp, lead_bp=l$lead_bp, message=msg)
                }
            } else {
                region_sum[[length(region_sum)+1L]] <- data.table(trait=tr, lead_chr=l$lead_chr, lead_snp=l$lead_snp, lead_bp=l$lead_bp, start=st, end=en, n_gwas=nrow(reg), n_asnp_hits=0L, n_significant_asnp=0L, best_p=NA_real_, best_asnp=NA_character_, best_haplotype=NA_character_, status="no_asnp_overlap")
            }
        }
    }
    hits <- rbindlist(all_hits, fill = TRUE)
    if (!nrow(hits)) {
        write_empty(report0)
    } else {
        fwrite(hits, file.path(report0, "viome_aSNP_hits.tsv"), sep = "\t")
        fwrite(hits[pass_viome == TRUE], file.path(report0, "viome_inherited_segments.tsv"), sep = "\t")
    }
    rs <- rbindlist(region_sum, fill = TRUE)
    if (nrow(rs) && "max_ld_r2_with_lead" %in% names(rs)) rs[!is.finite(max_ld_r2_with_lead), max_ld_r2_with_lead := NA_real_]
    fwrite(rs, file.path(report0, "viome_region_summary.tsv"), sep = "\t")
    nl <- rbindlist(network_log, fill = TRUE)
    fwrite(nl, file.path(report0, "viome_network_log.tsv"), sep = "\t")
    message("viome done: ", nrow(hits), " aSNP overlap rows; ", if (nrow(hits)) nrow(hits[pass_viome == TRUE]) else 0L, " rows pass viome filters")
    invisible(NULL)
}

args <- commandArgs(TRUE)
if (length(args) == 0L) usage()
cmd <- args[1]
cmd_args <- args[-1]
switch(cmd,
    prep_input = run_prep_input(cmd_args),
    prep_input_local = run_prep_input_local(cmd_args),
    add_positive_loci = run_add_positive_loci(cmd_args),
    prep_archaic = run_prep_archaic(cmd_args),
    prep_input_archaic = run_prep_archaic(cmd_args),
    make_hap = run_make_hap(cmd_args),
    filter_hap = run_filter_hap(cmd_args),
    make_phy = run_make_phy(cmd_args),
    make_tree = run_make_tree(cmd_args),
    viome = run_viome(cmd_args),
    hap_sample_map = run_hap_sample_map(cmd_args),
    positive_loci_fate = run_positive_loci_fate(cmd_args),
    usage()
)
