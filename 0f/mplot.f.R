# Reusable self Manhattan plotting core.
#
# Public API:
#   mh_plot(gwas_file, ..., mirror = magma_data)
#
# gwas_file must contain CHR/POS/P.  The optional mirror data frame must contain
# GENE/CHR/START/STOP/P and is drawn below Y=0 using the same cumulative genomic
# x coordinates.  cis/trans locus annotation remains available through cis_bed
# and mh_plot_bed.

library(lattice)
.mh_plot_vectors<-function(chr, pos, pvalue, 
	sig.level=NA, annotate=NULL, ann.default=list(),
	mirror=NULL, mirror.sig.level=NA, mirror.label.n=10L,
	mirror.col=c("gray", "darkgray"), mirror.sig.col="green",
	mirror.cis.col="red", threshold.col="skyblue", mirror.threshold.col="skyblue",
	chr.label.col="skyblue",
	should.thin=T, thin.pos.places=2, thin.logp.places=2, 
	xlab="Chromosome", ylab=expression(-log[10](p-value)),
	col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {

	if (length(chr)==0) stop("chromosome vector is empty")
	if (length(pos)==0) stop("position vector is empty")
	if (length(pvalue)==0) stop("pvalue vector is empty")

	# Normalise the optional lower, mirrored data set.  Its genomic regions are
	# projected into the same cumulative chromosome coordinates as the SNPs.
	mirror.plot <- NULL
	mirror.bottom <- 0
	if (!is.null(mirror)) {
		mirror.plot <- as.data.frame(mirror, stringsAsFactors=FALSE)
		need.mirror <- c("GENE", "CHR", "START", "STOP", "P")
		if (!all(need.mirror %in% names(mirror.plot))) stop("mirror must contain GENE, CHR, START, STOP, and P")
		mirror.plot$CHR <- suppressWarnings(as.integer(mirror.plot$CHR))
		mirror.plot$START <- suppressWarnings(as.numeric(mirror.plot$START))
		mirror.plot$STOP <- suppressWarnings(as.numeric(mirror.plot$STOP))
		mirror.plot$P <- suppressWarnings(as.numeric(mirror.plot$P))
		mirror.plot <- mirror.plot[is.finite(mirror.plot$CHR) & is.finite(mirror.plot$START) &
			is.finite(mirror.plot$STOP) & is.finite(mirror.plot$P) & mirror.plot$P > 0 & mirror.plot$P <= 1, , drop=FALSE]
		if (!nrow(mirror.plot)) mirror.plot <- NULL
	}

	# Use the union of chromosomes so both halves share exactly the same x scale.
	chr.values <- suppressWarnings(as.integer(as.character(chr)))
	chr.levels <- sort(unique(c(chr.values, if (!is.null(mirror.plot)) mirror.plot$CHR else integer())))
	chr.levels <- as.character(chr.levels[is.finite(chr.levels)])
	chr <- ordered(as.character(chr.values), levels=chr.levels)

	#make sure positions are in kbp
	if (any(pos>1e6)) pos<-pos/1e6;

	#calculate absolute genomic position
	#from relative chromosomal positions
	posmin <- tapply(pos,chr, min);
	posmax <- tapply(pos,chr, max);
	if (!is.null(mirror.plot)) {
		if (any(c(mirror.plot$START, mirror.plot$STOP)>1e6)) {
			mirror.plot$START <- mirror.plot$START/1e6
			mirror.plot$STOP <- mirror.plot$STOP/1e6
		}
		mirror.min <- tapply(mirror.plot$START, factor(mirror.plot$CHR, levels=as.integer(chr.levels)), min)
		mirror.max <- tapply(mirror.plot$STOP, factor(mirror.plot$CHR, levels=as.integer(chr.levels)), max)
		posmin <- pmin(posmin, mirror.min, na.rm=TRUE)
		posmax <- pmax(posmax, mirror.max, na.rm=TRUE)
	}
	posmin[!is.finite(posmin)] <- 0
	posmax[!is.finite(posmax)] <- posmin[!is.finite(posmax)]
	posshift <- head(c(0,cumsum(posmax)),-1);
	names(posshift) <- chr.levels
	genpos <- pos + posshift[chr];
	getGenPos<-function(cchr, cpos) {
		p<-posshift[as.character(cchr)]+cpos
		return(p)
	}
	if (!is.null(mirror.plot)) {
		mirror.plot$X1 <- mirror.plot$START + posshift[as.character(mirror.plot$CHR)]
		mirror.plot$X2 <- mirror.plot$STOP + posshift[as.character(mirror.plot$CHR)]
		mirror.plot$XMID <- (mirror.plot$X1 + mirror.plot$X2)/2
		mirror.plot$Y <- -log10(mirror.plot$P)
		if (length(mirror.sig.level) != 1 || !is.finite(mirror.sig.level) || mirror.sig.level <= 0 || mirror.sig.level > 1)
			mirror.sig.level <- 2.5e-6
		if (!"CLASS" %in% names(mirror.plot))
			mirror.plot$CLASS <- ifelse(mirror.plot$P <= mirror.sig.level, "trans", "background")
		mirror.plot$CLASS <- tolower(as.character(mirror.plot$CLASS))
		mirror.plot$CLASS[!mirror.plot$CLASS %in% c("background", "trans", "cis")] <- "background"
		mirror.plot$CLASS[mirror.plot$P > mirror.sig.level] <- "background"
		mirror.bottom <- -(ceiling(max(mirror.plot$Y, -log10(mirror.sig.level), na.rm=TRUE)) + 4)
	}

	#parse annotations
	grp <- NULL
	ann.settings <- list()
	label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5, 
		col=NULL, fontface=NULL, fontsize=NULL, show=F)
	parse.label<-function(rawval, groupname) {
		r<-list(text=groupname)
		if(is.logical(rawval)) {
			if(!rawval) {r$show <- F}
		} else if (is.character(rawval) || is.expression(rawval)) {
			if(nchar(rawval)>=1) {
				r$text <- rawval
			}
		} else if (is.list(rawval)) {
			r <- modifyList(r, rawval)
		}
		return(r)
	}

	if(!is.null(annotate)) {
		if (is.list(annotate)) {
			grp <- annotate[[1]]
		} else {
			grp <- annotate
		} 
		if (!is.factor(grp)) {
			grp <- factor(grp)
		}
	} else {
		grp <- factor(rep(1, times=length(pvalue)))
	}
  
	ann.settings<-vector("list", length(levels(grp)))
	ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)

	if (length(ann.settings)>1) { 
		lcols<-trellis.par.get("superpose.symbol")$col 
		lfills<-trellis.par.get("superpose.symbol")$fill
		for(i in 2:length(levels(grp))) {
			ann.settings[[i]]<-list(pch=pch, 
				col=lcols[(i-2) %% length(lcols) +1 ], 
				fill=lfills[(i-2) %% length(lfills) +1 ], 
				cex=cex, label=label.default);
			ann.settings[[i]]$label$show <- T
		}
	}
	names(ann.settings)<-levels(grp)
	for(i in 1:length(ann.settings)) {
		if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
		ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
			parse.label(ann.settings[[i]]$label, levels(grp)[i]))
	}
	if(is.list(annotate) && length(annotate)>1) {
		user.cols <- 2:length(annotate)
		ann.cols <- c()
		if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
			ann.cols<-match(names(annotate)[-1], names(ann.settings))
		} else {
			ann.cols<-user.cols-1
		}
		for(i in seq_along(user.cols)) {
			if(!is.null(annotate[[user.cols[i]]]$label)) {
				annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
					levels(grp)[ann.cols[i]])
			}
			ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
				annotate[[user.cols[i]]])
		}
	}
 	rm(annotate)

	#reduce number of points plotted
	if(should.thin) {
		thinned <- unique(data.frame(
			logp=round(-log10(pvalue),thin.logp.places), 
			pos=round(genpos,thin.pos.places), 
			chr=chr,
			grp=grp)
		)
		logp <- thinned$logp
		genpos <- thinned$pos
		chr <- thinned$chr
		grp <- thinned$grp
		rm(thinned)
	} else {
		logp <- -log10(pvalue)
	}
	rm(pos, pvalue)
	gc()

	# Print chromosome names on the ordinary bottom axis only when there is no
	# lower panel.  With a mirrored panel, labels are drawn at the shared Y=0
	# baseline so they stay visually attached to the GWAS panel.
	axis.chr <- function(side,...) {
		if(side=="bottom" && is.null(mirror.plot)) {
			chr.labels <- chr.levels
			chr.labels[chr.labels == "23"] <- "X"
			chr.labels[chr.labels == "24"] <- "Y"
			panel.axis(side=side, outside=T,
				at=((posmax+posmin)/2+posshift),
				labels=chr.labels, 
				ticks=F, rot=0,
				check.overlap=F
			)
		} else if (side=="bottom" || side=="top" || side=="right") {
			panel.axis(side=side, draw.labels=F, ticks=F);
		}
		else {
			axis.default(side=side,...);
		}
	 }

	#make sure the y-lim covers the range (plus a bit more to look nice)
	prepanel.chr<-function(x,y,...) { 
		A<-list();
		label.y <- unlist(lapply(ann.settings, function(z) {
			y0 <- z$label$y
			if (is.numeric(y0)) y0 else numeric()
		}), use.names=FALSE)
		maxy<-ceiling(max(y, label.y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+1;
		miny <- 0
		if (!is.null(mirror.plot)) {
			miny <- mirror.bottom
		}
		A$ylim=c(miny,maxy);
		A;
	}

	xyplot(logp~genpos, chr=chr, groups=grp,
		axis=axis.chr, ann.settings=ann.settings, 
		prepanel=prepanel.chr, scales=list(axs="i"),
		panel=function(x, y, ..., getgenpos) {
			if (!is.null(mirror.plot)) {
				total.width <- sum(posmax)
				visible.half <- total.width/(2*4000)
				half <- pmax((mirror.plot$X2-mirror.plot$X1)/2, visible.half)
				draw1 <- mirror.plot$XMID-half
				draw2 <- mirror.plot$XMID+half
				chr.index <- match(as.character(mirror.plot$CHR), chr.levels)
				seg.col <- mirror.col[(chr.index-1) %% length(mirror.col)+1]
				significant <- mirror.plot$P <= mirror.sig.level
				cis.significant <- significant & mirror.plot$CLASS == "cis"
				trans.significant <- significant & !cis.significant
				background <- !significant
				if (any(background)) panel.segments(draw1[background], -mirror.plot$Y[background],
					draw2[background], -mirror.plot$Y[background],
					col=adjustcolor(seg.col[background], alpha.f=0.82), lwd=1.25)
				if (any(trans.significant)) panel.segments(draw1[trans.significant], -mirror.plot$Y[trans.significant],
					draw2[trans.significant], -mirror.plot$Y[trans.significant],
					col=mirror.sig.col, lwd=2.1)
				if (any(cis.significant)) panel.segments(draw1[cis.significant], -mirror.plot$Y[cis.significant],
					draw2[cis.significant], -mirror.plot$Y[cis.significant],
					col=mirror.cis.col, lwd=2.1)
				# Mirror values are drawn at log10(P), so their threshold must also
				# be negative.  Keep both significance guides sky blue and unlabeled.
				panel.abline(h=log10(mirror.sig.level), lty=2, col=mirror.threshold.col)
				panel.abline(h=0, col="#4B5563", lwd=1.1)
				label.n <- min(as.integer(mirror.label.n), sum(significant))
				if (is.finite(label.n) && label.n > 0) {
					ii <- which(significant)
					ii <- head(ii[order(mirror.plot$P[ii])], label.n)
					label.col <- ifelse(mirror.plot$CLASS[ii] == "cis", mirror.cis.col, mirror.sig.col)
					panel.text(mirror.plot$XMID[ii], -mirror.plot$Y[ii]-0.65,
						labels=mirror.plot$GENE[ii], srt=90, adj=c(1,0.5),
						col=label.col, cex=0.58, font=2)
				}
			}
			if(!is.na(sig.level)) {
				#add significance line (if requested)
				panel.abline(h=-log10(sig.level), lty=2, col=threshold.col);
			}
			panel.superpose(x, y, ..., getgenpos=getgenpos);
			if (!is.null(mirror.plot)) {
				chr.labels <- chr.levels
				chr.labels[chr.labels == "23"] <- "X"
				chr.labels[chr.labels == "24"] <- "Y"
				panel.text((posmax+posmin)/2+posshift, 0, labels=chr.labels,
					adj=c(0.5, 1.25), col=chr.label.col, cex=0.76, font=2)
			}
			if(!is.null(panel.extra)) {
				panel.extra(x,y, getgenpos, ...)
			}
		},
		panel.groups = function(x,y,..., subscripts, group.number) {
			A<-list(...)
			#allow for different annotation settings
			gs <- ann.settings[[group.number]]
			A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
			A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
			A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
			A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
			A$x <- x
			A$y <- y
			do.call("panel.xyplot", A)
			#draw labels (if requested)
			if(gs$label$show) {
				gt<-gs$label
				names(gt)[which(names(gt)=="text")]<-"labels"
				gt$show<-NULL
				if(is.character(gt$x) | is.character(gt$y)) {
					peak = which.max(y)
					center = mean(range(x))
					if (is.character(gt$x)) {
						if(gt$x=="peak") {gt$x<-x[peak]}
						if(gt$x=="center") {gt$x<-center}
					}
					if (is.character(gt$y)) {
						if(gt$y=="peak") {gt$y<-y[peak]}
					}
				}
				if(is.list(gt$x)) {
					gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
				}
				do.call("panel.text", gt)
			}
		},
		xlab=xlab, ylab=ylab, 
		panel.extra=panel.extra, getgenpos=getGenPos, ...
	);
}

# Draw a Manhattan plot from a standardized small GWAS file. cis_bed supplies
# CHR, START, END, and GENE; cis_gene identifies the red target gene.
mh_plot <- function(smalled_gwas_file, col=c("gray", "darkgray"), cis_gene=NULL, cis_bed=NULL, mh_plot_bed=NULL,
		cis_color="red", other_top_color="green", other_top_max_per_chr=2,
		locus_size=1e6, genomewide_p=5e-8, main=NULL, mirror=NULL,
		mirror.sig.level=2.5e-6, ...) {
	if (!file.exists(smalled_gwas_file)) stop("GWAS file does not exist: ", smalled_gwas_file)
	if (requireNamespace("data.table", quietly=TRUE)) {
		dat <- data.table::fread(smalled_gwas_file, showProgress=FALSE, data.table=FALSE)
	} else {
		con <- if (grepl("\\.gz$", smalled_gwas_file, ignore.case=TRUE)) gzfile(smalled_gwas_file) else smalled_gwas_file
		dat <- read.table(con, header=TRUE, sep="\t", quote="", comment.char="", check.names=FALSE)
	}
	need <- c("CHR", "POS", "P")
	if (!all(need %in% names(dat))) stop("GWAS file must contain columns: ", paste(need, collapse=", "))
	dat$CHR <- sub("^chr", "", as.character(dat$CHR), ignore.case=TRUE)
	dat$CHR[dat$CHR == "X"] <- "23"; dat$CHR[dat$CHR == "Y"] <- "24"
	dat$CHR <- suppressWarnings(as.integer(dat$CHR)); dat$POS <- suppressWarnings(as.numeric(dat$POS)); dat$P <- suppressWarnings(as.numeric(dat$P))
	positive <- dat$P[is.finite(dat$P) & dat$P > 0]
	p_floor <- if (length(positive)) max(.Machine$double.xmin, min(positive) * 0.1) else .Machine$double.xmin
	dat$P[is.finite(dat$P) & dat$P <= 0] <- p_floor
	dat <- dat[is.finite(dat$CHR) & dat$CHR >= 1 & dat$CHR <= 24 & is.finite(dat$POS) & dat$POS > 0 & is.finite(dat$P) & dat$P > 0 & dat$P <= 1, , drop=FALSE]
	if (!nrow(dat)) stop("No valid CHR/POS/P rows in ", smalled_gwas_file)
	group <- rep("Background", nrow(dat))
	locus_size <- as.numeric(locus_size)
	if (!is.finite(locus_size) || locus_size <= 0) stop("locus_size must be positive")
	read_gene_bed <- function(x, arg_name) {
		if (is.null(x)) return(NULL)
		if (is.character(x) && length(x) == 1) {
			if (!file.exists(x)) stop(arg_name, " does not exist: ", x)
			b <- read.table(x, header=FALSE, comment.char="#", stringsAsFactors=FALSE)
		} else b <- as.data.frame(x, stringsAsFactors=FALSE)
		if (ncol(b) < 4) stop(arg_name, " must contain CHR, START, END, and GENE")
		b <- b[, 1:4, drop=FALSE]; names(b) <- c("CHR", "START", "END", "GENE")
		b$CHR <- sub("^chr", "", as.character(b$CHR), ignore.case=TRUE)
		b$CHR[b$CHR == "X"] <- "23"; b$CHR[b$CHR == "Y"] <- "24"
		b$CHR <- suppressWarnings(as.integer(b$CHR)); b$START <- suppressWarnings(as.numeric(b$START)); b$END <- suppressWarnings(as.numeric(b$END)); b$GENE <- as.character(b$GENE)
		b[is.finite(b$CHR) & is.finite(b$START) & is.finite(b$END) & nzchar(b$GENE), , drop=FALSE]
	}
	bed <- read_gene_bed(cis_bed, "cis_bed")
	primary_bed <- read_gene_bed(mh_plot_bed, "mh_plot_bed")
	cis <- NULL
	if (!is.null(cis_gene)) {
		if (is.null(bed)) stop("cis_bed is required when cis_gene is specified")
		cis <- bed[bed$GENE == as.character(cis_gene), , drop=FALSE]
		if (!nrow(cis)) stop("No cis_bed entry for cis_gene: ", cis_gene)
		cis_center <- (min(cis$START) + max(cis$END)) / 2
		cis_start <- cis_center - locus_size / 2; cis_end <- cis_center + locus_size / 2
		group[dat$CHR == cis$CHR[1] & dat$POS >= cis_start & dat$POS <= cis_end] <- as.character(cis_gene)
	}
	if (!is.null(mirror)) {
		mirror <- as.data.frame(mirror, stringsAsFactors=FALSE)
		if (!all(c("CHR", "START", "STOP", "P") %in% names(mirror)))
			stop("mirror must contain CHR, START, STOP, and P")
		if (length(mirror.sig.level) != 1 || !is.finite(mirror.sig.level) ||
			mirror.sig.level <= 0 || mirror.sig.level > 1) mirror.sig.level <- 2.5e-6
		mirror$CHR <- suppressWarnings(as.integer(mirror$CHR))
		mirror$START <- suppressWarnings(as.numeric(mirror$START))
		mirror$STOP <- suppressWarnings(as.numeric(mirror$STOP))
		mirror$P <- suppressWarnings(as.numeric(mirror$P))
		mirror$CLASS <- ifelse(is.finite(mirror$P) & mirror$P <= mirror.sig.level, "trans", "background")
		if (!is.null(cis) && nrow(cis)) {
			cis.overlap <- mirror$CHR == cis$CHR[1] & mirror$STOP >= cis_start & mirror$START <= cis_end
			mirror$CLASS[cis.overlap & mirror$P <= mirror.sig.level] <- "cis"
		}
	}
	max_per_chr <- as.integer(other_top_max_per_chr)
	if (!is.finite(max_per_chr) || max_per_chr < 0) stop("other_top_max_per_chr must be a non-negative integer")
	eligible <- which(group == "Background" & dat$P <= genomewide_p)
	other_groups <- other_labels <- character()
	other_label_y <- numeric()
	if (length(eligible) && max_per_chr > 0) {
		# Define significant loci from SNP coordinates first, independently of gene coverage.
		loci <- do.call(rbind, lapply(split(eligible, dat$CHR[eligible]), function(ii) {
			ii <- ii[order(dat$POS[ii])]
			loc <- cumsum(c(TRUE, diff(dat$POS[ii]) > locus_size))
			do.call(rbind, lapply(split(ii, loc), function(jj) {
				lead <- jj[which.min(dat$P[jj])]
				data.frame(CHR=dat$CHR[lead], POS=dat$POS[lead], P=dat$P[lead],
					LOCUS_START=min(dat$POS[jj]) - locus_size/2,
					LOCUS_END=max(dat$POS[jj]) + locus_size/2)
			}))
		}))
		loci <- loci[order(loci$CHR, loci$P, loci$POS), , drop=FALSE]
		keep <- unlist(lapply(split(seq_len(nrow(loci)), loci$CHR), head, max_per_chr), use.names=FALSE)
		loci <- loci[keep, , drop=FALSE]

		near_genes <- function(b, chr, pos, radius=1e6) {
			if (is.null(b) || !nrow(b)) return(NULL)
			x <- b[b$CHR == chr, , drop=FALSE]
			if (!nrow(x)) return(NULL)
			x$DIST <- pmax(x$START - pos, pos - x$END, 0)
			x <- x[x$DIST <= radius, , drop=FALSE]
			x[order(x$DIST, x$START, x$GENE), , drop=FALSE]
		}
		gene_label <- function(chr, pos) {
			# Exact overlap, primary first and cis BED only as fallback.
			x <- near_genes(primary_bed, chr, pos, 0)
			if (is.null(x) || !nrow(x)) x <- near_genes(bed, chr, pos, 0)
			if (!is.null(x) && nrow(x)) return(paste(head(unique(x$GENE), 3), collapse=","))
			# No overlap: genes within +/-1 Mb, again primary then fallback.
			x <- near_genes(primary_bed, chr, pos, 1e6)
			if (is.null(x) || !nrow(x)) x <- near_genes(bed, chr, pos, 1e6)
			if (!is.null(x) && nrow(x)) return(paste(head(unique(x$GENE), 3), collapse=","))
			# No nearby gene: report nearest upstream and downstream genes with
			# signed boundary distance (upstream negative, downstream positive).
			x <- if (!is.null(primary_bed) && any(primary_bed$CHR == chr)) primary_bed[primary_bed$CHR == chr, , drop=FALSE] else bed[bed$CHR == chr, , drop=FALSE]
			if (is.null(x) || !nrow(x)) return(paste0("chr", chr, ":", format(pos, scientific=FALSE)))
			up <- x[x$END < pos, , drop=FALSE]; down <- x[x$START > pos, , drop=FALSE]
			labs <- character()
			if (nrow(up)) { u <- up[which.max(up$END), ]; labs <- c(labs, sprintf("%s (%+.1fMb)", u$GENE, -(pos-u$END)/1e6)) }
			if (nrow(down)) { d <- down[which.min(down$START), ]; labs <- c(labs, sprintf("%s (%+.1fMb)", d$GENE, (d$START-pos)/1e6)) }
			paste(labs, collapse=",")
		}
		for (i in seq_len(nrow(loci))) {
			key <- paste0("Other locus ", i)
			inside <- group == "Background" & dat$CHR == loci$CHR[i] & dat$POS >= loci$LOCUS_START[i] & dat$POS <= loci$LOCUS_END[i]
			group[inside] <- key
			other_groups <- c(other_groups, key)
			other_labels <- c(other_labels, gene_label(loci$CHR[i], loci$POS[i]))
			lead_y <- -log10(loci$P[i])
			stagger_y <- -log10(genomewide_p) + 2.5 + ((i - 1) %% 4) * 4
			other_label_y <- c(other_label_y, max(lead_y + 1.5, stagger_y))
		}
	}
	group <- droplevels(factor(group, levels=unique(c("Background", as.character(cis_gene), other_groups))))
	ann <- list(group)
	if (!is.null(cis_gene) && as.character(cis_gene) %in% levels(group)) ann[[as.character(cis_gene)]] <- list(
		col=cis_color, fill=cis_color, cex=1.05,
		label=list(show=TRUE, text=as.character(cis_gene), pos=1, offset=0.8, cex=0.85, fontface=2))
	for (i in seq_along(other_groups)) if (other_groups[i] %in% levels(group)) ann[[other_groups[i]]] <- list(
		col=other_top_color, fill=other_top_color, cex=1.05,
		label=list(show=TRUE, text=other_labels[i], y=other_label_y[i], pos=3, offset=0.35, cex=0.72))
	if (is.null(main)) main <- sub("\\.(gz|bgz|tsv|txt)$", "", basename(smalled_gwas_file), ignore.case=TRUE)
	.mh_plot_vectors(dat$CHR, dat$POS, dat$P, sig.level=genomewide_p, annotate=ann,
		mirror=mirror, mirror.sig.level=mirror.sig.level, mirror.col=col,
		mirror.sig.col=other_top_color, mirror.cis.col=cis_color,
		col=col, main=main, ...)
}
