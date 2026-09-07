#!/usr/bin/env Rscript
# Multi-GWAS QC comparison of standardized summary statistics.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
if (!length(args) || any(args %in% c('-h', '--help'))) {
  cat('Usage: gwas_compare.sh compare --gwas-files A.gz,B.gz[,C.gz] [options]\n',
      '--output-dir DIR (default ./gwas_compare) --labels CSV --grch 37|38\n',
      '--mplot TRUE --compare-beta TRUE --compare-EAF TRUE\n',
      '--p-threshold 5e-8 --significant first|either|both (default first)\n',
      'Input: SNP CHR POS EA NEA P; BETA/EAF required when compared. Same build required.\n',
      'First GWAS is compared with each follower. Palindromic SNPs and duplicate\n',
      'allele/position keys are excluded from scatter plots and counted in QC.\n',
      '\nExamples:\n',
      '  gwas_compare.sh compare --gwas-files A.gz,B.gz,C.gz --labels A,B,C --output-dir qc\n',
      '  gwas_compare.sh compare --gwas-files A.gz,B.gz --significant either --p-threshold 5e-8\n',
      '  gwas_compare.sh compare --gwas-files A.gz,B.gz --mplot FALSE --compare-beta TRUE --compare-EAF TRUE\n',
      'Use gwas_compare.sh -h for full paths and output descriptions.\n')
  quit(status=0)
}
opt <- list('output-dir'='gwas_compare', mplot='TRUE', 'compare-beta'='TRUE',
            'compare-eaf'='TRUE', 'p-threshold'='5e-8', significant='first')
allowed <- c(names(opt), 'gwas-files', 'labels', 'grch')
if (length(args) %% 2L) stop('Options require values')
for (i in seq(1, length(args), 2)) {
  k <- tolower(sub('^--', '', args[i]))
  if (!k %in% allowed) stop('Unknown option: ', args[i])
  opt[[k]] <- args[i+1]
}
flag <- function(x) { if (!toupper(x) %in% c('TRUE','FALSE')) stop('Expected TRUE/FALSE'); toupper(x)=='TRUE' }
mp <- flag(opt$mplot); cb <- flag(opt[['compare-beta']]); ce <- flag(opt[['compare-eaf']])
if (is.null(opt[['gwas-files']])) stop('--gwas-files required')
files <- trimws(strsplit(opt[['gwas-files']], ',', fixed=TRUE)[[1]])
if (length(files)<2 || any(!file.exists(files))) stop('Provide at least two existing GWAS files')
labels <- if (is.null(opt$labels)) basename(files) else strsplit(opt$labels, ',', fixed=TRUE)[[1]]
if (length(labels)!=length(files) || anyDuplicated(labels)) stop('Labels must be unique and match files')
threshold <- as.numeric(opt[['p-threshold']])
if (!is.finite(threshold) || threshold<=0 || threshold>1) stop('Invalid P threshold')
if (!opt$significant %in% c('first','either','both')) stop('Invalid significant selector')
if (!is.null(opt$grch) && !opt$grch %in% c('37','38')) stop('Invalid GRCh')
builds <- unique(unlist(lapply(files,function(f) if(file.exists(paste0(f,'.grch'))) trimws(readLines(paste0(f,'.grch'))) else NULL)))
if(length(unique(c(builds,opt$grch)))>1L) stop('Mixed GRCh builds: liftover before comparison')
out <- opt[['output-dir']]; dir.create(out, recursive=TRUE, showWarnings=FALSE)
comp <- function(x) chartr('ACGT', 'TGCA', x)
read_gwas <- function(f) {
  # fread reads one genome-wide file at a time. gzip shell command is quoted.
  d <- if (grepl('\\.(gz|bgz)$', f)) fread(cmd=paste('gzip -cd --', shQuote(f)), showProgress=FALSE) else fread(f, showProgress=FALSE)
  setnames(d, toupper(sub('^#','',names(d))))
  required <- c('SNP','CHR','POS','EA','NEA','P', if(cb) 'BETA', if(ce) 'EAF')
  if (length(setdiff(required,names(d)))) stop(f, ': missing ', paste(setdiff(required,names(d)),collapse=','))
  d <- d[, ..required]
  for (k in intersect(c('POS','P','EAF','BETA'),names(d))) set(d,j=k,value=suppressWarnings(as.numeric(d[[k]])))
  d[, CHR := toupper(sub('^CHR','',toupper(as.character(CHR))))]
  d[CHR=='X',CHR:='23']; d[CHR=='Y',CHR:='24']; d[CHR %in% c('MT','M'),CHR:='25']
  d[, CHR := suppressWarnings(as.integer(CHR))]
  d[, `:=`(EA=toupper(EA),NEA=toupper(NEA))]
  d[, valid := is.finite(P) & P>=0 & P<=1 & is.finite(POS) & POS>0 & POS==floor(POS) & CHR %in% 1:25 & !is.na(EA) & !is.na(NEA) & EA!=NEA]
  invalid <- d[valid==FALSE,.N]; d <- d[valid==TRUE]; d[,valid:=NULL]
  if(!nrow(d)) stop(f, ': no valid variants')
  # Canonicalize SNP complements; indels require literal matching (no strand inference).
  d[, pair := paste(pmin(EA,NEA),pmax(EA,NEA),sep='/')]
  d[, pal := pair %in% c('A/T','C/G')]
  d[, canonical := pair]
  d[nchar(EA)==1 & nchar(NEA)==1 & grepl('^[ACGT]$',EA) & grepl('^[ACGT]$',NEA),
    canonical := pmin(pair,paste(pmin(comp(EA),comp(NEA)),pmax(comp(EA),comp(NEA)),sep='/'))]
  d[, key := paste(CHR,POS,canonical,sep=':')]
  d[, duplicate := duplicated(key) | duplicated(key,fromLast=TRUE)]
  attr(d,'invalid') <- invalid
  d
}
# Keep significant identities first; follower-only significant sites are included on demand.
first <- read_gwas(files[1])
base <- first[!duplicate & !pal]
if(opt$significant %in% c('first','both')) base <- base[P<=threshold]
audit <- list(); summary <- list(); tracks <- list(); maxima <- numeric(25)
track <- function(d) {
  # Deterministic 20kb / 0.1-logP bins keep all significant points and visible background.
  x <- d[,.(CHR,POS,logP=-log10(pmax(P,1e-300)),P)]
  x[, bin := paste(CHR,POS %/% 20000, floor(logP*10),sep=':')]
  x <- x[P<=threshold | !duplicated(bin)]
  x[,c('bin','P'):=NULL]; x
}
for (i in seq_along(files)) {
  d <- if (i==1L) first else read_gwas(files[i])
  audit[[i]] <- data.table(file=files[i],label=labels[i],invalid_rows=attr(d,'invalid'),valid_rows=nrow(d),duplicate_rows=sum(d$duplicate),palindromic_rows=sum(d$pal))
  if (mp) tracks[[i]] <- track(d)
  z <- d[,.(end=max(POS)),by=CHR]; maxima[z$CHR] <- pmax(maxima[z$CHR],z$end)
  if (i==1L || !(cb || ce)) next
  if(opt$significant %in% c('first','both')) d <- d[key %in% base$key]
  x <- merge(base,d[!duplicate & !pal],by='key',suffixes=c('.first','.other'))
  x <- x[switch(opt$significant,first=P.first<=threshold,either=P.first<=threshold|P.other<=threshold,both=P.first<=threshold&P.other<=threshold)]
  x[, flip := EA.first==NEA.other & NEA.first==EA.other |
      (nchar(EA.first)==1 & nchar(NEA.first)==1 & EA.first==comp(NEA.other) & NEA.first==comp(EA.other))]
  if (cb) x[flip==TRUE,BETA.other := -BETA.other]
  if (ce) x[flip==TRUE,EAF.other := 1-EAF.other]
  tag <- sprintf('01_vs_%02d',i)
  fwrite(x,file.path(out,paste0(tag,'.harmonized.tsv.gz')),sep='\t')
  for (v in c(if(cb) 'BETA',if(ce) 'EAF')) {
    a <- x[[paste0(v,'.first')]]; b <- x[[paste0(v,'.other')]]
    ok <- is.finite(a)&is.finite(b)
    if (v=='EAF') ok <- ok&a>=0&a<=1&b>=0&b<=1
    a <- a[ok]; b <- b[ok]
    r <- if(length(a)>1 && sd(a)>0 && sd(b)>0) cor(a,b) else NA_real_
    summary[[length(summary)+1]] <- data.table(first=labels[1],other=labels[i],metric=v,matched_significant=nrow(x),n=length(a),flipped=sum(x$flip),r=r,mean_difference=if(length(a)) mean(b-a) else NA_real_)
    png(file.path(out,paste0(tag,'.',v,'.png')),width=1400,height=1400,res=180)
    if (length(a)) {
      lim <- range(c(a,b)); if(diff(lim)==0) lim <- lim+c(-.01,.01)
      plot(a,b,pch=16,cex=.45,col=adjustcolor('#286b9e',alpha.f=.35),xlim=lim,ylim=lim,
           xlab=paste(labels[1],v),ylab=paste(labels[i],v),main=sprintf('n=%d; r=%.3f',length(a),r)); abline(0,1,col='firebrick')
    } else { plot.new(); title(main=paste(v,': no eligible matched significant variants')) }
    dev.off()
  }
}
if (mp) {
  ch <- which(maxima>0); offsets <- c(0,head(cumsum(maxima[ch]+1e6),-1)); names(offsets)<-ch
  ymax <- max(1,-log10(threshold),unlist(lapply(tracks,function(x) x$logP)))
  png(file.path(out,'manhattan.compare.png'),width=2600,height=max(1000,600*length(files)),res=180)
  par(mfrow=c(length(files),1),mar=c(3,4,2,1))
  for(i in seq_along(files)) {
    x <- tracks[[i]]
    plot(x$POS+offsets[as.character(x$CHR)],x$logP,pch=16,cex=.25,
         col=c('#315780','#d49740')[1+x$CHR%%2],xlim=c(0,sum(maxima[ch]+1e6)),ylim=c(0,ymax),xaxt='n',xlab='',ylab='-log10(P)',main=labels[i])
    axis(1,at=offsets+maxima[ch]/2,labels=ifelse(ch==23,'X',ifelse(ch==24,'Y',ifelse(ch==25,'MT',ch)))); abline(h=-log10(threshold),lty=2,col='firebrick')
  }
  dev.off()
}
fwrite(rbindlist(audit),file.path(out,'input_qc.tsv'),sep='\t')
if(length(summary)) fwrite(rbindlist(summary),file.path(out,'comparison_qc.tsv'),sep='\t')
writeLines(c(paste('GRCh:',if(is.null(opt$grch)) 'unspecified; caller must ensure same build' else opt$grch),
             paste('Significance:',opt$significant,'P <=',threshold),capture.output(sessionInfo())),file.path(out,'compare.log'))
cat('Comparison written to ',normalizePath(out),'\n',sep='')
