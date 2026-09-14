# Coverage retains the genome-wide definition at every zoom level: merge all
# reference calls per individual before averaging over tested individuals.
gu_density_region <- function(d, r, Q, dataset, build, lineage="Neanderthal") {
  ch <- as.character(r$chr[[1]])
  edges <- unique(round(seq(r$start[[1]],r$end[[1]],length.out=81)))
  bins <- data.frame(chr=ch,start=head(edges,-1),end=tail(edges,-1))
  z <- matrix(0,nrow(d$samples),nrow(bins))
  tested <- which(as.character(d$bins$chr)==ch)
  eligible <- if(length(tested)) rowSums(is.finite(d$z[,tested,drop=FALSE]))>0 else rep(FALSE,nrow(z))
  calls <- Q("SELECT sample_id,start,end FROM segments INDEXED BY idx_segments_region WHERE dataset_id=? AND genome_build=? AND chr=? AND start<? AND end>? AND method='ibdmix' AND source_class=?",list(dataset,build,ch,r$end[[1]],r$start[[1]],lineage))
  if(nrow(calls)) for(s in split(calls,calls$sample_id)) {
    i <- match(s$sample_id[[1]],d$samples$sample_id); if(is.na(i)) next
    intervals <- reduce_intervals(s)
    for(j in seq_len(nrow(intervals))) z[i,] <- z[i,]+pmax(0,pmin(bins$end,intervals$end[j])-pmax(bins$start,intervals$start[j]))
  }
  ploidy <- if(ch=="X" && isTRUE(d$manifest$x_male_only))1 else 2
  z <- sweep(z,2,ploidy*(bins$end-bins$start),"/")*100
  z[!eligible,] <- NA_real_
  list(samples=d$samples,bins=bins,z=z)
}

gu_density_geometry <- function(bins) {
  chromosomes <- intersect(c(as.character(1:22),"X"),unique(as.character(bins$chr)))
  lengths <- setNames(vapply(chromosomes,function(ch)max(bins$end[bins$chr==ch]),numeric(1)),chromosomes)
  list(chromosomes=chromosomes,lengths=lengths,offsets=setNames(c(0,head(cumsum(lengths),-1)),chromosomes))
}
