# FINAL C3 robustness synthesis.
# This is deliberately orthogonal to the native C3 calling rule. It summarizes
# prior robustness, H3-vs-H4 separation, credible-set resolution, SuSiE support,
# LD diagnostics and long-range-LD regions without choosing whichever statistic
# gives the strongest result.

le8_c3_final_robustness <- function(res, variants=tibble::tibble(), susie=list(), layer=c('protein','metabolite'), rawdir, outdir=NULL) {
  layer <- match.arg(layer)
  `%||%` <- function(x,y) if(is.null(x)) y else x
  finite_max <- function(x) if(any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
  finite_min <- function(x) if(any(is.finite(x))) min(x[is.finite(x)]) else NA_real_
  z <- dplyr::as_tibble(res)
  if(!nrow(z)) {
    out <- tibble::tibble(status='unavailable', reason='No C3 loci')
    write_raw_csv(out,'c3.robustness_consensus.csv',rawdir); return(out)
  }
  neednum <- c('PP.H3','PP.H4','PP.H4_p12_conservative','PP.H4_p12_default','PP.H4_p12_liberal','PP.H4_robust_min','credible_set_n','lead_shared_pp')
  for(nm in neednum) if(!nm %in% names(z)) z[[nm]] <- NA_real_
  if(!'status' %in% names(z)) z$status <- 'unknown'
  if(!'locus_class' %in% names(z)) z$locus_class <- NA_character_
  if(!'chr' %in% names(z)) z$chr <- sub('^chr([^:]+):.*$','\\1',as.character(z$locus))
  if(!'start' %in% names(z)) z$start <- suppressWarnings(as.numeric(sub('^.*:([0-9.]+)-.*$','\\1',as.character(z$locus))))
  if(!'end' %in% names(z)) z$end <- suppressWarnings(as.numeric(sub('^.*-([0-9.]+)$','\\1',as.character(z$locus))))

  # SuSiE summaries are signal-pair specific; collapse only after requiring
  # sensitivity across all p12 values for the same signal pair.
  ss <- list(); dd <- list(); kk <- 0L; kd <- 0L
  if(length(susie)) for(s in susie) {
    sm <- s$summary %||% tibble::tibble()
    pp <- s$pip %||% s$variants %||% tibble::tibble()
    dg <- s$diagnostics %||% tibble::tibble()
    if(is.data.frame(sm) && nrow(sm) && all(c('feature','locus') %in% names(sm))) {
      if('PP.H4.abf' %in% names(sm)) {
        grp <- intersect(c('feature','locus','idx1','idx2'),names(sm))
        a <- sm |> dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
          dplyr::summarise(susie_H4_min=finite_min(.data$PP.H4.abf),susie_H4_max=finite_max(.data$PP.H4.abf),.groups='drop') |>
          dplyr::group_by(.data$feature,.data$locus) |>
          dplyr::summarise(susie_H4_robust=finite_max(.data$susie_H4_min),susie_signal_pairs=dplyr::n(),.groups='drop')
        kk<-kk+1L; ss[[kk]]<-a
      }
    }
    if(is.data.frame(pp) && nrow(pp) && all(c('feature','locus','QTL_PIP','outcome_PIP') %in% names(pp))) {
      a <- pp |> dplyr::mutate(CLPP=.data$QTL_PIP*.data$outcome_PIP) |>
        dplyr::group_by(.data$feature,.data$locus) |>
        dplyr::summarise(max_CLPP=finite_max(.data$CLPP),max_QTL_PIP=finite_max(.data$QTL_PIP),max_outcome_PIP=finite_max(.data$outcome_PIP),.groups='drop')
      kk<-kk+1L; ss[[kk]]<-a
    }
    if(is.data.frame(dg) && nrow(dg) && all(c('locus','ld_summary_mismatch') %in% names(dg))) {
      if(!'feature'%in%names(dg) && is.data.frame(sm) && nrow(sm) && 'feature'%in%names(sm)) dg$feature <- sm$feature[1]
      if('feature'%in%names(dg)) {
        kd<-kd+1L; dd[[kd]] <- dg |> dplyr::group_by(.data$feature,.data$locus) |>
          dplyr::summarise(ld_summary_mismatch=finite_max(.data$ld_summary_mismatch),.groups='drop')
      }
    }
  }
  if(length(ss)) {
    ssum <- dplyr::bind_rows(ss) |>
      dplyr::group_by(.data$feature,.data$locus) |>
      dplyr::summarise(dplyr::across(dplyr::everything(),finite_max),.groups='drop')
    z <- dplyr::left_join(z,ssum,by=c('feature','locus'))
  }
  if(length(dd)) {
    dsum <- dplyr::bind_rows(dd) |> dplyr::group_by(.data$feature,.data$locus) |>
      dplyr::summarise(ld_summary_mismatch=finite_max(.data$ld_summary_mismatch),.groups='drop')
    z <- dplyr::left_join(z,dsum,by=c('feature','locus'))
  }
  for(nm in c('susie_H4_robust','susie_signal_pairs','max_CLPP','max_QTL_PIP','max_outcome_PIP','ld_summary_mismatch')) if(!nm%in%names(z)) z[[nm]]<-NA_real_

  post_thr <- le8_evidence_policy()$posterior
  z <- z |> dplyr::mutate(
    H34_sum=.data$PP.H3+.data$PP.H4,
    H4_H3_ratio=.data$PP.H4/pmax(.data$PP.H3,1e-12),
    coloc_H34_separation=is.finite(.data$H34_sum)&is.finite(.data$H4_H3_ratio)&.data$H34_sum>.99&.data$H4_H3_ratio>5,
    prior_robust=is.finite(.data$PP.H4_robust_min)&.data$PP.H4_robust_min>=post_thr,
    default_H4=is.finite(.data$PP.H4)&.data$PP.H4>=post_thr,
    susie_robust=is.finite(.data$susie_H4_robust)&.data$susie_H4_robust>=post_thr,
    clpp_high=is.finite(.data$max_CLPP)&.data$max_CLPP>=.01,
    one_snp_concentration=is.finite(.data$credible_set_n)&.data$credible_set_n==1&dplyr::coalesce(.data$lead_shared_pp,0)>.999,
    long_range_LD=as.character(.data$chr)%in%c('6','chr6') & dplyr::coalesce(.data$end,0)>=25000000 & dplyr::coalesce(.data$start,1e12)<=34000000,
    ld_mismatch_warning=is.finite(.data$ld_summary_mismatch)&.data$ld_summary_mismatch>.05,
    robustness_score=2L*as.integer(.data$prior_robust)+2L*as.integer(.data$coloc_H34_separation)+
      as.integer(dplyr::coalesce(.data$susie_robust,FALSE))+as.integer(dplyr::coalesce(.data$clpp_high,FALSE))-
      as.integer(.data$one_snp_concentration)-as.integer(.data$long_range_LD)-as.integer(.data$ld_mismatch_warning),
    robustness_tier=dplyr::case_when(
      .data$status!='ok' ~ 'U: unavailable/failed',
      .data$prior_robust & .data$coloc_H34_separation & !.data$long_range_LD & !.data$ld_mismatch_warning & (dplyr::coalesce(.data$susie_robust,FALSE)|!is.finite(.data$susie_H4_robust)) ~ 'A: robust shared-signal evidence',
      .data$prior_robust & .data$coloc_H34_separation ~ 'B: robust ABF with sensitivity warning',
      .data$default_H4 ~ 'C: default-prior H4 only',
      .data$PP.H4>.5 ~ 'D: suggestive colocalization', TRUE ~ 'E: little shared-signal support'),
    interpretation='ABF, H3/H4 separation, SuSiE and LD diagnostics are complementary; tier is not a causal posterior'
  ) |> dplyr::arrange(factor(.data$robustness_tier,levels=c('A: robust shared-signal evidence','B: robust ABF with sensitivity warning','C: default-prior H4 only','D: suggestive colocalization','E: little shared-signal support','U: unavailable/failed')),dplyr::desc(.data$robustness_score),dplyr::desc(.data$PP.H4))

  write_raw_csv(z,'c3.robustness_consensus.csv',rawdir)
  cnt <- z |> dplyr::count(.data$robustness_tier,name='loci')
  write_raw_csv(cnt,'c3.robustness_tiers.csv',rawdir)
  if(!is.null(outdir)&&requireNamespace('ggplot2',quietly=TRUE)) {
    top<-z |> dplyr::filter(.data$status=='ok') |> dplyr::slice_head(n=40) |>
      dplyr::mutate(label=factor(paste(.data$feature,.data$locus,sep=' | '),levels=rev(paste(.data$feature,.data$locus,sep=' | '))))
    if(nrow(top)) {
      p<-ggplot2::ggplot(top,ggplot2::aes(.data$PP.H4_robust_min,.data$label,color=.data$robustness_tier))+
        ggplot2::geom_vline(xintercept=post_thr,linetype=2,color='grey60')+ggplot2::geom_point(size=2.3)+
        ggplot2::scale_x_continuous(limits=c(0,1))+ggplot2::labs(title='C3 colocalization robustness',x='Minimum H4 across p12 sensitivity',y=NULL,color=NULL,
          caption='Long-range LD, LD mismatch and one-SNP concentration are explicit warnings; SuSiE is used when valid LD is available.')+
        theme_5c(8.5)+ggplot2::theme(legend.position='bottom')
      save_plot(p,'c3.Fig6.robustness_consensus.png',16,11,outdir=outdir)
    }
  }
  z
}
