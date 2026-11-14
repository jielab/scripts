# FINAL C2 robustness synthesis.
# Does not replace native MR/MR-link-2/DANDELION analyses. It prevents a single
# small P value from becoming the whole causal story by summarizing independent
# diagnostics under a prespecified evidence rule.

le8_c2_final_robustness <- function(mr, reverse_mr=tibble::tibble(), mrlink2=list(), assoc=tibble::tibble(),
                                    layer=c('protein','metabolite'), rawdir, outdir=NULL) {
  layer <- match.arg(layer)
  `%||%` <- function(x,y) if(is.null(x)) y else x
  safe <- function(d,n,default=NA_real_) if(n %in% names(d)) d[[n]] else rep(default,nrow(d))
  chr <- function(d,n,default=NA_character_) if(n %in% names(d)) as.character(d[[n]]) else rep(default,nrow(d))
  primary <- if(layer=='protein') 'cis' else 'local'
  if(!nrow(mr)) {
    out<-tibble::tibble(status='unavailable',reason=paste('No',primary,'MR rows'))
    write_raw_csv(out,'c2.robustness_consensus.csv',rawdir);return(out)
  }
  m0 <- dplyr::as_tibble(mr) |> dplyr::filter(.data$analysis==primary)
  z <- m0 |>
    dplyr::mutate(
      feature=.data$exposure,
      MR_beta=.data$b, MR_se=.data$se, MR_p=.data$pval,
      MR_FDR=if('FDR_all'%in%names(m0)) .data$FDR_all else p.adjust(.data$pval,'BH'),
      n_IV=safe(m0,'n_IV',0),
      median_F=safe(m0,'median_F'),
      Q_p=safe(m0,'Q_p'),
      egger_intercept_p=safe(m0,'egger_intercept_p'),
      weighted_median_beta=safe(m0,'tsmr_weighted_median_b'),
      steiger_support_fraction=safe(m0,'steiger_support_fraction'),
      forward_support=is.finite(.data$MR_FDR)&.data$MR_FDR<le8_evidence_policy()$mr_fdr,
      strong_F=is.finite(.data$median_F)&.data$median_F>=10,
      very_strong_F=is.finite(.data$median_F)&.data$median_F>=30,
      heterogeneity_warning=ifelse(is.finite(.data$Q_p)&.data$n_IV>1,.data$Q_p<.05,NA),
      egger_warning=ifelse(is.finite(.data$egger_intercept_p)&.data$n_IV>=3,.data$egger_intercept_p<.05,NA),
      weighted_median_concordant=ifelse(.data$n_IV>=3 & is.finite(.data$weighted_median_beta),
                                       sign(.data$weighted_median_beta)==sign(.data$MR_beta),NA),
      steiger_support=ifelse(is.finite(.data$steiger_support_fraction),.data$steiger_support_fraction>=.5,NA),
      steiger_strong=ifelse(is.finite(.data$steiger_support_fraction),.data$steiger_support_fraction>=.8,NA)
    )
  if(!nrow(z)) {
    out<-tibble::tibble(status='unavailable',reason=paste('No',primary,'MR rows'))
    write_raw_csv(out,'c2.robustness_consensus.csv',rawdir);return(out)
  }

  # Reverse direction: disease liability -> omic. This is a warning about direction,
  # not proof that the forward effect is false.
  if(nrow(reverse_mr)) {
    rv<-dplyr::as_tibble(reverse_mr)
    id <- if('feature'%in%names(rv)) 'feature' else if('exposure'%in%names(rv)) 'exposure' else NA_character_
    if(!is.na(id)) {
      rv<-rv |> dplyr::transmute(feature=.data[[id]],reverse_beta=safe(rv,'b'),reverse_p=safe(rv,'pval'),
                                reverse_FDR=if('FDR_reverse'%in%names(rv)) .data$FDR_reverse else p.adjust(safe(rv,'pval'),'BH')) |>
        dplyr::group_by(.data$feature) |> dplyr::slice_min(.data$reverse_p,n=1,with_ties=FALSE) |> dplyr::ungroup()
      z<-dplyr::left_join(z,rv,by='feature')
    }
  }
  if(!'reverse_FDR'%in%names(z)) z$reverse_FDR<-NA_real_
  if(!'reverse_beta'%in%names(z)) z$reverse_beta<-NA_real_
  z$reverse_support <- is.finite(z$reverse_FDR)&z$reverse_FDR<le8_evidence_policy()$mr_fdr

  # Observational direction is triangulation only.
  if(nrow(assoc)&&all(c('term','beta')%in%names(assoc))) {
    ob<-dplyr::as_tibble(assoc) |> dplyr::transmute(feature=.data$term,observed_beta=.data$beta,
          observed_p=if('p.value'%in%names(assoc)) .data$p.value else NA_real_) |>
      dplyr::group_by(.data$feature) |> dplyr::slice_min(dplyr::coalesce(.data$observed_p,Inf),n=1,with_ties=FALSE) |> dplyr::ungroup()
    z<-dplyr::left_join(z,ob,by='feature')
  }
  if(!'observed_beta'%in%names(z)) z$observed_beta<-NA_real_
  z$observed_concordant <- ifelse(is.finite(z$observed_beta),sign(z$observed_beta)==sign(z$MR_beta),NA)

  # MR-link-2 is an orthogonal regional LD-aware sensitivity. Column names vary
  # across upstream builds, so map conservatively and never manufacture support.
  ml <- mrlink2$results %||% tibble::tibble()
  if(nrow(ml)) {
    pick <- function(pattern) {x<-grep(pattern,names(ml),ignore.case=TRUE,value=TRUE);if(length(x))x[1]else NA_character_}
    fc<-pick('^(trait|exposure)$'); bc<-pick('^(alpha|effect|beta)$'); pc<-pick('p.*alpha|alpha.*p|^pval$|^p_value$')
    if(!is.na(fc)) {
      mm<-tibble::tibble(feature=as.character(ml[[fc]]),mrlink2_beta=if(is.na(bc))NA_real_ else suppressWarnings(as.numeric(ml[[bc]])),
                         mrlink2_p=if(is.na(pc))NA_real_ else suppressWarnings(as.numeric(ml[[pc]]))) |>
        dplyr::group_by(.data$feature) |> dplyr::slice_min(dplyr::coalesce(.data$mrlink2_p,Inf),n=1,with_ties=FALSE) |> dplyr::ungroup()
      z<-dplyr::left_join(z,mm,by='feature')
    }
  }
  if(!'mrlink2_beta'%in%names(z)) z$mrlink2_beta<-NA_real_
  if(!'mrlink2_p'%in%names(z)) z$mrlink2_p<-NA_real_
  z$mrlink2_concordant<-ifelse(is.finite(z$mrlink2_beta),sign(z$mrlink2_beta)==sign(z$MR_beta),NA)

  z <- z |>
    dplyr::mutate(
      diagnostics_available = rowSums(cbind(is.finite(.data$Q_p),is.finite(.data$egger_intercept_p),
        !is.na(.data$weighted_median_concordant),!is.na(.data$steiger_support),is.finite(.data$reverse_FDR),is.finite(.data$mrlink2_p))),
      warning_count = rowSums(cbind(dplyr::coalesce(.data$heterogeneity_warning,FALSE),dplyr::coalesce(.data$egger_warning,FALSE),
        dplyr::coalesce(!.data$weighted_median_concordant,FALSE),dplyr::coalesce(!.data$steiger_support,FALSE),
        dplyr::coalesce(.data$reverse_support,FALSE),dplyr::coalesce(!.data$mrlink2_concordant,FALSE))),
      robustness_score = 2*as.integer(.data$forward_support)+as.integer(.data$strong_F)+as.integer(.data$very_strong_F)+
        as.integer(dplyr::coalesce(!.data$heterogeneity_warning,FALSE))+as.integer(dplyr::coalesce(!.data$egger_warning,FALSE))+
        as.integer(dplyr::coalesce(.data$weighted_median_concordant,FALSE))+as.integer(dplyr::coalesce(.data$steiger_support,FALSE))+
        as.integer(dplyr::coalesce(.data$mrlink2_concordant,FALSE))-2*as.integer(dplyr::coalesce(.data$reverse_support,FALSE)),
      evidence_tier=dplyr::case_when(
        .data$forward_support & .data$strong_F & .data$warning_count==0 & .data$diagnostics_available>=2 ~ 'A: forward MR + compatible diagnostics',
        .data$forward_support & .data$strong_F ~ 'B: forward MR with sensitivity warning/incomplete diagnostics',
        .data$forward_support ~ 'C: MR support but weak/incomplete instrument evidence',
        is.finite(.data$MR_p)&.data$MR_p<.05 ~ 'D: nominal MR only', TRUE ~ 'E: no forward MR support'),
      interpretation='Tier is a robustness synthesis, not a posterior probability of causality'
    ) |>
    dplyr::arrange(factor(.data$evidence_tier,levels=c('A: forward MR + compatible diagnostics','B: forward MR with sensitivity warning/incomplete diagnostics','C: MR support but weak/incomplete instrument evidence','D: nominal MR only','E: no forward MR support')),dplyr::desc(.data$robustness_score),.data$MR_p)

  write_raw_csv(z,'c2.robustness_consensus.csv',rawdir)
  counts<-z |> dplyr::count(.data$evidence_tier,name='features')
  write_raw_csv(counts,'c2.robustness_tiers.csv',rawdir)
  if(!is.null(outdir)&&requireNamespace('ggplot2',quietly=TRUE)) {
    top<-z |> dplyr::slice_head(n=40) |> dplyr::mutate(feature=factor(.data$feature,levels=rev(.data$feature)))
    p<-ggplot2::ggplot(top,ggplot2::aes(.data$robustness_score,.data$feature,color=.data$evidence_tier))+
      ggplot2::geom_point(size=2.4)+ggplot2::geom_vline(xintercept=0,linetype=2,color='grey70')+
      ggplot2::labs(title='C2 causal robustness synthesis',x='Prespecified robustness score',y=NULL,color=NULL,
                    caption='Forward cis/local MR remains primary; sensitivity methods are not combined by minimum P.')+
      theme_5c(9)+ggplot2::theme(legend.position='bottom')
    save_plot(p,'c2.Fig15.robustness_consensus.png',14,10,outdir=outdir)
  }
  z
}
