# Extended S2a: non-linear LE8 and out-of-fold omics-score associations with incident disease.
suppressPackageStartupMessages({
  fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
  source(file.path(fdir, "comm.f.R"))
})
LE8_JOB <- "s2_nonlin"
outdir <- file.path(out.base, LE8_JOB); setwd2(outdir); rawdir <- le8_job_dir(outdir, LE8_JOB); dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
nonlin_cache <- file.path(rawdir,"s2_nonlin.res.rds")
penalty_cache <- file.path(rawdir,"s2.penalty.res.rds")
if(cache_valid(nonlin_cache) && cache_valid(penalty_cache)){
  cache_message("S2/non-linear",nonlin_cache)
  cache_message("S2/penalty",penalty_cache)
} else {
need <- unique(c("eid", "ethnic.c", vars.basic, vars.le8, "sleep_duration", "sleep.hours", "sleep",
                 "birth_date", "date_attend", "date_lost", "date_death", paste0("fod_icd10_", Y)))
dat <- read_all(need) |> filter_analysis_cohort() |> make_outcome(Y)
tvar <- paste0(Y, ".t2e"); evar <- paste0(Y, ".Yt2e"); covs <- intersect(vars.basic, names(dat))

load_oof_score <- function(layer) {
  score_name <- paste0(layer, "_oof_score")
  empty_score <- function() tibble::tibble(eid = dat$eid[0])
  od <- if (layer == "prot") out.prot else out.met
  f <- file.path(le8_job_dir(od, "c5_consolidate"), "c5.res.rds")
  if (!file.exists(f)) return(empty_score())
  z <- readRDS(f); pr <- z$scores$rows %||% tibble()
  if (!all(c("method", "eid", "score_z", "split") %in% names(pr))) return(empty_score())
  preferred <- c("C4 YSplus", "Parsimonious", "Pradeep / glmnet")
  best <- preferred[preferred %in% unique(pr$method)][1]
  if (!length(best) || is.na(best)) return(empty_score())
  # Training predictions are out-of-fold and validation predictions are fully
  # held out; prevalent rows are excluded from this incident-risk analysis.
  sc <- pr |> filter(method == best, split %in% c("training", "validation"), is.finite(score_z)) |>
    transmute(eid, score = score_z) |> distinct(eid, .keep_all = TRUE)
  names(sc)[2] <- score_name
  sc
}
for (layer in c(if (prot_DO) "prot" else NULL, if (met_DO) "met" else NULL)) dat <- left_join(dat, load_oof_score(layer), by = "eid")

vars0 <- unique(c(intersect(c("sleep_duration", "sleep.hours", "sleep"), names(dat)), intersect(vars.le8, names(dat)), grep("_oof_score$", names(dat), value = TRUE)))
fit_nonlin <- function(v, return_curve = FALSE) {
  cv <- setdiff(covs, v); d <- dat[, unique(c(tvar, evar, v, cv)), drop = FALSE]; d <- d[complete.cases(d), , drop = FALSE]
  d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
  if (nrow(d) < 1000 || sum(d[[evar]] == 1) < 50 || !is.finite(sd(d[[v]])) || sd(d[[v]]) == 0) return(tibble())
  rhs_lin <- c(paste0("scale(", bt(v), ")"), bt(cv)); rhs_spl <- c(paste0("splines::ns(", bt(v), ",df=4)"), bt(cv))
  fl <- tryCatch(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(rhs_lin, collapse = " + "))), d), error = function(e) NULL)
  fs <- tryCatch(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(rhs_spl, collapse = " + "))), d, x = TRUE), error = function(e) NULL)
  if (is.null(fl) || is.null(fs)) return(tibble())
  lr <- tryCatch(anova(fl, fs, test = "LRT"), error = function(e) NULL)
  pnon <- if (!is.null(lr) && nrow(lr) >= 2) lr$`Pr(>|Chi|)`[2] else NA_real_
  if (!return_curve) return(tibble(variable = v, n = nrow(d), events = sum(d[[evar]] == 1), p_nonlin = pnon, AIC_linear = AIC(fl), AIC_spline = AIC(fs)))
  xg <- seq(quantile(d[[v]], .02), quantile(d[[v]], .98), length.out = 120)
  nd <- data.frame(xg); names(nd) <- v
  for (cc in cv) {
    if (is.numeric(d[[cc]])) nd[[cc]] <- median(d[[cc]], na.rm = TRUE)
    else nd[[cc]] <- factor(names(sort(table(d[[cc]]), decreasing = TRUE))[1], levels = levels(d[[cc]]))
  }
  nd_ref <- nd[1, , drop = FALSE]; nd_ref[[v]] <- median(d[[v]], na.rm = TRUE)
  tt <- delete.response(terms(fs)); X <- model.matrix(tt, data = nd, xlev = fs$xlevels); X0 <- model.matrix(tt, data = nd_ref, xlev = fs$xlevels)
  cn <- names(coef(fs)); X <- X[, intersect(cn, colnames(X)), drop = FALSE]; X0 <- X0[, colnames(X), drop = FALSE]
  b <- coef(fs)[colnames(X)]; V <- vcov(fs)[colnames(X), colnames(X), drop = FALSE]
  XD <- sweep(X, 2, as.numeric(X0[1, ]), "-"); lp <- as.numeric(XD %*% b); se <- sqrt(pmax(0, rowSums((XD %*% V) * XD)))
  tibble(variable = v, x = xg, HR = exp(lp), lo = exp(lp - 1.96 * se), hi = exp(lp + 1.96 * se), p_nonlin = pnon)
}
res <- map_dfr(vars0, fit_nonlin, return_curve = FALSE)
if (!nrow(res)) res <- tibble(variable = character(), n = integer(), events = integer(), p_nonlin = numeric(), AIC_linear = numeric(), AIC_spline = numeric())
res <- res |> mutate(FDR_nonlin = p.adjust(p_nonlin, "BH")) |> arrange(p_nonlin)
curves <- map_dfr(vars0, fit_nonlin, return_curve = TRUE)
if (!nrow(curves)) curves <- tibble(variable = character(), x = numeric(), HR = numeric(), lo = numeric(), hi = numeric(), p_nonlin = numeric())
curves <- curves |> left_join(res |> select(variable, FDR_nonlin), by = "variable")
nadir <- if(!nrow(curves))tibble() else curves |> group_by(variable) |> slice_min(HR,n=1,with_ties=FALSE) |> ungroup() |>
  transmute(variable,nadir_x=x,nadir_HR=HR,p_nonlin,FDR_nonlin)
res <- res |> left_join(nadir |> select(variable,nadir_x,nadir_HR),by="variable")
write_raw_csv(res, "s2.nonlin_tests.csv", rawdir); write_raw_csv(curves, "s2.nonlin_curves.csv", rawdir)
pA <- if(!nrow(res))blank_plot("A. Evidence for non-linearity") else ggplot(res,aes(-log10(pmax(p_nonlin,1e-300)),fct_reorder(variable,p_nonlin,.desc=TRUE)))+
  geom_vline(xintercept=-log10(.05/nrow(res)),linetype=2,color="grey55")+geom_segment(aes(x=0,xend=-log10(pmax(p_nonlin,1e-300)),yend=variable),color="grey75")+
  geom_point(aes(color=FDR_nonlin<.05),size=2)+scale_color_manual(values=c(`TRUE`="#D95F02",`FALSE`="#4C78A8"),guide="none")+
  labs(title="A. Likelihood-ratio evidence for non-linearity",subtitle="Spline (4 df) versus linear Cox term; dashed line is Bonferroni 0.05/N",x=expression(-log[10](P[nonlinear])),y=NULL)+theme_5c(9)
pB <- if (!nrow(curves)) blank_plot("B. Non-linear risk chart", "No variable met the sample/event requirements") else
  ggplot(curves, aes(x, HR)) + geom_hline(yintercept = 1, linetype = 2, color = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80", alpha = .55) + geom_line(color = "#2C7FB8", linewidth = 1) +
  geom_point(data=nadir,aes(nadir_x,nadir_HR),inherit.aes=FALSE,color="#D95F02",size=1.7)+
  facet_wrap(~variable, scales = "free_x", ncol = 3) + scale_y_log10() +
  labs(title = "B. Non-linear LE8 and omics-score risk chart", subtitle = "Natural-spline HRs (95% CI), median reference; orange point is the fitted minimum within the central 96%", x = NULL, y = "Hazard ratio (log scale)") + theme_5c(9)
p<-pA|pB
save_plot(p, "s2_nonlin.Fig01.spline_patterns.png", 12, 9, outdir = outdir)
write_xlsx2(list(nonlin_tests = res, spline_curves = curves,nadir=nadir), "s2_nonlin.Fig01.spline_patterns.out.xlsx")
saveRDS(list(tests=res,curves=curves,nadir=nadir),nonlin_cache,compress="xz")

# Nested-CV alternatives to the simple LE8 sum, including bottleneck and
# cross-validated penalty rules. This is part of the same S2 job and output tree.
penalty_dat<-dat;comps<-intersect(vars.le8,names(penalty_dat));covs<-intersect(vars.basic,names(penalty_dat))
if(length(comps)<2)stop("S2 penalty analysis requires at least two LE8 component variables.",call.=FALSE)
penalty_dat<-penalty_dat[complete.cases(penalty_dat[,unique(c(tvar,evar,comps,covs)),drop=FALSE])&penalty_dat[[tvar]]>0,,drop=FALSE]
K<-as.integer(Sys.getenv("S2_PENALTY_FOLDS",unset="5"));K_INNER<-as.integer(Sys.getenv("S2_PENALTY_INNER_FOLDS",unset="4"));gammas<-c(0,.02,.04,.06,.08,.10,.15,.20,.30)
gamma_name<-function(g)paste0("penalty_",gsub("\\.","p",format(g,scientific=FALSE,trim=TRUE)))
fit_scaler<-function(x){lo<-vapply(x,min,numeric(1),na.rm=TRUE);hi<-vapply(x,max,numeric(1),na.rm=TRUE);list(lo=lo,span=pmax(hi-lo,1e-8))}
apply_scaler<-function(x,scaler){z<-sweep(as.matrix(x),2,scaler$lo,"-");z<-sweep(z,2,scaler$span,"/");pmin(pmax(z,0),1)}
add_scores<-function(d,z,gamma_values=gammas){d$LE8_mean<-rowMeans(z);d$LE8_min<-apply(z,1,min);d$LE8_geomean<-exp(rowMeans(log(pmax(z,.01))));d$n_fail<-rowSums(z<.4);d$LE8_passfail<-d$LE8_mean-.08*pmax(d$n_fail-1,0);for(g in gamma_values)d[[gamma_name(g)]]<-d$LE8_mean-g*pmax(d$n_fail-1,0)^2;d}
prepare_pair<-function(tr,te){sc<-fit_scaler(tr[,comps,drop=FALSE]);list(tr=add_scores(tr,apply_scaler(tr[,comps,drop=FALSE],sc)),te=add_scores(te,apply_scaler(te[,comps,drop=FALSE],sc)))}
fit_test_score<-function(tr,te,score){
  m<-mean(tr[[score]],na.rm=TRUE);ss<-sd(tr[[score]],na.rm=TRUE);if(!is.finite(ss)||ss==0)ss<-1;tr[[score]]<-(tr[[score]]-m)/ss;te[[score]]<-(te[[score]]-m)/ss
  fit<-tryCatch(coxph(as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",paste(bt(c(score,covs)),collapse=" + "))),tr),error=function(e)NULL)
  if(is.null(fit))return(tibble(score=score,beta=NA_real_,HR=NA_real_,p=NA_real_,cindex=NA_real_));lp<-tryCatch(as.numeric(predict(fit,newdata=te,type="lp")),error=function(e)rep(NA_real_,nrow(te)));sm<-coef(summary(fit));tibble(score=score,beta=sm[score,"coef"],HR=exp(sm[score,"coef"]),p=sm[score,"Pr(>|z|)"],cindex=fcidx(Surv(te[[tvar]],te[[evar]]),lp))
}
select_gamma_inner<-function(tr,seed){
  folds<-make_folds(tr,evar,K_INNER,seed);rows<-lapply(seq_len(K_INNER),function(fd){a<-tr[folds!=fd,,drop=FALSE];b<-tr[folds==fd,,drop=FALSE];pp<-prepare_pair(a,b);map_dfr(gammas,function(g)fit_test_score(pp$tr,pp$te,gamma_name(g))|>mutate(gamma=g,inner_fold=fd))})
  perf<-bind_rows(rows)|>group_by(gamma)|>summarise(mean_cindex=mean(cindex,na.rm=TRUE),valid_folds=sum(is.finite(cindex)),.groups="drop")|>arrange(desc(mean_cindex),gamma);if(!nrow(perf)||!any(is.finite(perf$mean_cindex)))return(list(gamma=0,performance=perf));list(gamma=perf$gamma[which.max(perf$mean_cindex)],performance=perf)
}
fold<-make_folds(penalty_dat,evar,K,SEED);rows<-list();chosen<-list();inner_all<-list();rr<-0L;base_rules<-c("LE8_mean","LE8_min","LE8_geomean","LE8_passfail")
for(fd in seq_len(K)){
  tr0<-penalty_dat[fold!=fd,,drop=FALSE];te0<-penalty_dat[fold==fd,,drop=FALSE];sel<-select_gamma_inner(tr0,SEED+fd);inner_best<-if(nrow(sel$performance)&&any(is.finite(sel$performance$mean_cindex)))max(sel$performance$mean_cindex,na.rm=TRUE)else NA_real_
  chosen[[fd]]<-tibble(fold=fd,chosen_gamma=sel$gamma,chosen_penalty=gamma_name(sel$gamma),inner_mean_cindex=inner_best);inner_all[[fd]]<-sel$performance|>mutate(outer_fold=fd);pp<-prepare_pair(tr0,te0)
  for(score in c(base_rules,gamma_name(sel$gamma))){rr<-rr+1L;rows[[rr]]<-fit_test_score(pp$tr,pp$te,score)|>mutate(fold=fd,gamma=ifelse(str_detect(score,"^penalty_"),sel$gamma,NA_real_),rule=ifelse(str_detect(score,"^penalty_"),"CV-selected big penalty",score))}
}
perf<-bind_rows(rows);choose<-bind_rows(chosen);inner_perf<-bind_rows(inner_all);summ<-perf|>group_by(rule)|>summarise(mean_cindex=mean(cindex,na.rm=TRUE),sd_cindex=sd(cindex,na.rm=TRUE),mean_HR=mean(HR,na.rm=TRUE),.groups="drop")|>arrange(desc(mean_cindex))
write_raw_csv(perf,"s2.penalty_CV_by_fold.csv",rawdir);write_raw_csv(choose,"s2.penalty_selected_gamma.csv",rawdir);write_raw_csv(inner_perf,"s2.penalty_inner_CV.csv",rawdir);write_raw_csv(summ,"s2.penalty_summary.csv",rawdir)
pA<-if(!nrow(perf))blank_plot("Cross-validated LE8 aggregation algorithms")else ggplot(perf,aes(cindex,fct_reorder(rule,cindex,mean),color=rule))+geom_point(position=position_jitter(height=.08),alpha=.6)+stat_summary(fun=mean,geom="point",size=3)+labs(title="A. Nested-CV LE8 aggregation algorithms",x="Held-out C-index",y=NULL,color=NULL)+theme_5c(10)+theme(legend.position="none")
best_gamma<-if(nrow(choose))as.numeric(names(sort(table(choose$chosen_gamma),decreasing=TRUE))[1])else 0;full<-prepare_pair(penalty_dat,penalty_dat)$tr;bestp<-gamma_name(best_gamma);full$best_penalty<-full[[bestp]];full$risk_group<-factor(ntile(-full$best_penalty,5),levels=1:5,labels=paste0("Risk Q",1:5))
fit<-tryCatch(coxph(as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ risk_group + ",paste(bt(covs),collapse=" + "))),full),error=function(e)NULL)
hr<-if(is.null(fit))tibble()else{sm<-coef(summary(fit));tibble(term=rownames(sm),HR=exp(sm[,"coef"]),lo=exp(sm[,"coef"]-1.96*sm[,"se(coef)"]),hi=exp(sm[,"coef"]+1.96*sm[,"se(coef)"]),p=sm[,"Pr(>|z|)"])|>mutate(group=str_extract(term,"[2-5]"),group=factor(group,levels=2:5))}
pB<-ggplot(full,aes(n_fail,best_penalty))+geom_boxplot(aes(group=factor(n_fail)),fill="grey88",outlier.alpha=.08)+geom_smooth(method="loess",se=FALSE,color="#D95F02")+labs(title=paste0("B. Most often selected gamma = ",best_gamma),x="Number of LE8 pillars below 40%",y="Penalized LE8 score")+theme_5c(10)
pC<-if(!nrow(hr))blank_plot("C. Disease risk across penalized-score quintiles")else ggplot(hr,aes(HR,group))+geom_vline(xintercept=1,color="grey60")+geom_errorbarh(aes(xmin=lo,xmax=hi),height=.15)+geom_point(size=2,color="#4C78A8")+scale_x_log10()+labs(title="C. Disease risk across penalized-score quintiles",x="Hazard ratio versus lowest-risk quintile",y="Higher-risk quintile")+theme_5c(10)
save_plot(pA/(pB|pC)+plot_layout(heights=c(.8,1.1)),"s2_nonlin.Fig02.pass_fail_penalty.png",13,9,outdir=outdir);write_xlsx2(list(CV_by_fold=perf,selected_penalty=choose,inner_CV=inner_perf,summary=summ,risk_group_HR=hr),"s2_nonlin.Fig02.pass_fail_penalty.out.xlsx")
saveRDS(list(performance=perf,selected=choose,inner=inner_perf,summary=summ,risk_HR=hr),penalty_cache,compress="xz");finalize_outputs(LE8_JOB,outdir)
}
