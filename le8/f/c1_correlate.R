# C1: Correlation — incident/prospective omics associations, baseline-prevalent
# associations, Yin-Yang temporal patterns, data-driven trajectory clusters,
# adjusted risk gradients, and functional enrichment.
#
# Key changes in this version
#   1) prevalent disease is analysed by baseline logistic regression, not Cox.
#   2) baseline = 0 on all YY plots: b2e < 0 prevalent; b2e > 0 incident.
#   3) basic vs adj2 attenuation is recomputed on the SAME complete-case sample.
#   4) YY spline annotations use a global trajectory P value, not one spline basis beta.
#   5) C1_CLUSTER_TOP defaults to 100; k is selected using silhouette + gap + stability.
#   6) quintile HRs use the same adj2 covariates as the primary C1 association model.
#   7) old trajectory Fig7 is removed; enrichment is now Fig7.

suppressPackageStartupMessages({
  source(file.path(Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f")), "comm.f.R"))
})
LE8_JOB <- "c1_correlate"

TOP_N              <- as.integer(Sys.getenv("C1_TOP_N", unset = "30"))
YY_TOP              <- as.integer(Sys.getenv("C1_YY_TOP", unset = "6"))
GRADIENT_TOP        <- as.integer(Sys.getenv("C1_GRADIENT_TOP", unset = "20"))
CLUSTER_TOP         <- as.integer(Sys.getenv("C1_CLUSTER_TOP", unset = "100"))
ATTENUATION_TOP     <- as.integer(Sys.getenv("C1_ATTENUATION_TOP", unset = "500"))
YY_BINS             <- as.integer(Sys.getenv("C1_YY_BINS", unset = "28"))
YY_MAX_YEAR         <- as.numeric(Sys.getenv("C1_YY_MAX_YEAR", unset = "20"))
GRADIENT_STEP       <- as.numeric(Sys.getenv("C1_GRADIENT_STEP", unset = "0.25"))
CLUSTER_STEP        <- as.numeric(Sys.getenv("C1_CLUSTER_STEP", unset = "0.5"))
MIN_BIN_N           <- as.integer(Sys.getenv("C1_MIN_BIN_N", unset = "20"))
CLUSTER_K_MAX       <- as.integer(Sys.getenv("C1_CLUSTER_K_MAX", unset = "8"))
CLUSTER_GAP_B       <- as.integer(Sys.getenv("C1_CLUSTER_GAP_B", unset = "30"))
CLUSTER_STABILITY_B <- as.integer(Sys.getenv("C1_CLUSTER_STABILITY_B", unset = "40"))
C1_MIN_EVENT        <- 20L

# Reproduce the validated backup style used by the former Fig3/Fig4: binned
# means are smoothed only inside the observed range, the gradient is shown as
# a colour/size bubble matrix, and clusters retain faint individual curves.
make_backup_gradient_style <- function(dat, features, bvar, layer, year_max=14.4,
                                       year_step=.2, k=3L, min_bin_n=MIN_BIN_N) {
  features<-intersect(features,names(dat)); if(length(features)<2)return(NULL)
  d<-dat|>filter(is.finite(.data[[bvar]]),.data[[bvar]]<0)|>mutate(year=-.data[[bvar]])|>
    filter(year>=1,year<=year_max)
  max_bin<-floor(year_max/year_step+sqrt(.Machine$double.eps));min_bin<-ceiling(1/year_step-sqrt(.Machine$double.eps))
  raw<-map_dfr(features,function(nm){
    x<-suppressWarnings(as.numeric(d[[nm]]));z<-as.numeric(scale(x));bid<-floor(d$year/year_step+sqrt(.Machine$double.eps))
    tibble(bin=bid,z=z)|>filter(is.finite(z),bin>=min_bin,bin<=max_bin)|>group_by(bin)|>
      summarise(z=mean(z),N=n(),.groups="drop")|>filter(N>=min_bin_n)|>mutate(feature=nm,year=bin*year_step)
  })
  smooth<-raw|>group_by(feature)|>group_modify(~{
    z<-.x|>arrange(year); if(nrow(z)<4)return(tibble())
    fit<-tryCatch(loess(z~year,data=z,span=.55,degree=2),error=function(e)NULL)
    pred<-if(is.null(fit))z$z else predict(fit,newdata=z$year)
    z|>mutate(estimate=ifelse(is.finite(pred),pred,.data$z))|>select(year,N,z,estimate)
  })|>ungroup()
  if(!nrow(smooth)||n_distinct(smooth$feature)<2)return(NULL)
  mat<-smooth|>select(feature,year,estimate)|>pivot_wider(names_from=year,values_from=estimate)|>as.data.frame()
  rownames(mat)<-mat$feature;mat<-as.matrix(mat[,-1,drop=FALSE])
  for(i in seq_len(nrow(mat))){ok<-which(is.finite(mat[i,]));if(length(ok)>=2)mat[i,!is.finite(mat[i,])]<-approx(ok,mat[i,ok],xout=which(!is.finite(mat[i,])),rule=2)$y}
  mat<-t(scale(t(mat)));mat[!is.finite(mat)]<-0
  kk<-min(as.integer(k),nrow(mat));cl<-tibble(feature=rownames(mat),cluster=cutree(hclust(dist(mat),method="ward.D2"),k=kk))
  bubble<-smooth|>mutate(feature=factor(feature,levels=rev(features)),years_before=-year,size=pmin(abs(estimate),.8))
  p1<-ggplot(bubble,aes(years_before,feature))+geom_point(aes(size=size,color=estimate),alpha=.72)+
    scale_size_continuous(range=c(.7,6),name="|mean z|")+scale_color_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0,name="mean z")+
    labs(title=paste0("Pre-diagnostic ",layer," gradients"),x=paste0("Years before ",Y," diagnosis"),y=NULL)+theme_5c(12)
  lines<-smooth|>left_join(cl,by="feature")|>mutate(years_before=-year)
  cnt<-cl|>count(cluster,name="n");lines<-lines|>left_join(cnt,by="cluster")|>mutate(panel=paste0("Cluster ",cluster," (N=",n,")"))
  p2<-ggplot(lines,aes(years_before,estimate,group=feature))+geom_line(alpha=.32,linewidth=.55,color="grey45")+
    stat_summary(aes(group=cluster,color=factor(cluster)),fun=mean,geom="line",linewidth=1.5)+facet_wrap(~panel,nrow=1,scales="free_y")+
    scale_color_brewer(palette="Dark2",guide="none")+labs(title=paste0("Clustered pre-diagnostic ",layer," patterns"),x=paste0("Years before ",Y," diagnosis"),y="Mean z-score")+theme_5c(11)
  list(bubble=p1,clusters=p2,smooth=smooth,cluster=cl)
}

# Publication typography. This overrides the shared theme only inside C1.
theme_5c <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size * 1.14, hjust = 0),
      plot.subtitle = element_text(face = "bold", size = base_size * .94, color = "grey30"),
      axis.title = element_text(face = "bold", size = base_size * 1.08),
      axis.text = element_text(face = "bold", size = base_size, color = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size * 1.02),
      panel.grid.major.y = element_line(color = "grey91", linewidth = .25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(9, 13, 9, 13)
    )
}
forest_theme <- function(base_size = 10) theme_5c(base_size) + theme(panel.grid.major.y = element_blank())

assoc_has_results <- function(x) is.data.frame(x) && nrow(x) > 0 && any(is.finite(x$p.value))

assoc_empty_message <- function(x, case_label, min_event = C1_MIN_EVENT) {
  nt <- suppressWarnings(max(x$N_total, na.rm = TRUE)); ne <- suppressWarnings(max(x$N_event, na.rm = TRUE))
  if (!is.finite(nt)) nt <- NA_integer_; if (!is.finite(ne)) ne <- NA_integer_
  paste0("Complete-case N = ", format(nt, big.mark = ","), "\n",
         case_label, " cases = ", format(ne, big.mark = ","),
         if (is.finite(ne) && ne < min_event) paste0(" (<", min_event, " required)") else "")
}

assoc_blank_plot <- function(title, message) {
  ggplot() + annotate("text", x = 0, y = .16, label = title, fontface = "bold", size = 7) +
    annotate("text", x = 0, y = -.12, label = message, fontface = "bold", color = "grey25", size = 5) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void(base_size = 16)
}

# ----------------------------------------------------------------------------
# Baseline prevalent association: logistic regression
# ----------------------------------------------------------------------------
make_prevalent_status <- function(dat, outcome = Y) {
  ydate <- paste0("fod_icd10_", outcome)
  if (!all(c(ydate, "date_attend") %in% names(dat))) stop("Cannot construct baseline prevalent status.", call. = FALSE)
  yd <- as.Date(dat[[ydate]]); ba <- as.Date(dat$date_attend)
  # Same-day diagnoses are deliberately excluded from the prevalent comparison.
  case_when(!is.na(yd) & !is.na(ba) & yd < ba ~ 1,
            is.na(yd) ~ 0,
            !is.na(yd) & !is.na(ba) & yd > ba ~ 0,
            TRUE ~ NA_real_)
}

logistic_scan <- function(dat, xs, covars, y, scale_x = TRUE, min_n = 500, min_case = 20) {
  xs <- intersect(xs, names(dat)); covars <- intersect(covars, names(dat))
  bind_rows(parallel_map(xs, function(x) {
    d <- dat[, unique(c(y, x, covars)), drop = FALSE]
    d <- d[complete.cases(d), , drop = FALSE]
    nc <- sum(d[[y]] == 1, na.rm = TRUE)
    empty <- tibble(term = x, estimate = NA_real_, beta = NA_real_, std.error = NA_real_,
                    conf.low = NA_real_, conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_,
                    N_total = nrow(d), N_event = nc)
    if (nrow(d) < min_n || nc < min_case || length(unique(d[[y]])) < 2) return(empty)
    d[[x]] <- suppressWarnings(as.numeric(d[[x]])); sx <- sd(d[[x]], na.rm = TRUE)
    if (!is.finite(sx) || sx <= 0) return(empty)
    if (scale_x) d[[x]] <- as.numeric(scale(d[[x]]))
    f <- reformulate(c(x, covars), response = y)
    fit <- tryCatch(glm(f, data = d, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(empty)
    sm <- coef(summary(fit)); if (!x %in% rownames(sm)) return(empty)
    b <- sm[x, "Estimate"]; se <- sm[x, "Std. Error"]
    tibble(term = x, estimate = exp(b), beta = b, std.error = se,
           conf.low = exp(b - 1.96 * se), conf.high = exp(b + 1.96 * se),
           statistic = b / se, p.value = 2 * pnorm(abs(b / se), lower.tail = FALSE),
           N_total = nrow(d), N_event = nc)
  })) |> mutate(FDR = p.adjust(p.value, "BH")) |> arrange(p.value)
}

# Same-N attenuation: both Cox models use the complete adj2 sample for each biomarker.
same_sample_attenuation <- function(dat, features, tvar, evar, covs_basic, covs_adj2) {
  features <- intersect(features, names(dat))
  bind_rows(parallel_map(features, function(x) {
    need <- unique(c(tvar, evar, x, covs_adj2)); d <- dat[, need, drop = FALSE]
    d <- d[complete.cases(d), , drop = FALSE]
    ne <- sum(d[[evar]] == 1, na.rm = TRUE)
    if (nrow(d) < 500 || ne < C1_MIN_EVENT) return(tibble(term=x,N=nrow(d),events=ne,
      beta_basic_sameN=NA_real_,p_basic_sameN=NA_real_,beta_adj2_sameN=NA_real_,p_adj2_sameN=NA_real_,attenuation_pct=NA_real_))
    d[[x]] <- as.numeric(scale(suppressWarnings(as.numeric(d[[x]]))))
    fit_one <- function(covars) {
      f <- as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(c(x, covars)), collapse=" + ")))
      fit <- tryCatch(coxph(f, d, ties="efron"), error=function(e) NULL)
      if (is.null(fit)) return(c(beta=NA_real_,p=NA_real_))
      sm <- coef(summary(fit)); if (!x %in% rownames(sm)) return(c(beta=NA_real_,p=NA_real_))
      c(beta=sm[x,"coef"], p=sm[x,"Pr(>|z|)"])
    }
    a <- fit_one(covs_basic); b <- fit_one(covs_adj2)
    tibble(term=x,N=nrow(d),events=ne,beta_basic_sameN=a[["beta"]],p_basic_sameN=a[["p"]],
           beta_adj2_sameN=b[["beta"]],p_adj2_sameN=b[["p"]],
           attenuation_pct=100*(abs(a[["beta"]])-abs(b[["beta"]]))/pmax(abs(a[["beta"]]),1e-12))
  })) |> arrange(p_basic_sameN)
}

# ----------------------------------------------------------------------------
# Yin-Yang: baseline = 0, prevalent < 0, incident > 0
# ----------------------------------------------------------------------------
residualize_z <- function(d, x, covars) {
  covars <- intersect(covars, names(d)); z <- suppressWarnings(as.numeric(d[[x]]))
  if (!length(covars)) return(as.numeric(scale(z)))
  dd <- d[, c(x, covars), drop=FALSE]; dd[[x]] <- z
  f <- reformulate(covars, response=x)
  fit <- tryCatch(lm(f, dd, na.action=na.exclude), error=function(e) NULL)
  if (is.null(fit)) return(as.numeric(scale(z)))
  as.numeric(scale(residuals(fit)))
}

global_trajectory_p <- function(time, z, df = 3) {
  ok <- is.finite(time) & is.finite(z)
  if (sum(ok) < 80 || length(unique(time[ok])) < 8) return(NA_real_)
  d <- data.frame(time=time[ok], z=z[ok])
  f0 <- tryCatch(lm(z ~ 1, d), error=function(e) NULL)
  f1 <- tryCatch(lm(z ~ splines::ns(time, df=df), d), error=function(e) NULL)
  if (is.null(f0) || is.null(f1)) return(NA_real_)
  a <- tryCatch(anova(f0,f1), error=function(e) NULL)
  if (is.null(a) || nrow(a) < 2) NA_real_ else as.numeric(a$`Pr(>F)`[2])
}

make_yy_panel <- function(dat, features, bvar, covars, title, ylab="Mean biomarker z-score",
                          bins=YY_BINS, max_year=YY_MAX_YEAR, min_bin_n=MIN_BIN_N) {
  features <- intersect(features, names(dat)); if (!length(features)) return(blank_plot(title,"No feature available"))
  d0 <- dat |> filter(is.finite(.data[[bvar]]))
  lim <- min(max_year, max(abs(d0[[bvar]]), na.rm=TRUE)); if (!is.finite(lim) || lim <= 0) lim <- max_year
  br <- seq(-lim, lim, length.out=bins+1); mids <- (head(br,-1)+tail(br,-1))/2
  hist0 <- hist(d0[[bvar]], breaks=br, plot=FALSE)
  h <- tibble(mid=mids, xmin=head(br,-1), xmax=tail(br,-1), n=hist0$counts,
              side=ifelse(mids<0,"Prevalent (Yang)","Incident (Yin)"))
  rows <- map_dfr(features, function(x) {
    need <- unique(c(bvar,x,covars)); dd <- d0[, need, drop=FALSE]; dd <- dd[complete.cases(dd),,drop=FALSE]
    if (nrow(dd)<80) return(tibble())
    dd$z <- residualize_z(dd,x,covars)
    dd$bin <- cut(dd[[bvar]], br, include.lowest=TRUE, labels=FALSE)
    pg <- global_trajectory_p(dd[[bvar]],dd$z)
    sm <- dd |> filter(is.finite(z),!is.na(bin)) |> group_by(bin) |>
      summarise(mean=mean(z),sd=sd(z),se=sd/sqrt(n()),N=n(),.groups="drop") |>
      filter(N>=min_bin_n) |> mutate(time=mids[bin],feature=x,p_global=pg,N_total=nrow(dd))
    sm
  })
  if (!nrow(rows)) return(blank_plot(title,"No bin reached the minimum sample size"))
  labs0 <- rows |> group_by(feature) |> summarise(p_global=first(p_global),N_total=first(N_total),.groups="drop") |>
    mutate(label=sprintf("%s (N=%s, Pglobal=%s)",feature,format(N_total,big.mark=","),
                         ifelse(is.finite(p_global),format.pval(p_global,digits=2,eps=1e-99),"NA")))
  rows <- rows |> left_join(labs0 |> select(feature,label),by="feature")
  yr <- range(rows$mean,na.rm=TRUE); if(!all(is.finite(yr))||diff(yr)<.2) yr <- c(-1,1)
  # Put the actual case-count histogram on the primary scale and biomarker
  # trajectories on the secondary scale, matching yy.R Fig1 panel a.
  pad <- .14*diff(yr); yr <- yr+c(-pad,pad); hmax <- max(h$n,1)
  plot_max <- hmax*1.12
  scale_factor <- plot_max/diff(yr)
  z_to_count <- function(z)(z-yr[1])*scale_factor
  rows <- rows |> mutate(y_plot=z_to_count(mean))
  ggplot() +
    geom_col(data=h,aes(mid,n,fill=side),width=diff(br)[1]*.98,color="white",alpha=.55,inherit.aes=FALSE) +
    geom_vline(xintercept=0,linetype=2,color="grey38",linewidth=.75) +
    geom_hline(yintercept=z_to_count(0),linetype=3,color="grey62") +
    geom_line(data=rows,aes(time,y_plot,color=label,group=label),linewidth=.95) +
    geom_point(data=rows,aes(time,y_plot,color=label),size=2.45) +
    scale_fill_manual(values=c(`Prevalent (Yang)`="grey58",`Incident (Yin)`="#78B7E5"),guide="none") +
    scale_y_continuous(name="Patient count",limits=c(0,plot_max),
      sec.axis=sec_axis(~./scale_factor+yr[1],name=ylab)) +
    labs(title=title,x="Years relative to baseline",color=NULL) +
    coord_cartesian(xlim=c(-lim,lim),clip="off") + theme_5c(11) +
    guides(color=guide_legend(ncol=2,byrow=TRUE,override.aes=list(linewidth=1.1,size=2.8))) +
    theme(legend.position=c(.985,.985),legend.justification=c(1,1),legend.direction="horizontal",legend.box="horizontal",
          legend.background=element_rect(fill=scales::alpha("white",.82),color="grey75",linewidth=.25),
          legend.text=element_text(size=8.2,face="bold"),legend.margin=margin(3,4,3,4))
}

# ----------------------------------------------------------------------------
# Fine temporal gradient (no extrapolated smoothed curve)
# ----------------------------------------------------------------------------
make_gradient_panel <- function(dat, features, bvar, side=c("incident","prevalent"), step=.25,
                                min_bin_n=20, title=NULL) {
  side <- match.arg(side); features <- intersect(features,names(dat))
  dd <- dat |> filter(is.finite(.data[[bvar]]), if(side=="incident") .data[[bvar]]>0 else .data[[bvar]]<0)
  rows <- map_dfr(features,function(x){
    x0 <- suppressWarnings(as.numeric(dd[[x]])); if(sum(is.finite(x0))<50)return(tibble())
    z <- as.numeric(scale(x0)); tb <- floor(dd[[bvar]]/step)*step + step/2
    tibble(time=tb,z=z) |> filter(is.finite(time),is.finite(z)) |> group_by(time) |>
      summarise(mean=mean(z),N=n(),.groups="drop") |> filter(N>=min_bin_n) |> mutate(feature=x)
  })
  if(!nrow(rows)) return(blank_plot(title %||% "Temporal gradient","No sufficiently populated time bins"))
  ord <- rev(features[features %in% unique(rows$feature)]); rows$feature <- factor(rows$feature,levels=ord)
  ggplot(rows,aes(time,feature)) +
    geom_vline(xintercept=0,linetype=2,color="grey45") +
    geom_point(aes(size=pmin(N,200),fill=mean),shape=21,color="grey55",stroke=.15,alpha=.88) +
    scale_fill_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0,name="Mean z") +
    scale_size_continuous(range=c(.6,4),name="Bin N") +
    labs(title=title,x="Years relative to baseline",y=NULL) + theme_5c(10) +
    theme(legend.position="right",axis.text.y=element_text(size=8.5,face="bold"))
}

# ----------------------------------------------------------------------------
# Publication-style YY circle matrix used as the first row of Fig4.
# No line joins adjacent circles: each dot is an observed time-bin summary.
# ----------------------------------------------------------------------------
make_gradient_circle_panel <- function(dat,features,bvar,step=GRADIENT_STEP,min_bin_n=MIN_BIN_N,
                                       max_year=YY_MAX_YEAR,title="a. Temporal biomarker gradients") {
  features<-intersect(features,names(dat)); d<-dat|>filter(is.finite(.data[[bvar]]),abs(.data[[bvar]])<=max_year)
  rows<-map_dfr(features,function(x){
    x0<-suppressWarnings(as.numeric(d[[x]])); if(sum(is.finite(x0))<50)return(tibble())
    z<-as.numeric(scale(x0)); tb<-floor(d[[bvar]]/step)*step+step/2
    tibble(time=tb,z=z)|>filter(is.finite(time),is.finite(z))|>group_by(time)|>
      summarise(mean=mean(z),N=n(),.groups="drop")|>filter(N>=min_bin_n)|>mutate(feature=x)
  })
  if(!nrow(rows))return(blank_plot(title,"No sufficiently populated time bins"))
  # Keep strongest proteins at the top, matching the conventional circle-matrix style.
  ord<-rev(features[features%in%unique(rows$feature)]);rows$feature<-factor(rows$feature,levels=ord)
  ggplot(rows,aes(time,feature))+
    geom_vline(xintercept=0,linetype=2,color="grey38",linewidth=.7)+
    geom_point(aes(size=N,fill=mean),shape=21,color="grey72",stroke=.22,alpha=.95)+
    scale_fill_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0,name="Mean z")+
    scale_size_continuous(range=c(.8,4.6),name="Bin N")+
    scale_x_continuous(breaks=pretty(c(-max_year,max_year),n=7))+
    labs(title=title,x="Years relative to baseline",y=NULL)+theme_5c(10)+
    theme(legend.position="right",axis.text.y=element_text(size=8.6,face="bold"))
}

# ----------------------------------------------------------------------------
# Data-driven clustering: silhouette + gap + temporal-bin bootstrap stability
# ----------------------------------------------------------------------------
adjusted_rand <- function(a,b) {
  tab <- table(a,b); n <- sum(tab); if(n<2)return(NA_real_)
  c2 <- function(x) x*(x-1)/2
  nij <- sum(c2(tab)); ai <- sum(c2(rowSums(tab))); bj <- sum(c2(colSums(tab))); tot <- c2(n)
  expected <- ai*bj/tot; denom <- .5*(ai+bj)-expected
  if(!is.finite(denom)||denom==0) return(NA_real_)
  (nij-expected)/denom
}

mean_silhouette <- function(distmat, cl) {
  n <- length(cl); if(n<3 || length(unique(cl))<2)return(NA_real_)
  s <- vapply(seq_len(n),function(i){
    same <- which(cl==cl[i] & seq_len(n)!=i); a <- if(length(same))mean(distmat[i,same]) else 0
    oth <- setdiff(unique(cl),cl[i]); b <- min(vapply(oth,function(g)mean(distmat[i,cl==g]),numeric(1)))
    if(max(a,b)==0)0 else (b-a)/max(a,b)
  },numeric(1))
  mean(s,na.rm=TRUE)
}

within_ss <- function(x,cl) {
  sum(vapply(unique(cl),function(g){z<-x[cl==g,,drop=FALSE]; if(nrow(z)<2)return(0);sum(rowSums((z-matrix(colMeans(z),nrow(z),ncol(z),byrow=TRUE))^2))},numeric(1)))
}

trajectory_matrix <- function(dat,features,bvar,step=.5,min_bin_n=15,max_year=YY_MAX_YEAR) {
  features<-intersect(features,names(dat)); d<-dat|>filter(is.finite(.data[[bvar]]),abs(.data[[bvar]])<=max_year)
  grid<-seq(-max_year+step/2,max_year-step/2,by=step)
  raw<-map_dfr(features,function(x){
    z<-as.numeric(scale(suppressWarnings(as.numeric(d[[x]])))); tb<-floor((d[[bvar]]+max_year)/step)*step-max_year+step/2
    tibble(time=tb,z=z)|>filter(is.finite(z),is.finite(time))|>group_by(time)|>summarise(mean=mean(z),N=n(),.groups="drop")|>
      filter(N>=min_bin_n)|>mutate(feature=x)
  })
  if(!nrow(raw))return(NULL)
  wide<-raw|>select(feature,time,mean)|>pivot_wider(names_from=time,values_from=mean)
  rn<-wide$feature; X<-as.matrix(wide[,-1,drop=FALSE]);rownames(X)<-rn
  # Internal interpolation only for clustering. Plotted raw trajectories are never extrapolated.
  for(i in seq_len(nrow(X))){
    ok<-which(is.finite(X[i,])); if(length(ok)>=2){
      miss<-which(!is.finite(X[i,]) & seq_len(ncol(X))>=min(ok) & seq_len(ncol(X))<=max(ok))
      if(length(miss))X[i,miss]<-approx(ok,X[i,ok],xout=miss,rule=1)$y
    }
    m<-mean(X[i,],na.rm=TRUE); if(!is.finite(m))m<-0; X[i,!is.finite(X[i,])]<-m
  }
  X<-t(scale(t(X)));X[!is.finite(X)]<-0
  list(X=X,raw=raw,grid=grid)
}

choose_cluster_k <- function(X,kmax=8,gap_B=30,stability_B=40,seed=SEED) {
  n<-nrow(X); ks<-2:min(kmax,n-1); if(!length(ks))return(list(k=1,metrics=tibble()))
  D<-as.matrix(dist(X)); hc<-hclust(as.dist(D),method="ward.D2")
  sil<-vapply(ks,function(k)mean_silhouette(D,cutree(hc,k)),numeric(1))
  # Gap statistic on a uniform reference box with the same feature dimensions.
  set.seed(seed); obsW<-vapply(ks,function(k)within_ss(X,cutree(hc,k)),numeric(1)); obsW<-pmax(obsW,1e-12)
  ref_logW<-matrix(NA_real_,nrow=gap_B,ncol=length(ks))
  lo<-apply(X,2,min);hi<-apply(X,2,max)
  for(b in seq_len(gap_B)){
    XR<-sweep(matrix(runif(n*ncol(X)),nrow=n),2,hi-lo,"*");XR<-sweep(XR,2,lo,"+")
    hR<-hclust(dist(XR),method="ward.D2")
    ref_logW[b,]<-vapply(ks,function(k)log(pmax(within_ss(XR,cutree(hR,k)),1e-12)),numeric(1))
  }
  gap<-colMeans(ref_logW,na.rm=TRUE)-log(obsW); gap_se<-apply(ref_logW,2,sd,na.rm=TRUE)*sqrt(1+1/gap_B)
  # Stability: resample time dimensions, cluster again, compare with original labels by ARI.
  stab<-vapply(ks,function(k){
    orig<-cutree(hc,k); z<-rep(NA_real_,stability_B)
    for(b in seq_len(stability_B)){
      # Sparse outcomes can leave only one or two populated temporal bins.
      # Never request more columns than exist when bootstrapping stability.
      n_take<-min(ncol(X),max(1L,ceiling(.8*ncol(X))))
      cols<-sort(sample(seq_len(ncol(X)),n_take,replace=FALSE))
      cb<-cutree(hclust(dist(X[,cols,drop=FALSE]),method="ward.D2"),k)
      z[b]<-adjusted_rand(orig,cb)
    }
    median(z,na.rm=TRUE)
  },numeric(1))
  met<-tibble(k=ks,silhouette=sil,gap=gap,gap_se=gap_se,stability=stab)
  rank_metric<-function(x)rank(-ifelse(is.finite(x),x,-Inf),ties.method="min")
  met<-met|>mutate(rank_sum=rank_metric(silhouette)+rank_metric(gap)+rank_metric(stability))
  best<-met|>arrange(rank_sum,desc(silhouette),desc(stability),desc(gap))|>slice(1)|>pull(k)
  list(k=best,metrics=met,hc=hc)
}

make_cluster_figure <- function(dat,features,bvar,step=CLUSTER_STEP,max_year=YY_MAX_YEAR,
                                side=c("incident","prevalent"),panel_label="b",fixed_k=NULL) {
  side<-match.arg(side)
  dat<-dat|>filter(if(side=="incident") .data[[bvar]]>0 else .data[[bvar]]<0)
  tm<-trajectory_matrix(dat,features,bvar,step=step,min_bin_n=max(10L,floor(MIN_BIN_N*.75)),max_year=max_year)
  if(is.null(tm)||nrow(tm$X)<3){
    p0<-blank_plot(paste0(panel_label,". ",ifelse(side=="incident","Incident","Baseline-prevalent")," trajectory clusters"),
                   "Insufficient trajectories")
    return(list(plot=p0,cluster_plot=p0,diagnostic_plot=p0,cluster=tibble(),metrics=tibble(),raw=tibble(),k=1L))
  }
  ck<-choose_cluster_k(tm$X,CLUSTER_K_MAX,CLUSTER_GAP_B,CLUSTER_STABILITY_B)
  if(!is.null(fixed_k))ck$k<-min(max(1L,as.integer(fixed_k)),nrow(tm$X))
  cl<-tibble(feature=rownames(tm$X),cluster=cutree(ck$hc,ck$k))
  lines<-tm$raw|>inner_join(cl,by="feature")
  cnt<-cl|>count(cluster,name="proteins")
  lines<-lines|>left_join(cnt,by="cluster")|>mutate(panel=factor(cluster,levels=sort(unique(cluster)),labels=paste0("Cluster ",sort(unique(cluster))," (N=",cnt$proteins[match(sort(unique(cluster)),cnt$cluster)],")")))
  pA<-ggplot(lines,aes(time,mean,group=feature))+
    geom_vline(xintercept=0,linetype=2,color="grey45")+geom_line(color="grey55",alpha=.24,linewidth=.42)+
    stat_summary(aes(group=cluster,color=factor(cluster)),fun=mean,geom="line",linewidth=1.45)+
    facet_wrap(~panel,nrow=1,scales="free_y")+scale_color_brewer(palette="Dark2",guide="none")+
    labs(title=paste0(panel_label,". ",ifelse(side=="incident","Incident","Baseline-prevalent"),
                      " trajectory clusters (top ",nrow(tm$X)," biomarkers; k = ",ck$k,")"),
         x="Years relative to baseline",y="Mean biomarker z-score")+theme_5c(10)
  # Do NOT min-max rescale the three criteria together. A constant stability
  # series can legitimately equal 1, but rescaling formerly collapsed the
  # display into horizontal 0/1 lines and hid the silhouette/gap structure.
  md<-ck$metrics|>select(k,silhouette,gap,stability)|>pivot_longer(-k,names_to="criterion",values_to="value")|>
    mutate(criterion=factor(criterion,levels=c("silhouette","gap","stability"),
      labels=c("Silhouette","Gap statistic","Bootstrap stability (ARI)")))
  bestpts<-md|>filter(k==ck$k)
  pB<-ggplot(md,aes(k,value,group=criterion))+geom_vline(xintercept=ck$k,linetype=2,color="grey40")+
    geom_line(linewidth=.85,color="#3F6F9F")+geom_point(size=2,color="#3F6F9F")+
    geom_point(data=bestpts,size=3.1,shape=21,fill="white",stroke=.9,color="#B2182B")+
    facet_wrap(~criterion,nrow=1,scales="free_y")+scale_x_continuous(breaks=unique(md$k))+
    labs(title="c. Cluster-number diagnostics",x="Number of clusters (k)",y="Native criterion value")+
    theme_5c(9)+theme(legend.position="none")
  list(plot=pA/pB+plot_layout(heights=c(3.2,1)),cluster_plot=pA,diagnostic_plot=pB,
       cluster=cl,metrics=ck$metrics,raw=tm$raw,k=ck$k)
}

# ----------------------------------------------------------------------------
# Adjusted quintile analysis
# ----------------------------------------------------------------------------
plot_quantile_top <- function(dat, features, tvar, evar, covars, title_prefix) {
  features <- intersect(features,names(dat)); covars<-intersect(covars,names(dat))
  map(features,function(x){
    need<-unique(c(tvar,evar,x,covars)); d<-dat[,need,drop=FALSE];d<-d[complete.cases(d),,drop=FALSE]
    if(nrow(d)<500||sum(d[[evar]]==1)<20)return(blank_plot(x,"Insufficient complete-case follow-up"))
    d$value<-suppressWarnings(as.numeric(d[[x]]));d$quintile<-factor(ntile(d$value,5),levels=1:5,labels=paste0("Q",1:5))
    sf<-tryCatch(survfit(Surv(d[[tvar]],d[[evar]])~d$quintile),error=function(e)NULL)
    if(is.null(sf))return(blank_plot(x,"Survival model failed"))
    ss<-summary(sf); sdf<-tibble(time=ss$time,surv=ss$surv,lo=ss$lower,hi=ss$upper,strata=as.character(ss$strata))|>
      mutate(quintile=str_extract(strata,"Q[1-5]"))
    f<-as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ relevel(quintile,'Q1') + ",paste(bt(covars),collapse=" + ")))
    fit<-tryCatch(coxph(f,d,ties="efron"),error=function(e)NULL);sm<-if(is.null(fit))NULL else coef(summary(fit))
    rn<-if(is.null(sm))character() else rownames(sm);i5<-grep("Q5",rn,fixed=TRUE)[1]
    ann<-if(is.na(i5)||!length(i5))"" else sprintf("Adjusted Q5 vs Q1: HR %.2f (%.2f–%.2f), P=%s",
      exp(sm[i5,"coef"]),exp(sm[i5,"coef"]-1.96*sm[i5,"se(coef)"]),exp(sm[i5,"coef"]+1.96*sm[i5,"se(coef)"]),
      format.pval(sm[i5,"Pr(>|z|)"],digits=2,eps=1e-99))
    ggplot(sdf,aes(time,1-surv,color=quintile,fill=quintile))+geom_step(linewidth=.7)+
      geom_ribbon(aes(ymin=1-hi,ymax=1-lo),alpha=.07,color=NA)+annotate("text",x=Inf,y=Inf,label=ann,hjust=1.03,vjust=1.5,size=3,fontface="bold")+
      scale_color_brewer(palette="Spectral",direction=-1)+scale_fill_brewer(palette="Spectral",direction=-1)+
      labs(title=x,x="Years after baseline",y="Cumulative incidence",color=NULL,fill=NULL)+
      theme_5c(10)+theme(legend.position="top")
  })
}

# ----------------------------------------------------------------------------
# Functional enrichment (same logic as previous C1; all assayed proteins as bg)
# ----------------------------------------------------------------------------
run_functional_enrichment <- function(assoc, layer, background=assoc$term) {
  sig<-assoc|>filter(is.finite(p.value),p.value*nrow(assoc)<.05)|>pull(term)|>unique()
  background<-intersect(unique(as.character(background)),unique(as.character(assoc$term)))
  empty<-tibble(source=character(),term_id=character(),term_name=character(),adjusted_p=numeric(),intersection_size=integer(),significant_n=integer(),background_n=integer())
  if(!length(sig)||layer!="protein")return(empty)
  ans<-NULL;used_custom_bg<-TRUE
  if(requireNamespace("gprofiler2",quietly=TRUE)){
    gp<-tryCatch(gprofiler2::gost(sig,organism="hsapiens",correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),custom_bg=background,domain_scope="custom",user_threshold=1),error=function(e)NULL)
    if(is.null(gp)||is.null(gp$result)||!nrow(gp$result)){
      used_custom_bg<-FALSE
      gp<-tryCatch(gprofiler2::gost(sig,organism="hsapiens",correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),domain_scope="annotated",user_threshold=1),error=function(e)NULL)
    }
    if(!is.null(gp)&&!is.null(gp$result)&&nrow(gp$result))ans<-as_tibble(gp$result)|>transmute(source,term_id=native,term_name=name,adjusted_p=p_value,intersection_size=as.integer(intersection_size))
  }
  if(is.null(ans)){
    py<-Sys.which("python3");script<-file.path(Sys.getenv("LE8_FDIR"),"c1_enrich.py")
    if(!nzchar(py)||!file.exists(script))return(empty)
    sf<-tempfile(fileext=".txt");bf<-tempfile(fileext=".txt");of<-tempfile(fileext=".csv");on.exit(unlink(c(sf,bf,of)),add=TRUE)
    writeLines(sig,sf);writeLines(background,bf);status<-suppressWarnings(system2(py,c(script,sf,bf,of),stdout=FALSE,stderr=FALSE))
    if(status==0&&file.exists(of))ans<-as_tibble(data.table::fread(of,showProgress=FALSE))
  }
  if(is.null(ans)||!nrow(ans))return(empty)
  ans|>mutate(significant_n=length(sig),background_n=if(used_custom_bg)length(background) else NA_integer_)|>arrange(adjusted_p)
}

plot_functional_enrichment <- function(enrich,n_sig,assoc=NULL) {
  if(!nrow(enrich)){
    d<-assoc|>filter(is.finite(p.value),p.value*nrow(assoc)<.05)|>slice_min(p.value,n=25,with_ties=FALSE)|>
      mutate(score=-log10(pmax(p.adjust(p.value,"BH"),.Machine$double.xmin)),direction=ifelse(beta>=0,"Positive","Negative"),term=factor(term,levels=rev(term)))
    if(!nrow(d))return(blank_plot(paste0("Functional enrichment of ",n_sig," significant biomarkers"),"No significant enrichment input was available"))
    return(ggplot(d,aes(score,term,fill=direction))+geom_col(width=.76)+
      scale_fill_manual(values=c(Positive="#D7301F",Negative="#2C7FB8"))+
      labs(title=paste0("Functional-enrichment input: ",n_sig," significant biomarkers"),
           subtitle="Offline fallback: GO/KEGG/Reactome service unavailable; bars show biomarker FDR, not pathway enrichment",
           x=expression(-log[10](biomarker~FDR)),y=NULL,fill="Association")+
      theme_5c(10)+theme(legend.position="top"))
  }
  levels0<-c("GO:BP","GO:MF","GO:CC","KEGG","REAC");labs0<-c(`GO:BP`="BP",`GO:MF`="MF",`GO:CC`="CC",KEGG="KEGG",REAC="Reactome")
  d<-enrich|>filter(is.finite(adjusted_p))|>group_by(source)|>slice_min(adjusted_p,n=5,with_ties=FALSE)|>ungroup()|>
    mutate(source=factor(source,levels=levels0,labels=unname(labs0[levels0])),score=-log10(pmax(adjusted_p,.Machine$double.xmin)),key=paste(source,term_id,sep="___"))|>
    arrange(source,score)|>mutate(key=factor(key,levels=unique(key)))
  ggplot(d,aes(key,score,fill=source))+geom_col(width=.78)+geom_text(aes(label=intersection_size),vjust=-.25,size=2.8,fontface="bold")+
    facet_grid(~source,scales="free_x",space="free_x")+scale_x_discrete(labels=setNames(str_wrap(d$term_name,18),as.character(d$key)))+
    labs(title=paste0("Functional enrichment of ",n_sig," Bonferroni-significant biomarkers"),subtitle=paste0(if(any(is.na(d$background_n)))"Standard annotated-genome background (custom assay background returned no terms)" else "Custom background = all assayed proteins","; ",sum(d$adjusted_p<.05)," displayed terms have FDR < 0.05"),x=NULL,y=expression(-log[10](FDR)),fill="Ontology")+
    theme_5c(10)+theme(legend.position="top",axis.text.x=element_text(angle=58,hjust=1,vjust=1,face="bold",size=8.5))
}

# ----------------------------------------------------------------------------
# Main C1
# ----------------------------------------------------------------------------
run_c1_layer <- function(layer=c("protein","metabolite")) {
  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir)
  rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
  # Remove obsolete numbering from previous C1 versions so output_index is unambiguous.
  unlink(file.path(outdir,c("c1.Fig3.gradient.png","c1.Fig4.gradient.png","c1.Fig5.gradient_top.png","c1.Fig5.gradient_cluster.png","c1.Fig4.gradient.png",
                            "c1.Fig7.trajectory_top.png","c1.Fig8.enrich_sig.png","c1.Fig01.mh.png","c1.Fig02.vc.png")),force=TRUE)
  scan_cache<-file.path(rawdir,"c1.scan.rds");selected_cache<-file.path(rawdir,"c1.res.rds")

  biom0<-if(layer=="protein")read_prot() else read_met();features_all<-setdiff(names(biom0),"eid")
  scan<-if(cache_valid(scan_cache))tryCatch(readRDS(scan_cache),error=function(e)NULL) else NULL
  reuse<-!is.null(scan)

  # If the expensive scans exist, load only the features needed for figures + same-N table.
  if(reuse){
    a0<-scan$association_adj2 %||% scan$association_basic; p0<-scan$prevalent_adj2 %||% scan$prevalent_basic
    fig_features<-unique(c(
      a0|>filter(is.finite(p.value))|>slice_min(p.value,n=max(CLUSTER_TOP,ATTENUATION_TOP,TOP_N),with_ties=FALSE)|>pull(term),
      p0|>filter(is.finite(p.value))|>slice_min(p.value,n=max(GRADIENT_TOP,TOP_N),with_ties=FALSE)|>pull(term)))
    biom<-biom0[,intersect(c("eid",fig_features),names(biom0)),drop=FALSE]
  } else biom<-biom0

  need<-unique(c("eid","ethnic.c",vars.basic,vars.le8,"birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y)))
  # IMPORTANT: inner_join restricts C1 to participants with the omics assay.
  dat<-read_all(need)|>inner_join(biom,by="eid")|>filter_analysis_cohort()|>make_outcome(Y)
  rm(biom,biom0);invisible(gc())
  features<-intersect(features_all,names(dat));covs_basic<-intersect(vars.basic,names(dat));covs_adj2<-intersect(vars.adj2,names(dat))
  tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e");bvar<-paste0(Y,".b2e")
  dat$prevalent_status<-make_prevalent_status(dat,Y)

  if(!reuse){
    message("C1/",layer,": incident Cox scans: basic + adj2")
    assoc_basic<-cox_scan(dat,features,covs_basic,Y,time_var=tvar,event_var=evar)
    assoc_adj2<-cox_scan(dat,features,covs_adj2,Y,time_var=tvar,event_var=evar)
    message("C1/",layer,": baseline-prevalent logistic scans: basic + adj2")
    prevalent_basic<-logistic_scan(dat,features,covs_basic,"prevalent_status")
    prevalent_adj2<-logistic_scan(dat,features,covs_adj2,"prevalent_status")
    # Same-N attenuation only needs the strongest features; set C1_ATTENUATION_TOP=0 for all.
    att_features<-assoc_basic|>filter(is.finite(p.value))|>arrange(p.value)|>pull(term)
    if(ATTENUATION_TOP>0)att_features<-head(att_features,ATTENUATION_TOP)
    attenuation_sameN<-same_sample_attenuation(dat,att_features,tvar,evar,covs_basic,covs_adj2)
    scan<-list(association_basic=assoc_basic,association_adj2=assoc_adj2,
               prevalent_basic=prevalent_basic,prevalent_adj2=prevalent_adj2,attenuation_sameN=attenuation_sameN)
    saveRDS(scan,scan_cache,compress="xz")
  } else {
    assoc_basic<-scan$association_basic;assoc_adj2<-scan$association_adj2
    prevalent_basic<-scan$prevalent_basic;prevalent_adj2<-scan$prevalent_adj2
    attenuation_sameN<-scan$attenuation_sameN %||% tibble()
    message("C1/",layer,": reuse association scans and regenerate figures")
  }

  assoc<-if(covs_use_name=="adj2")assoc_adj2 else assoc_basic
  assoc_prevalent<-if(covs_use_name=="adj2")prevalent_adj2 else prevalent_basic
  prefix<-ifelse(layer=="protein","pwas","mwas")
  write_raw_csv(assoc_basic,paste0(prefix,"_incident_basic.csv"),rawdir);write_raw_csv(assoc_adj2,paste0(prefix,"_incident_adj2.csv"),rawdir)
  write_raw_csv(prevalent_basic,paste0(prefix,"_prevalent_basic.csv"),rawdir);write_raw_csv(prevalent_adj2,paste0(prefix,"_prevalent_adj2.csv"),rawdir)
  write_raw_csv(attenuation_sameN,"adj2_attenuation_sameN.csv",rawdir)

  cohort<-tibble(layer,N_omics=nrow(dat),incident_events=sum(dat[[evar]]==1,na.rm=TRUE),
                 prevalent_cases=sum(dat$prevalent_status==1,na.rm=TRUE),features=length(features_all),
                 LE8_components=length(intersect(vars.le8,names(dat))),covs_use=covs_use_name)
  write_raw_csv(cohort,"c1.cohort.csv",rawdir)

  top<-assoc|>filter(is.finite(p.value))|>slice_min(p.value,n=TOP_N,with_ties=FALSE)|>pull(term)
  top6<-head(top,YY_TOP);top_inc<-assoc|>filter(is.finite(p.value))|>slice_min(p.value,n=GRADIENT_TOP,with_ties=FALSE)|>pull(term)
  top_prev<-assoc_prevalent|>filter(is.finite(p.value))|>slice_min(p.value,n=GRADIENT_TOP,with_ties=FALSE)|>pull(term)
  cluster_features<-assoc|>filter(is.finite(p.value))|>slice_min(p.value,n=CLUSTER_TOP,with_ties=FALSE)|>pull(term)

  # Fig1 + Fig2: four association views; prevalent panels are logistic OR analyses.
  if(layer=="protein"){
    p1a<-plot_pwas_manhattan(assoc_basic,features_all,paste0("a. Incident ",Y," — basic Cox"))
    p1b<-plot_pwas_manhattan(prevalent_basic,features_all,paste0("b. Baseline-prevalent ",Y," — basic logistic"))
    p1c<-plot_pwas_manhattan(assoc_adj2,features_all,paste0("c. Incident ",Y," — adj2 Cox"))
    p1d<-plot_pwas_manhattan(prevalent_adj2,features_all,paste0("d. Baseline-prevalent ",Y," — adj2 logistic"))
    save_plot((p1a|p1c)/(p1b|p1d),"c1.Fig1.mh.png",24,14,outdir=outdir)
    p2a<-plot_volcano(assoc_basic,"beta","term","a. Incident ProtWAS — basic",.05/max(1,nrow(assoc_basic)),TOP_N)
    p2b<-plot_volcano(prevalent_basic,"beta","term","b. Prevalent ProtWAS — basic logistic",.05/max(1,nrow(prevalent_basic)),TOP_N)
    p2c<-plot_volcano(assoc_adj2,"beta","term","c. Incident ProtWAS — adj2",.05/max(1,nrow(assoc_adj2)),TOP_N)
    p2d<-plot_volcano(prevalent_adj2,"beta","term","d. Prevalent ProtWAS — adj2 logistic",.05/max(1,nrow(prevalent_adj2)),TOP_N)
    save_plot((p2a|p2c)/(p2b|p2d),"c1.Fig2.vc.png",19,14,outdir=outdir)
  } else {
    ann<-read_met_annotation(features_all);annot<-function(z)z|>mutate(trait=term)|>left_join(ann,by="trait")|>mutate(label=coalesce(label,term),group=coalesce(group,"Other"))
    p1a<-plot_met_circle(annot(assoc_basic),paste0("a. Incident ",Y," MWAS — basic"),label_n=40)
    p1b<-plot_met_circle(annot(prevalent_basic),paste0("b. Prevalent ",Y," MWAS — basic logistic"),label_n=40)
    p1c<-plot_met_circle(annot(assoc_adj2),paste0("c. Incident ",Y," MWAS — adj2"),label_n=40)
    p1d<-plot_met_circle(annot(prevalent_adj2),paste0("d. Prevalent ",Y," MWAS — adj2 logistic"),label_n=40)
    save_plot((p1a|p1c)/(p1b|p1d),"c1.Fig1.circular.png",24,16,outdir=outdir)
    p2a<-plot_volcano(assoc_basic,"beta","term","a. Incident MWAS — basic",.05/max(1,nrow(assoc_basic)),TOP_N)
    p2b<-plot_volcano(prevalent_basic,"beta","term","b. Prevalent MWAS — basic logistic",.05/max(1,nrow(prevalent_basic)),TOP_N)
    p2c<-plot_volcano(assoc_adj2,"beta","term","c. Incident MWAS — adj2",.05/max(1,nrow(assoc_adj2)),TOP_N)
    p2d<-plot_volcano(prevalent_adj2,"beta","term","d. Prevalent MWAS — adj2 logistic",.05/max(1,nrow(prevalent_adj2)),TOP_N)
    save_plot((p2a|p2c)/(p2b|p2d),"c1.Fig2.vc.png",19,14,outdir=outdir)
  }

  # Fig3: same participants in basic and adj2 residualized YY panels.
  yy_need<-unique(c(bvar,top6,covs_adj2));yy_dat<-dat[complete.cases(dat[,intersect(yy_need,names(dat)),drop=FALSE]) & is.finite(dat[[bvar]]),,drop=FALSE]
  p3a<-make_yy_panel(yy_dat,top6,bvar,covs_basic,"a. Yin–Yang protein trajectories — basic adjustment",
                     ifelse(layer=="protein","Mean protein z-score","Mean metabolite z-score"))
  p3b<-make_yy_panel(yy_dat,top6,bvar,covs_adj2,"b. Yin–Yang trajectories — basic + LE8 adjustment",
                     ifelse(layer=="protein","Mean protein z-score","Mean metabolite z-score"))
  save_plot(p3a/p3b+plot_layout(heights=c(1,1)),"c1.Fig3.yy_top.png",16,13.5,outdir=outdir)

  # Fig4: old publication-style pre-diagnostic dot matrix. Dots are independent
  # bin summaries; there are deliberately no connecting lines.
  grad4<-make_backup_gradient_style(dat,top_inc,bvar,layer,k=3L)
  p4<-if(is.null(grad4))blank_plot(paste0("Pre-diagnostic ",layer," gradients"))else grad4$bubble
  save_plot(p4,"c1.Fig4.gradient.png",10.8,8.2,outdir=outdir)

  # Fig5a: exactly one row containing three clusters of the top 100 proteins,
  # reproducing the validated backup figure structure.
  grad5<-make_backup_gradient_style(dat,cluster_features,bvar,layer,k=3L)
  cl<-if(is.null(grad5))make_cluster_figure(dat,cluster_features,bvar,CLUSTER_STEP,YY_MAX_YEAR,"prevalent","a",fixed_k=3L) else
    list(cluster_plot=grad5$clusters,cluster=grad5$cluster,metrics=tibble())
  save_plot(cl$cluster_plot,"c1.Fig5.gradient_cluster.png",15,5.8,outdir=outdir)
  write_raw_csv(cl$cluster,"c1.cluster_membership.csv",rawdir);write_raw_csv(cl$metrics,"c1.cluster_selection.csv",rawdir)

  # Fig6: adjusted Q5-vs-Q1 HRs on the same adj2 complete-case sample; no subtitle.
  qplots<-plot_quantile_top(dat,top6,tvar,evar,covs_adj2,paste0("Adjusted for ",paste(covs_adj2,collapse=", ")))
  save_plot(wrap_plots(qplots,ncol=3),"c1.Fig6.quantile_top.png",16,10,outdir=outdir)

  # Fig7: functional coherence of Bonferroni-significant incident associations.
  enrich<-run_functional_enrichment(assoc,layer,features_all);n_sig<-sum(is.finite(assoc$p.value)&assoc$p.value*nrow(assoc)<.05)
  write_raw_csv(enrich,"c1.enrichment_sig.csv",rawdir);save_plot(plot_functional_enrichment(enrich,n_sig,assoc),"c1.Fig7.enrich_sig.png",16,9,outdir=outdir)

  out<-list(meta=module_meta(layer,extra=list(N=nrow(dat),events=sum(dat[[evar]]==1,na.rm=TRUE),covs_use=covs_use_name)),
            association=assoc,association_basic=assoc_basic,association_adj2=assoc_adj2,association_LE8=assoc_adj2,
            prevalent=assoc_prevalent,prevalent_basic=prevalent_basic,prevalent_adj2=prevalent_adj2,
            attenuation=attenuation_sameN,enrichment=enrich,top_features=top,clusters=cl$cluster,cluster_selection=cl$metrics)
  if(layer=="protein"){out$pwas_incident<-assoc;out$pwas_prevalent<-assoc_prevalent;out$top_proteins<-top}
  else {out$MWAS<-assoc;out$MWAS_prevalent<-assoc_prevalent}
  write_xlsx2(list(cohort=cohort,association=assoc,prevalent=assoc_prevalent,incident_basic=assoc_basic,incident_adj2=assoc_adj2,
                   prevalent_basic=prevalent_basic,prevalent_adj2=prevalent_adj2,attenuation_sameN=attenuation_sameN,
                   cluster_membership=cl$cluster,cluster_selection=cl$metrics,enrichment=enrich),"c1.out.xlsx")
  finalize_outputs(LE8_JOB,outdir);saveRDS(out,selected_cache,compress="xz");out
}

if(prot_DO){invisible(run_c1_layer("protein"));gc(full=TRUE)}
if(met_DO){invisible(run_c1_layer("metabolite"));gc(full=TRUE)}
