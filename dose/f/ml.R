#!/usr/bin/env Rscript
# ml: cox, plot
# Launch with dose.sh; source this file to register definitions for reuse.
if (!exists("dose_step_code", mode="function")) {
    source(file.path(Sys.getenv("DOSE_SCRIPT_DIR", if (.Platform$OS.type == "windows")
        "D:/scripts/dose/f" else "/mnt/d/scripts/dose/f"), "0.common.R"), local=TRUE)
}

# ===========================================================================
# ml / cox
# ===========================================================================
dose_register_step("ml", "cox", expression({
# ======================================================================
# ml.R
# MYTH-like cumulative Cox path for the TotalGain-ranked food groups.
#
# ML residualized variables are used only for ranking.
# Cumulative scores are rebuilt from RAW food-group intake categories:
# proper/lowest-risk category = 1, otherwise = 0.
# M2 primary model; no alcohol; no M3; no strata(region).
# ======================================================================

pacman::p_load(dplyr, tidyr, survival, broom, ggplot2, patchwork)

DIR_PROJECT <- dose_output_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_FOOD <- file.path(DIR_PROJECT, "ukb/02foodwas")
DIR_GROUP <- file.path(DIR_FOOD, "group")
DIR_ML <- file.path(DIR_PROJECT, "ukb/03ml")
DIR_FIG <- file.path(DIR_ML, "figures")
dir.create(DIR_FIG, recursive=TRUE, showWarnings=FALSE)

cfg <- readRDS(file.path(DIR_DATA, "analysis_config.rds"))
dat <- readRDS(file.path(DIR_FOOD, "dat.group.rds")) %>% mutate(eid=as.character(eid))
imp <- read.csv(file.path(DIR_ML, "Food_Importance.csv"), check.names=FALSE) %>% arrange(desc(TotalGain_cv))
proper <- read.csv(file.path(DIR_GROUP, "Group_proper_reference.csv"), check.names=FALSE)

imp_lst <- as.character(imp$Food); P <- length(imp_lst)
proper_use <- proper %>% filter(exposure %in% imp_lst) %>% distinct(exposure,.keep_all=TRUE) %>% select(exposure,best_ref)
if(nrow(proper_use)!=P) stop("Not every ranked food has one proper reference.")
proper_use <- tibble(Food=imp_lst,Rank=seq_along(imp_lst),TotalGain_cv=imp$TotalGain_cv) %>%
	left_join(rename(proper_use,Food=exposure),by="Food")

func_quin <- function(x){
	x <- as.numeric(x); ok <- is.finite(x); out <- rep(NA_integer_,length(x))
	if(sum(ok)<2) return(factor(out))
	br0 <- quantile(x[ok],probs=seq(0,1,1/5),na.rm=TRUE,names=FALSE); bru <- unique(br0)
	if(length(bru)==2) out[ok] <- ifelse(x[ok]>0,2L,1L) else {
		br <- if(length(bru)<6) c(min(x[ok])-1,bru) else bru
		out[ok] <- as.integer(cut(x[ok],breaks=br,include.lowest=TRUE,labels=FALSE))
	}
	lv <- sort(unique(out[!is.na(out)])); out <- match(out,lv)
	factor(out,levels=seq_along(lv),labels=paste0("C",seq_along(lv)),ordered=FALSE)
}

food_cat <- dat %>% select(eid); qc <- list()
for(food in imp_lst){
	food_cat[[food]] <- func_quin(dat[[food]])
	ref <- proper_use$best_ref[match(food,proper_use$Food)]; lv <- levels(food_cat[[food]])
	if(is.na(ref)||!ref %in% lv) stop("Proper category not found for: ",food)
	qc[[food]] <- tibble(Food=food,best_ref=ref,n_cat=length(lv),levels=paste(lv,collapse=";"),
		n_best=sum(as.character(food_cat[[food]])==ref,na.rm=TRUE),
		pct_best=mean(as.character(food_cat[[food]])==ref,na.rm=TRUE)*100,
		n_missing=sum(is.na(food_cat[[food]])))
}
qc <- bind_rows(qc)
if(any(qc$n_missing>0)) stop("Food category missingness found; fixed-cohort score path required.")
write.csv(proper_use %>% left_join(qc,by=c("Food","best_ref")),file.path(DIR_ML,"Imp_food_proper_reference.csv"),row.names=FALSE)

diet_bin <- dat %>% select(eid)
for(food in imp_lst){
	ref <- proper_use$best_ref[match(food,proper_use$Food)]
	diet_bin[[food]] <- as.integer(as.character(food_cat[[food]])==ref)
}
dose_save_rds(diet_bin,file.path(DIR_ML,"Imp_food_binary_components.rds"))

score_path <- dat %>% select(eid)
for(k in seq_len(P)) score_path[[paste0("Imp",k)]] <- rowSums(diet_bin[,imp_lst[seq_len(k)],drop=FALSE])
std_z <- function(x){s <- sd(x,na.rm=TRUE); if(!is.finite(s)||s==0) rep(NA_real_,length(x)) else (x-mean(x,na.rm=TRUE))/s}
for(k in seq_len(P)) score_path[[paste0("Imp",k,"_z")]] <- std_z(score_path[[paste0("Imp",k)]])
dose_save_rds(score_path,file.path(DIR_ML,"Imp_score_path.rds"))

# Prespecified primary development outcome; do not switch based on post hoc results.
Y <- "cvd_stroke_i"; y <- paste0(Y,".Yt2e"); t <- paste0(Y,".t2e")
covs.m2 <- unique(c(cfg$covs.m1,cfg$covs.m2.add)); stopifnot(!"alcohol.c" %in% covs.m2)
d <- dat %>% select(eid,all_of(c(t,y,covs.m2))) %>% left_join(score_path,by="eid")

fit_one <- function(k){
	raw <- paste0("Imp",k); z <- paste0(raw,"_z")
	fz <- coxph(reformulate(c(z,covs.m2),response=paste0("Surv(",t,",",y,")")),data=d)
	tz <- broom::tidy(fz,exponentiate=TRUE,conf.int=TRUE) %>% filter(term==z)
	fr <- coxph(reformulate(c(raw,covs.m2),response=paste0("Surv(",t,",",y,")")),data=d)
	tr <- broom::tidy(fr,exponentiate=TRUE,conf.int=TRUE) %>% filter(term==raw)
	ph <- tryCatch(cox.zph(fz)$table[z,"p"],error=function(e) NA_real_); x <- d[[raw]]
	tibble(k=k,Food_added=imp_lst[k],Foods=paste(imp_lst[seq_len(k)],collapse=";"),N=fz$n,events=fz$nevent,
		score_min=min(x),score_max=max(x),score_mean=mean(x),score_sd=sd(x),
		HR_perSD=tz$estimate,LCI_perSD=tz$conf.low,HCI_perSD=tz$conf.high,P_perSD=tz$p.value,
		neglog10P=-log10(pmax(tz$p.value,1e-300)),PH_P=ph,
		HR_perPoint=tr$estimate,LCI_perPoint=tr$conf.low,HCI_perPoint=tr$conf.high,P_perPoint=tr$p.value)
}
res <- purrr::map_dfr(seq_len(P),fit_one)
write.csv(res,file.path(DIR_ML,"Cox_Imp_food_pattern.csv"),row.names=FALSE)

sfs <- read.csv(file.path(DIR_ML,"SequentialSelection_TotalGain.csv"),check.names=FALSE)
combined <- left_join(res,sfs,by="k",suffix=c("_cox","_ml"))
write.csv(combined,file.path(DIR_ML,"Selection_path_combined.csv"),row.names=FALSE)
sel <- read.csv(file.path(DIR_ML,"Selection_Summary.csv")); k_primary <- sel$k[sel$criterion=="DeltaLogLoss005"][1]

p1 <- ggplot(combined,aes(k,loss_cv_mean))+geom_line()+geom_point()+geom_vline(xintercept=k_primary,linetype=2)+
	labs(x="Number of food components (k)",y="CV log loss",title="LightGBM prefix-selection path")+theme_bw()
p2 <- ggplot(res,aes(k,HR_perSD))+geom_hline(yintercept=1,linetype=2)+geom_errorbar(aes(ymin=LCI_perSD,ymax=HCI_perSD),width=.15)+
	geom_line()+geom_point()+geom_vline(xintercept=k_primary,linetype=2)+labs(x="Number of food components (k)",y="HR per 1-SD higher cumulative score",title="Cox association path")+theme_bw()
p3 <- ggplot(res,aes(k,neglog10P))+geom_hline(yintercept=-log10(.05),linetype=2)+geom_line()+geom_point()+geom_vline(xintercept=k_primary,linetype=2)+
	labs(x="Number of food components (k)",y=expression(-log[10](italic(P))),title="Statistical-significance path")+theme_bw()
fig <- p1/(p2|p3)+plot_annotation(tag_levels="A")
ggsave(file.path(DIR_FIG,"Imp_Cox_path.pdf"),fig,width=12,height=8)

bestp <- res %>% slice_min(P_perSD,n=1,with_ties=FALSE)
cat("MYTH-like cumulative Cox analysis complete\n")
cat(" - ranked food groups:",P,"\n")
cat(" - primary Delta<0.005 k:",k_primary,"\n")
cat(" - cohort:",nrow(d),"\n")
cat(" - events:",sum(d[[y]]==1),"\n")
cat(" - best P at k =",bestp$k,"; HR/SD =",sprintf("%.3f",bestp$HR_perSD),"; P =",format(bestp$P_perSD,scientific=TRUE,digits=3),"\n")
cat(" - Cox path:",file.path(DIR_ML,"Cox_Imp_food_pattern.csv"),"\n")
}))

# ===========================================================================
# ml / plot
# ===========================================================================
dose_register_step("ml", "plot", expression({
# ======================================================================
# ml.R
# Final DOSE machine-learning figure
#
# Visual style adapted from the user's 02_ML.plot.R.
#
# Panel A+B:
#   LightGBM TotalGain importance + sequential CV log-loss curve
#   + cross-method feature-selection tracking
#
# Additional MYTH-like cumulative association panel:
#   bars = -log10(P) for the cumulative score after each added food
#   line/ribbon = HR and 95% CI per 1-SD higher cumulative score
#   primary cutoff = TotalGain ranking + first Delta log-loss < 0.005
#
# Panel C:
#   Pooled out-of-fold ROC curves under LightGBM for the four feature sets
#
# Panel D (inset):
#   paired geographic-fold AUC comparison:
#   all candidate foods vs primary Delta-log-loss-selected set
#
# The primary k is always read from Selection_Summary.csv.
# No k is hard-coded.
# ======================================================================

pacman::p_load(
	data.table, dplyr, tidyr, ggplot2, patchwork,
	pROC, forcats, scales
)

DIR_PROJECT <- dose_output_dir()
DIR_ML <- file.path(DIR_PROJECT, "ukb/03ml")
DIR_FIG <- file.path(DIR_ML, "figures")
dir.create(DIR_FIG, recursive=TRUE, showWarnings=FALSE)

rd <- function(f){
	p <- file.path(DIR_ML, f)
	if(!file.exists(p)) stop("Missing ML result file: ", p)
	as.data.frame(data.table::fread(p))
}

# ======================================================================
# 1. Inputs
# ======================================================================

sel <- rd("Selection_Summary.csv")
imp <- rd("Food_Importance.csv")
sfs <- rd("SequentialSelection_TotalGain.csv")
cox_path <- rd("Cox_Imp_food_pattern.csv")

xgb_imp <- rd("XGBoost_Importance.csv")
rf_imp <- rd("RandomForest_Importance.csv")
shap_rk <- rd("SHAP_Summary.csv")
rfe_rk <- rd("RFE_Ranks.csv")

val_fold <- rd("SelectionValidation_ByFold.csv")
val_pred <- rd("SelectionValidation_Predictions.csv")

# Lightweight post-ML evaluation outputs from ml.py
f_boot <- file.path(DIR_ML, "Bootstrap_Performance_Summary.csv")
f_cm <- file.path(DIR_ML, "Confusion_Matrix_Primary.csv")
f_featcmp <- file.path(DIR_ML, "Paired_FeatureSet_Comparison_LightGBM.csv")
f_modelcmp <- file.path(DIR_ML, "Paired_Model_Comparison_Primary.csv")

primary <- sel %>%
	filter(criterion == "DeltaLogLoss005") %>%
	slice(1)

if(!nrow(primary)) stop("DeltaLogLoss005 primary selection is missing.")

k <- as.integer(primary$k)
stop_k <- as.integer(primary$stop_at_k)

food_lab <- function(x){
	x <- trimws(gsub("_", " ", x))
	x <- gsub("Apples Pears", "Apples/Pears", x, fixed=TRUE)
	x <- gsub("High fat dairy", "High-fat dairy", x, fixed=TRUE)
	x <- gsub("Low fat dairy", "Low-fat dairy", x, fixed=TRUE)
	x <- gsub("Green leafy vegetables", "Green leafy vegetables", x, fixed=TRUE)
	x
}

thm <- theme_bw(base_size=10) +
	theme(
		axis.title=element_text(size=11, face="bold"),
		axis.text=element_text(size=9),
		legend.title=element_blank(),
		panel.grid.minor=element_blank(),
		panel.grid.major.x=element_blank(),
		plot.margin=margin(10,10,10,10)
	)

# Colours retained from 02_ML.plot.R
C <- list(
	hi="#2F5D8A",
	lo="#D9E1EA",
	line="#3B668F",
	rib="#B8C9D9",
	cut="#7B8FA6",
	ora="#D55E00",
	dot="#4A7EBB",
	bg=scales::alpha("#EAF1F7", 0.35)
)

col_frame <- "grey35"
lw_frame <- 0.55
left_margin_size <- 40

# ======================================================================
# 2. Panel A+B
#    Importance + sequential loss + 10-method selection tracking
# ======================================================================

# Merge the TotalGain ranking with the prefix-selection curve.
dfAB <- imp %>%
	arrange(desc(TotalGain_cv)) %>%
	transmute(
		Food,
		Importance=TotalGain_cv,
		food_idx=row_number()
	) %>%
	left_join(
		sfs %>%
			select(
				Food=Food_added,
				loss_cv_mean,
				loss_cv_sd,
				loss_cv_se,
				Delta_loss,
				auc_cv_mean,
				auprc_cv_mean
			),
		by="Food"
	) %>%
	mutate(
		F_lab=factor(food_lab(Food), levels=food_lab(Food)),
		topk=food_idx <= k,
		loss_lower=loss_cv_mean-loss_cv_sd,
		loss_upper=loss_cv_mean+loss_cv_sd
	)

if(anyNA(dfAB$loss_cv_mean))
	stop("Food_Importance and SequentialSelection_TotalGain do not align.")

# ----------------------------------------------------------------------
# 2a. Define the feature sets shown in the method-tracking panel.
#
# There are 10 validation methods:
#   3 model-based: LGB / XGB / RF
#   3 SHAP-based : LGB / XGB / RF
#   4 RFE-based  : LGB / XGB / RF / Logistic
#
# LGB model-based uses the actual primary Delta-log-loss selected set (k).
# All alternative methods use the primary selection size k;
# RFE uses its actual Selected==TRUE indicators.
# ----------------------------------------------------------------------

top_n <- function(x, n){
	x <- unique(na.omit(as.character(x)))
	x[seq_len(min(n, length(x)))]
}

top_by_rank <- function(dat, rank_col, n){
	dat %>%
		arrange(.data[[rank_col]], Food) %>%
		pull(Food) %>%
		top_n(n)
}

shap_top <- function(method, n=k){
	shap_rk %>%
		filter(Method == method) %>%
		arrange(Rank, Food) %>%
		pull(Food) %>%
		top_n(n)
}

rfe_selected <- function(method){
	rfe_rk %>%
		filter(Method == method, Selected) %>%
		arrange(RFE_rank, Food) %>%
		pull(Food) %>%
		unique() -> selected
	if(length(selected) != k)
		stop("RFE selection size differs from primary k; rerun ml.py.")
	selected
}

fs <- list(
	LGB_MB  = dfAB$Food[dfAB$food_idx <= k],
	XGB_MB  = top_by_rank(xgb_imp, "Rank", k),
	RF_MB   = top_by_rank(rf_imp, "Rank", k),

	LGB_SHAP = shap_top("LightGBM", k),
	XGB_SHAP = shap_top("XGBoost", k),
	RF_SHAP  = shap_top("RandomForest", k),

	RFE_LGB = rfe_selected("LightGBM"),
	RFE_XGB = rfe_selected("XGBoost"),
	RFE_RF  = rfe_selected("RandomForest"),
	RFE_GLM = rfe_selected("Logistic")
)

info <- data.frame(
	Method=c(
		"LGB_MB","XGB_MB","RF_MB",
		"LGB_SHAP","XGB_SHAP","RF_SHAP",
		"RFE_LGB","RFE_XGB","RFE_RF","RFE_GLM"
	),
	lab=c(
		"LightGBM","XGBoost","Random forest",
		"LightGBM","XGBoost","Random forest",
		"LightGBM","XGBoost","Random forest",
		"Logistic regression"
	),
	grp=c(
		rep("Model-based approach",3),
		rep("SHAP value-based approach",3),
		rep("RFE-based approach",4)
	),
	y=c(12,11,10, 7.5,6.5,5.5, 3,2,1,0),
	stringsAsFactors=FALSE
)

grp_header_df <- info %>%
	group_by(grp) %>%
	summarise(y=max(y)+0.8, .groups="drop")

sel_long <- bind_rows(lapply(names(fs), function(nm){
	data.frame(
		Method=nm,
		Food=fs[[nm]],
		sel=1L,
		stringsAsFactors=FALSE
	)
}))

track_df <- tidyr::expand_grid(
	Method=info$Method,
	Food=dfAB$Food
) %>%
	left_join(sel_long, by=c("Method","Food")) %>%
	left_join(info, by="Method") %>%
	mutate(
		sel=tidyr::replace_na(sel,0L),
		F_lab=factor(food_lab(Food),levels=levels(dfAB$F_lab)),
		primary=match(Food,dfAB$Food) <= k
	)

# ----------------------------------------------------------------------
# 2b. Top: TotalGain bars + transformed CV log-loss curve
# ----------------------------------------------------------------------

n_food <- nlevels(dfAB$F_lab)
imp_max <- max(dfAB$Importance,na.rm=TRUE)*1.05

loss_range <- max(dfAB$loss_upper,na.rm=TRUE)-min(dfAB$loss_lower,na.rm=TRUE)
if(!is.finite(loss_range) || loss_range<=0) stop("Invalid loss range.")

a <- imp_max/loss_range
b <- -a*min(dfAB$loss_lower,na.rm=TRUE)

ptAB <- dfAB[dfAB$food_idx==k,,drop=FALSE]

pAB_top <- ggplot(dfAB,aes(food_idx)) +
	annotate(
		"rect", xmin=.5, xmax=k+.5, ymin=-Inf, ymax=Inf,
		fill=C$bg
	) +
	geom_col(aes(y=Importance,fill=topk),width=.64) +
	geom_ribbon(
		aes(
			x=food_idx,
			ymin=a*loss_lower+b,
			ymax=a*loss_upper+b,
			group=1
		),
		fill=C$rib,alpha=.18,inherit.aes=FALSE
	) +
	geom_line(
		aes(x=food_idx,y=a*loss_cv_mean+b,group=1),
		colour=C$line,linewidth=1,inherit.aes=FALSE
	) +
	geom_point(
		aes(x=food_idx,y=a*loss_cv_mean+b,colour=topk),
		size=2.2,inherit.aes=FALSE
	) +
	geom_vline(
		xintercept=k+.5,
		linetype="22",colour=C$cut,linewidth=.7
	) +
	geom_vline(
		xintercept=c(.5,n_food+.5),
		colour=col_frame,linewidth=lw_frame
	) +
	annotate(
		"text",
		x=n_food+.2,
		y=imp_max*.98,
		label=paste0(
			"n = ",k,
			"; first Delta loss < 0.005 at k = ",stop_k
		),
		colour=C$ora,fontface="bold",size=3.4,hjust=1
	) +
	scale_fill_manual(
		values=c(`TRUE`=C$hi,`FALSE`=C$lo),
		guide="none"
	) +
	scale_colour_manual(
		values=c(`TRUE`=C$ora,`FALSE`=C$line),
		guide="none"
	) +
	scale_x_continuous(limits=c(.5,n_food+.5),breaks=seq_len(n_food),
		labels=food_lab(dfAB$Food),expand=expansion(mult=0)) +
	scale_y_continuous(
		name="Food importance",
		expand=expansion(c(0,.05)),
		sec.axis=sec_axis(~(. - b)/a,name="CV log loss")
	) +
	thm +
	theme(
		panel.border=element_blank(),
		axis.line.x.top=element_line(colour=col_frame,linewidth=lw_frame),
		axis.line.x.bottom=element_blank(),
		axis.line.y=element_blank(),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.title.y=element_text(size=12,face="bold"),
		plot.margin=margin(6,10,-1,left_margin_size)
	)

# ----------------------------------------------------------------------
# 2c. Bottom: 10-method selection track
# ----------------------------------------------------------------------

count_df <- track_df %>%
	group_by(Food,F_lab) %>%
	summarise(SelectionCount=sum(sel),.groups="drop") %>%
	mutate(y=-1.5)

stopifnot(length(fs)==10L, !"RFE_SVM" %in% names(fs),
          all(count_df$SelectionCount <= length(fs)),
          sum(count_df$SelectionCount)==sum(lengths(fs)))
write.csv(count_df %>% transmute(Food,SelectionCount,MethodsCompared=length(fs)),
          file.path(DIR_FIG,"Selection_Count_displayed_methods.csv"),row.names=FALSE)

line_df <- track_df %>%
	filter(sel==1) %>%
	group_by(F_lab) %>%
	summarise(
		y_min=min(y),
		y_max=max(y),
		count=n(),
		.groups="drop"
	) %>%
	filter(count>1)

pAB_track <- ggplot(track_df,aes(as.integer(F_lab),y)) +
	annotate(
		"rect", xmin=.5, xmax=k+.5, ymin=-Inf, ymax=Inf,
		fill=C$bg
	) +
	geom_point(colour="grey92",size=1.5) +
	geom_segment(
		data=line_df,
		aes(x=as.integer(F_lab),xend=as.integer(F_lab),y=y_min,yend=y_max),
		colour=C$dot,linewidth=.7
	) +
	geom_point(
		data=filter(track_df,sel==1),
		shape=21,size=3,fill=C$dot,colour="white",stroke=.4
	) +
	geom_text(
		data=count_df,
		aes(x=as.integer(F_lab),y=-1.5,label=SelectionCount),
		size=3.5,fontface="bold",colour=C$hi
	) +
	annotate(
		"text",x=.5,y=-1.5,label="Selection count",
		hjust=1,fontface="bold",size=3.2,colour="grey35"
	) +
	geom_hline(
		yintercept=c(8.75,4.25),
		colour="grey88",linewidth=.4
	) +
	geom_vline(
		xintercept=c(.5,n_food+.5),
		colour=col_frame,linewidth=lw_frame
	) +
	geom_text(
		data=grp_header_df,
		aes(x=.5,y=y,label=grp),
		hjust=1,fontface="bold",size=3.2,
		colour="grey35",inherit.aes=FALSE
	) +
	scale_x_continuous(limits=c(.5,n_food+.5),breaks=seq_len(n_food),
		labels=food_lab(dfAB$Food),expand=expansion(mult=0)) +
	scale_y_continuous(
		breaks=c(-1.5,info$y),
		labels=c("",info$lab),
		expand=expansion(add=c(1.0,1.2))
	) +
	coord_cartesian(clip="off") +
	theme_minimal(base_size=8.8) +
	theme(
		panel.grid=element_blank(),
		axis.line.x.bottom=element_line(
			colour=col_frame,linewidth=lw_frame
		),
		axis.text.x=element_text(
			angle=90,hjust=1,vjust=.5,colour="black",
			size=9,face="bold"
		),
		axis.text.y=element_text(size=8.8,colour="black"),
		plot.margin=margin(-1,10,8,left_margin_size)
	)

pAB <- pAB_top/pAB_track +
	plot_layout(heights=c(2.2,1.8)) &
	labs(x=NULL,y=NULL)

ggsave(
	file.path(DIR_FIG,"Fig_ML_pAB.pdf"),
	plot=pAB,width=8.5,height=7.2,device=cairo_pdf
)

# ======================================================================
# 3. MYTH-like cumulative association path
#
# Similar information content to the reference figure:
#   bars: -log10(P) for Imp1 ... ImpP
#   line: HR per 1-SD higher cumulative score
#   ribbon: 95% CI
#
# Styling follows Fig_ML_pAB rather than the rainbow/reference palette.
# The primary selected set is the first k foods defined by:
#   TotalGain ranking + first Delta log-loss < 0.005.
# ======================================================================

cox_plot_df <- cox_path %>%
	arrange(k) %>%
	mutate(
		food_idx=row_number(),
		F_lab=factor(food_lab(Food_added), levels=food_lab(Food_added)),
		topk=k <= !!k,
		logP=-log10(pmax(P_perSD, 1e-300))
	)

if(nrow(cox_plot_df) != nrow(dfAB))
	stop("Cox path and ML ranking have different numbers of foods.")

if(!all(cox_plot_df$Food_added == dfAB$Food))
	stop("Cox cumulative path is not aligned with the TotalGain ranking.")

# Transform the HR scale onto the -log10(P) plotting scale.
# Use the observed confidence interval range with modest padding so the ribbon
# remains visible and the inverse secondary axis is exact.
hr_lo <- min(cox_plot_df$LCI_perSD, na.rm=TRUE)
hr_hi <- max(cox_plot_df$HCI_perSD, na.rm=TRUE)
hr_pad <- max((hr_hi-hr_lo)*.08, .005)
hr_min <- hr_lo-hr_pad
hr_max <- hr_hi+hr_pad

p_max <- max(cox_plot_df$logP, na.rm=TRUE)*1.08
if(!is.finite(p_max) || p_max<=0) p_max <- 1
if(!is.finite(hr_min) || !is.finite(hr_max) || hr_max<=hr_min)
	stop("Invalid HR range in Cox_Imp_food_pattern.csv.")

hr_to_y <- function(x) (x-hr_min)/(hr_max-hr_min)*p_max
y_to_hr <- function(y) hr_min + y/p_max*(hr_max-hr_min)

# Split the HR curve at the prespecified primary cutoff so that the selected
# prefix is emphasized while later additions remain visible as sensitivity.
cox_sel <- cox_plot_df %>% filter(k <= !!k)
cox_post <- cox_plot_df %>% filter(k >= !!k)

pCox <- ggplot(cox_plot_df, aes(F_lab)) +
	annotate(
		"rect", xmin=.5, xmax=k+.5, ymin=-Inf, ymax=Inf,
		fill=C$bg
	) +
	geom_col(
		aes(y=logP, fill=topk),
		width=.64
	) +
	geom_ribbon(
		aes(
			x=food_idx,
			ymin=hr_to_y(LCI_perSD),
			ymax=hr_to_y(HCI_perSD),
			group=1
		),
		fill=C$rib, alpha=.18, inherit.aes=FALSE
	) +
	geom_line(
		data=cox_post,
		aes(x=food_idx, y=hr_to_y(HR_perSD), group=1),
		colour=C$cut, linewidth=.9, inherit.aes=FALSE
	) +
	geom_point(
		data=cox_post,
		aes(x=food_idx, y=hr_to_y(HR_perSD)),
		shape=21, size=2.2, fill="white", colour=C$cut,
		stroke=.7, inherit.aes=FALSE
	) +
	geom_line(
		data=cox_sel,
		aes(x=food_idx, y=hr_to_y(HR_perSD), group=1),
		colour=C$line, linewidth=1.05, inherit.aes=FALSE
	) +
	geom_point(
		data=cox_sel,
		aes(x=food_idx, y=hr_to_y(HR_perSD)),
		shape=21, size=2.4, fill=C$ora, colour=C$line,
		stroke=.65, inherit.aes=FALSE
	) +
	geom_vline(
		xintercept=k+.5,
		linetype="22", colour=C$cut, linewidth=.7
	) +
	geom_vline(
		xintercept=c(.5,n_food+.5),
		colour=col_frame, linewidth=lw_frame
	) +
	annotate(
		"text",
		x=k+.58,
		y=p_max*.97,
		label=paste0(
			"Primary n = ", k,
			"\nTotalGain + Delta log loss < 0.005",
			"\nstop before ", food_lab(primary$stop_food)
		),
		hjust=0, vjust=1,
		colour=C$ora, fontface="bold", size=3.3
	) +
	scale_fill_manual(
		values=c(`TRUE`=C$hi, `FALSE`=C$lo),
		guide="none"
	) +
	scale_y_continuous(
		name=expression(-log[10](italic(P))),
		limits=c(0,p_max),
		expand=expansion(mult=c(0,.02)),
		sec.axis=sec_axis(
			~ y_to_hr(.),
			name="HR per 1-SD higher cumulative score"
		)
	) +
	labs(x=NULL) +
	thm +
	theme(
		panel.border=element_blank(),
		axis.line.x=element_line(colour=col_frame, linewidth=lw_frame),
		axis.line.y.left=element_line(colour=col_frame, linewidth=lw_frame),
		axis.line.y.right=element_line(colour=col_frame, linewidth=lw_frame),
		axis.text.x=element_text(
			angle=45, hjust=1, colour="black", size=9.5, face="bold"
		),
		axis.title.y.left=element_text(size=11, face="bold"),
		axis.title.y.right=element_text(size=11, face="bold"),
		plot.margin=margin(8,12,12,12)
	)

ggsave(
	file.path(DIR_FIG,"Fig_ML_Cox_path.pdf"),
	pCox,width=11.5,height=5.8,device=cairo_pdf
)

# Mirror layout: ML extends upward and Cox evidence extends downward.
# Both panels use exactly the same food coordinates and selection boundary.
# Negative drawing coordinates are only a layout device: Cox P values and
# hazard ratios are displayed on their original scales via exact inverses.
mirror_x <- function(){
    scale_x_continuous(limits=c(.5,n_food+.5), breaks=seq_len(n_food),
                       expand=expansion(mult=0))
}
mirror_theme <- thm + theme(
    panel.border=element_blank(), panel.grid.major.x=element_blank(),
    axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.title.x=element_blank(),
    axis.line.y=element_line(colour=col_frame, linewidth=.45),
    plot.title=element_text(face="bold", size=12),
    plot.subtitle=element_text(size=9, colour="grey35"),
    plot.margin=margin(6,12,6,12)
)

pMirrorML <- ggplot(dfAB, aes(food_idx)) +
    annotate("rect", xmin=.5,xmax=k+.5,ymin=-Inf,ymax=Inf,fill=C$bg) +
    geom_col(aes(y=Importance,fill=topk),width=.64) +
    geom_ribbon(aes(ymin=a*loss_lower+b,ymax=a*loss_upper+b),
                fill=C$rib,alpha=.3) +
    geom_line(aes(y=a*loss_cv_mean+b),colour=C$line,linewidth=1) +
    geom_point(aes(y=a*loss_cv_mean+b,colour=topk),size=2.3) +
    geom_vline(xintercept=k+.5,linetype="22",colour=C$cut,linewidth=.7) +
    geom_hline(yintercept=0,colour=col_frame,linewidth=.5) +
    scale_fill_manual(values=c(`TRUE`=C$hi,`FALSE`=C$lo),guide="none") +
    scale_colour_manual(values=c(`TRUE`=C$ora,`FALSE`=C$line),guide="none") +
    mirror_x() +
    scale_y_continuous(name="Food importance (TotalGain)",
        limits=c(0,imp_max*1.04),expand=expansion(mult=0),
        sec.axis=sec_axis(~(. - b)/a,name="CV log loss")) +
    labs(title="A. Machine-learning selection",
        subtitle=paste0("LightGBM | Bars: importance; line and band: mean CV log loss +/- 1 SD | Primary: ",
                        k," foods; first improvement < 0.005 at step ",stop_k)) +
    mirror_theme

# One shared label strip at the center of the mirrored panels.
pMirrorLabels <- ggplot(dfAB,aes(food_idx,.5)) +
    annotate("rect",xmin=.5,xmax=k+.5,ymin=0,ymax=1,fill=C$bg) +
    geom_vline(xintercept=k+.5,linetype="22",colour=C$cut,linewidth=.7) +
    geom_text(aes(label=food_lab(Food)),angle=38,size=3.1,
              fontface="bold",hjust=.5,vjust=.5) +
    mirror_x() + scale_y_continuous(limits=c(0,1),expand=expansion(mult=0)) +
    coord_cartesian(clip="off") + theme_void() +
    theme(plot.margin=margin(0,12,0,12))

pMirrorCox <- ggplot(cox_plot_df,aes(food_idx)) +
    annotate("rect",xmin=.5,xmax=k+.5,ymin=-Inf,ymax=Inf,fill=C$bg) +
    geom_col(aes(y=-logP,fill=topk),width=.64) +
    geom_ribbon(aes(ymin=-hr_to_y(HCI_perSD),ymax=-hr_to_y(LCI_perSD)),
                fill=C$rib,alpha=.3) +
    geom_line(aes(y=-hr_to_y(HR_perSD)),colour=C$line,linewidth=1) +
    geom_point(aes(y=-hr_to_y(HR_perSD),colour=topk),size=2.3) +
    geom_hline(yintercept=0,colour=col_frame,linewidth=.5) +
    geom_vline(xintercept=k+.5,linetype="22",colour=C$cut,linewidth=.7) +
    scale_fill_manual(values=c(`TRUE`=C$hi,`FALSE`=C$lo),guide="none") +
    scale_colour_manual(values=c(`TRUE`=C$ora,`FALSE`=C$line),guide="none") +
    mirror_x() +
    scale_y_continuous(name=expression("Cox association: "*-log[10](italic(P))),
        limits=c(-p_max,0),labels=function(x) format(abs(x),trim=TRUE),
        expand=expansion(mult=0),
        sec.axis=sec_axis(~y_to_hr(-.),name="HR per 1-SD higher cumulative score")) +
    labs(caption="B. Cumulative Cox association | Bars: -log10(P); line and band: HR and 95% CI\nEach position adds that food to all preceding foods; Cox axes increase downward.") +
    mirror_theme + theme(plot.caption=element_text(hjust=0,size=9,colour="grey30"))

# Equal plot widths are enforced by one flat patchwork, including the track.
# Use numeric coordinates here too: annotations must not extend a discrete
# scale and shift the selection dots relative to the mirrored bars.
mirror_track_df <- track_df %>% mutate(food_idx=match(Food,dfAB$Food))
mirror_line_df <- line_df %>% mutate(food_idx=as.integer(F_lab))
mirror_count_df <- count_df %>% mutate(food_idx=match(Food,dfAB$Food))
pMirrorTrack <- ggplot(mirror_track_df,aes(food_idx,y)) +
    annotate("rect",xmin=.5,xmax=k+.5,ymin=-Inf,ymax=Inf,fill=C$bg) +
    geom_point(colour="grey92",size=1.5) +
    geom_segment(data=mirror_line_df,
        aes(x=food_idx,xend=food_idx,y=y_min,yend=y_max),
        colour=C$dot,linewidth=.7) +
    geom_point(data=filter(mirror_track_df,sel==1),
        shape=21,size=3,fill=C$dot,colour="white",stroke=.4) +
    geom_text(data=mirror_count_df,aes(y=-1.5,label=SelectionCount),
        size=3.5,fontface="bold",colour=C$hi) +
    geom_hline(yintercept=c(8.75,4.25),colour="grey88",linewidth=.4) +
    geom_vline(xintercept=k+.5,linetype="22",colour=C$cut,linewidth=.7) +
    geom_text(data=grp_header_df,aes(x=.5,y=y,label=grp),
        hjust=1,fontface="bold",size=3.2,colour="grey35",inherit.aes=FALSE) +
    mirror_x() +
    scale_y_continuous(breaks=c(-1.5,info$y),labels=c("Selection count",info$lab),
        expand=expansion(add=c(1,1.2))) +
    coord_cartesian(clip="off") +
    labs(x=NULL,y=NULL,title="C. Selection agreement across 10 methods") +
    theme_minimal(base_size=8.8) +
    theme(panel.grid=element_blank(),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),axis.text.y=element_text(colour="black"),
          plot.title=element_text(face="bold",size=12),
          plot.margin=margin(12,12,8,40))

# Preserve the original numeric scales despite the mirrored display.
stopifnot(isTRUE(all.equal(y_to_hr(hr_to_y(cox_plot_df$HR_perSD)),
                           cox_plot_df$HR_perSD)))

pMirror <- wrap_plots(pMirrorML,pMirrorLabels,pMirrorCox,
                     ncol=1,heights=c(3,1.25,3))
# Preserve the original pAB order: ML, selection track with food labels,
# then append the Cox mirror below. A flat stack aligns all food columns.
pAB_Cox <- wrap_plots(
    pAB_top + labs(x=NULL),
    pAB_track + labs(x=NULL,y=NULL),
    pMirrorCox + theme(plot.margin=margin(0,10,12,left_margin_size)),
    ncol=1,heights=c(2.5,2.3,2.5)
)
ggsave(file.path(DIR_FIG,"Fig_ML_Cox_mirror.pdf"),pMirror,
       width=14,height=8,device=cairo_pdf)
ggsave(file.path(DIR_FIG,"Fig_ML_pAB_Cox.pdf"),pAB_Cox,
       width=8.5,height=9.5,device=cairo_pdf)

# ======================================================================
# 4. Panel C+D
#    LightGBM ROC across the four validation feature sets
#    + paired fold-AUC inset for all candidates vs primary selected set
# ======================================================================

primary_set_name <- unique(
	val_fold$FeatureSet[
		grepl("^PrimaryDelta005_",val_fold$FeatureSet)
	]
)[1]

if(is.na(primary_set_name))
	stop("PrimaryDelta005 feature set missing from validation output.")

all_set_name <- unique(
	val_fold$FeatureSet[
		grepl("^All[0-9]+$", val_fold$FeatureSet)
	]
)[1]

if(is.na(all_set_name))
	stop("All-candidate feature set is missing from validation output.")

all_n <- as.integer(sub("^All","",all_set_name))
shap_set_name <- paste0("SHAP_Consensus", k)
rfe_set_name <- paste0("RFE_Consensus", k)
if(primary_set_name != paste0("PrimaryDelta005_", k) ||
   !all(c(shap_set_name, rfe_set_name) %in% val_fold$FeatureSet) ||
   !all(c(shap_set_name, rfe_set_name) %in% val_pred$FeatureSet))
	stop("Validation sets do not match primary k; rerun ml.py.")

feature_order <- c(
	all_set_name,
	primary_set_name,
	shap_set_name,
	rfe_set_name
)
feature_order <- feature_order[
	feature_order %in% unique(val_pred$FeatureSet)
]

feature_labs <- setNames(
	c(paste0("SHAP consensus ", k), paste0("RFE consensus ", k)),
	c(shap_set_name, rfe_set_name)
)
feature_labs[all_set_name] <- paste0("All ",all_n)
feature_labs[primary_set_name] <- paste0("Primary selected ",k)

# Pooled OOF ROC under the primary ML algorithm (LightGBM).
roc_dat <- val_pred %>%
	filter(
		Method=="LightGBM",
		FeatureSet %in% feature_order
	) %>%
	group_split(FeatureSet) %>%
	lapply(function(x){
		ro <- pROC::roc(
			x$y_true,x$score,
			quiet=TRUE,direction="<"
		)
		data.frame(
			FPR=1-ro$specificities,
			TPR=ro$sensitivities,
			FeatureSet=x$FeatureSet[1],
			AUC=as.numeric(pROC::auc(ro))
		)
	}) %>%
	bind_rows()

auc_labels <- roc_dat %>%
	group_by(FeatureSet) %>%
	summarise(AUC=first(AUC),.groups="drop")

legend_labels <- setNames(
	paste0(
		feature_labs[auc_labels$FeatureSet],
		" (AUC = ",
		sprintf("%.3f",auc_labels$AUC),
		")"
	),
	auc_labels$FeatureSet
)

roc_cols <- setNames(c("#4A7EBB", "#9AA4AD"), c(shap_set_name, rfe_set_name))
roc_cols[all_set_name] <- "#7B8FA6"
roc_cols[primary_set_name] <- C$ora
roc_cols <- roc_cols[feature_order]

pC <- ggplot(
	roc_dat,
	aes(
		FPR,TPR,
		colour=factor(FeatureSet,levels=feature_order)
	)
) +
	geom_path(linewidth=1.1) +
	geom_abline(
		slope=1,linetype=2,
		colour="grey80",linewidth=.7
	) +
	coord_equal() +
	scale_colour_manual(
		values=roc_cols,
		labels=legend_labels[feature_order]
	) +
	labs(
		x="False positive rate",
		y="True positive rate",
		title="D. LightGBM out-of-fold ROC"
	) +
	thm +
	theme(
		panel.grid=element_blank(),
		legend.position=c(.37,.79),
		legend.background=element_blank(),
		legend.key=element_blank(),
		legend.text=element_text(size=7.7,face="bold"),
		plot.margin=margin(10,20,10,10)
	)

# ----------------------------------------------------------------------
# Inset D: paired geographic-fold AUC, All17 vs Primary selected set.
# ----------------------------------------------------------------------

dfD <- val_fold %>%
	filter(
		Method=="LightGBM",
		FeatureSet %in% c(all_set_name,primary_set_name)
	) %>%
	mutate(
		Block=ifelse(
			FeatureSet==all_set_name,
			paste0("All ",all_n),
			paste0("Primary ",k)
		),
		Block=factor(
			Block,
			levels=c(paste0("All ",all_n),paste0("Primary ",k))
		)
	)

wideD <- dfD %>%
	select(Fold,FeatureSet,AUC) %>%
	pivot_wider(names_from=FeatureSet,values_from=AUC)

if(all(c(all_set_name,primary_set_name) %in% names(wideD))){
	d_auc <- wideD[[primary_set_name]]-wideD[[all_set_name]]
	md <- mean(d_auc,na.rm=TRUE)
	pv <- tryCatch(
		stats::wilcox.test(
			wideD[[primary_set_name]],
			wideD[[all_set_name]],
			paired=TRUE,
			exact=FALSE
		)$p.value,
		error=function(e) NA_real_
	)
} else {
	md <- NA_real_
	pv <- NA_real_
}

ym <- max(dfD$AUC,na.rm=TRUE)

labD <- if(is.finite(pv)){
	sprintf("Delta AUC = %.3f; P = %.3g",md,pv)
} else {
	sprintf("Delta AUC = %.3f",md)
}

pD <- ggplot(dfD,aes(Block,AUC)) +
	geom_line(
		aes(group=Fold),
		colour="grey88",alpha=.55,linewidth=.3
	) +
	geom_boxplot(
		outlier.shape=NA,width=.36,
		fill="grey98",colour="grey40",
		linewidth=.3,alpha=.5
	) +
	geom_point(
		aes(colour=Block),
		position=position_jitter(width=.05,seed=12345),
		size=1.7,alpha=.65
	) +
	stat_summary(
		fun=mean,geom="point",
		shape=23,size=2.2,
		fill=C$ora,stroke=.4
	) +
	annotate(
		"segment",
		x=c(1,2,1),xend=c(1,2,2),
		y=c(ym+.006,ym+.006,ym+.009),
		yend=c(ym+.009,ym+.009,ym+.009),
		linewidth=.3
	) +
	annotate(
		"text",
		x=1.5,y=ym+.017,
		label=labD,
		size=2.2,fontface="bold"
	) +
	scale_colour_manual(
		values=setNames(
			c("#666666",C$ora),
			c(paste0("All ",all_n),paste0("Primary ",k))
		),
		guide="none"
	) +
	labs(x=NULL,y="AUC") +
	theme_bw(base_size=8) +
	theme(
		axis.text=element_text(size=7),
		axis.title.y=element_text(size=8,face="bold"),
		panel.grid=element_blank(),
		aspect.ratio=1.1,
		plot.background=element_rect(
			fill="white",colour="grey85",linewidth=.3
		),
		plot.margin=margin(4,4,4,4)
	)

pCD <- pC +
	inset_element(
		pD,
		left=.50,bottom=.045,right=.985,top=.50
	)

ggsave(
	file.path(DIR_FIG,"Fig_ML_pCD.pdf"),
	pCD,width=5.2,height=5.2,
	device=cairo_pdf
)

# ======================================================================
# 5. Main figure: mirrored ML/Cox selection and predictive validation
# ======================================================================

# Keep the full aligned mirror stack intact next to the predictive ROC panel.
f_main <- wrap_plots(
    wrap_elements(full=patchworkGrob(pAB_Cox)),
    pCD,
    nrow=1, widths=c(8.5,5.0)
)

ggsave(
	file.path(DIR_FIG,"Fig_ML.pdf"),
	f_main,width=13.5,height=9.5,
	device=cairo_pdf
)

# ======================================================================
# 6. Optional supplementary validation figure
# ======================================================================

val_sum <- rd("SelectionValidation_Summary.csv")

pS <- ggplot(
	val_sum,
	aes(
		x=FeatureSet,
		y=AUC_mean,
		shape=Method,
		group=Method
	)
) +
	geom_line(
		position=position_dodge(width=.25),
		colour="grey75",linewidth=.4
	) +
	geom_point(
		position=position_dodge(width=.25),
		size=2.3
	) +
	geom_errorbar(
		aes(
			ymin=AUC_mean-AUC_sd,
			ymax=AUC_mean+AUC_sd
		),
		position=position_dodge(width=.25),
		width=.10,linewidth=.35
	) +
	labs(
		x=NULL,y="Mean geographic-CV AUC",
		title="Predictive validation across feature sets and algorithms"
	) +
	thm +
	theme(
		axis.text.x=element_text(angle=25,hjust=1),
		legend.position="right"
	)

ggsave(
	file.path(DIR_FIG,"Fig_ML_validation.pdf"),
	pS,width=10,height=5.8,
	device=cairo_pdf
)

# ======================================================================
# 7. MYTH-like ML performance figure with bootstrap 95% CIs
#
# A: pooled OOF ROC curves for four core models using the PRIMARY food set
# B: four confusion matrices using fixed OOF Youden cutoffs
# C: performance metrics with bootstrap 95% CIs
#
# This section is run only when ml.py outputs exist.
# ======================================================================

if(all(file.exists(c(f_boot, f_cm)))){
	boot <- as.data.frame(fread(f_boot))
	cm <- as.data.frame(fread(f_cm))

	core_models <- c(
		"LightGBM","XGBoost","RandomForest","Logistic"
	)

	primary_set_eval <- unique(
		boot$FeatureSet[
			grepl("^PrimaryDelta005_", boot$FeatureSet)
		]
	)[1]

	if(!is.na(primary_set_eval)){

		model_labs <- c(
			LightGBM="LightGBM",
			XGBoost="XGBoost",
			RandomForest="Random forest",
			Logistic="Logistic regression"
		)

		model_cols <- c(
			LightGBM=C$hi,
			XGBoost=C$ora,
			RandomForest="#6A994E",
			Logistic="#9AA4AD"
		)

		# --------------------------------------------------------------
		# A. ROC
		# --------------------------------------------------------------
		roc_primary <- val_pred %>%
			filter(
				FeatureSet==primary_set_eval,
				Method %in% core_models
			) %>%
			group_split(Method) %>%
			lapply(function(x){
				ro <- pROC::roc(
					x$y_true, x$score,
					quiet=TRUE, direction="<"
				)
				data.frame(
					FPR=1-ro$specificities,
					TPR=ro$sensitivities,
					Method=x$Method[1],
					AUC=as.numeric(pROC::auc(ro))
				)
			}) %>%
			bind_rows()

		auc_ci <- boot %>%
			filter(
				FeatureSet==primary_set_eval,
				Method %in% core_models,
				Metric=="AUC"
			) %>%
			mutate(
				label=paste0(
					model_labs[Method],
					" (AUC = ",
					sprintf("%.3f",Estimate),
					" [",
					sprintf("%.3f",Lower95),
					"-",
					sprintf("%.3f",Upper95),
					"])"
				)
			)

		roc_labels <- setNames(
			auc_ci$label,
			auc_ci$Method
		)

		pPerfA <- ggplot(
			roc_primary,
			aes(
				FPR,TPR,
				colour=factor(Method,levels=core_models)
			)
		) +
			geom_path(linewidth=1.1) +
			geom_abline(
				slope=1,linetype=2,
				colour="grey75",linewidth=.6
			) +
			coord_equal() +
			scale_colour_manual(
				values=model_cols,
				labels=roc_labels[core_models],
				drop=FALSE
			) +
			labs(
				x="False positive rate",
				y="True positive rate",
				title="A."
			) +
			thm +
			theme(
				panel.grid=element_blank(),
				legend.position=c(.43,.19),
				legend.background=element_rect(
					fill=scales::alpha("white",.85),
					colour=NA
				),
				legend.key=element_blank(),
				legend.text=element_text(
					size=7.8,face="bold"
				),
				plot.margin=margin(6,8,6,6)
			)

		# --------------------------------------------------------------
		# B. Confusion matrices
		# --------------------------------------------------------------
		cm <- cm %>%
			filter(
				FeatureSet==primary_set_eval,
				Method %in% core_models
			) %>%
			mutate(
				True_lab=factor(
					ifelse(True==1,"Stroke","No stroke"),
					levels=c("No stroke","Stroke")
				),
				Pred_lab=factor(
					ifelse(Predicted==1,"Stroke","No stroke"),
					levels=c("No stroke","Stroke")
				)
			)

		cm_plot_one <- function(method){
			d <- cm %>% filter(Method==method)
			if(!nrow(d)) return(ggplot()+theme_void())

			cap <- sprintf(
				"Spec=%.3f | Sens=%.3f | Acc=%.3f",
				d$Specificity[1],
				d$Sensitivity[1],
				d$Accuracy[1]
			)

			ggplot(
				d,
				aes(Pred_lab,True_lab,fill=N)
			) +
				geom_tile(
					colour="white",linewidth=.8
				) +
				geom_text(
					aes(label=scales::comma(N)),
					size=3.1,fontface="bold",
					colour="grey15"
				) +
				scale_fill_gradient(
					low="#EEF3F8",
					high=C$hi,
					guide="none"
				) +
				labs(
					x="Predicted label",
					y="True label",
					title=model_labs[method],
					caption=cap
				) +
				coord_equal() +
				theme_minimal(base_size=8) +
				theme(
					panel.grid=element_blank(),
					axis.text.x=element_text(
						size=7.2,angle=0
					),
					axis.text.y=element_text(
						size=7.2
					),
					axis.title=element_text(
						size=7.5,face="bold"
					),
					plot.title=element_text(
						size=9,face="bold",
						hjust=.5
					),
					plot.caption=element_text(
						size=6.8,hjust=.5,
						face="bold"
					),
					plot.margin=margin(3,3,3,3)
				)
		}

		p_cm_list <- lapply(
			core_models,
			cm_plot_one
		)

		pPerfB0 <- wrap_plots(
			p_cm_list,
			ncol=2
		) +
			plot_annotation(
				title="B.",
				theme=theme(
					plot.title=element_text(
						size=13,face="plain",hjust=0
					)
				)
			)

		# Protect the internal 2x2 layout from outer patchwork flattening.
		pPerfB <- wrap_elements(
			full=patchworkGrob(pPerfB0)
		)

		# --------------------------------------------------------------
		# C. Metric bars + bootstrap 95% CIs
		# --------------------------------------------------------------
		metric_keep <- c(
			"AUC","Accuracy","Precision","Recall","F1"
		)

		metric_labs <- c(
			AUC="AUC",
			Accuracy="Accuracy",
			Precision="Precision",
			Recall="Recall",
			F1="F1"
		)

		perf_bar <- boot %>%
			filter(
				FeatureSet==primary_set_eval,
				Method %in% core_models,
				Metric %in% metric_keep
			) %>%
			mutate(
				Method=factor(Method,levels=core_models),
				Metric=factor(
					Metric,
					levels=metric_keep,
					labels=metric_labs[metric_keep]
				)
			)

		pPerfC <- ggplot(
			perf_bar,
			aes(
				Method,Estimate,
				fill=Metric
			)
		) +
			geom_col(
				position=position_dodge(width=.78),
				width=.68
			) +
			geom_errorbar(
				aes(
					ymin=Lower95,
					ymax=Upper95
				),
				position=position_dodge(width=.78),
				width=.13,
				linewidth=.35
			) +
			scale_x_discrete(
				labels=model_labs
			) +
			scale_fill_manual(
				values=c(
					AUC=C$hi,
					Accuracy=C$ora,
					Precision="#6A994E",
					Recall="#B8A23A",
					F1="#4A7EBB"
				)
			) +
			scale_y_continuous(
				limits=c(0,1),
				expand=expansion(mult=c(0,.03))
			) +
			labs(
				x=NULL,
				y="Performance metric",
				title="C."
			) +
			thm +
			theme(
				panel.grid.minor=element_blank(),
				axis.text.x=element_text(
					angle=18,hjust=1,
					size=8.2
				),
				legend.position="top",
				legend.text=element_text(size=7.5),
				plot.margin=margin(6,8,6,6)
			)

		# --------------------------------------------------------------
		# MYTH-like composition:
		# A left-top; C left-bottom; B occupies right side.
		# --------------------------------------------------------------
		pLeft <- wrap_plots(
			pPerfA,pPerfC,
			ncol=1,
			heights=c(1.15,1)
		)

		pLeftGrob <- wrap_elements(
			full=patchworkGrob(pLeft)
		)

		pPerf <- wrap_plots(
			pLeftGrob,pPerfB,
			ncol=2,
			widths=c(1.03,1.25)
		)

		ggsave(
			file.path(
				DIR_FIG,
				"Fig_ML_performance_bootstrap.pdf"
			),
			pPerf,
			width=13.2,height=7.4,
			device=cairo_pdf
		)
	}
}

# ----------------------------------------------------------------------
# Optional compact feature-set comparison plot:
# LightGBM Primary vs All/SHAP/RFE using paired bootstrap CIs.
# ----------------------------------------------------------------------
if(file.exists(f_featcmp)){
	fcmp <- as.data.frame(fread(f_featcmp))

	if(nrow(fcmp)){
		fcmp_long <- bind_rows(
			fcmp %>%
				transmute(
					Comparator,
					Metric="Delta AUC",
					Estimate=Delta_AUC,
					Lower=Delta_AUC_Lower95_boot,
					Upper=Delta_AUC_Upper95_boot,
					P=DeLong_P
				),
			fcmp %>%
				transmute(
					Comparator,
					Metric="Delta AUPRC",
					Estimate=Delta_AUPRC,
					Lower=Delta_AUPRC_Lower95_boot,
					Upper=Delta_AUPRC_Upper95_boot,
					P=Bootstrap_P_DeltaAUPRC
				)
		) %>%
			mutate(
				Comparator=gsub(
					"_"," ",Comparator
				)
			)

		pCmp <- ggplot(
			fcmp_long,
			aes(
				Estimate,
				forcats::fct_rev(Comparator)
			)
		) +
			geom_vline(
				xintercept=0,
				linetype=2,
				colour="grey65"
			) +
			geom_errorbar(
				aes(xmin=Lower,xmax=Upper),
				orientation="y", width=.15,
				linewidth=.45
			) +
			geom_point(
				size=2.4,
				colour=C$ora
			) +
			facet_wrap(
				~Metric,
				scales="free_x",
				ncol=1
			) +
			labs(
				x="Primary - comparator",
				y=NULL
			) +
			thm +
			theme(
				panel.grid.minor=element_blank(),
				panel.grid.major.y=element_blank(),
				strip.text=element_text(face="bold")
			)

		ggsave(
			file.path(
				DIR_FIG,
				"Fig_ML_paired_comparison.pdf"
			),
			pCmp,
			width=6.6,height=5.2,
			device=cairo_pdf
		)
	}
}

cat("Final ML plots complete\n")
cat(" - primary criterion: DeltaLogLoss005\n")
cat(" - selected k:",k,"\n")
cat(" - feature-selection validation methods: 10\n")
cat("   3 model-based + 3 SHAP + 4 RFE\n")
cat(" - predictive validation: 4 feature sets x",
	length(unique(val_fold$Method)),
	"algorithms =",4*length(unique(val_fold$Method)),
	"model-feature-set combinations\n")
cat(" - main figure:",file.path(DIR_FIG,"Fig_ML.pdf"),"\n")
cat(" - cumulative Cox-style figure:",file.path(DIR_FIG,"Fig_ML_Cox_path.pdf"),"\n")
cat(" - mirrored ML/Cox figure:",file.path(DIR_FIG,"Fig_ML_Cox_mirror.pdf"),"\n")
cat(" - combined pAB + Cox figure:",file.path(DIR_FIG,"Fig_ML_pAB_Cox.pdf"),"\n")
if(file.exists(f_boot)){
	cat(" - bootstrap performance figure:",file.path(DIR_FIG,"Fig_ML_performance_bootstrap.pdf"),"\n")
	cat(" - paired-comparison figure:",file.path(DIR_FIG,"Fig_ML_paired_comparison.pdf"),"\n")
}
}))

if (sys.nframe() == 0L) dose_run_module("ml")
