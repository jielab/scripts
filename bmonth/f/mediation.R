# Linear-by-Cox coefficient products; exploratory associational mediation.
# Use internal .x/.m names so biological feature names cannot break formulas.
mediation_data <- function(d, x, m, y=NULL) {
	vars <- unique(c(x,m,setdiff(analysis_covars_for(x),x),if(!is.null(y)) paste0(y,c(".t2e",".Yt2e"))))
	if (!is.null(y)) d <- ensure_outcome(d,y)
	z <- complete_model_data(d,vars)
	if(nrow(z)<1000 || length(unique(z[[x]]))<2 || length(unique(z[[m]]))<2) return(NULL)
	if(!is.null(y) && sum(z[[paste0(y,".Yt2e")]]==1)<50) return(NULL)
	z$.x <- if(x %in% c("bmonth_winter_summer","bmonth_spring_autumn","border12","border1n")) z[[x]] else as.numeric(scale(z[[x]]))
	z$.m <- as.numeric(scale(z[[m]]))
	if(!is.null(y)) { z$.time<-z[[paste0(y,".t2e")]]; z$.event<-z[[paste0(y,".Yt2e")]] }
	z
}
a_formula <- function(x) reformulate(c(".x",setdiff(analysis_covars_for(x),x)),response=".m")
total_formula <- function(x) reformulate(c(".x",setdiff(analysis_covars_for(x),x)),response="Surv(.time,.event)")
direct_formula <- function(x) reformulate(c(".x",".m",setdiff(analysis_covars_for(x),x)),response="Surv(.time,.event)")
model_coef <- function(fit, term, type) {
	if(is.null(fit)) return(NULL)
	sm <- coef(summary(fit))
	cols <- if(type=="lm") c("Estimate","Std. Error","Pr(>|t|)") else c("coef","se(coef)","Pr(>|z|)")
	if(!term %in% rownames(sm) || any(!is.finite(sm[term,cols]))) return(NULL)
	unname(sm[term,cols])
}
screen_a <- function(d,x,m) {
	z<-mediation_data(d,x,m)
	r<-tibble(mediator=m,N=if(is.null(z)) NA_integer_ else nrow(z),a=NA_real_,p_a=NA_real_,note="insufficient data/variation")
	if(is.null(z)) return(r)
	fit<-model_try(lm(a_formula(x),z)); a<-model_coef(fit$fit,".x","lm");r$note<-fit$note
	if(!is.null(a)) {r$a<-a[1];r$p_a<-a[3]} else if(!nzchar(r$note)) r$note<-"Exposure coefficient not estimable"
	r
}
screen_b <- function(d,x,m,y) {
	z<-mediation_data(d,x,m,y)
	r<-tibble(mediator=m,N_b=if(is.null(z)) NA_integer_ else nrow(z),b=NA_real_,p_b=NA_real_,note_b="insufficient data/events/variation")
	if(is.null(z)) return(r)
	fit<-model_try(coxph(direct_formula(x),z));b<-model_coef(fit$fit,".m","cox");r$note_b<-fit$note
	if(!is.null(b)) {r$b<-b[1];r$p_b<-b[3]} else if(!nzchar(r$note_b)) r$note_b<-"Mediator coefficient not estimable"
	r
}
empty_mediation <- function() tibble(exposure=character(),mediator=character(),outcome=character(),
	analysis_track=character(),covariate_model=character(),N=integer(),events=integer(),
	a_beta=double(),b_beta=double(),total_beta=double(),direct_beta=double(),indirect_beta=double(),indirect_se=double(),
	p_sobel=double(),p_boot=double(),indirect_lo=double(),indirect_hi=double(),prop_mediated=double(),
	bootstrap_requested=integer(),bootstrap_valid=integer(),bootstrap_sufficient=logical(),exposure_unit=character(),
	status=character(),note=character())
mediate_feature <- function(d,x,m,y,boot=100) {
	z<-mediation_data(d,x,m,y)
	row<-tibble(exposure=x,mediator=m,outcome=y,
		analysis_track=analysis_track(x),covariate_model=analysis_model_label(x),N=if(is.null(z)) NA_integer_ else nrow(z),
		events=if(is.null(z)) NA_integer_ else sum(z$.event==1),a_beta=NA_real_,b_beta=NA_real_,total_beta=NA_real_,direct_beta=NA_real_,
		indirect_beta=NA_real_,indirect_se=NA_real_,p_sobel=NA_real_,p_boot=NA_real_,indirect_lo=NA_real_,indirect_hi=NA_real_,
		prop_mediated=NA_real_,bootstrap_requested=as.integer(boot),bootstrap_valid=0L,bootstrap_sufficient=FALSE,
		exposure_unit=if(x %in% c("bmonth_winter_summer","bmonth_spring_autumn","border12","border1n")) "0 to 1 contrast" else "per SD in complete-case mediation sample",
		status="skipped",note="insufficient data/events/variation")
	if(is.null(z)) return(row)
	fa<-model_try(lm(a_formula(x),z));ft<-model_try(coxph(total_formula(x),z));fd<-model_try(coxph(direct_formula(x),z))
	a<-model_coef(fa$fit,".x","lm");t<-model_coef(ft$fit,".x","cox")
	b<-model_coef(fd$fit,".m","cox");direct<-model_coef(fd$fit,".x","cox")
	row$note<-paste(Filter(nzchar,c(fa$note,ft$note,fd$note)),collapse="; ")
	if(any(vapply(list(a,b,t,direct),is.null,logical(1)))) {
		row$status<-"failed"; if(!nzchar(row$note))row$note<-"Coefficient not estimable"; return(row)
	}
	row$a_beta<-a[1];row$b_beta<-b[1];row$total_beta<-t[1];row$direct_beta<-direct[1]
	row$indirect_beta<-a[1]*b[1]; row$indirect_se<-sqrt(b[1]^2*a[2]^2+a[1]^2*b[2]^2)
	if(row$indirect_se>0) row$p_sobel<-2*pnorm(abs(row$indirect_beta/row$indirect_se),lower.tail=FALSE)
	if(abs(t[1])>1e-6) row$prop_mediated<-row$indirect_beta/t[1]
	bb<-rep(NA_real_,boot)
	for(i in seq_len(boot)) {
		q<-z[sample.int(nrow(z),replace=TRUE),,drop=FALSE]
		aa<-model_coef(model_try(lm(a_formula(x),q))$fit,".x","lm")
		bt<-model_coef(model_try(coxph(direct_formula(x),q))$fit,".m","cox")
		if(!is.null(aa)&&!is.null(bt)) bb[i]<-aa[1]*bt[1]
	}
	bb<-bb[is.finite(bb)];row$bootstrap_valid<-length(bb)
	row$bootstrap_sufficient<-boot>0 && length(bb)>=max(50,ceiling(.8*boot))
	if(row$bootstrap_sufficient) {
		row$p_boot<-min(1,2*min((sum(bb<=0)+1)/(length(bb)+1),(sum(bb>=0)+1)/(length(bb)+1)))
		row$indirect_lo<-quantile(bb,.025,names=FALSE);row$indirect_hi<-quantile(bb,.975,names=FALSE)
	}
	row$status<-"ok";row
}
