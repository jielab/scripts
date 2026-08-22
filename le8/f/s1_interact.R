# S1: sex-specific omics associations plus pairwise LE8 interactions.
suppressPackageStartupMessages({source(file.path(Sys.getenv("LE8_FDIR",unset=file.path(Sys.getenv("DIRSCRIPT"),"f")),"comm.f.R"))})
LE8_JOB<-"s1_interact";MAX_FEATURES<-as.integer(Sys.getenv("S1_SEX_MAX",unset="120"))

run_sex_layer<-function(layer=c("protein","metabolite")){
 layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir);rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
 c1f<-file.path(le8_job_dir(outdir,"c1_correlate"),"c1.res.rds");if(!file.exists(c1f))stop("Run C1 first.",call.=FALSE);c1<-readRDS(c1f);a<-(c1$association%||%c1$pwas_incident%||%c1$MWAS)|>as_tibble();features<-a|>arrange(p.value)|>slice_head(n=MAX_FEATURES)|>pull(term)
 biom<-if(layer=="protein")read_prot() else read_met();features<-intersect(features,setdiff(names(biom),"eid"));biom<-biom[,c("eid",features),drop=FALSE];need<-unique(c("eid","ethnic.c",vars.basic,"sex","birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y)))
 d<-read_all(need)|>left_join(biom,by="eid")|>filter_analysis_cohort()|>make_outcome(Y);tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e");covs<-setdiff(intersect(vars.basic,names(d)),"sex");d$sex_f<-factor(d$sex);features<-intersect(features,names(d))
 res<-map_dfr(features,function(f){
   z<-d[,unique(c(tvar,evar,f,"sex_f",covs)),drop=FALSE];z<-z[complete.cases(z),,drop=FALSE];z$sex_f<-droplevels(z$sex_f)
   if(nrow(z)<500||sum(z[[evar]]==1)<30||nlevels(z$sex_f)<2)return(tibble());z[[f]]<-std_num(z[[f]])
   rhs<-c(paste0(bt(f),"*sex_f"),bt(covs));fit<-tryCatch(coxph(as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",paste(rhs,collapse=" + "))),z),error=function(e)NULL)
   if(is.null(fit))return(tibble());sm<-coef(summary(fit));ir<-grep(":",rownames(sm),value=TRUE)[1];if(!length(ir))ir<-NA_character_
   st<-map_dfr(levels(z$sex_f),function(sx){
     zz<-droplevels(z[z$sex_f==sx,,drop=FALSE]);rhs0<-c(bt(f),bt(covs));ff<-tryCatch(coxph(as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",paste(rhs0,collapse=" + "))),zz),error=function(e)NULL)
     if(is.null(ff))return(tibble(sex=sx,beta=NA_real_,se=NA_real_,p=NA_real_));ss<-coef(summary(ff));if(!f%in%rownames(ss))return(tibble(sex=sx,beta=NA_real_,se=NA_real_,p=NA_real_))
     tibble(sex=sx,beta=ss[f,"coef"],se=ss[f,"se(coef)"],p=ss[f,"Pr(>|z|)"])
   })|>pivot_wider(names_from=sex,values_from=c(beta,se,p),names_glue="{.value}_sex{sex}")
   bind_cols(tibble(feature=f,n=nrow(z),events=sum(z[[evar]]==1),beta_interaction=if(!is.na(ir))sm[ir,"coef"]else NA_real_,se_interaction=if(!is.na(ir))sm[ir,"se(coef)"]else NA_real_,p_interaction=if(!is.na(ir))sm[ir,"Pr(>|z|)"]else NA_real_),st)
 })
 if(!nrow(res))res<-tibble(feature=character(),n=integer(),events=integer(),beta_interaction=numeric(),se_interaction=numeric(),p_interaction=numeric())
 res<-res|>mutate(FDR_interaction=p.adjust(p_interaction,"BH"))|>arrange(p_interaction)
 write_raw_csv(res,"s1.sex_interaction.csv",rawdir)
 bcols<-grep("^beta_sex",names(res),value=TRUE);if(length(bcols)>=2){pcols<-sub("beta_","p_",bcols);x<-bcols[1];y<-bcols[2];res$label<-ifelse(min_rank(res$p_interaction)<=20,res$feature,NA_character_)
  pA<-ggplot(res,aes(.data[[x]],.data[[y]]))+geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+geom_hline(yintercept=0,color="grey75")+geom_vline(xintercept=0,color="grey75")+geom_point(aes(color=FDR_interaction<.05),size=2,alpha=.8)+ggrepel::geom_text_repel(aes(label=label),size=2.7,seed=20,max.overlaps=25,na.rm=TRUE)+scale_color_manual(values=c(`TRUE`="#D95F02",`FALSE`="grey70"),guide="none")+labs(title=paste0("A. Sex-specific ",layer," associations"),x=paste0("Effect in sex group ",sub("beta_sex","",x)),y=paste0("Effect in sex group ",sub("beta_sex","",y)))+theme_5c(11)
  top<-res|>slice_head(n=40)|>mutate(feature=factor(feature,levels=rev(feature)),lo=beta_interaction-1.96*se_interaction,hi=beta_interaction+1.96*se_interaction)
  pB<-ggplot(top,aes(beta_interaction,feature))+geom_vline(xintercept=0,color="grey60")+geom_errorbarh(aes(xmin=lo,xmax=hi),height=.12)+geom_point(aes(color=FDR_interaction<.05),size=2)+scale_color_manual(values=c(`TRUE`="#D95F02",`FALSE`="#4C78A8"),guide="none")+labs(title="B. Omic-by-sex interaction",x="Interaction log-hazard coefficient",y=NULL)+forest_theme(9)
  save_plot(pA/pB+plot_layout(heights=c(.8,1.2)),"s1_int.Fig01.sex_interaction.png",10,11,outdir=outdir)
 }else save_plot(blank_plot("Sex interaction","Sex-specific estimates could not be formed"),"s1_int.Fig01.sex_interaction.png",9,5,outdir=outdir)
 write_xlsx2(list(sex_interaction=res),"s1_int.Fig01.sex_interaction.out.xlsx");finalize_outputs(LE8_JOB,outdir);res
}
run_le8_interactions <- function() {
  outdir <- file.path(out.base, LE8_JOB); setwd2(outdir)
  rawdir <- le8_job_dir(outdir, LE8_JOB); dir.create(rawdir, recursive=TRUE, showWarnings=FALSE)
  need <- unique(c("eid","ethnic.c",vars.basic,vars.le8,"birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y)))
  dat <- read_all(need) |> filter_analysis_cohort() |> make_outcome(Y)
  tvar <- paste0(Y,".t2e"); evar <- paste0(Y,".Yt2e")
  comps <- intersect(vars.le8,names(dat)); covs <- intersect(vars.basic,names(dat))
  if(length(comps)<2)stop("S1 requires at least two LE8 component variables.",call.=FALSE)
  interaction_formula <- function(x,z,cv)as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",paste(c(paste0(bt(x),"*",bt(z)),bt(cv)),collapse=" + ")))
  res <- map_dfr(combn(comps,2,simplify=FALSE),function(v){
    x<-v[[1]];z<-v[[2]];cv<-setdiff(covs,c(x,z));d<-dat[,unique(c(tvar,evar,x,z,cv)),drop=FALSE];d<-d[complete.cases(d),,drop=FALSE]
    if(nrow(d)<1000||sum(d[[evar]]==1)<50)return(tibble());d[[x]]<-std_num(d[[x]]);d[[z]]<-std_num(d[[z]])
    fit<-tryCatch(coxph(interaction_formula(x,z,cv),d),error=function(e)NULL);if(is.null(fit))return(tibble())
    sm<-coef(summary(fit));ir<-grep(":",rownames(sm),value=TRUE)[1];if(is.na(ir)||!ir%in%rownames(sm))return(tibble())
    tibble(xvar=x,zvar=z,x=sub("\\.pts$","",x),z=sub("\\.pts$","",z),n=nrow(d),events=sum(d[[evar]]==1),beta_interaction=sm[ir,"coef"],se_interaction=sm[ir,"se(coef)"],p_interaction=sm[ir,"Pr(>|z|)"])
  })
  if(!nrow(res))res<-tibble(x=character(),z=character(),xvar=character(),zvar=character(),n=integer(),events=integer(),beta_interaction=numeric(),se_interaction=numeric(),p_interaction=numeric())
  res<-res|>mutate(FDR_interaction=p.adjust(p_interaction,"BH"),pair=paste(coalesce(LE8_LABS[x],x),"×",coalesce(LE8_LABS[z],z)))|>arrange(p_interaction)
  write_raw_csv(res,"s1.LE8_pairwise_interactions.csv",rawdir)
  mat<-bind_rows(res|>transmute(row=x,col=z,beta=beta_interaction,p=p_interaction,FDR=FDR_interaction),res|>transmute(row=z,col=x,beta=beta_interaction,p=p_interaction,FDR=FDR_interaction))|>
    mutate(row=factor(row,levels=rev(names.le8)),col=factor(col,levels=names.le8),label=ifelse(FDR<.05,"*",ifelse(p<.05,"·","")))
  pA<-if(!nrow(mat))blank_plot("A. All pairwise LE8 interactions","No interaction model met the sample/event requirements")else ggplot(mat,aes(col,row,fill=cap(beta,.25)))+geom_tile(color="white")+geom_text(aes(label=label),fontface="bold",size=5)+scale_fill_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0,name="Interaction β")+labs(title="A. All pairwise LE8 interactions",subtitle="* FDR < 0.05; · nominal P < 0.05",x=NULL,y=NULL)+theme_5c(10)+theme(axis.text.x=element_text(angle=35,hjust=1))
  top<-res|>slice_head(n=4)
  grid<-if(!nrow(top))tibble()else map_dfr(seq_len(nrow(top)),function(i){
    rw<-top[i,,drop=FALSE];x<-rw$x[[1]];z<-rw$z[[1]];xvar<-rw$xvar[[1]];zvar<-rw$zvar[[1]];cv<-setdiff(covs,c(xvar,zvar));d<-dat[,unique(c(tvar,evar,xvar,zvar,cv)),drop=FALSE];d<-d[complete.cases(d),,drop=FALSE]
    if(nrow(d)<1000||sum(d[[evar]]==1)<50)return(tibble());d[[xvar]]<-std_num(d[[xvar]]);d[[zvar]]<-std_num(d[[zvar]]);fit<-tryCatch(coxph(interaction_formula(xvar,zvar,cv),d),error=function(e)NULL);if(is.null(fit))return(tibble())
    g<-expand_grid(xvalue=seq(-2,2,length.out=45),zvalue=seq(-2,2,length.out=45));nd<-g;names(nd)[1:2]<-c(xvar,zvar)
    for(v in cv)if(is.numeric(d[[v]]))nd[[v]]<-median(d[[v]],na.rm=TRUE)else{mode_value<-names(sort(table(d[[v]]),decreasing=TRUE))[1];nd[[v]]<-if(is.factor(d[[v]]))factor(mode_value,levels=levels(d[[v]]))else mode_value}
    lp<-tryCatch(as.numeric(predict(fit,newdata=nd,type="lp")),error=function(e)rep(NA_real_,nrow(nd)));g|>mutate(pair=rw$pair[[1]],relative_hazard=exp(lp-median(lp,na.rm=TRUE)))
  })|>filter(is.finite(relative_hazard))
  pB<-if(!nrow(grid))blank_plot("B. Interaction surfaces","No finite adjusted predictions were available")else ggplot(grid,aes(xvalue,zvalue,fill=pmin(relative_hazard,3)))+geom_raster()+geom_contour(aes(z=relative_hazard),color="white",alpha=.55,bins=6)+facet_wrap(~pair,nrow=1)+scale_fill_viridis_c(name="Relative hazard",option="C")+labs(title="B. Joint LE8 risk surfaces",subtitle="Axes are standardized LE8 component scores; surfaces are adjusted associations",x="First LE8 component (SD)",y="Second LE8 component (SD)")+theme_5c(9)
  save_plot(pA/pB+plot_layout(heights=c(.9,1.1)),"s1_int.Fig02.LE8_component_interactions.png",13,10,outdir=outdir)
  write_xlsx2(list(interaction_tests=res,prediction_surfaces=grid),"s1_int.Fig02.LE8_component_interactions.out.xlsx")
  saveRDS(list(results=res,surfaces=grid),file.path(rawdir,"s1.interactions.res.rds"),compress="xz");finalize_outputs(LE8_JOB,outdir)
}

if(prot_DO)invisible(run_sex_layer("protein"))
if(met_DO)invisible(run_sex_layer("metabolite"))
invisible(run_le8_interactions())
