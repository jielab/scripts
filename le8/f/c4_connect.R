# C4: Connection — LE8-supervised proxy discovery, module assignment, and associational mediation.
suppressPackageStartupMessages({
  source(file.path(Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f")), "comm.f.R"))
})
LE8_JOB <- "c4_connect"
MAX_N <- as.integer(Sys.getenv("C4_MAX_N", unset = "60000"))
BLOCK <- as.integer(Sys.getenv("C4_BLOCK", unset = "80"))
FDR_CUT <- as.numeric(Sys.getenv("C4_FDR", unset = "0.05"))
SPEC_CUT <- as.numeric(Sys.getenv("C4_SPECIFICITY", unset = "0.35"))
YSP_PLUS_N <- as.integer(Sys.getenv("C4_YSP_PLUS", unset = "80"))
NS_N <- as.integer(Sys.getenv("C4_NS_TOP", unset = "160"))
MED_TOP <- as.integer(Sys.getenv("C4_MED_TOP", unset = "30"))
MED_BOOT <- as.integer(Sys.getenv("C4_MED_BOOT", unset = "100"))
MED_MAX_N <- as.integer(Sys.getenv("C4_MED_MAX_N", unset = "30000"))

stratified_sample <- function(d, event, max_n) {
  if (max_n <= 0 || nrow(d) <= max_n) return(d)
  if (!event %in% names(d)) return(slice_sample(d, n = max_n))
  n_all <- nrow(d)
  z <- d |> group_by(.data[[event]]) |> mutate(.rand = runif(n()), .rank = rank(.rand, ties.method = "first"),
                                                # Use doubles for the proportional allocation: max_n and
                                                # n() are integers, and their product can exceed 2^31 - 1.
                                                .take = max(1, ceiling(as.double(max_n) * n() / as.double(n_all)))) |>
    filter(.rank <= .take) |> ungroup() |> select(-.rand, -.rank, -.take)
  if(nrow(z) > max_n) z <- slice_sample(z, n = max_n)
  z
}

stratified_split <- function(d, event) {
  out <- rep(NA_character_, nrow(d)); set.seed(SEED)
  if (!event %in% names(d)) return(sample(rep(c("discovery", "replication"), length.out = nrow(d))))
  strata <- ifelse(is.na(d[[event]]), "missing", as.character(d[[event]]))
  for (v in unique(strata)) {
    id <- which(strata == v)
    out[id] <- sample(rep(c("discovery", "replication"), length.out = length(id)))
  }
  out
}

proxy_scan <- function(d, features, components, covars, split_name) {
  d <- d[complete.cases(d[, unique(c(components,covars)),drop=FALSE]),,drop=FALSE]
  features <- intersect(features,names(d)); components<-intersect(components,names(d));covars<-intersect(covars,names(d))
  if(!length(features)||!length(components)||nrow(d)<500)return(tibble())
  blocks<-split(features,ceiling(seq_along(features)/BLOCK));out<-list();k<-0L
  for(cmp in components){
    other<-setdiff(components,cmp);mm<-model.matrix(reformulate(unique(c(covars,other))),d);q<-qr(mm)
    y<-std_num(d[[cmp]]);yr<-qr.resid(q,y);sy<-sqrt(sum(yr^2));df<-max(3,nrow(d)-ncol(mm)-2)
    for(bb in blocks){
      x<-as.matrix(d[,bb,drop=FALSE]);storage.mode(x)<-"double"
      for(j in seq_len(ncol(x))){m<-median(x[,j],na.rm=TRUE);if(!is.finite(m))m<-0;x[!is.finite(x[,j]),j]<-m;x[,j]<-std_num(x[,j]);x[!is.finite(x[,j]),j]<-0}
      xr<-qr.resid(q,x);den<-sqrt(colSums(xr^2))*sy;r<-as.numeric(crossprod(xr,yr)/pmax(den,1e-12));r<-cap(r,.999)
      z<-r*sqrt(df/pmax(1-r^2,1e-9));p<-2*pt(abs(z),df=df,lower.tail=FALSE);k<-k+1L
      out[[k]]<-tibble(feature=bb,component=sub("\\.pts$","",cmp),component_var=cmp,r=r,z=z,p.value=p,n=nrow(d),split=split_name)
    }
  }
  bind_rows(out)|>group_by(component)|>mutate(FDR_component=p.adjust(p.value,"BH"))|>ungroup()|>mutate(FDR_global=p.adjust(p.value,"BH"))
}

make_proxy_sets <- function(disc,rep,disease_assoc,annotation) {
  x<-disc|>select(feature,component,component_var,r_disc=r,z_disc=z,p_disc=p.value,FDR_disc=FDR_component)|>
    left_join(rep|>select(feature,component,r_rep=r,z_rep=z,p_rep=p.value,FDR_rep=FDR_component),by=c("feature","component"))
  spec<-disc|>group_by(feature)|>mutate(absz=abs(z),specificity=absz/sum(absz,na.rm=TRUE))|>slice_max(absz,n=1,with_ties=FALSE)|>
    select(feature,primary_component=component,primary_component_var=component_var,specificity)
  x_primary<-x|>left_join(spec,by="feature");primary<-x_primary|>filter(component==primary_component)|>
    mutate(same_direction=is.finite(r_rep)&sign(r_disc)==sign(r_rep),
           strict_YS=coalesce(FDR_disc<FDR_CUT&FDR_rep<FDR_CUT&same_direction&specificity>=SPEC_CUT,FALSE),
           YS_model=strict_YS,selection_rule=ifelse(strict_YS,"YS_strict: replicated and pillar-specific","not supervised"))
  if(sum(primary$strict_YS,na.rm=TRUE)<8){
    fallback<-primary|>filter(!strict_YS)|>group_by(primary_component)|>arrange(FDR_disc,FDR_rep,desc(specificity))|>slice_head(n=5)|>ungroup()|>pull(feature)
    primary<-primary|>mutate(YS_model=strict_YS|feature%in%fallback,
      selection_rule=case_when(strict_YS~"YS_strict: replicated and pillar-specific",feature%in%fallback~"YS_exploratory: LE8-ranked fallback",TRUE~"not supervised"))
  }
  primary<-primary|>mutate(YS=YS_model,proxy_level=case_when(strict_YS~"YS_strict",YS_model~"YS_exploratory",TRUE~"not supervised"))
  ys_strict<-primary|>filter(strict_YS)|>pull(feature)|>unique();ys<-primary|>filter(YS_model)|>pull(feature)|>unique()
  dis <- disease_assoc
  if (!"beta" %in% names(dis)) dis$beta <- if ("estimate" %in% names(dis)) safe_log(dis$estimate) else NA_real_
  if (!"FDR" %in% names(dis)) dis$FDR <- p.adjust(dis$p.value, "BH")
  dis <- dis |> arrange(p.value)
  plus<-dis|>filter(!term%in%ys)|>slice_head(n=YSP_PLUS_N)|>pull(term)|>unique()
  ns<-dis|>slice_head(n=NS_N)|>pull(term)|>unique()
  membership<-tibble(feature=unique(c(ys,plus,ns)))|>mutate(in_YS=feature%in%ys,in_YSP_plus=feature%in%plus,in_NS=feature%in%ns,
    set=case_when(in_YS~"YS",in_YSP_plus~"YSP_plus",TRUE~"NS"))|>
    left_join(primary|>select(feature,primary_component,specificity,strict_YS,YS_model,proxy_level,selection_rule,r_disc,r_rep,FDR_disc,FDR_rep),by="feature")|>
    left_join(dis|>select(feature=term,disease_beta=beta,disease_p=p.value,disease_FDR=FDR),by="feature")|>left_join(annotation,by="feature")
  list(primary=primary,YS=ys,YS_strict=ys_strict,YSP_plus=plus,YSP=unique(c(ys,plus)),NS=ns,membership=membership,
       YS_edges=x_primary|>filter(feature%in%ys,component==primary_component)|>
         left_join(primary|>select(feature,strict_YS,YS_model,proxy_level,selection_rule),by="feature")|>left_join(annotation,by="feature"))
}

mediation_one <- function(dat, component_var, feature, covars, tvar, evar, B=100) {
  need<-unique(c(component_var,feature,covars,tvar,evar));d<-dat[,intersect(need,names(dat)),drop=FALSE];d<-d[complete.cases(d),,drop=FALSE]
  if(nrow(d)<1000||sum(d[[evar]]==1)<50)return(tibble())
  d[[component_var]]<-std_num(d[[component_var]]);d[[feature]]<-std_num(d[[feature]])
  fa<-tryCatch(lm(reformulate(c(component_var,covars),response=feature),d),error=function(e)NULL)
  ft<-tryCatch(coxph(as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",paste(bt(c(component_var,covars)),collapse=" + "))),d),error=function(e)NULL)
  fd<-tryCatch(coxph(as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",paste(bt(c(component_var,feature,covars)),collapse=" + "))),d),error=function(e)NULL)
  if(any(vapply(list(fa,ft,fd),is.null,logical(1))))return(tibble())
  sma <- coef(summary(fa)); smt <- coef(summary(ft)); smd <- coef(summary(fd))
  if (!component_var %in% rownames(sma) || !component_var %in% rownames(smt) ||
      !component_var %in% rownames(smd) || !feature %in% rownames(smd)) return(tibble())
  a<-sma[component_var,"Estimate"];ase<-sma[component_var,"Std. Error"];b<-smd[feature,"coef"];bse<-smd[feature,"se(coef)"]
  total<-smt[component_var,"coef"];direct<-smd[component_var,"coef"];ind<-a*b;seind<-sqrt(b^2*ase^2+a^2*bse^2);p<-ifelse(is.finite(seind)&&seind>0,2*pnorm(abs(ind/seind),lower.tail=FALSE),NA_real_)
  boot<-numeric();if(B>0){
    for(i in seq_len(B)){id<-sample.int(nrow(d),replace=TRUE);db<-d[id,,drop=FALSE];za<-tryCatch(coef(lm(reformulate(c(component_var,covars),response=feature),db))[[component_var]],error=function(e)NA_real_);zb<-tryCatch(coef(coxph(as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ ",paste(bt(c(component_var,feature,covars)),collapse=" + "))),db))[[feature]],error=function(e)NA_real_);boot[i]<-za*zb}
  }
  boot<-boot[is.finite(boot)];blo<-if(length(boot)>=20)as.numeric(quantile(boot,.025,names=FALSE))else NA_real_;bhi<-if(length(boot)>=20)as.numeric(quantile(boot,.975,names=FALSE))else NA_real_
  tibble(component=sub("\\.pts$","",component_var),component_var,feature,n=nrow(d),events=sum(d[[evar]]==1),a_beta=a,b_beta=b,total_beta=total,direct_beta=direct,indirect_beta=ind,indirect_se=seind,indirect_p=p,
         indirect_lo=blo,indirect_hi=bhi,bootstrap_valid=length(boot),prop_mediated=ifelse(total!=0,ind/total,NA_real_))
}

plot_c4 <- function(scan,sets,med,layer,outdir) {
  topf<-sets$primary|>arrange(desc(YS),FDR_disc)|>slice_head(n=80)|>pull(feature)
  hm<-scan|>filter(split=="discovery",feature%in%topf)|>mutate(component=factor(component,levels=names.le8),feature=factor(feature,levels=rev(topf)),zcap=cap(z,8))
  p1<-if(!nrow(hm)) blank_plot(paste0("LE8-supervised ",layer," proxy map"),"No discovery association was available") else ggplot(hm,aes(component,feature,fill=zcap))+geom_tile(color="white",linewidth=.15)+scale_fill_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0,name="z")+
    labs(title=paste0("LE8-supervised ",layer," proxy map"),subtitle="Discovery associations; proxy status requires replication and pillar specificity",x=NULL,y=NULL)+theme_5c(8)+theme(axis.text.x=element_text(angle=35,hjust=1))
  save_plot(p1,"c4.Fig01.proxy_heatmap.png",10,11,outdir=outdir)

  flow<-sets$membership|>filter(set%in%c("YS","YSP_plus"))|>count(group=coalesce(group,"Other"),primary_component=coalesce(primary_component,"Disease-selected"),set)
  if(requireNamespace("ggalluvial",quietly=TRUE)&&nrow(flow)){
    p2<-ggplot(flow,aes(axis1=group,axis2=primary_component,axis3=set,y=n))+ggalluvial::geom_alluvium(aes(fill=primary_component),width=1/12,alpha=.75)+ggalluvial::geom_stratum(width=1/8,fill="grey94",color="grey55")+ggalluvial::stat_stratum(geom="text",aes(label=after_stat(stratum)),size=2.6)+scale_x_discrete(limits=c("Omic group","LE8 pillar","Proxy set"),expand=c(.08,.08))+scale_fill_manual(values=cols_le8,na.value="grey65",guide="none")+labs(title="Omic group → LE8 pillar → supervised set",y="Number of features",x=NULL)+theme_5c(10)
  } else if(nrow(flow)) {
    p2<-ggplot(flow,aes(primary_component,fct_reorder(group,n,sum),size=n,color=set))+geom_point(alpha=.8)+scale_size_continuous(range=c(2,10))+scale_color_manual(values=c(YS="#D95F02",YSP_plus="#1B9E77"))+labs(title="Omic groups assigned to LE8 pillars",x=NULL,y=NULL,color=NULL,size="Features")+theme_5c(10)+theme(axis.text.x=element_text(angle=35,hjust=1))
  } else p2<-blank_plot("Omic groups assigned to LE8 pillars","No proxy-set flow could be formed")
  save_plot(p2,"c4.Fig02.group_pillar_flow.png",12,7,outdir=outdir)

  edges<-sets$YS_edges|>arrange(FDR_disc)|>slice_head(n=100)
  if(nrow(edges)){
    pillars<-tibble(name=unique(edges$component),type="pillar",angle=seq(0,2*pi,length.out=n_distinct(edges$component)+1)[-(n_distinct(edges$component)+1)],r=.35)
    feats<-tibble(name=unique(edges$feature),type="feature",angle=seq(0,2*pi,length.out=n_distinct(edges$feature)+1)[-(n_distinct(edges$feature)+1)]+.12,r=1)
    nodes<-bind_rows(pillars,feats)|>mutate(x=r*cos(angle),y=r*sin(angle));ee<-edges|>left_join(nodes|>select(component=name,x1=x,y1=y),by="component")|>left_join(nodes|>select(feature=name,x2=x,y2=y),by="feature")
    p3<-ggplot()+geom_curve(data=ee,aes(x=x1,y=y1,xend=x2,yend=y2,color=component,alpha=pmin(abs(z_disc)/8,1)),curvature=.12,linewidth=.55)+geom_point(data=nodes,aes(x,y,shape=type),size=2.5)+ggrepel::geom_text_repel(data=nodes,aes(x,y,label=name),size=2.3,max.overlaps=35,seed=11)+scale_color_manual(values=cols_le8)+scale_alpha(range=c(.2,.9),guide="none")+coord_equal()+theme_void()+labs(title="LE8 pillar–proxy module network",color="Pillar")+theme(plot.title=element_text(face="bold"),legend.position="bottom")
  } else p3<-blank_plot("LE8 pillar–proxy module network")
  save_plot(p3,"c4.Fig03.module_network.png",11,9,outdir=outdir)

  # Globe/orbit view: the radius encodes pillar specificity and the angle encodes the primary LE8 component.
  orbit<-sets$membership|>filter(set%in%c("YS","YSP_plus"))|>mutate(comp=factor(coalesce(primary_component,"disease"),levels=c(names.le8,"disease")),angle=2*pi*(as.numeric(comp)-1)/nlevels(comp)+runif(n(),-.18,.18),radius=ifelse(set=="YS",.65+.35*coalesce(specificity,0),1.18),x=radius*cos(angle),y=radius*sin(angle),lab=ifelse(min_rank(disease_p)<=25|set=="YS"&min_rank(FDR_disc)<=20,feature,NA_character_))
  p4<-if(!nrow(orbit)) blank_plot("LE8–omics relationship globe","No supervised or plus feature was available") else ggplot(orbit,aes(x,y))+annotate("path",x=cos(seq(0,2*pi,length.out=300)),y=sin(seq(0,2*pi,length.out=300)),color="grey80")+geom_segment(aes(x=0,y=0,xend=x,yend=y,color=primary_component),alpha=.12)+geom_point(aes(color=primary_component,shape=set,size=-log10(pmax(disease_p,1e-30))),alpha=.8)+ggrepel::geom_text_repel(aes(label=lab),size=2.4,max.overlaps=30,seed=12,na.rm=TRUE)+scale_color_manual(values=c(cols_le8,disease="grey35"),na.value="grey60")+scale_size_continuous(range=c(1.5,5),name="Disease evidence")+coord_equal()+theme_void()+labs(title="LE8–omics relationship globe",color="Primary pillar",shape="Set")+theme(plot.title=element_text(face="bold"),legend.position="bottom")
  save_plot(p4,"c4.Fig04.relationship_globe.png",11,9,outdir=outdir)

  if(nrow(med)){
    wheel<-med|>filter(is.finite(prop_mediated))|>arrange(indirect_p)|>slice_head(n=60)|>mutate(wheel_label=paste(feature,component,sep=" | "),wheel_label=factor(wheel_label,levels=wheel_label),id=row_number(),angle=90-360*(id-.5)/n(),hjust=ifelse(angle< -90,1,0),angle=ifelse(angle< -90,angle+180,angle),value=100*cap(prop_mediated,1))
    p5<-if(!nrow(wheel)) blank_plot("Associational mediation wheel","Mediation effects were estimable, but no finite proportion mediated was available") else ggplot(wheel,aes(id,abs(value),fill=component))+geom_col(width=.8)+geom_hline(yintercept=0)+geom_text(aes(y=abs(value)+3,label=wheel_label,angle=angle,hjust=hjust),size=2,na.rm=TRUE)+coord_polar(clip="off")+scale_fill_manual(values=cols_le8)+ylim(-10,max(abs(wheel$value),na.rm=TRUE)+15)+labs(title="Associational mediation wheel",subtitle="Percent mediated is shown for screened LE8 → omic → disease paths",fill="LE8 pillar")+theme_void(10)+theme(plot.title=element_text(face="bold",hjust=.5),plot.subtitle=element_text(hjust=.5),legend.position="bottom",plot.margin=margin(30,55,30,55))
    save_plot(p5,"c4.Fig05.mediation_wheel.png",11,10,outdir=outdir)
    md<-med|>filter(is.finite(indirect_beta),is.finite(total_beta),is.finite(direct_beta))|>arrange(indirect_p)|>slice_head(n=50)|>mutate(label=paste(feature,component,sep=" | "),label=factor(label,levels=rev(label)))
    if(!nrow(md)) {
      save_plot(blank_plot("Mediation diagnostics","No finite mediation diagnostic row was available"),"c4.Fig06.mediation_diagnostics.png",9,5,outdir=outdir)
    } else {
      pa<-ggplot(md,aes(indirect_beta,label))+geom_vline(xintercept=0,color="grey60")+geom_errorbarh(aes(xmin=coalesce(indirect_lo,indirect_beta-1.96*indirect_se),xmax=coalesce(indirect_hi,indirect_beta+1.96*indirect_se)),height=0)+geom_point(aes(color=component),size=2)+scale_color_manual(values=cols_le8)+labs(title="A. Indirect effects",x="Indirect log-hazard effect",y=NULL,color=NULL)+forest_theme(9)
      pb<-ggplot(md,aes(total_beta,direct_beta,color=component))+geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+geom_point(aes(size=abs(prop_mediated)),alpha=.8)+scale_color_manual(values=cols_le8)+labs(title="B. Total versus direct effects",x="Total effect",y="Direct effect",color=NULL,size="|Proportion mediated|")+theme_5c(10)
      pc<-ggplot(md,aes(100*cap(prop_mediated,1),-log10(pmax(indirect_p,1e-30)),color=component))+geom_vline(xintercept=0,color="grey60")+geom_point(alpha=.8)+scale_color_manual(values=cols_le8)+labs(title="C. Proportion and evidence",x="Percent mediated (plot capped)",y=expression(-log[10](P[indirect])),color=NULL)+theme_5c(10)
      save_plot(pa | (pb/pc),"c4.Fig06.mediation_diagnostics.png",14,10,outdir=outdir)
    }
  } else {
    save_plot(blank_plot("Associational mediation wheel","No screened YS path could be estimated"),"c4.Fig05.mediation_wheel.png",9,5,outdir=outdir)
    save_plot(blank_plot("Mediation diagnostics","No screened YS path could be estimated"),"c4.Fig06.mediation_diagnostics.png",9,5,outdir=outdir)
  }
}

run_c4_layer <- function(layer=c("protein","metabolite")) {
  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir);rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE);cache<-file.path(rawdir,"c4.res.rds")
  if(cache_valid(cache)){cache_message(paste0("C4/",layer),cache);return(readRDS(cache))}
  c1f<-file.path(le8_job_dir(outdir,"c1_correlate"),"c1.res.rds");if(!file.exists(c1f))stop("Run C1 first.",call.=FALSE);c1<-readRDS(c1f);disease<-(c1$association%||%c1$pwas_incident%||%c1$MWAS)|>as_tibble();if(!"beta"%in%names(disease))disease<-disease|>mutate(beta=safe_log(estimate))
  biom<-if(layer=="protein")read_prot() else read_met();features<-setdiff(names(biom),"eid");ann<-layer_annotation(layer,features)
  need<-unique(c("eid","ethnic.c",vars.basic,vars.le8,"birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y)))
  dat0<-read_all(need)|>filter_analysis_cohort()|>make_outcome(Y);tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e")
  comps<-intersect(vars.le8,names(dat0));covs<-intersect(vars.basic,names(dat0))
  dat0<-dat0[complete.cases(dat0[,comps,drop=FALSE]),,drop=FALSE]|>stratified_sample(evar,MAX_N)
  # Sample participants before joining the wide omics matrix; this materially reduces C4 peak memory.
  dat<-inner_join(dat0,biom|>filter(eid%in%dat0$eid),by="eid")
  fold<-stratified_split(dat,evar);disc<-proxy_scan(dat[fold=="discovery",],features,comps,covs,"discovery");rep<-proxy_scan(dat[fold=="replication",],features,comps,covs,"replication")
  if(!nrow(disc)||!nrow(rep)) stop("C4/",layer,": discovery or replication proxy scan returned no result; check LE8 variables and sample size.",call.=FALSE)
  scan<-bind_rows(disc,rep);sets<-make_proxy_sets(disc,rep,disease,ann)
  write_raw_csv(scan,"c4.LE8_feature_associations.csv",rawdir);write_raw_csv(sets$primary,"c4.primary_pillar_assignment.csv",rawdir);write_raw_csv(sets$membership,"c4.proxy_membership_YS_YSP_NS.csv",rawdir);write_raw_csv(sets$YS_edges,"c4.YS_edges.csv",rawdir)
  pairs<-sets$YS_edges|>left_join(disease|>select(feature=term,disease_p=p.value),by="feature")|>arrange(desc(strict_YS),FDR_disc,disease_p)|>slice_head(n=MED_TOP)
  med_dat<-stratified_sample(dat,evar,MED_MAX_N)
  med<-if(!nrow(pairs))tibble() else map_dfr(seq_len(nrow(pairs)),function(i){
    rw<-pairs[i,,drop=FALSE]
    mediation_one(med_dat,rw$component_var[[1]],rw$feature[[1]],unique(c(covs,setdiff(comps,rw$component_var[[1]]))),tvar,evar,MED_BOOT)
  })
  if(nrow(med))med<-med|>mutate(FDR_indirect=p.adjust(indirect_p,"BH"))|>arrange(indirect_p)
  write_raw_csv(med,"c4.mediation_all.csv",rawdir);plot_c4(scan,sets,med,layer,outdir)
  lists<-list(YS_by_component=split(sets$YS_edges$feature,sets$YS_edges$component),YS_all=sets$YS,YS_strict=sets$YS_strict,YSP_plus=sets$YSP_plus,YSP_all=sets$YSP,NS=sets$NS,index=sets$membership)
  saveRDS(lists,file.path(rawdir,paste0("c4.",layer,"_lists.rds")),compress="xz");if(layer=="protein")saveRDS(lists,file.path(rawdir,"c4.protein_lists.rds"),compress="xz")
  out<-list(meta=module_meta(layer,extra=list(scan_N=nrow(dat),discovery_N=sum(fold=="discovery"),replication_N=sum(fold=="replication"),mediation_N=nrow(med_dat),YS_strict_n=length(sets$YS_strict),YS_model_n=length(sets$YS))),scan=scan,primary=sets$primary,membership=sets$membership,YS_edges=sets$YS_edges,mediation=med,lists=lists)
  saveRDS(out,cache,compress="xz");write_xlsx2(list(primary_assignment=sets$primary,proxy_membership=sets$membership,YS_edges=sets$YS_edges,all_LE8_associations=scan,mediation=med),"c4.out.xlsx");finalize_outputs(LE8_JOB,outdir);out
}

if(prot_DO)run_c4_layer("protein")
if(met_DO)run_c4_layer("metabolite")
