# C4: Connection — LE8-supervised proxy discovery, module assignment, and associational mediation.
suppressPackageStartupMessages({
  fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
  source(file.path(fdir, "comm.f.R"))
  source(file.path(fdir, "c4_state_network.R"))
})
LE8_JOB <- "c4_connect"
C4_CODE_VERSION <- "2026-09-05.5c-audit-v1"
MAX_N <- as.integer(Sys.getenv("C4_MAX_N", unset = "60000"))
BLOCK <- as.integer(Sys.getenv("C4_BLOCK", unset = "80"))
FDR_CUT <- as.numeric(Sys.getenv("C4_FDR", unset = "0.05"))
SPEC_CUT <- as.numeric(Sys.getenv("C4_SPECIFICITY", unset = "0.35"))
YSP_PLUS_N <- as.integer(Sys.getenv("C4_YSP_PLUS", unset = "80"))
NS_N <- as.integer(Sys.getenv("C4_NS_TOP", unset = "160"))
MED_TOP <- as.integer(Sys.getenv("C4_MED_TOP", unset = "30"))
MED_BOOT <- as.integer(Sys.getenv("C4_MED_BOOT", unset = "100"))
MED_MAX_N <- as.integer(Sys.getenv("C4_MED_MAX_N", unset = "30000"))
C4_MODULE_MAX <- as.integer(Sys.getenv("C4_MODULE_MAX",unset="1200"))
C4_MODULE_BOOT <- as.integer(Sys.getenv("C4_MODULE_BOOT",unset="100"))
C4_MODULE_STABILITY <- as.numeric(Sys.getenv("C4_MODULE_STABILITY",unset="0.70"))
C4_MODULE_K_MAX <- as.integer(Sys.getenv("C4_MODULE_K_MAX",unset="8"))

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

make_genetic_edges <- function(disc,rep,disease,membership){
  # Preserve the edge schema when PRS inputs are absent or a scan is not estimable.
  # Plotting and list export must handle unavailable evidence without inventing edges.
  if(!nrow(disc)||!nrow(rep)) {
    message("C4: genetic bridges unavailable; discovery rows=", nrow(disc),
      "; replication rows=", nrow(rep))
    return(tibble(feature=character(),prs=character(),prs_var=character(),
      r_disc=double(),z_disc=double(),p_disc=double(),FDR_disc=double(),
      r_rep=double(),z_rep=double(),p_rep=double(),FDR_rep=double(),
      same_direction=logical(),replicated=logical(),disease_beta=double(),
      disease_p=double(),set=character(),primary_component=character(),bridge_score=double()))
  }
  dis<-disease;if(!"beta"%in%names(dis))dis$beta<-safe_log(dis$estimate)
  disc|>select(feature,prs=component,prs_var=component_var,r_disc=r,z_disc=z,p_disc=p.value,FDR_disc=FDR_component)|>
    left_join(rep|>select(feature,prs=component,r_rep=r,z_rep=z,p_rep=p.value,FDR_rep=FDR_component),by=c("feature","prs"))|>
    mutate(same_direction=is.finite(r_rep)&sign(r_disc)==sign(r_rep),replicated=coalesce(FDR_disc<FDR_CUT&FDR_rep<FDR_CUT&same_direction,FALSE))|>
    left_join(dis|>select(feature=term,disease_beta=beta,disease_p=p.value),by="feature")|>
    left_join(membership|>select(feature,set,primary_component),by="feature")|>
    mutate(bridge_score=-log10(pmax(FDR_disc,1e-300))+-log10(pmax(FDR_rep,1e-300))+-log10(pmax(disease_p,1e-300)))|>
    arrange(desc(replicated),desc(bridge_score))
}

supervised_modules <- function(disc,rep,sets){
  p<-disc|>select(feature,component,z_disc=z)|>inner_join(rep|>select(feature,component,z_rep=z),by=c("feature","component"))|>
    mutate(same_direction=is.finite(z_rep)&sign(z_disc)==sign(z_rep),z_joint=ifelse(same_direction,(z_disc+z_rep)/sqrt(2),.25*(z_disc+z_rep)),strength=pmax(abs(z_disc),abs(z_rep),na.rm=TRUE))
  if(!nrow(p))return(list(membership=tibble(),metrics=tibble(),profile=tibble(),YS_core=character()))
  keep<-p|>group_by(feature)|>summarise(strength=max(strength,na.rm=TRUE),.groups="drop")|>arrange(desc(strength))|>slice_head(n=C4_MODULE_MAX)|>pull(feature)
  pw<-p|>filter(feature%in%keep)|>select(feature,component,z_joint)|>pivot_wider(names_from=component,values_from=z_joint,values_fill=0)
  rn<-pw$feature;X<-as.matrix(pw[,-1,drop=FALSE]);rownames(X)<-rn;X[!is.finite(X)]<-0
  X<-t(scale(t(X)));X[!is.finite(X)]<-0
  n<-nrow(X);ks<-if(n>=4)2:min(C4_MODULE_K_MAX,n-1)else integer();D<-as.matrix(dist(X))
  sil_one<-function(cl){
    mean(vapply(seq_len(n),function(i){same<-which(cl==cl[i]&seq_len(n)!=i);a<-if(length(same))mean(D[i,same])else 0;other<-setdiff(unique(cl),cl[i]);if(!length(other))return(0);b<-min(vapply(other,function(g)mean(D[i,cl==g]),numeric(1)));if(max(a,b)==0)0 else (b-a)/max(a,b)},numeric(1)),na.rm=TRUE)
  }
  hc<-hclust(as.dist(D),method="ward.D2")
  metrics<-if(length(ks))map_dfr(ks,function(k){cl<-cutree(hc,k);tibble(k,silhouette=sil_one(cl),min_module_n=min(table(cl)))})else tibble(k=1L,silhouette=NA_real_,min_module_n=n)
  eligible<-metrics|>filter(min_module_n>=max(3L,ceiling(.03*n)));if(!nrow(eligible))eligible<-metrics|>slice(1)
  k<-eligible|>arrange(desc(silhouette),k)|>slice(1)|>pull(k);orig<-if(k==1)rep(1L,n)else cutree(hc,k)
  stable<-rep(1,n)
  if(k>1&&C4_MODULE_BOOT>0&&ncol(X)>1){
    set.seed(SEED+404);hit<-numeric(n)
    for(b in seq_len(C4_MODULE_BOOT)){
      cols<-sample(seq_len(ncol(X)),ncol(X),replace=TRUE);cb<-cutree(hclust(dist(X[,cols,drop=FALSE]),method="ward.D2"),k)
      tab<-table(cb,orig);mp<-apply(tab,1,function(v)as.integer(colnames(tab)[which.max(v)]));aligned<-mp[as.character(cb)]
      hit<-hit+(aligned==orig)
    }
    stable<-hit/C4_MODULE_BOOT
  }
  membership<-tibble(feature=rn,module=orig,module_stability=stable)|>
    left_join(sets$primary|>select(feature,primary_component,strict_YS,YS_model,FDR_disc,FDR_rep),by="feature")|>
    mutate(YS_core=coalesce(YS_model,FALSE)&module_stability>=C4_MODULE_STABILITY)
  prof<-p|>filter(feature%in%rn)|>left_join(membership|>select(feature,module,module_stability,YS_core,primary_component),by="feature")
  list(membership=membership,metrics=metrics|>mutate(selected_k=k),profile=prof,YS_core=membership|>filter(YS_core)|>pull(feature))
}

plot_supervised_atlas <- function(modules,sets,med,layer,outdir){
  mem<-modules$membership;pr<-modules$profile
  if(!nrow(mem)||!nrow(pr)){
    save_plot(blank_plot("LE8-supervised biomarker atlas","No replicated profile was clusterable"),"c4.Fig7.supervised_atlas.png",10,6,outdir=outdir)
    return(invisible(NULL))
  }
  top<-mem|>arrange(desc(YS_core),desc(module_stability),FDR_disc)|>slice_head(n=100)|>pull(feature)
  hm<-pr|>filter(feature%in%top)|>mutate(feature=factor(feature,levels=rev(top)),component=factor(component,levels=names.le8),z=cap(z_joint,8))
  pa<-ggplot(hm,aes(component,feature,fill=z))+geom_tile(color="white",linewidth=.15)+
    facet_grid(module~.,scales="free_y",space="free_y")+scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,name="replicated z")+
    labs(title=paste0("a. LE8-supervised ",layer," modules"),subtitle="Top replicated eight-pillar profiles; clustering is limited to the strongest 1,200 features",x=NULL,y=NULL)+theme_5c(7)+theme(axis.text.x=element_text(angle=35,hjust=1))
  pb<-ggplot(mem,aes(module_stability,factor(module),color=YS_core))+geom_jitter(height=.16,width=0,alpha=.45,size=1.4)+
    geom_vline(xintercept=C4_MODULE_STABILITY,linetype=2)+scale_color_manual(values=c(`FALSE`="grey65",`TRUE`="#D73027"))+
    labs(title="b. Perturbation-bootstrap module stability",x="Feature stability",y="Module",color="YS core")+theme_5c(9)
  funnel<-tibble(stage=factor(c("All profiled","YS model","YS strict","YS core"),levels=rev(c("All profiled","YS model","YS strict","YS core"))),
                 n=c(nrow(mem),sum(mem$YS_model,na.rm=TRUE),sum(mem$strict_YS,na.rm=TRUE),sum(mem$YS_core,na.rm=TRUE)))
  pc<-ggplot(funnel,aes(n,stage,fill=stage))+geom_col(width=.72)+geom_text(aes(label=n),hjust=-.15,fontface="bold")+
    scale_x_continuous(expand=expansion(mult=c(0,.18)))+guides(fill="none")+labs(title="c. Supervised selection funnel",x="Features",y=NULL)+theme_5c(9)
  md<-med|>filter(is.finite(indirect_beta))|>arrange(indirect_p)|>slice_head(n=20)|>mutate(label=factor(paste(feature,component,sep=" | "),levels=rev(paste(feature,component,sep=" | "))))
  pd<-if(!nrow(md))blank_plot("d. LE8–omic–disease paths","No estimable path")else ggplot(md,aes(indirect_beta,label,color=component))+
    geom_vline(xintercept=0,color="grey60")+geom_errorbar(aes(xmin=coalesce(indirect_lo,indirect_beta-1.96*indirect_se),xmax=coalesce(indirect_hi,indirect_beta+1.96*indirect_se)),orientation="y",width=.05)+geom_point()+
    scale_color_manual(values=cols_le8)+labs(title="d. Associational bridge estimates",x="Indirect log-hazard effect",y=NULL,color=NULL)+forest_theme(8)
  save_plot((pa|pb)/(pc|pd)+plot_layout(widths=c(1.35,1)),"c4.Fig7.supervised_atlas.png",17,14,outdir=outdir)
  save_plot(pc|pd,"c4.Fig9.selection_mediation.png",14,8,outdir=outdir)
}

plot_module_globe <- function(modules,outdir){
  d0<-as_tibble(modules$profile%||%tibble())
  if(!nrow(d0)||!all(c("YS_core","feature","module","z_joint","module_stability","primary_component")%in%names(d0))){save_plot(blank_plot("Stable LE8 module globe","No YS-core feature passed stability"),"c4.Fig8.network_globe.png",9,6,outdir=outdir);return(invisible(NULL))}
  d<-d0|>filter(YS_core)|>group_by(feature,module)|>slice_max(abs(z_joint),n=1,with_ties=FALSE)|>ungroup()
  if(!nrow(d)){save_plot(blank_plot("Stable LE8 module globe","No YS-core feature passed stability"),"c4.Fig8.network_globe.png",9,6,outdir=outdir);return(invisible(NULL))}
  mods<-sort(unique(d$module));cent<-tibble(module=mods,a=head(seq(0,2*pi,length.out=length(mods)+1),-1),x=.28*cos(a),y=.28*sin(a))
  d<-d|>group_by(module)|>mutate(a=head(seq(0,2*pi,length.out=n()+1),-1)+first(module)*.31,x=cos(a),y=sin(a))|>ungroup()|>left_join(cent|>select(module,mx=x,my=y),by="module")
  p<-ggplot(d)+geom_curve(aes(x=mx,y=my,xend=x,yend=y,color=primary_component,alpha=module_stability),curvature=.08,linewidth=.5)+
    geom_point(aes(x,y,size=abs(z_joint),fill=primary_component),shape=21,color="grey30")+geom_point(data=cent,aes(x,y),shape=23,size=5,fill="#FFD92F")+
    geom_text(data=cent,aes(x,y,label=paste0("M",module)),fontface="bold")+geom_text_repel(aes(x,y,label=feature),size=2.3,max.overlaps=35,seed=14)+
    scale_color_manual(values=cols_le8)+scale_fill_manual(values=cols_le8)+coord_equal(clip="off")+theme_void()+
    labs(title="Stable LE8-supervised biomarker globe",subtitle="Yellow diamonds are learned modules; outer nodes are stable YS-core biomarkers",color="Primary pillar",fill="Primary pillar",size="Profile strength",alpha="Stability")+
    theme(plot.title=element_text(face="bold",size=15),legend.position="bottom",plot.margin=margin(20,45,20,45))
  save_plot(p,"c4.Fig8.network_globe.png",13,11,outdir=outdir)
}

plot_c4 <- function(scan,sets,med,genetics,layer,outdir){
  topf<-sets$primary|>arrange(desc(strict_YS),desc(YS_model),FDR_disc)|>slice_head(n=55)|>pull(feature)
  hm<-scan|>filter(split=="discovery",feature%in%topf)|>mutate(component=factor(component,levels=names.le8),
    feature=factor(feature,levels=rev(topf)),zcap=cap(z,8))
  p1<-if(!nrow(hm))blank_plot(paste0("LE8-supervised ",layer," proxy map"))else ggplot(hm,aes(component,feature,fill=zcap))+
    geom_tile(color="white",linewidth=.18)+scale_fill_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0,name="Partial z")+
    labs(title=paste0("LE8-supervised ",layer," proxy map"),subtitle="Discovery associations; YS requires directionally concordant replication and pillar specificity",
      x=NULL,y=NULL)+theme_5c(8)+theme(axis.text.x=element_text(angle=35,hjust=1))
  save_plot(p1,"c4.Fig1.proxy_heatmap.png",10.5,11,outdir=outdir)

  flow<-sets$membership|>filter(set%in%c("YS","YSP_plus"))|>count(group=coalesce(group,"Other"),primary_component=coalesce(primary_component,"Disease-selected"),set)
  if(requireNamespace("ggalluvial",quietly=TRUE)&&nrow(flow)){
    p2<-ggplot(flow,aes(axis1=group,axis2=primary_component,axis3=set,y=n))+ggalluvial::geom_alluvium(aes(fill=primary_component),width=1/12,alpha=.75)+
      ggalluvial::geom_stratum(width=1/8,fill="grey96",color="grey55")+ggalluvial::stat_stratum(geom="text",aes(label=after_stat(stratum)),size=2.5)+
      scale_x_discrete(limits=c("Omic group","LE8 pillar","Proxy set"),expand=c(.08,.08))+scale_fill_manual(values=cols_le8,na.value="grey65",guide="none")+
      labs(title="Omic group → LE8 pillar → supervised set",subtitle="YS is LE8-supervised; YSP_plus adds disease-ranked features",y="Features",x=NULL)+theme_5c(10)
  }else if(nrow(flow))p2<-ggplot(flow,aes(primary_component,fct_reorder(group,n,sum),size=n,color=set))+geom_point(alpha=.8)+
    scale_size_continuous(range=c(2,10))+scale_color_manual(values=c(YS="#D95F02",YSP_plus="#1B9E77"))+
    labs(title="Omic groups assigned to LE8 pillars",x=NULL,y=NULL,color=NULL,size="Features")+theme_5c(10)
  else p2<-blank_plot("Omic groups assigned to LE8 pillars")
  save_plot(p2,"c4.Fig2.group_pillar_flow.png",13,7.5,outdir=outdir)

  le<-sets$YS_edges|>arrange(FDR_disc)|>slice_head(n=18)|>mutate(side="LE8")
  ge<-genetics|>filter(replicated)|>arrange(desc(bridge_score))|>slice_head(n=18)|>mutate(side="Genetics")
  if(nrow(le)||nrow(ge)){
    lfeat<-unique(le$feature);rfeat<-unique(ge$feature);lpil<-unique(le$component);rprs<-unique(ge$prs)
    make_nodes<-function(nms,type,x,lo=.05,hi=.95){
      nms<-unique(na.omit(as.character(nms)))
      if(!length(nms))return(tibble(name=character(),type=character(),x=numeric(),y=numeric()))
      tibble(name=nms,type=type,x=x,y=seq(lo,hi,length.out=length(nms)))
    }
    node<-bind_rows(make_nodes(lpil,"LE8",-2,.08,.92),make_nodes(lfeat,"omic_left",-1),
      tibble(name=Y,type="disease",x=0,y=.5),make_nodes(rfeat,"omic_right",1),make_nodes(rprs,"PRS",2,.08,.92))
    el<-le|>left_join(node|>filter(type=="LE8")|>select(component=name,x1=x,y1=y),by="component")|>
      left_join(node|>filter(type=="omic_left")|>select(feature=name,x2=x,y2=y),by="feature")
    er<-ge|>left_join(node|>filter(type=="PRS")|>select(prs=name,x1=x,y1=y),by="prs")|>
      left_join(node|>filter(type=="omic_right")|>select(feature=name,x2=x,y2=y),by="feature")
    toY<-bind_rows(node|>filter(type=="omic_left")|>transmute(x1=x,y1=y,x2=0,y2=.5),node|>filter(type=="omic_right")|>transmute(x1=x,y1=y,x2=0,y2=.5))
    p3<-ggplot()+geom_curve(data=el,aes(x=x1,y=y1,xend=x2,yend=y2,color=component),curvature=.08,linewidth=.65,
        arrow=grid::arrow(length=grid::unit(1.6,"mm")))+
      geom_curve(data=er,aes(x=x1,y=y1,xend=x2,yend=y2),curvature=-.08,color="#6A3D9A",linewidth=.65,
        arrow=grid::arrow(length=grid::unit(1.6,"mm")))+
      geom_curve(data=toY,aes(x=x1,y=y1,xend=x2,yend=y2),curvature=.06,color="grey68",alpha=.55,
        arrow=grid::arrow(length=grid::unit(1.5,"mm")))+
      geom_point(data=node,aes(x,y,shape=type),size=2.5)+geom_text(data=node,aes(x,y,label=name),size=2.4,nudge_y=.022,check_overlap=TRUE)+
      scale_color_manual(values=cols_le8,guide="none")+coord_cartesian(xlim=c(-2.25,2.25),ylim=c(0,1),clip="off")+theme_void()+
      labs(title="LE8 → omics → disease ← omics ← genetics",subtitle="Top replicated bridges; PRS associations are genetic anchors, not proof of mediation")+
      theme(plot.title=element_text(face="bold",size=15),plot.subtitle=element_text(color="grey35"),plot.margin=margin(12,25,12,25))
  }else p3<-blank_plot("LE8 → omics → disease ← omics ← genetics","No replicated bridge was available")
  save_plot(p3,"c4.Fig3.connection_bridge.png",15,8.5,outdir=outdir)

  evf<-unique(c(head(le$feature,25),head(ge$feature,25)))
  mat1<-sets$primary|>filter(feature%in%evf)|>transmute(feature,metric="LE8 proxy",value=z_disc,set=ifelse(strict_YS,"YS strict",ifelse(YS_model,"YS exploratory","Other")))
  mat2<-sets$membership|>filter(feature%in%evf)|>transmute(feature,metric="Disease",value=sign(disease_beta)*-log10(pmax(disease_p,1e-300)),set=coalesce(set,"Other"))
  mat3<-genetics|>filter(feature%in%evf,replicated)|>group_by(feature)|>slice_max(abs(z_disc),n=1,with_ties=FALSE)|>ungroup()|>transmute(feature,metric="Genetic anchor",value=z_disc,set=coalesce(set,"Other"))
  em<-bind_rows(mat1,mat2,mat3)|>mutate(value=cap(value,8),feature=factor(feature,levels=rev(evf)),metric=factor(metric,levels=c("LE8 proxy","Disease","Genetic anchor")))
  p4<-if(!nrow(em))blank_plot("Connection evidence matrix")else ggplot(em,aes(metric,feature))+
    geom_point(aes(size=abs(value),fill=value,shape=set),color="grey40")+scale_fill_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0)+
    scale_size_continuous(range=c(1.5,7))+labs(title="Connection evidence matrix",subtitle="Signed evidence across modifiable, proximal disease, and genetic dimensions",
      x=NULL,y=NULL,fill="Signed evidence",size="|Evidence|",shape="Set")+theme_5c(9)+theme(legend.position="bottom")
  save_plot(p4,"c4.Fig4.connection_evidence.png",11.5,9.5,outdir=outdir)

  md<-med|>filter(is.finite(indirect_beta))|>arrange(indirect_p)|>slice_head(n=35)|>mutate(label=paste(feature,component,sep=" | "),label=factor(label,levels=rev(label)))
  p5<-if(!nrow(md))blank_plot("Associational LE8–omic–disease paths")else ggplot(md,aes(indirect_beta,label,color=component))+
    geom_vline(xintercept=0,color="grey60")+geom_errorbar(aes(xmin=coalesce(indirect_lo,indirect_beta-1.96*indirect_se),xmax=coalesce(indirect_hi,indirect_beta+1.96*indirect_se)),orientation="y",width=.08)+
    geom_point(aes(size=pmin(abs(prop_mediated),1)),alpha=.85)+scale_color_manual(values=cols_le8)+
    labs(title="Associational LE8–omic–disease paths",subtitle="Same-baseline LE8 and omics: indirect effects do not establish temporal mediation",
      x="Indirect log-hazard effect",y=NULL,color="LE8 pillar",size="|Proportion mediated|")+forest_theme(8)+theme(legend.position="bottom")
  save_plot(p5,"c4.Fig5.mediation_forest.png",12,10,outdir=outdir)

  if(nrow(md)){
    pb<-ggplot(md,aes(total_beta,direct_beta,color=component))+geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+
      geom_point(aes(size=abs(prop_mediated)),alpha=.8)+scale_color_manual(values=cols_le8)+
      labs(title="a. Total versus direct effects",x="Total effect",y="Direct effect",color=NULL,size="|Proportion mediated|")+theme_5c(9)
    pc<-ggplot(md,aes(100*cap(prop_mediated,1),-log10(pmax(indirect_p,1e-30)),color=component))+geom_vline(xintercept=0,color="grey60")+
      geom_point(alpha=.8)+scale_color_manual(values=cols_le8)+labs(title="b. Proportion and indirect-effect evidence",x="Percent mediated (capped)",y=expression(-log[10](P[indirect])),color=NULL)+theme_5c(9)
    p6<-pb|pc
  }else p6<-blank_plot("Mediation diagnostics")
  save_plot(p6,"c4.Fig6.mediation_diagnostics.png",14,7.5,outdir=outdir)
}

# Retained only to reproduce the pre-genetics C4 figure set.  Keeping the same
# name used to override the six-argument publication function above and made
# every current C4 run fail with an unused `genetics` argument.
plot_c4_legacy <- function(scan,sets,med,layer,outdir) {
  topf<-sets$primary|>arrange(desc(YS),FDR_disc)|>slice_head(n=80)|>pull(feature)
  hm<-scan|>filter(split=="discovery",feature%in%topf)|>mutate(component=factor(component,levels=names.le8),feature=factor(feature,levels=rev(topf)),zcap=cap(z,8))
  p1<-if(!nrow(hm)) blank_plot(paste0("LE8-supervised ",layer," proxy map"),"No discovery association was available") else ggplot(hm,aes(component,feature,fill=zcap))+geom_tile(color="white",linewidth=.15)+scale_fill_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0,name="z")+
    labs(title=paste0("LE8-supervised ",layer," proxy map"),subtitle="Discovery associations; proxy status requires replication and pillar specificity",x=NULL,y=NULL)+theme_5c(8)+theme(axis.text.x=element_text(angle=35,hjust=1))
  save_plot(p1,"c4.Fig1.proxy_heatmap.png",10,11,outdir=outdir)

  flow<-sets$membership|>filter(set%in%c("YS","YSP_plus"))|>count(group=coalesce(group,"Other"),primary_component=coalesce(primary_component,"Disease-selected"),set)
  if(requireNamespace("ggalluvial",quietly=TRUE)&&nrow(flow)){
    p2<-ggplot(flow,aes(axis1=group,axis2=primary_component,axis3=set,y=n))+ggalluvial::geom_alluvium(aes(fill=primary_component),width=1/12,alpha=.75)+ggalluvial::geom_stratum(width=1/8,fill="grey94",color="grey55")+ggalluvial::stat_stratum(geom="text",aes(label=after_stat(stratum)),size=2.6)+scale_x_discrete(limits=c("Omic group","LE8 pillar","Proxy set"),expand=c(.08,.08))+scale_fill_manual(values=cols_le8,na.value="grey65",guide="none")+labs(title="Omic group → LE8 pillar → supervised set",y="Number of features",x=NULL)+theme_5c(10)
  } else if(nrow(flow)) {
    p2<-ggplot(flow,aes(primary_component,fct_reorder(group,n,sum),size=n,color=set))+geom_point(alpha=.8)+scale_size_continuous(range=c(2,10))+scale_color_manual(values=c(YS="#D95F02",YSP_plus="#1B9E77"))+labs(title="Omic groups assigned to LE8 pillars",x=NULL,y=NULL,color=NULL,size="Features")+theme_5c(10)+theme(axis.text.x=element_text(angle=35,hjust=1))
  } else p2<-blank_plot("Omic groups assigned to LE8 pillars","No proxy-set flow could be formed")
  save_plot(p2,"c4.Fig2.group_pillar_flow.png",12,7,outdir=outdir)

  edges<-sets$YS_edges|>arrange(FDR_disc)|>slice_head(n=100)
  if(nrow(edges)){
    pillars<-tibble(name=unique(edges$component),type="pillar",angle=seq(0,2*pi,length.out=n_distinct(edges$component)+1)[-(n_distinct(edges$component)+1)],r=.35)
    feats<-tibble(name=unique(edges$feature),type="feature",angle=seq(0,2*pi,length.out=n_distinct(edges$feature)+1)[-(n_distinct(edges$feature)+1)]+.12,r=1)
    nodes<-bind_rows(pillars,feats)|>mutate(x=r*cos(angle),y=r*sin(angle));ee<-edges|>left_join(nodes|>select(component=name,x1=x,y1=y),by="component")|>left_join(nodes|>select(feature=name,x2=x,y2=y),by="feature")
    p3<-ggplot()+geom_curve(data=ee,aes(x=x1,y=y1,xend=x2,yend=y2,color=component,alpha=pmin(abs(z_disc)/8,1)),curvature=.12,linewidth=.55)+geom_point(data=nodes,aes(x,y,shape=type),size=2.5)+ggrepel::geom_text_repel(data=nodes,aes(x,y,label=name),size=2.3,max.overlaps=35,seed=11)+scale_color_manual(values=cols_le8)+scale_alpha(range=c(.2,.9),guide="none")+coord_equal()+theme_void()+labs(title="LE8 pillar–proxy module network",color="Pillar")+theme(plot.title=element_text(face="bold"),legend.position="bottom")
  } else p3<-blank_plot("LE8 pillar–proxy module network")
  save_plot(p3,"c4.Fig3.module_network.png",11,9,outdir=outdir)

  # Globe/orbit view: the radius encodes pillar specificity and the angle encodes the primary LE8 component.
  orbit<-sets$membership|>filter(set%in%c("YS","YSP_plus"))|>mutate(comp=factor(coalesce(primary_component,"disease"),levels=c(names.le8,"disease")),angle=2*pi*(as.numeric(comp)-1)/nlevels(comp)+runif(n(),-.18,.18),radius=ifelse(set=="YS",.65+.35*coalesce(specificity,0),1.18),x=radius*cos(angle),y=radius*sin(angle),lab=ifelse(min_rank(disease_p)<=25|set=="YS"&min_rank(FDR_disc)<=20,feature,NA_character_))
  p4<-if(!nrow(orbit)) blank_plot("LE8–omics relationship globe","No supervised or plus feature was available") else ggplot(orbit,aes(x,y))+annotate("path",x=cos(seq(0,2*pi,length.out=300)),y=sin(seq(0,2*pi,length.out=300)),color="grey80")+geom_segment(aes(x=0,y=0,xend=x,yend=y,color=primary_component),alpha=.12)+geom_point(aes(color=primary_component,shape=set,size=-log10(pmax(disease_p,1e-30))),alpha=.8)+ggrepel::geom_text_repel(aes(label=lab),size=2.4,max.overlaps=30,seed=12,na.rm=TRUE)+scale_color_manual(values=c(cols_le8,disease="grey35"),na.value="grey60")+scale_size_continuous(range=c(1.5,5),name="Disease evidence")+coord_equal()+theme_void()+labs(title="LE8–omics relationship globe",color="Primary pillar",shape="Set")+theme(plot.title=element_text(face="bold"),legend.position="bottom")
  save_plot(p4,"c4.Fig4.relationship_globe.png",11,9,outdir=outdir)

  if(nrow(med)){
    wheel<-med|>filter(is.finite(prop_mediated))|>arrange(indirect_p)|>slice_head(n=60)|>mutate(wheel_label=paste(feature,component,sep=" | "),wheel_label=factor(wheel_label,levels=wheel_label),id=row_number(),angle=90-360*(id-.5)/n(),hjust=ifelse(angle< -90,1,0),angle=ifelse(angle< -90,angle+180,angle),value=100*cap(prop_mediated,1))
    p5<-if(!nrow(wheel)) blank_plot("Associational mediation wheel","Mediation effects were estimable, but no finite proportion mediated was available") else ggplot(wheel,aes(id,abs(value),fill=component))+geom_col(width=.8)+geom_hline(yintercept=0)+geom_text(aes(y=abs(value)+3,label=wheel_label,angle=angle,hjust=hjust),size=2,na.rm=TRUE)+coord_polar(clip="off")+scale_fill_manual(values=cols_le8)+ylim(-10,max(abs(wheel$value),na.rm=TRUE)+15)+labs(title="Associational mediation wheel",subtitle="Percent mediated is shown for screened LE8 → omic → disease paths",fill="LE8 pillar")+theme_void(10)+theme(plot.title=element_text(face="bold",hjust=.5),plot.subtitle=element_text(hjust=.5),legend.position="bottom",plot.margin=margin(30,55,30,55))
    save_plot(p5,"c4.Fig5.mediation_wheel.png",11,10,outdir=outdir)
    md<-med|>filter(is.finite(indirect_beta),is.finite(total_beta),is.finite(direct_beta))|>arrange(indirect_p)|>slice_head(n=50)|>mutate(label=paste(feature,component,sep=" | "),label=factor(label,levels=rev(label)))
    if(!nrow(md)) {
      save_plot(blank_plot("Mediation diagnostics","No finite mediation diagnostic row was available"),"c4.Fig6.mediation_diagnostics.png",9,5,outdir=outdir)
    } else {
      pa<-ggplot(md,aes(indirect_beta,label))+geom_vline(xintercept=0,color="grey60")+geom_errorbar(aes(xmin=coalesce(indirect_lo,indirect_beta-1.96*indirect_se),xmax=coalesce(indirect_hi,indirect_beta+1.96*indirect_se)),orientation="y",width=0)+geom_point(aes(color=component),size=2)+scale_color_manual(values=cols_le8)+labs(title="A. Indirect effects",x="Indirect log-hazard effect",y=NULL,color=NULL)+forest_theme(9)
      pb<-ggplot(md,aes(total_beta,direct_beta,color=component))+geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+geom_point(aes(size=abs(prop_mediated)),alpha=.8)+scale_color_manual(values=cols_le8)+labs(title="B. Total versus direct effects",x="Total effect",y="Direct effect",color=NULL,size="|Proportion mediated|")+theme_5c(10)
      pc<-ggplot(md,aes(100*cap(prop_mediated,1),-log10(pmax(indirect_p,1e-30)),color=component))+geom_vline(xintercept=0,color="grey60")+geom_point(alpha=.8)+scale_color_manual(values=cols_le8)+labs(title="C. Proportion and evidence",x="Percent mediated (plot capped)",y=expression(-log[10](P[indirect])),color=NULL)+theme_5c(10)
      save_plot(pa | (pb/pc),"c4.Fig6.mediation_diagnostics.png",14,10,outdir=outdir)
    }
  } else {
    save_plot(blank_plot("Associational mediation wheel","No screened YS path could be estimated"),"c4.Fig5.mediation_wheel.png",9,5,outdir=outdir)
    save_plot(blank_plot("Mediation diagnostics","No screened YS path could be estimated"),"c4.Fig6.mediation_diagnostics.png",9,5,outdir=outdir)
  }
}

run_c4_layer <- function(layer=c("protein","metabolite")) {
  # LE8_REVISION_UPDATES: guarded caches, definitions and additions.
  layer <- match.arg(layer)
  le8_load_revision(layer, "c4_connect")
  .le8_review_env <- environment()
  on.exit(le8_finish_revision(layer, "c4_connect", .le8_review_env), add=TRUE)

  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir);rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE);cache<-file.path(rawdir,"c4.res.rds")
  if(cache_valid(cache)){
    old<-tryCatch(readRDS(cache),error=function(e)NULL)
    if(!is.null(old)&&identical(old$meta$code_version%||%NA_character_,C4_CODE_VERSION)){
      cache_message(paste0("C4/",layer),cache);return(old)
    }
    message("C4/",layer,": cached result predates ",C4_CODE_VERSION,"; recomputing")
  }
  c1f<-file.path(le8_job_dir(outdir,"c1_correlate"),"c1.res.rds");if(!file.exists(c1f))stop("Run C1 first.",call.=FALSE);c1<-readRDS(c1f);disease<-(c1$association%||%c1$pwas_incident%||%c1$MWAS)|>as_tibble();if(!"beta"%in%names(disease))disease<-disease|>mutate(beta=safe_log(estimate))
  biom<-if(layer=="protein")read_prot() else read_met();features<-setdiff(names(biom),"eid");ann<-layer_annotation(layer,features)
  all0<-read_all();prs_vars<-find_prs_vars(all0,Y,4)
  need<-unique(c("eid","ethnic.c",vars.basic,vars.le8,prs_vars,"birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y)))
  dat0<-all0[,intersect(need,names(all0)),drop=FALSE]|>filter_analysis_cohort()|>make_outcome(Y);rm(all0);invisible(gc())
  tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e");bvar<-paste0(Y,".b2e")
  comps<-intersect(vars.le8,names(dat0));covs<-intersect(vars.basic,names(dat0))
  dat0<-dat0[complete.cases(dat0[,comps,drop=FALSE]),,drop=FALSE]|>stratified_sample(evar,MAX_N)
  # Sample participants before joining the wide omics matrix; this materially reduces C4 peak memory.
  dat<-inner_join(dat0,biom|>filter(eid%in%dat0$eid),by="eid")
  fold<-stratified_split(dat,evar)
  scan_cache<-file.path(rawdir,"c4.LE8_proxy_scan.rds");scan_stage<-read_stage_cache(scan_cache)
  if(is.null(scan_stage)){
    disc<-proxy_scan(dat[fold=="discovery",],features,comps,covs,"discovery")
    rep<-proxy_scan(dat[fold=="replication",],features,comps,covs,"replication")
    scan_stage<-list(discovery=disc,replication=rep);write_stage_cache(scan_stage,scan_cache)
  }else{disc<-scan_stage$discovery;rep<-scan_stage$replication;message("C4/",layer,": reuse LE8 proxy-scan stage")}
  if(!nrow(disc)||!nrow(rep)){
    reason<-paste0("No discovery/replication proxy scan was estimable (discovery rows=",nrow(disc),
      "; replication rows=",nrow(rep),"). Check LE8 variables and complete-case sample size.")
    message("C4/",layer,": ",reason," Writing blank panels and continuing.")
    suffix<-if(layer=="protein")c("proxy_heatmap","group_pillar_flow","connection_bridge",
      "connection_evidence","mediation_forest","mediation_diagnostics","supervised_atlas",
      "network_globe","selection_mediation","state_network_remodeling","state_network_edges") else c("proxy_heatmap","group_pillar_flow",
      "module_network","relationship_globe","mediation_wheel","mediation_diagnostics",
      "supervised_atlas","network_globe","selection_mediation","state_network_remodeling","state_network_edges")
    for(i in seq_along(suffix))save_plot(blank_plot(paste0("C4 Figure ",i),reason),
      paste0("c4.Fig",i,".",suffix[[i]],".png"),10,6,outdir=outdir)
    lists<-list(YS_by_component=list(),YS_all=character(),YS_strict=character(),YS_core=character(),
      YSP_plus=character(),YSP_all=character(),NS=character(),index=tibble(),
      genetic_anchors=character(),PRS_variables=prs_vars)
    status<-tibble(status="unavailable",detail=reason,discovery_rows=nrow(disc),replication_rows=nrow(rep))
    out<-list(meta=module_meta(layer,extra=list(code_version=C4_CODE_VERSION,status="unavailable")),
      status=status,scan=bind_rows(disc,rep),primary=tibble(),membership=tibble(),YS_edges=tibble(),
      modules=list(membership=tibble(),metrics=tibble(),YS_core=character()),genetic_scan=tibble(),
      genetic_edges=tibble(),PRS_variables=prs_vars,mediation=tibble(),state_network=list(),lists=lists)
    write_raw_csv(status,"c4.status.csv",rawdir);saveRDS(out,cache,compress="xz")
    write_xlsx2(list(status=status,all_LE8_associations=out$scan),"c4.out.xlsx")
    finalize_outputs(LE8_JOB,outdir);return(out)
  }
  scan<-bind_rows(disc,rep);sets<-make_proxy_sets(disc,rep,disease,ann)
  modules<-supervised_modules(disc,rep,sets)
  write_raw_csv(scan,"c4.LE8_feature_associations.csv",rawdir);write_raw_csv(sets$primary,"c4.primary_pillar_assignment.csv",rawdir);write_raw_csv(sets$membership,"c4.proxy_membership_YS_YSP_NS.csv",rawdir);write_raw_csv(sets$YS_edges,"c4.YS_edges.csv",rawdir)
  write_raw_csv(modules$membership,"c4.supervised_module_membership.csv",rawdir);write_raw_csv(modules$metrics,"c4.supervised_module_selection.csv",rawdir)
  genetic_cache<-file.path(rawdir,"c4.PRS_proxy_scan.rds");genetic_stage<-read_stage_cache(genetic_cache)
  if(is.null(genetic_stage)){
    gdisc<-if(length(prs_vars))proxy_scan(dat[fold=="discovery",],features,prs_vars,covs,"discovery")else tibble()
    grep0<-if(length(prs_vars))proxy_scan(dat[fold=="replication",],features,prs_vars,covs,"replication")else tibble()
    genetic_stage<-list(discovery=gdisc,replication=grep0);write_stage_cache(genetic_stage,genetic_cache)
  }else{gdisc<-genetic_stage$discovery;grep0<-genetic_stage$replication;message("C4/",layer,": reuse PRS proxy-scan stage")}
  gscan<-bind_rows(gdisc,grep0);gedges<-make_genetic_edges(gdisc,grep0,disease,sets$membership)
  write_raw_csv(gscan,"c4.PRS_feature_associations.csv",rawdir);write_raw_csv(gedges,"c4.genetic_omic_disease_bridges.csv",rawdir)
  pairs<-sets$YS_edges|>left_join(disease|>select(feature=term,disease_p=p.value),by="feature")|>arrange(desc(strict_YS),FDR_disc,disease_p)|>slice_head(n=MED_TOP)
  med_dat<-stratified_sample(dat,evar,MED_MAX_N)
  mediation_cache<-file.path(rawdir,"c4.mediation_stage.rds");med<-read_stage_cache(mediation_cache)
  if(is.null(med)){
    med<-if(!nrow(pairs))tibble() else map_dfr(seq_len(nrow(pairs)),function(i){
      rw<-pairs[i,,drop=FALSE]
      mediation_one(med_dat,rw$component_var[[1]],rw$feature[[1]],unique(c(covs,setdiff(comps,rw$component_var[[1]]))),tvar,evar,MED_BOOT)
    })
    if(nrow(med))med<-med|>mutate(FDR_indirect=p.adjust(indirect_p,"BH"))|>arrange(indirect_p)
    write_stage_cache(med,mediation_cache)
  }else message("C4/",layer,": reuse mediation stage")
  write_raw_csv(med,"c4.mediation_all.csv",rawdir);plot_c4(scan,sets,med,gedges,layer,outdir)
  plot_supervised_atlas(modules,sets,med,layer,outdir);plot_module_globe(modules,outdir)
  state_network<-run_c4_state_network(dat,disease,sets,modules,features,covs,tvar,evar,bvar,
    layer,rawdir,outdir)
  write_raw_csv(state_network$status%||%tibble(),"c4.state_network_status.csv",rawdir)
  write_raw_csv(state_network$state_counts%||%tibble(),"c4.state_network_counts.csv",rawdir)
  write_raw_csv(state_network$edges%||%tibble(),"c4.state_network_edges.csv",rawdir)
  write_raw_csv(state_network$hubs%||%tibble(),"c4.state_network_hubs.csv",rawdir)
  lists<-list(YS_by_component=split(sets$YS_edges$feature,sets$YS_edges$component),YS_all=sets$YS,YS_strict=sets$YS_strict,YS_core=modules$YS_core,YSP_plus=sets$YSP_plus,YSP_all=sets$YSP,NS=sets$NS,index=sets$membership,
              genetic_anchors=gedges|>filter(replicated)|>pull(feature)|>unique(),PRS_variables=prs_vars)
  saveRDS(lists,file.path(rawdir,paste0("c4.",layer,"_lists.rds")),compress="xz");if(layer=="protein")saveRDS(lists,file.path(rawdir,"c4.protein_lists.rds"),compress="xz")
  out<-list(meta=module_meta(layer,extra=list(code_version=C4_CODE_VERSION,status="ok",scan_N=nrow(dat),discovery_N=sum(fold=="discovery"),replication_N=sum(fold=="replication"),mediation_N=nrow(med_dat),YS_strict_n=length(sets$YS_strict),YS_model_n=length(sets$YS),PRS_n=length(prs_vars))),
            scan=scan,primary=sets$primary,membership=sets$membership,YS_edges=sets$YS_edges,modules=modules,
            genetic_scan=gscan,genetic_edges=gedges,PRS_variables=prs_vars,mediation=med,
            state_network=state_network,lists=lists)
  saveRDS(out,cache,compress="xz");write_xlsx2(list(primary_assignment=sets$primary,proxy_membership=sets$membership,YS_edges=sets$YS_edges,
    supervised_modules=modules$membership,module_selection=modules$metrics,all_LE8_associations=scan,PRS_associations=gscan,genetic_omic_bridges=gedges,mediation=med,
    state_network_status=state_network$status%||%tibble(),state_counts=state_network$state_counts%||%tibble(),
    state_network_edges=state_network$edges%||%tibble(),state_network_hubs=state_network$hubs%||%tibble()),
    "c4.out.xlsx");finalize_outputs(LE8_JOB,outdir);out
}

if(prot_DO)run_c4_layer("protein")
if(met_DO)run_c4_layer("metabolite")
