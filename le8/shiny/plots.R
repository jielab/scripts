# One reactive figure for the selected table / assay. No 3000-file plot loop.
deviation_task_label <- function(d) {
  labels<-c(objective_plaque='客观冠脉病变',coronary_calcium='冠脉钙化',known_cad_state='已知CAD状态（探索性）')
  task<-if('task'%in%names(d))unique(as.character(d$task))else character()
  if(length(task)==1L&&!is.na(task)&&task%in%names(labels))unname(labels[task])else'记录状态/结局状态'
}
make_plot <- function(d,kind,title,x=NULL,y=NULL,group=NULL) {
  if(!nrow(d))return(empty_plot())
  base<-theme_bw(base_size=13)+theme(plot.title=element_text(size=13),legend.position='bottom')
  if(kind=='temporal') {
    tc<-if('time'%in%names(d))'time' else 'landmark_years'
    fc<-intersect(c('feature','term','protein'),names(d));d$feature_label<-if(length(fc))as.character(d[[fc[1]]])else'result'
    p<-ggplot(d,aes(x=.data[[tc]],y=beta,colour=feature_label))+geom_point()+geom_line()+geom_hline(yintercept=0,linetype=2)+labs(x=tc,y='Reported beta',colour='分子',title=title)+base
    if('effect_measure'%in%names(d))p<-p+facet_wrap(~effect_measure,scales='free')
    return(p)
  }
  if(kind=='proxy_reconstruction')return(ggplot(d,aes(component,delta_R2_vs_NS,colour=model))+geom_point()+geom_hline(yintercept=0,linetype=2)+coord_flip()+labs(x=NULL,y='Held-out ΔR² versus NS',title=title)+base)
  if(kind=='coloc') {
    pc<-intersect(c('PP_H4','PP.H4.abf','PP.H4'),names(d))[1];fc<-intersect(c('feature','assay','protein','gene','locus'),names(d))
    d$label<-if(length(fc))as.character(d[[fc[1]]]) else as.character(seq_len(nrow(d)))
    if(nrow(d)>80)return(ggplot(d,aes(x=.data[[pc]]))+geom_histogram(bins=30)+labs(x='Regional PP.H4 (not SNP posterior conditional on H4)',title=title)+base)
    return(ggplot(d,aes(x=.data[[pc]],y=label))+geom_point()+labs(x='Regional PP.H4',y=NULL,title=title)+base)
  }
  if(kind=='confidence'&&all(c('mean_confidence','Brier_IPCW','model')%in%names(d)))
    return(ggplot(d,aes(mean_confidence,Brier_IPCW,colour=model))+geom_point()+geom_line()+labs(x='事先 confidence（参考支持度，不是正确概率）',y='留出集 IPCW Brier',title=title)+base)
  if(kind=='roc')return(ggplot(d,aes(FPR,TPR,colour=model))+geom_line()+geom_abline(slope=1,intercept=0,linetype=2)+labs(title=title)+base)
  if(kind=='decision')return(ggplot(d,aes(threshold,net_benefit,colour=model))+geom_line()+geom_hline(yintercept=0,linetype=2)+labs(x='风险阈值',y='Net benefit',title=title)+base)
  if(kind=='deviation_state') {
    if(!'feature'%in%names(d))return(empty_plot())
    return(ggplot(d,aes(feature,mean_deviation,colour=factor(state)))+geom_point(size=3)+geom_hline(yintercept=0,linetype=2)+coord_flip()+labs(colour=deviation_task_label(d),x=NULL,y='遗传校正偏离（参考残差 SD）',title=title)+base)
  }
  if(kind=='deviation_histogram') {
    d<-d[is.finite(d$lower)&is.finite(d$upper),,drop=FALSE];d$mid<-(d$lower+d$upper)/2
    return(ggplot(d,aes(mid,N,colour=factor(state),group=interaction(feature,state)))+geom_line(na.rm=TRUE)+geom_point(na.rm=TRUE)+labs(x='遗传校正偏离 Z',y='人数（小于10的单元格未显示）',colour=deviation_task_label(d),title=title)+base)
  }
  if(kind=='pgs_capture')return(ggplot(d,aes(own_PGS_partial_R2))+geom_histogram(bins=40)+geom_vline(xintercept=0,linetype=2)+labs(x='OOF own-PGS partial R²',y='Assays',title=title)+base)
  if(kind=='coverage') {
    if(!'subset'%in%names(d))return(empty_plot())
    d<-d[grepl('^support_quality:',d$subset),,drop=FALSE]
    return(ggplot(d,aes(coverage,AUC_IPCW,colour=model))+geom_line()+geom_point()+labs(x='实际覆盖率',y='IPCW AUC',title=title)+base)
  }
  if(kind=='calibration')return(ggplot(d,aes(predicted,observed,colour=model))+geom_point()+geom_line()+geom_abline(slope=1,intercept=0,linetype=2)+labs(title=title)+base)
  if(kind=='contrasts') {
    d$label<-paste(d$contrast,if('model'%in%names(d))d$model else '',sep=' | ')
    return(ggplot(d,aes(delta,label))+geom_segment(aes(x=lower,xend=upper,yend=label))+geom_point()+geom_vline(xintercept=0,linetype=2)+facet_wrap(~metric,scales='free_x')+labs(x='配对差值：模型 − 对照',y=NULL,title=title)+base)
  }
  if(kind=='concordance') {
    p<-ggplot(d,aes(pgs_full_beta,observed_beta))+geom_point(alpha=.55)+geom_vline(xintercept=0,linetype=2)+geom_hline(yintercept=0,linetype=2)+labs(x='分子 PGS 的疾病关联 β',y='实测分子的疾病关联 β',title=title)+base
    if('analysis'%in%names(d))p<-p+facet_wrap(~analysis,scales='free')
    return(p)
  }
  if(kind=='effect'&&'beta'%in%names(d)) {
    labelcols<-intersect(c('feature','analysis','series','scope','model','component','term'),names(d))
    if(!length(labelcols))labelcols<-names(d)[1]
    d$label<-apply(d[,labelcols,drop=FALSE],1,paste,collapse=' | ')
    se<-intersect(c('std.error','se','SE'),names(d))
    if(!all(c('lo','hi')%in%names(d))&&length(se)){d$lo<-d$beta-1.96*d[[se[1]]];d$hi<-d$beta+1.96*d[[se[1]]]}
    if(nrow(d)>80) {
      pc<-intersect(c('p.value','p'),names(d));if(!length(pc))return(empty_plot('Select one molecule for the effect plot.'))
      d$minus_logp<--log10(pmax(as.numeric(d[[pc[1]]]),1e-300))
      return(ggplot(d,aes(beta,minus_logp))+geom_point(alpha=.5)+labs(x='Reported β; filter analysis/units before comparing',y='−log10(P)',title=title)+base)
    }
    p<-ggplot(d,aes(beta,label))+geom_point()+geom_vline(xintercept=0,linetype=2)+labs(x='Reported β (not pooled across effect measures)',y=NULL,title=title)+base
    if(all(c('covariates_before','covariates_after')%in%names(d)))p<-p+labs(x='调整后平均偏离差 β（状态1 − 状态0）',subtitle=deviation_task_label(d))
    if(all(c('lo','hi')%in%names(d)))p<-p+geom_segment(aes(x=lo,xend=hi,yend=label))
    if('effect_measure'%in%names(d))p<-p+facet_wrap(~effect_measure,scales='free')
    return(p)
  }
  if(kind=='models')return(ggplot(d,aes(AUC,model))+geom_point(size=3)+labs(title=title)+base)
  if(!is.null(x)&&!is.null(y)&&x%in%names(d)&&y%in%names(d)) {
    p<-ggplot(d,aes(x=.data[[x]],y=.data[[y]]))
    if(!is.null(group)&&nzchar(group)&&group%in%names(d))p<-p+aes(colour=factor(.data[[group]]))
    return(p+geom_point(alpha=.7)+labs(title=title)+base)
  }
  empty_plot('This table has no automatic plot. Choose numeric X/Y, or another result table.')
}
