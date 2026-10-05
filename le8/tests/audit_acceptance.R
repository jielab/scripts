# A01--A07 acceptance cases use the actual parsed production function bodies.
# Loading definitions avoids sourcing the full UKB pipeline during synthetic tests.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(purrr); library(data.table)})
script <- sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1])
root <- dirname(dirname(normalizePath(script)))
load_functions <- function(file, only=NULL) {
	walk <- function(x) {
		if (missing(x)) return(invisible(NULL))
		if (!is.call(x) && !is.expression(x)) return(invisible(NULL))
		if (is.call(x) && identical(x[[1]],as.name('<-')) && length(x)==3L && is.symbol(x[[2]]) &&
			is.call(x[[3]]) && identical(x[[3]][[1]],as.name('function'))) {
			if(is.null(only) || as.character(x[[2]]) %in% only) eval(x,envir=.GlobalEnv); return(invisible(NULL))
		}
		for (item in as.list(x)) walk(item)
	}
	walk(parse(file))
}
for (file in c('0.common.R','c1.correlate.R','c2.cause.R','c3.coloc.R','c4.connect.R')) load_functions(file.path(root,'f',file))
load_functions(file.path(root,"f/final.R"),"make_outer_split")
if (exists("split",envir=.GlobalEnv,inherits=FALSE)) rm(split,envir=.GlobalEnv)
.le8_identity_hashes <- new.env(parent=emptyenv())
.le8_sumstat_qc <- new.env(parent=emptyenv()); .le8_sumstat_qc$rows <- list()
args <- commandArgs(TRUE)
if (length(args) && args[1]=='--outcomes') {
	d <- as.data.frame(fread(args[2],colClasses='character',na.strings=''))
	z <- t2e(d,NA,'fod_icd10_cvd_cad','birth_date','date_attend','date_lost','date_death',le8_follow_end())
	cat(jsonlite::toJSON(list(prevalent=as.integer(as.Date(d$fod_icd10_cvd_cad)<=as.Date(d$date_attend) & !is.na(d$fod_icd10_cvd_cad)),event=z$Yt2e,time=z$t2e),na='null',digits=NA))
	quit(status=0)
}
Sys.unsetenv(c('LE8_GROUP_FILE','LE8_GROUP_COLUMN','PGS_GROUP_COLUMN','LE8_OUTER_ROSTER','LE8_GWAS_MANIFEST','LE8_EVIDENCE_LEVEL','C3_H4','C3_MR_FDR','C4_LE8_MEASURE_MAP'))
SEED <- 2026L; C4_MODULE_MAX <- 1200L; C4_MODULE_K_MAX <- 5L; C4_MODULE_BOOT <- 0L; C4_MODULE_STABILITY <- .7
scratch <- tempfile('le8-R-acceptance-',tmpdir='/tmp');dir.create(scratch)
passed <- 0L
check <- function(name,code) {
	tryCatch(force(code),error=function(e) stop(name,': ',conditionMessage(e),call.=FALSE))
	passed <<- passed+1L; cat('PASS ',name,'\n',sep='')
}
expect_error <- function(code) { e<-tryCatch({force(code);NULL},error=identity);stopifnot(inherits(e,'error'));invisible(e) }

check('A01 shared families, string IDs, conflict and outer roster',{
	d <- data.frame(eid=c('00001','00002','00003','00004'))
	f <- file.path(scratch,'groups.csv');fwrite(data.frame(eid=d$eid,group=c('001','001','002','002')),f)
	Sys.setenv(LE8_GROUP_FILE=f)
	stopifnot(identical(le8_participant_groups(d),c('001','001','002','002')))
	r <- file.path(scratch,'roster.csv');fwrite(data.frame(eid=d$eid,role=c('test','test','training','training')),r)
	stopifnot(identical(le8_outer_roles(d,r),c('validation','validation','training','training')))
	fwrite(data.frame(eid=d$eid,role=c('test','training','training','training')),r);expect_error(le8_outer_roles(d,r))
	Sys.setenv(LE8_GROUP_COLUMN='family');expect_error(le8_participant_groups(transform(d,family='wrong')))
	Sys.unsetenv(c('LE8_GROUP_FILE','LE8_GROUP_COLUMN'))
	d<-data.frame(eid=paste0('x',1:200),.group=rep(paste0('family',1:100),each=2))
	s<-make_outer_split(d,test_frac=.3,seed=12);stopifnot(sum(s=='validation')==60L,all(vapply(split(s,d$.group),function(x) length(unique(x)),integer(1))==1L))
})
check('A02 complete replacement excludes raw/score/derived design columns',{
	map<-c4_le8_measure_map(); measured<-c('smoke.pts','sbp','hba1c_ngsp','bmi','nonhdl',unique(map$component),'le8.pts','weight')
	r<-c4_replacement_background(c('age','sex'),measured,'disease_prs',unique(map$component),map)
	stopifnot(setequal(r$background,c('age','sex','disease_prs')))
	set.seed(1);d<-as.data.frame(setNames(replicate(3,rnorm(100),simplify=FALSE),r$background))
	x<-le8_prepare_prediction_matrix(d,d,r$background)
	c4_assert_replacement(colnames(x$train),r$forbidden)
	expect_error(c4_assert_replacement(c(colnames(x$train),'bmi__high'),r$forbidden))
})
check('A02 partial replacement preserves other domains and additive definition',{
	base<-c('age','sex');measured<-c('bmi','weight','bmi.pts','sbp','bp.pts','le8.pts')
	additive<-unique(c(base,measured,'disease_prs'))
	r<-c4_replacement_background(base,measured,'disease_prs','bmi.pts')
	stopifnot(!any(c('bmi','bmi.pts','weight','le8.pts') %in% r$background),all(c('sbp','bp.pts') %in% r$background),all(measured %in% additive))
})
check('A03 deployed scores change fidelity while frozen panel OLS is unchanged',{
	set.seed(2);d<-data.frame(eid=paste0('s',1:260),.group=paste0('g',rep(1:130,each=2)),age=rnorm(260),F1=rnorm(260),F2=rnorm(260))
	d$bmi.pts<-50+8*d$F1+2*d$F2+rnorm(260);tr<-d[1:180,];te<-d[181:260,]
	x<-c4_concept_transform(tr,te,c('F1','F2'),'bmi.pts')
	cz<-list(status='ok',test=cbind(te,x$test),features='concept_bmi.pts')
	a<-c4_deployed_concept_fidelity(cz,tr,te,'bmi.pts','YS',B=30L)
	cz$test$concept_bmi.pts<-cz$test$concept_bmi.pts+20
	b<-c4_deployed_concept_fidelity(cz,tr,te,'bmi.pts','YS',B=30L)
	stopifnot(a$estimate[a$metric=='RMSE']<b$estimate[b$metric=='RMSE'],all(a$bootstrap_valid==30L),all(is.finite(a$lower)))
	post1<-focus_proxy_accuracy(tr,te,c('F1','F2'),'bmi.pts','age','YS')
	post2<-focus_proxy_accuracy(tr,te,c('F1','F2'),'bmi.pts','age','YS')
	stopifnot(identical(post1,post2),all(post1$evaluation=='posthoc_panel_reconstruction'))
	te$bmi.pts<-rnorm(nrow(te))*10000
	x2<-c4_concept_transform(tr,te,c('F1','F2'),'bmi.pts')
	stopifnot(identical(x$test,x2$test))
})
check('A04 cis null/trans strong remain in separate layers and failed rows survive',{
	mr<-tibble(exposure=c('P','P','failed'),analysis=c('cis','trans','cis'),b=c(.1,.8,NA),pval=c(.3,1e-10,NA),FDR_all=c(.4,1e-9,NA))
	rv<-tibble(feature='P',b=.1,pval=.8,FDR_reverse=.8,n_IV=4L)
	e<-le8_bidirectional_evidence(mr,rv,'protein')
	stopifnot(!e$forward_support[e$feature=='P'&e$scope=='primary'],e$forward_support[e$feature=='P'&e$scope=='secondary'],e$forward_status[e$feature=='failed'&e$scope=='primary']=='failed_or_not_estimable')
	stopifnot(le8_select_mr(mr,'protein')$pval[1]==.3)
})
mr<-tibble(exposure='P',analysis='cis',FDR_all=.01,instrument_chr='1',instrument_pos_min=100,instrument_pos_max=100,instrument_positions='1:100',instrument_snps='rsA')
co<-tibble(feature='P',locus='chr1:90-110',chr='1',start=90,end=110,status='ok',PP.H4_robust_min=.75,locus_class='cis',aligned_MR_IVs='rsA',beta_scale_status='verified',prior_complete=TRUE)
check('A05 posterior 0.75 fails 0.8 everywhere; passes 0.7 as region evidence',{
	Sys.setenv(C3_H4='.8');p8<-le8_evidence_policy();stopifnot(!le8_same_locus_evidence(mr,co,'protein')$eligible,!le8_coloc_pass(co),!length(le8_same_locus_candidates(mr,co,'protein')))
	Sys.setenv(C3_H4='.7');p7<-le8_evidence_policy();stopifnot(le8_same_locus_evidence(mr,co,'protein')$eligible,le8_coloc_pass(co),identical(le8_same_locus_candidates(mr,co,'protein'),'P'),p8$hash!=p7$hash)
})
check('A05 incomplete priors, FDR and signal policy cannot be bypassed',{
	co2<-co;co2$prior_complete<-FALSE;stopifnot(!le8_same_locus_evidence(mr,co2,'protein')$eligible)
	Sys.setenv(C3_MR_FDR='.005');stopifnot(!le8_same_locus_evidence(mr,co,'protein')$eligible)
	Sys.setenv(C3_MR_FDR='.05',LE8_EVIDENCE_LEVEL='signal_only');stopifnot(!le8_same_locus_evidence(mr,co,'protein')$eligible)
	co$MR_signal_status<-'all_instrument_signals_supported';stopifnot(le8_same_locus_evidence(mr,co,'protein')$eligible)
	Sys.unsetenv('LE8_EVIDENCE_LEVEL')
})
check('A06 SNV swapping, build and ALT identity',{
	x<-standardize_sumstat(tibble(SNP='rsA',CHR='1',POS=5,EA='C',NEA='A',BETA=.2,SE=.1,P=.05,N=100,EAF=.2))
	y<-x;y$EA<-'A';y$NEA<-'C';y$BETA<--.2;y$EAF<-.8
	stopifnot(nrow(harmonize_sumstats(x,y))==1,harmonize_sumstats(x,y)$BETA_y==.2)
	x$BUILD<-'37';y$BUILD<-'38';stopifnot(nrow(harmonize_sumstats(x,y))==0)
	y$BUILD<-'37';y$EA<-'G';y$NEA<-'A';stopifnot(nrow(harmonize_sumstats(x,y))==0)
})
check('A06 unverified INDELs are excluded; reference norm matches shifted alleles',{
	fa<-file.path(scratch,'reference.fa');writeLines(c('>1','CAAAAAAGTCGATCGATCGA'),fa);stopifnot(system2('samtools',c('faidx',fa))==0)
	meta<-list(build='37',reference_fasta=fa,normalization_proof=NA_character_)
	make<-function(pos) { z<-standardize_sumstat(tibble(SNP=paste0('shift',pos),CHR='1',POS=pos,REF='AA',ALT='A',EA='A',NEA='AA',BETA=.2,SE=.1,P=.05,N=100,EAF=.2));le8_variant_identity(z,meta) }
	x<-make(2);y<-make(4);stopifnot(nrow(harmonize_sumstats(x,y))==0)
	x<-le8_normalize_variants(x,meta);y<-le8_normalize_variants(y,meta)
	stopifnot(identical(x$variant_id,y$variant_id),nrow(harmonize_sumstats(x,y))==1,all(x$normalization_status=='reference_checked_left_aligned_split'))
	bad<-make(2);bad$REF<-'GG';bad$NEA<-'GG';expect_error(le8_normalize_variants(bad,meta))
})
check('A06 conflicting duplicate effects fail, including flipped effects',{
	d<-tibble(SNP=c('a','b'),CHR='1',POS=8,EA='A',NEA='C',BETA=c(.2,.8),SE=.1,P=.05,N=100,EAF=.2)
	expect_error(standardize_sumstat(d));d$BETA[2]<-.2;stopifnot(nrow(standardize_sumstat(d))==1)
})
check('A07 singleton silhouette is zero; no feasible K and one-feature states explicit',{
	D<-dist(matrix(c(0,1,10),ncol=1));stopifnot(cluster::silhouette(c(1L,1L,2L),D)[3,'sil_width']==0)
	d<-expand_grid(feature=paste0('f',1:3),component=c('bmi','bp'));d$r<-c(.8,.1,.7,.2,.1,.8)
	f<-fit_le8_modules(d,d);stopifnot(f$k==1L,f$status=='no_feasible_K_single_module',!any(f$metrics$feasible))
	f<-fit_le8_modules(d[1:2,],d[1:2,]);stopifnot(f$k==1L,f$status=='too_few_features_single_module')
})
check('A07 original and every bootstrap use the same profile and size rules',{
	set.seed(12);n<-500L
	d<-data.frame(eid=paste0('id',1:n),.group=paste0('family',rep(1:(n/2),each=2)),age=rnorm(n),bmi.pts=rnorm(n),bp.pts=rnorm(n))
	fs<-paste0('F',1:12)
	for(j in seq_along(fs)) d[[fs[j]]]<-if(j<=6) d$bmi.pts+rnorm(n,sd=.15) else d$bp.pts+rnorm(n,sd=.15)
	d$.le8_proxy_half<-stratified_split(d,NULL)
	disc<-proxy_scan(d[d$.le8_proxy_half=='discovery',],fs,c('bmi.pts','bp.pts'),'age','discovery')
	rep<-proxy_scan(d[d$.le8_proxy_half=='replication',],fs,c('bmi.pts','bp.pts'),'age','replication')
	sets<-list(primary=tibble(feature=fs,primary_component=rep(c('bmi','bp'),each=6),strict_YS=TRUE,YS_model=TRUE,FDR_disc=.001,FDR_rep=.001))
	fit<-fit_le8_modules(disc,rep);stopifnot(fit$status=='ok',all(fit$metrics$min_module_n[fit$metrics$k==fit$k]>=fit$metrics$minimum_required[1]))
	C4_MODULE_BOOT<-0L;z<-supervised_modules(disc,rep,sets,d,c('bmi.pts','bp.pts'),'age');stopifnot(all(is.na(z$membership$module_stability)),all(z$membership$stability_status=='not_requested_B0'))
	C4_MODULE_BOOT<-20L;z<-supervised_modules(disc,rep,sets,d,c('bmi.pts','bp.pts'),'age')
	stopifnot(nrow(z$bootstrap)==20L,all(z$bootstrap$status=='ok'),all(z$bootstrap$k==fit$k),all(is.finite(z$membership$module_stability)))
})
check('A01 selected diagnosis alias survives phenotype column selection',{
	Y<-'cad';le8_custom_covars<-character();Sys.setenv(LE8_Y_DATE='custom_date')
	d<-data.frame(eid='x',custom_date='2020-01-01');d<-le8_select_phenotypes(d)
	stopifnot(identical(d$fod_icd10_cad,d$custom_date))
	Sys.unsetenv('LE8_Y_DATE')
	Sys.setenv(DATE_FOLLOW_END='2020-12-31');stopifnot(le8_follow_end()==as.Date('2020-12-31'));Sys.setenv(DATE_FOLLOW_END='2020-02-31');expect_error(le8_follow_end());Sys.unsetenv('DATE_FOLLOW_END')
})
check('A02 actual fitted replacement coefficients exclude all replaced measurements',{
	set.seed(73);d<-data.frame(eid=paste0('x',1:400),age=rnorm(400),sex=rbinom(400,1,.5),bmi=rnorm(400),bmi.pts=rnorm(400),sbp=rnorm(400),disease_prs=rnorm(400),concept_bmi.pts=rnorm(400),time=rexp(400)+.1,event=rbinom(400,1,.4))
	r<-c4_replacement_background(c('age','sex'),c('bmi','bmi.pts','sbp'),'disease_prs',unique(c4_le8_measure_map()$component))
	fit<-le8_fit_budget_model(d[1:300,],d[301:400,],r$background,'concept_bmi.pts','time','event')
	stopifnot(fit$status=='ok');c4_assert_replacement(fit$coefficient$variable,r$forbidden)
	stopifnot(all(c('age','sex','disease_prs','concept_bmi.pts') %in% fit$coefficient$variable),max(abs(rowSums(fit$contributions)-fit$lp-fit$lp_center))<1e-10)
})
check('A03 real nested concept fit is unchanged by poisoning held-out LE8 labels',{
	set.seed(74);n<-800L;d<-data.frame(eid=paste0('n',1:n),.group=paste0('f',rep(1:(n/2),each=2)),age=rnorm(n),F1=rnorm(n),F2=rnorm(n))
	d$bmi.pts<-50+8*d$F1+6*d$F2+rnorm(n)
	tr<-d[1:600,];te<-d[601:800,];yang<-d[FALSE,]
	z1<-c4_fit_concepts(tr,te,yang,c('F1','F2'),'bmi.pts','age',c('F1','F2'),'Yin',2L)
	te$bmi.pts<-rnorm(nrow(te))*1e6
	z2<-c4_fit_concepts(tr,te,yang,c('F1','F2'),'bmi.pts','age',c('F1','F2'),'Yin',2L)
	stopifnot(z1$status=='ok',z2$status=='ok',identical(z1$test[z1$features],z2$test[z2$features]),identical(z1$fold_panels,z2$fold_panels),identical(z1$coefficients,z2$coefficients))
})
check('A05 policy invalidates affected cache keys and is available to early completion check',{
	Sys.setenv(C3_H4='.7');a<-le8_module_policy('c3_coloc');c1<-le8_module_policy('c1_correlate')
	Sys.setenv(C3_H4='.8');expect_error(le8_coloc_pass(transform(co,policy_hash=a$evidence$hash)));stopifnot(!identical(a,le8_module_policy('c3_coloc')),identical(c1,le8_module_policy('c1_correlate')))
	Sys.setenv(C3_H4='.7')
	r<-file.path(scratch,'completed');dir.create(file.path(r,'cad/prot/c2_cause'),recursive=TRUE)
	saveRDS(list(meta=list(trait='cad',layer='protein',module='c2_cause',generated='synthetic',seed=2026L,analysis_options=le8_analysis_options('cad'),module_policy=le8_module_policy('c2_cause'))),file.path(r,'cad/prot/c2_cause/c2.res.rds'))
	status<-system2(file.path(R.home('bin'),'Rscript'),c(shQuote(file.path(root,'f/0.common.R')),'--check-completed',r,'cad','prot','c2_cause',file.path(scratch,'completed.csv')),stdout=file.path(scratch,'check.log'),stderr=file.path(scratch,'check.err'))
	if(status!=0) stop(paste(readLines(file.path(scratch,'check.err')),collapse='\n'))
	stopifnot(file.exists(file.path(scratch,'completed.csv')))
})
check('A06 normalized proof reuse rejects stale file or incompatible effect alleles',{
	f<-file.path(scratch,'proven.csv');proof<-file.path(scratch,'proven.json')
	d<-tibble(SNP='indel',CHR='1',POS=1,REF='CA',ALT='C',EA='C',NEA='CA',BETA=.2,SE=.1,P=.05,N=100,EAF=.2);fwrite(d,f)
	p<-list(file_sha256=le8_file_sha256(f),reference_sha256=le8_file_sha256(fa),tool='bcftools',tool_version=system2('bcftools','--version-only',stdout=TRUE),build='37',reference_checked=TRUE,left_aligned=TRUE,multiallelic_split=TRUE)
	jsonlite::write_json(p,proof,auto_unbox=TRUE)
	z<-standardize_sumstat(d,source_file=f);meta<-list(build='37',normalization_proof=proof,reference_fasta=NA_character_)
	stopifnot(nrow(le8_variant_identity(z[FALSE,],meta))==0L)
	z<-le8_variant_identity(z,meta);stopifnot(z$variant_id=='37:1:1:CA:C',nrow(harmonize_sumstats(z,z))==1)
	z$EA<-'T';expect_error(le8_variant_identity(z,meta));cat('\n',file=f,append=TRUE);expect_error(le8_variant_identity(z,meta))
})
check('A06 real prepared pair normalizes INDELs and preserves CPU BF definition',{
	qf<-file.path(scratch,'qtl.csv');yf<-file.path(scratch,'outcome.csv');mf<-file.path(scratch,'gwas.csv')
	q<-tibble(SNP=c('iv_indel','snv1','snv2'),CHR='1',POS=c(4,8,9),REF=c('AA','G','T'),ALT=c('A','A','C'),EA=c('A','A','C'),NEA=c('AA','G','T'),BETA=c(.6,.2,-.1),SE=.1,P=c(1e-9,.05,.3),N=1000,EAF=.2)
	y<-q;y$SNP[1]<-'outcome_indel';y$POS[1]<-2;y$BETA<-c(.5,.1,-.1)
	fwrite(q,qf);fwrite(y,yf);fwrite(data.frame(file=c(qf,yf),build='37',sdY=1,beta_scale=c('SD','log_odds'),ancestry='EUR',reference_fasta=fa),mf)
	Sys.setenv(LE8_GWAS_MANIFEST=mf)
	p<-c3_prepare_locus_pair(qf,yf,'1',1,19,'cc',.3)
	stopifnot(nrow(p$pair)==3L,all(p$pair$normalization_status=='reference_checked_left_aligned_split'),all(p$lead_coverage),!anyNA(p$pair$variant_id))
	x<-c3_reference_posterior(p$pair$lbf1,p$pair$lbf2)
	fit<-suppressWarnings(coloc::coloc.abf(p$x,p$y,p1=1e-4,p2=1e-4,p12=1e-5))
	stopifnot(max(abs(unname(x)-unname(fit$summary[paste0('PP.H',0:4,'.abf')])))<1e-10)
	iv<-standardize_sumstat(tibble(SNP='iv_indel',CHR='1',POS=100,refA='A',bJ=.6,bJ_se=.1,pJ=1e-9),joint=TRUE)
	iv<-recover_qtl_alleles(iv,qf);iv<-le8_normalize_variants(iv,le8_gwas_metadata(qf))
	stopifnot(iv$BUILD=='37',iv$variant_key %in% p$pair$variant,!is.na(iv$variant_id))
	Sys.unsetenv('LE8_GWAS_MANIFEST')
})

unlink(scratch,recursive=TRUE)
cat(passed,'R acceptance cases passed\n')
