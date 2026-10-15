source(file.path(Sys.getenv('BMONTH_FDIR'),'0.common.R'))
dat<-read_all()%>%prep_covars()%>%make_birth_month()
for(y in OUTCOMES)dat<-ensure_outcome(dat,y)
keep<-unique(c('eid','birth_month','birth_north','birth_east',ALL_EXPOSURES,
               'bmonth_cos_x_north','bmonth_sin_x_north','bmonth_cos_x_east','bmonth_sin_x_east',covars,
               unlist(lapply(OUTCOMES,function(y)paste0(y,c('.t2e','.Yt2e'))))))
cohort<-dat[,intersect(keep,names(dat)),drop=FALSE]
attr(cohort,'bmonth_metadata')<-list(input=file.path(indir,'Rdata','all.rds'),outcomes=OUTCOMES,follow_end=as.character(follow_end()))
saveRDS(cohort,file.path(rawdir,'bmonth.cohort.rds'))
qc<-bind_rows(tibble(metric='all.rds N',value=nrow(dat)),
              tibble(metric='birth_month nonmissing',value=sum(!is.na(dat$birth_month))),
              tibble(metric='birth_north nonmissing',value=sum(!is.na(dat$birth_north))),
              tibble(metric='birth_east nonmissing',value=sum(!is.na(dat$birth_east))),
              bind_rows(lapply(ALL_EXPOSURES,function(x)tibble(metric=paste0(x,' nonmissing'),value=sum(!is.na(dat[[x]]))))))
write_csv2(qc,'cohort_qc.csv')
write_csv2(count(dat,birth_month,name='N'),'birth_month_counts.csv')
cat('Saved ',file.path(rawdir,'bmonth.cohort.rds'),'\n',sep='')
print(qc)
