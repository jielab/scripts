source(file.path(Sys.getenv("BORDER_FDIR"), "0.common.R"))
dat <- load_cohort()
res <- bind_rows(lapply(OUTCOMES, function(y) {
  message("[assoc] START ",y,": ",length(EXPOSURES)," birth-order models")
  r<-bind_rows(lapply(EXPOSURES, function(x) fit_cox(dat,x,y)))
  message("[assoc] DONE ",y,": ",sum(r$status=="ok"),"/",nrow(r)," successful fits")
  r
}))
if (nrow(res)) res <- res %>% group_by(exposure) %>% mutate(FDR=p.adjust(p,"BH")) %>% ungroup()
write_csv2(res,"main_associations.csv")
saveRDS(res,file.path(rawdir,"main_associations.rds"))
print(res, n=Inf)
