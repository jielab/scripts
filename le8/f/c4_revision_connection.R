# Connection is association, not intervention or automatically causal mediation.
if(!exists(".le8_original_mediation_one",envir=.GlobalEnv,inherits=FALSE))
  .le8_original_mediation_one<-mediation_one
mediation_one <- function(...) {
  z<-.le8_original_mediation_one(...)
  if(nrow(z)) {
    z$associational_product_ratio<-z$prop_mediated
    z$causal_mediation_identified<-FALSE
    z$modifiability_demonstrated<-FALSE
    z$interpretation<-"Same-baseline exposure and omic; coefficient-product ratio is descriptive, not an identified causal mediation proportion"
  }
  z
}
le8_c4_additions <- function(out,outdir) {
  scan<-as_tibble(out$scan%||%tibble())
  if(!nrow(scan)||!"component"%in%names(scan))return(list(status=tibble(status="connection scan unavailable")))
  scan<-scan|>mutate(connection_domain=ifelse(component%in%c("diet","pa","smoke","sleep"),
    "Behavioral LE4","Biological LE4"),causal_interpretation="not established")
  count<-scan|>group_by(connection_domain,component,split)|>
    summarise(tested=n(),FDR_supported=sum(is.finite(FDR_component)&FDR_component<.05),.groups="drop")
  rd<-le8_job_dir(outdir,"c4_connect")
  write_raw_csv(scan,"c4.behavior_biology_associations.csv",rd)
  write_raw_csv(count,"c4.behavior_biology_summary.csv",rd)
  list(behavior_biology=count,association_scope=scan,
    interpretation=tibble(statement=c("LE8 association does not prove that changing behavior will change this protein",
      "A non-PGS residual still includes uncaptured inherited variation",
      "Same-baseline path products are not an identified interventional mediation proportion")))
}
