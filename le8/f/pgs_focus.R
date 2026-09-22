# Focused rerun: no MR/ML dependency and no modification of existing PGS files.
fdir<-Sys.getenv("LE8_FDIR",unset=".")
source(file.path(fdir,"comm.f.R"))
source(file.path(fdir,"c1_pgs.R"))
source(file.path(fdir,"c3_pgs.R"))
for(layer in c(if(prot_DO)"protein",if(met_DO)"metabolite")) {
  outdir<-if(layer=="protein")out.prot else out.met
  x<-run_c1_pgs_focus(layer,outdir)
  cf<-file.path(outdir,"c3_coloc","c3.coloc_summary.csv")
  co<-if(file.exists(cf))as.data.frame(data.table::fread(cf,showProgress=FALSE))else data.frame()
  tri<-read_c3_pgs_integration(layer,outdir,co)
  pgs_plot_focus(x,outdir,tri,attr(tri,"loci"))
  plot_c3_pgs_integration(tri,outdir)
  rm(x,tri);invisible(gc())
}
