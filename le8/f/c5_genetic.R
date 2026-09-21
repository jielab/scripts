# Compatibility alias: joint validation for any Y, not separate layer cohorts.
fdir<-Sys.getenv("LE8_FDIR",file.path(Sys.getenv("DIRSCRIPT"),"f"))
source(file.path(fdir,"c5_consolidate.R"))
