# Publication assembly from aggregate module outputs only. No participant data,
# model fitting, validation-set selection, or internet access is required.
suppressPackageStartupMessages({library(data.table);library(dplyr);library(tidyr)
  library(ggplot2);library(patchwork);library(openxlsx);library(purrr)})
`%||%`<-function(x,y)if(is.null(x))y else x
source(file.path(Sys.getenv('LE8_FDIR', unset='.'),'figure_policy.R'))
source(file.path(Sys.getenv('LE8_FDIR', unset='.'),'le8_budget_config.R'))
source(file.path(Sys.getenv('LE8_FDIR', unset='.'),'publication_panels.R'))
traits<-Sys.getenv('Y',unset='cvd_cad')
layers<-strsplit(Sys.getenv('BIOM',unset='prot,met'),',',fixed=TRUE)[[1]]
for(layer in layers) publication_run(traits,layer,Sys.getenv('LE8_ANALYSIS_ROOT',unset='/mnt/d/analysis/le8'))
le8_layer_figure_correspondence(file.path(Sys.getenv('LE8_ANALYSIS_ROOT',unset='/mnt/d/analysis/le8'),traits))
