source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
if(!requireNamespace("terra",quietly=TRUE))
  stop("Package `terra` is required for the optional HadUK-Grid temperature linkage.")

tempdir<-Sys.getenv("BMONTH_TEMP_DIR",file.path(dir0,"data","ukb","temperature"))
resolution<-Sys.getenv("BMONTH_TEMP_RES","1km")
if(!nzchar(tempdir)||!dir.exists(tempdir))
  stop("Set --temp-dir to a directory containing HadUK-Grid monthly mean-temperature NetCDF files.")

cohort<-load_cohort()
req<-c("eid","birth_month","birth_year_bm","birth_east","birth_north")
miss<-setdiff(req,names(cohort));if(length(miss))stop("Cohort missing: ",paste(miss,collapse=", "))
z<-cohort[,req,drop=FALSE]
z$birth_year_bm<-as.integer(as_num(z$birth_year_bm))
z$birth_month<-as.integer(as_num(z$birth_month))

# Estimated conception month = 9 calendar months before birth month.
idx<-z$birth_year_bm*12L+(z$birth_month-1L)
cidx<-idx-9L
z$conception_year<-cidx%/%12L
z$conception_month<-(cidx%%12L)+1L

nc<-list.files(tempdir,pattern="\\.nc$",recursive=TRUE,full.names=TRUE)
if(!length(nc))stop("No NetCDF files under ",tempdir)
bn<-basename(nc)

find_year_file<-function(year){
  # Standard HadUK-Grid annual monthly file, e.g.
  # tas_hadukgrid_uk_1km_mon_195001-195012.nc
  ok<-grepl("^tas_hadukgrid_uk_",bn) &
      grepl(paste0("_",resolution,"_"),bn,fixed=TRUE) &
      grepl("_mon_",bn,fixed=TRUE) &
      grepl(sprintf("_mon_%04d01-%04d12[.]nc$",year,year),bn) &
      !grepl("30y|20y",bn)
  hits<-nc[ok]
  if(!length(hits))return(NA_character_)
  if(length(hits)>1L)stop("Ambiguous HadUK-Grid files for ",year,": ",paste(hits,collapse=", "))
  hits[1]
}
pick_month_layer<-function(r,year,month){
  tt<-try(terra::time(r),silent=TRUE)
  if(!inherits(tt,"try-error")&&(inherits(tt,"Date")||inherits(tt,"POSIXt"))&&
     length(tt)==terra::nlyr(r)&&any(!is.na(tt))){
    yy<-as.integer(format(as.Date(tt),"%Y"));mm<-as.integer(format(as.Date(tt),"%m"))
    hit<-which(yy==year&mm==month)
    if(length(hit)==1L)return(r[[hit]])
    stop("NetCDF time axis does not uniquely identify ",year,"-",month)
  }
  if(terra::nlyr(r)==12L)return(r[[month]])
  stop("Cannot identify month ",year,"-",sprintf("%02d",month)," in ",terra::sources(r)[1])
}
read_tas<-function(file){
  r<-terra::rast(file,subds="tas")
  if(is.na(terra::crs(r))||!nzchar(terra::crs(r)))stop("NetCDF missing coordinate reference system: ",file)
  u<-unique(tolower(terra::units(r)))
  if(length(u)!=1L)stop("Inconsistent temperature units: ",file)
  if(u%in%c("k","kelvin"))r<-r-273.15 else
    if(!u%in%c("c","degc","degree_celsius","degrees_celsius","celsius"))
      stop("Unsupported or missing temperature units in ",file,": ",u)
  r
}
eligible<-is.finite(z$birth_east)&is.finite(z$birth_north)&
          is.finite(z$birth_year_bm)&z$birth_month%in%1:12
years<-sort(unique(c(z$birth_year_bm[eligible],z$conception_year[eligible])))
if(!length(years))stop("No participants with valid birth year, month and coordinates")
year_files<-setNames(vapply(years,find_year_file,character(1)),years)
if(anyNA(year_files))stop("Missing HadUK-Grid annual monthly tas files for: ",
                         paste(names(year_files)[is.na(year_files)],collapse=", "))
if(isTRUE(attr(cohort,"bmonth_metadata")$birth_year_approx))
  warning("Birth year is approximate (assessment_year-age); historical temperature may be assigned to the wrong year.")
extract_monthly<-function(year,month,east,north){
  out<-rep(NA_real_,length(year))
  keys<-unique(data.frame(year=year,month=month))
  keys<-keys[is.finite(keys$year)&keys$month%in%1:12,,drop=FALSE]
  last_year<-NA_integer_;r<-NULL
  keys<-keys[order(keys$year,keys$month),,drop=FALSE]
  for(i in seq_len(nrow(keys))){
    yy<-keys$year[i];mm<-keys$month[i]
    take<-which(year==yy&month==mm&is.finite(east)&is.finite(north))
    if(!length(take))next
    if(is.na(last_year)||last_year!=yy){r<-read_tas(year_files[[as.character(yy)]]);last_year<-yy}
    lyr<-pick_month_layer(r,yy,mm)
    pts<-terra::vect(data.frame(x=east[take],y=north[take]),
                     geom=c("x","y"),crs="EPSG:27700")
    if(!is.na(terra::crs(lyr))&&nzchar(terra::crs(lyr)))
      pts<-terra::project(pts,terra::crs(lyr))
    v<-terra::extract(lyr,pts)
    out[take]<-as.numeric(v[[ncol(v)]])
    message("[temp] ",yy,"-",sprintf("%02d",mm),": ",length(take)," participants")
  }
  out
}

z$temp_birth_mean<-extract_monthly(z$birth_year_bm,z$birth_month,z$birth_east,z$birth_north)
z$temp_conception_mean<-extract_monthly(z$conception_year,z$conception_month,z$birth_east,z$birth_north)

# Optional 1961-1990 monthly climatology in the same download directory.
clim_hit<-nc[
  grepl("^tas_hadukgrid_uk_",bn)&grepl(paste0("_",resolution,"_"),bn,fixed=TRUE)&
  grepl("mon-30y|mon_30y",bn)&grepl("1961",bn)&grepl("1990",bn)
]
if(length(clim_hit)>1L)stop("Multiple 1961-1990 climatology files; keep one version in --temp-dir")
extract_clim<-function(month,east,north){
  out<-rep(NA_real_,length(month))
  if(!length(clim_hit))return(out)
  r<-read_tas(clim_hit[1])
  if(terra::nlyr(r)!=12L)stop("Monthly climatology must have 12 layers")
  for(mm in 1:12){
    take<-which(month==mm&is.finite(east)&is.finite(north));if(!length(take))next
    lyr<-if(terra::nlyr(r)==12L)r[[mm]] else pick_month_layer(r,1961,mm)
    pts<-terra::vect(data.frame(x=east[take],y=north[take]),
                     geom=c("x","y"),crs="EPSG:27700")
    if(!is.na(terra::crs(lyr))&&nzchar(terra::crs(lyr)))
      pts<-terra::project(pts,terra::crs(lyr))
    v<-terra::extract(lyr,pts);out[take]<-as.numeric(v[[ncol(v)]])
  }
  out
}
z$temp_birth_expected<-extract_clim(z$birth_month,z$birth_east,z$birth_north)
z$temp_conception_expected<-extract_clim(z$conception_month,z$birth_east,z$birth_north)
z$temp_birth_anom<-z$temp_birth_mean-z$temp_birth_expected
z$temp_conception_anom<-z$temp_conception_mean-z$temp_conception_expected

meta<-list(
  source="Met Office HadUK-Grid monthly mean air temperature (tas)",
  units="degrees Celsius",
  cohort_md5=unname(tools::md5sum(file.path(rawdir,"bmonth.cohort.rds"))),
  input_files=year_files,
  input_md5=tools::md5sum(c(unname(year_files),clim_hit)),
  directory=normalizePath(tempdir),
  resolution=resolution,
  conception_definition="9 calendar months before month of birth",
  climatology=if(length(clim_hit))normalizePath(clim_hit[1]) else NA_character_,
  birth_year_source=attr(cohort,"bmonth_metadata")$birth_year_source,
  birth_year_approx=attr(cohort,"bmonth_metadata")$birth_year_approx,
  warning="If birth year was derived as assessment_year-age, year can be off by one."
)
attr(z,"temperature_metadata")<-meta
saveRDS(z,file.path(rawdir,"bmonth_temperature.rds"))
write_csv2(
  tibble(metric=c("birth temperature N","conception temperature N",
                  "birth anomaly N","conception anomaly N"),
         value=c(sum(is.finite(z$temp_birth_mean)),sum(is.finite(z$temp_conception_mean)),
                 sum(is.finite(z$temp_birth_anom)),sum(is.finite(z$temp_conception_anom)))),
  "temperature_qc.csv")
writeLines(capture.output(str(meta)),file.path(rawdir,"temperature_metadata.txt"))
cat("Saved ",file.path(rawdir,"bmonth_temperature.rds"),"\n",sep="")
