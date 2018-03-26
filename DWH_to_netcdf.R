#!/usr/bin/env Rscript

# This script can be run from the command line to extract
# station data from DWH and save it in a netcdf

# packages
suppressMessages({
  library(dwhget)
  library(tidyverse)
  library(ncdf4)
  library(R.utils)
  library(lubridate)
})

# check command-line arguments
args <- R.utils::commandArgs(TRUE, asValue=TRUE)

if (length(args) < 2){
  cat("Usage info:

DWH_to_netcdf --param=<pname> --period=<period> --mstype=<mstype> --outfile=<outfile>
")
  q(save = 'no')
}

extractArgs <- function(args, argname, default){
  if (length(args) >= 1){
    out <- args[[argname]]
    if (is.null(out)){
      out <- args[1]
      args <<- args[-1]
    } else {
      args <<- args[-which(names(args) == argname)]
    }
  } else {
    out <- default
  }
  out
}


param <- extractArgs(args, 'param', 'tre200b0')
period <- extractArgs(args, 'period', '201601')
mstype <- extractArgs(args, 'mstype', '')
outfile <- extractArgs(args, 'outfile', '')

## warp command-line arguments
if (mstype == "") mstype <- NULL
if (is.null(mstype)) msstring <- "default"
fixperiod <- TRUE
if (nchar(period) == 6){
  p2 <- as.Date(paste0(period, '01'), format="%Y%m%d") + seq(0, 31)
  period <- format(p2[range(which(format(p2, "%Y%m") == period))+ c(0,1)], "%Y%m%d")
} else if (grepl("-", period)){
  period <- strsplit(period, "-")[[1]]
}
if (outfile == ""){
  outfile <- paste0("/store/msclim/bhendj/tmp/DWH/",
                    paste(param, paste(period, collapse ='-'), msstring, sep='_'),
                    ".nc")
}

print(paste("Parameter:       ", param))
print(paste("Period:          ", period[1], "to", period[2]))
print(paste("Measurement type:", mstype))
print(paste("Output file:     ", outfile))

## get all stations data
out <- dwhget(stats = 'all',
              mstypes = mstype,
              pars=param,
              t.beg=period[1],
              t.end=period[2],
              datform="yyyymmdd",
              filter = list(use.limitation.id = 40,
                iso.country.cd = "CH"),
              stat.id.type = 'station_id',
              additional.output = c("owner.id", "owner.name", "use.limitation.id")) %>%
  as_tibble() %>%
  mutate(time = as.POSIXct(as.character(time), format = '%Y%m%d%H%M'))

## hack to get rid of 00UTC from next month
if (fixperiod){
  out <- filter(out, format(time, "%Y%m%d") != period[2])
}

## metadata on parameter
par.meta <- try(dwh.para.info(pars=param), silent=TRUE)
if (class(par.meta) == "try-error"){
  par.meta <- data.frame(unit = ifelse(param == "tre200b0", "Celsius", '1'))
}


stations <- na.omit(unique(out$station))

## get station meta data
stat.meta <- dwh.stat.info(stats = stations,
                           mstypes = if (is.null(mstype)) {
                             iconv(dwh.mstype.info()[,2], from='latin1', to='UTF-8')
                           } else {
                             mstype
                           },
                           stat.id.type = 'station_id',
                           output = c("station_id", "nat_ind", "nat_abbr",
                             "wmo_ind", "lon", "lat", "height", "name"))
stat.meta$name <- iconv(stat.meta$name, from = 'latin1', to = 'UTF-8')
stat.meta <- as_tibble(stat.meta) %>%
  arrange(station_id) %>%
  full_join(out %>%
            select(station, owner.id:use.limitation.id) %>%
            unique() %>%
            rename(station_id = station) %>%
            mutate(owner.name = iconv(owner.name, from = 'latin1', to = 'UTF-8'))) %>%
  filter(lon > -1e11, lat > -1e11, height > -1e11)



## convert data frame
out.wide <- as_tibble(out) %>%
  select(-parameter, -owner.id:-use.limitation.id) %>%
  spread(station, value)
out.arr <- as.matrix(select(out.wide, -time))

## check that ordering of meta data is correct
stopifnot(as.character(stat.meta$station_id) == colnames(out.arr))

## Trying to follow CF Conventions 1.7

## NetCDF Dimensions
reftime <- min(out.wide$time)
tdiff <- as.numeric(difftime(out.wide$time, reftime, units = "hours"))
time.nc <- ncdim_def("time",
                     paste("hours since", format(reftime, "%Y-%m-%d %H:%M:%S")),
                     tdiff,
                     unlim=TRUE,
                     longname = 'time of measurement')
station.nc <- ncdim_def('station',
                        '',
                        1:nrow(stat.meta),
                        create_dimvar=FALSE)


## get number of characters necessary for the meta information
charnames <- stat.meta %>%
  select_if(is.character) %>%
  names()
numnames <- setdiff(names(stat.meta), charnames)
nchars <- sapply(stat.meta %>%
  select_if(is.character), function(x) max(nchar(x), na.rm=T))

nchar.nc <- lapply(nchars, function(x) ncdim_def(paste("string", x, sep="_"),
                                                 "", 1:x, create_dimvar=FALSE))

numunits <- sapply(numnames, function(x) "1", simplify=FALSE)
numunits[c("lon", "lat", "height")] <- list(lon = "degrees_east",
                                            lat = "degrees_north",
                                            height = "m")

## NetCDF Variables
vars.nc <- c(sapply(charnames, function(x){
               ncvar_def(x, "", list(nchar.nc[[x]], station.nc), prec = 'char')},
                  simplify=FALSE),
             sapply(numnames, function(x){
               ncvar_def(x, numunits[[x]], station.nc,
                         prec = ifelse(is.integer(stat.meta[[x]]), "integer", "float"),
                         missval = ifelse(is.integer(stat.meta[[x]]), -9999, -1e20))},
                    simplify=FALSE))
vars.nc[[param]] <- ncvar_def(param, par.meta$unit, list(station.nc, time.nc), missval=-1e20)


## open file and write to it
if (file.exists(outfile)){
  file.remove(outfile)
}
system(paste("mkdir -p", dirname(outfile))) 

ncout <- nc_create(outfile, 
          vars.nc, 
          force_v4 = TRUE)
for (nn in names(stat.meta)){
  ncvar_put(ncout, nn, stat.meta[[nn]])
}
ncvar_put(ncout, param, t(out.arr))

## additional attributes
standard_names <- c(lon = 'longitude', lat='latitude', height='altitude')
for (sn in names(standard_names)) ncatt_put(ncout, sn, "standard_name", standard_names[sn], prec = 'text')
ncatt_put(ncout, param, "coordinates",
          paste(names(stat.meta), collapse = " "), prec = 'text')
ncatt_put(ncout, 'station_id', "cf_role", "timeseries_id", prec='text')

## put global attributes
ncatt_put(ncout, 0, 'Conventions', "CF-1.7", prec='text')
ncatt_put(ncout, 0, "featureType", "timeSeries", prec='text')
ncatt_put(ncout, 0, "history", paste0(date(), ": Extracted from dwh using /user/bhendj/R/DWH_to_netcdf.R by ", Sys.getenv("USER")), prec='text')
##ncatt_put(ncout, 0, "history", Reduce(function(x,y) paste(x, "
##", y), as.list(c(
##  paste0(date(), ": Extracted from dwh using ~/R/extract_temperature_for_ML_project.R by ", Sys.getenv("USER")), 
##  capture.output(print(session_info()))))), prec='text')
nc_close(ncout)

q(save = 'no')
