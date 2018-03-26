
library(dwhget)
library(tidyverse)
library(ncdf4)

## temperature hourly values
param <- "tre200b0"
mstype <- "ANETZ"

## get all stations data
startdate <- as.Date("2015-09-30")
enddate <- Sys.Date() - 1
t2m <- dwhget(stats = 'all', mstypes=NULL, pars=param, 
              t.beg="2015.10.01", 
              t.end="2015.10.31")
stations <- unique(as.character(formatC(na.omit(t2m$station), width=4, flag=0)))

## get station meta data
stat.meta <- dwh.stat.info(stats = stations)
stat.meta$name <- iconv(stat.meta$name, from = 'latin1', to = 'UTF-8')
stat.meta <- stat.meta[order(stat.meta$nat_ind),]

## temperature data
t2m <- dwhget(stats=stations, 
              mstypes = NULL,
              pars=param, 
              t.beg=format(startdate, "%Y.%m.%d"), 
              t.end=format(enddate - 1, "%Y.%m.%d"))

t2m.wide <- as_tibble(t2m) %>%
  spread(station, value)
t2m.arr <- as.matrix(select(t2m.wide, -parameter, -time))

## check that ordering of meta data is correct
stopifnot(as.character(stat.meta$nat_ind) == names(t2m.wide)[-(1:2)])


## Trying to follow CF Conventions 1.7

## NetCDF Dimensions - fix timestamp to denote endtime
time.nc <- ncdim_def("time", 
                     paste0("hours since ", startdate, " 00:00:00"), 
                     1:nrow(t2m.wide) - 2, 
                     unlim=TRUE, longname = 'time of measurement')
station.nc <- ncdim_def('station', '1', 1:nrow(stat.meta))
nchar.nc <- ncdim_def("nchar", "", 1:max(nchar(stat.meta$name)), create_dimvar=FALSE)
nchar3.nc <- ncdim_def("nchar3", "", 1:3, create_dimvar=FALSE)
nchar4.nc <- ncdim_def("nchar4", "", 1:4, create_dimvar=FALSE)

## NetCDF Variables
vars.nc <- list()
for (nn in names(stat.meta)[c(1,2,6)]) {
  vars.nc[[nn]] <- ncvar_def(nn, "", list(list(nat_ind = nchar4.nc, nat_abbr = nchar3.nc, name = nchar.nc)[[nn]], station.nc),
                             prec = 'char')
}
for (nn in names(stat.meta)[3:5]){
  vars.nc[[nn]] <- ncvar_def(name = nn, 
                             units = c(lon = "degrees_east", lat = "degrees_north", height = 'm')[nn],
                             dim = station.nc,
                             longname = c(lon = "station longitude",
                                          lat = "station latitude", 
                                          height = "station altitude")[nn])
}
vars.nc$t2m <- ncvar_def("t2m", "deg. C", list(time.nc, station.nc), missval=-1e20)


## open file and write to it
ncfile <- paste0("/store/msclim/bhendj/verification/obs/t2m_SMNET_", 
                 format(startdate, "%Y%m%d"),
                 '-',
                 format(enddate, "%Y%m%d"), 
                 ".nc")
if (file.exists(ncfile)){
  file.remove(ncfile)
} else {
  dir.create(dirname(ncfile))
}
ncout <- nc_create(ncfile, 
          vars.nc, 
          force_v4 = TRUE)
for (nn in names(stat.meta)){
  ncvar_put(ncout, nn, stat.meta[[nn]])
}
ncvar_put(ncout, "t2m", t2m.arr)

## additional attributes
standard_names <- c(lon = 'longitude', lat='latitude', height='altitude', t2m = "air_temperature")
for (sn in names(standard_names)) ncatt_put(ncout, sn, "standard_name", standard_names[sn], prec = 'text')
ncatt_put(ncout, 't2m', "coordinates", "lat lon height name", prec = 'text')
ncatt_put(ncout, 'nat_ind', "cf_role", "timeseries_id", prec='text')

## put global attributes
ncatt_put(ncout, 0, 'Conventions', "CF-1.7", prec='text')
ncatt_put(ncout, 0, "featureType", "timeSeries", prec='text')
ncatt_put(ncout, 0, "history", paste0(date(), ": Extracted from dwh using ~/R/extract_temperature_for_ML_project.R by ", Sys.getenv("USER")), prec='text')
##ncatt_put(ncout, 0, "history", Reduce(function(x,y) paste(x, "
##", y), as.list(c(
##  paste0(date(), ": Extracted from dwh using ~/R/extract_temperature_for_ML_project.R by ", Sys.getenv("USER")), 
##  capture.output(print(session_info()))))), prec='text')
nc_close(ncout)
