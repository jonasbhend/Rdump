## compute seasonal percentiles from infile
library(ncdf4)
source("~/R/ncutils.R")

args <- commandArgs(TRUE)
infile <- args[1]
pctl <- args[2]
outfile <- args[3]

nc <- nc_open(infile)
dtmp <- ncvar_get(nc, names(nc$var)[1])
dtime <- nc_time(nc)
nc_close(nc)

## get seasonal index
seasi <- function(time){
  mon <- as.numeric(format(time, '%m'))
  seasstr <- c("DJF", 'MAM', 'JJA', 'SON')
  return(seasstr[(mon %/% 3) %% 4 + 1 ] )
}

dpctl <- apply(dtmp, 1:2, tapply, seasi(dtime), quantile, as.numeric(pctl)/100)

## find leap year
tyear <- names(which.max(table(format(dtime, '%Y'))))

## set up netcdf template
system(paste0("cdo -s seldate,", tyear, '-01-01,', tyear, '-12-31 ', infile, ' ', infile, '.tmp'))
nctmp <- nc_open(paste0(infile, '.tmp'))

## get time
tmptime <- nc_time(nctmp)

## write data
nc_write(paste0(infile, '.tmp'), outfile, names(nctmp$var)[1], data=aperm(dpctl[seasi(tmptime), , ], c(2,3,1)))

file.remove(paste0(infile, '.tmp'))

q(save='no')
