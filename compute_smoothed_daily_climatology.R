### compute_smoothed_daily_climatology.R --- 
## 
## Filename: compute_smoothed_daily_climatology.R
## Description: 
## Author: Jonas Bhend
## Maintainer: 
## Created: Wed Sep  3 14:33:12 2014 (+0200)
## Version: 
## Last-Updated: Wed Sep  3 17:09:26 2014 (+0200)
##           By: Jonas Bhend
##     Update #: 14
## URL: 
## Keywords: 
## Compatibility: 
## 
######################################################################
## 
### Commentary: 
## 
## 
## 
######################################################################
## 
### Change Log:
## 
## 
######################################################################
## 
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, write to
## the Free Software Foundation, Inc., 51 Franklin Street, Fifth
## Floor, Boston, MA 02110-1301, USA.
## 
######################################################################

## load packages and functions
library(ncdf)

## initialise
obspath <- '/store/msclim/Prediction/E-OBS/E-OBSv10.0'
varname <- 'tn'
period <- c(1981, 2012)
obsfile <- paste0(varname, '_0.22deg_rot_v10.0.nc')

## open netcdf file
nc <- open.ncdf(paste(obspath, obsfile, sep='/'))

## get NetCDF time
if (length(grep('days', nc$dim$time$units)) == 1){
  ntime <- as.Date(gsub('.*since ', '', nc$dim$time$units)) + nc$dim$time$vals
} else {
  stop('Do not know how to convert times')
}

alldat <- array(NA, c(nc$dim$x$len, nc$dim$y$len, 366, diff(period) + 1))
dimnames(alldat)[[4]] <- as.character(period[1]:period[2])
## loop through years
for (yy in as.character(seq(period[1], period[2]))){
  itime <- which(format(ntime, '%Y') == yy)
  dtmp <- get.var.ncdf(nc, varname, start=c(1,1,min(itime)), count=c(-1, -1, diff(range(itime)) + 1))
  ## replace with missing values
  dtmp[dtmp < -99] <- NA
  ## duplicate the Feb. 28 for non leap years
  if (dim(dtmp)[3] == 365) dtmp <- dtmp[,,sort(c(1:365, 59))]
  alldat[,,,yy] <- dtmp
}

## compute mean
datmn <- rowMeans(alldat, dims=3, na.rm=T)
## generate mask with at least 90% 
datna <- rowMeans(!is.na(alldat), dims=2)
datmn[rep(datna < 0.9, length=length(datmn))] <- NA

## fit Loess or other function to mean temperature
hf <- dnorm(seq(-3,3,length=91))
hf <- hf / sum(hf)
datfit2 <- aperm(apply(datmn, 1:2, filter, filter=hf, circular=TRUE), c(2,3,1))
## fit Loess function to repeated series
datfit <- aperm(apply(datmn, 1:2, function(x){
  if (any(is.na(x))){
    out <- rep(NA, length(x))
  } else {
    xdf <- data.frame(x=seq(1, 3*length(x)), y=rep(x, 3))
    xloess <- loess(y ~ x, xdf, span=151/nrow(xdf))
    ## cut out the middle part
    out <- xloess$fit[length(x) + seq(x)]
  }
  return(out)}), c(2,3,1))  
## is roughly 5 times slower than the filtering above



######################################################################
### compute_smoothed_daily_climatology.R ends here
