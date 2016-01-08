library(ncdf4)
library(maps)

setwd("/store/msclim/bhendj/tmp")

files <- list.files(pattern='20130501_global..nc')
tp <- list()
for (f in files){
  nc <- nc_open(f)
  tptmp <- ncvar_get(nc, 'tp')
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  nc_close(nc)
  tpdiff <- tptmp
  for (i in 2:214) tpdiff[,,,i] <- tptmp[,,,i] - tptmp[,,,i-1]

  tp[[f]] <- as.data.frame(which(tpdiff < 0, arr=T))
  tp[[f]]$diff <- tpdiff[as.matrix(tp[[f]])]
  tp[[f]]$lon <- lon[tp[[f]]$dim1]
  tp[[f]]$lat <- lat[tp[[f]]$dim2]
  tp[[f]] <- tp[[f]][order(tp[[f]]$diff), ]
  rm(tptmp, tpdiff)
  gc()
}

nc <- nc_open("20130501_leaptst_gb.nc")
tptmp <- ncvar_get(nc, 'tp')
nc_close(nc)
tpdiff <- tptmp ; for (i in 2:214) tpdiff[,,,i] <- tptmp[,,,i] - tptmp[,,,i-1]


## plot rounding errors in dependance of lead time
png('~/tmp/rounding_errors_by_lead_time.png', width=6, height=10, units='in', res=200)
par(mfrow=c(2,1), mar=c(5,5,3,1))
plot(table(tp[[2]]$dim4), xlab='lead time (days)', ylab='Number of rounding errors', xaxt='n')
axis(1)
plot(tapply(abs(tp[[2]]$diff), tp[[2]]$dim4, mean), type='l', lwd=2,
     xlab='lead time (days)', ylab='Magnitude of average rounding error',)
dev.off()
