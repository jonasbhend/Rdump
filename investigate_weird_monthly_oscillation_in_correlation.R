# this script is used to investigate monthly oscillations in 
# correlation of ECMWF forecasts with E-OBS (and ERA-INT to a lesser extent)
library(ncdf4)
library(geocors)
rlola <- geocors.trafo(x=13.9, y=61, from.type='lonlat', to.type='rotpol', 
                       to.pars=list(plon=-162, plat=39.25))

dpath <- '/store/msclim/bhendj/EUPORIAS/ecmwf-system4/eobs0.44/daily/tas/none'
files <- list.files(dpath, pattern="^....05.._.*nc$", full.names=TRUE)

nc <- nc_open(files[1])
rlon <- nc$dim$rlon$vals
rlat <- nc$dim$rlat$vals
loni <- which.min((rlon - rlola$rlon)**2)
lati <- which.min((rlat - rlola$rlat)**2)
nc_close(nc)

fcst <- array(NA, c(214, length(files), 51))
for (i in seq(files)){
  print(i)
  nc <- nc_open(files[i])
  fcst[,i,] <- t(ncvar_get(nc, "mean2t24", start=c(loni, lati,1,1), count=c(1,1,-1,-1)))
  nc_close(nc)
}
fcst.mn <- rowMeans(fcst, dims=2)

ofile <- "/store/msclim/bhendj/EUPORIAS/E-OBS/eobs0.44/daily/tas/tg_0.44deg_rot_v10.0.nc"
nc <- nc_open(ofile)
obs.tmp <- ncvar_get(nc, 'tg', start=c(loni, lati, 1), count=c(1,1,-1)) 
otime <- as.Date("1950-01-01") + nc$dim$time$vals
obs.i <- which(format(otime, '%m%d') == "0501" & as.numeric(format(otime, '%Y')) %in% 1981:2013)
obs <- array(obs.tmp[outer(0:213, obs.i, '+')], c(214, 33))
nc_close(nc)

f.corr <-diag(cor(t(obs), t(fcst.mn))) 
plot(f.corr, type='l', lwd=2, xlab='days after May 1', ylab='correlation')
abline(h=0, lty=2, lwd=2)
lines(filter(f.corr, rep(1/30,30)), lwd=2, col=2)
## abline(v=seq(16, 214, 30))

## is it the mean or standard deviation?
plot(apply(obs, 1, sd), type='l', lwd=2)
lines(sqrt(rowMeans(apply(fcst, c(1,3), sd)**2)), lwd=2, col=2)


a.obs <- pacf(as.vector(obs), plot=F)$acf
a.fcst <- array(NA, c(length(a.obs), 51))
for (i in 1:51) a.fcst[,i] <- pacf(as.vector(fcst[,,i]), plot=F)$acf
matplot(a.fcst, pch=15, col='grey')
lines(a.obs, type='h', lwd=3, lend=3, col=-sign(a.obs)+3)

## aggregate forecast and observations
ftime <- as.Date("1990-05-01") - 1 + seq(1:nrow(fcst))
fcst.mon <- apply(fcst, 2:3, tapply, format(ftime, '%m'), mean)
obs.mon <- apply(obs, 2, tapply, format(ftime, '%m'), mean)
fcst.mon.mn <- rowMeans(fcst.mon, dims=2)
fcst.mon.mn2 <- rowMeans(fcst.mon[,,1:15], dims=2)


## plot correlation for individual months (scatterplot)
library(ggplot2)
of.df <- data.frame(fcst=c(fcst.mon.mn - rowMeans(fcst.mon.mn)), obs=c(obs.mon - rowMeans(obs.mon)), lead=paste("lead time", 1:7), year=rep(1981:2013, each=7))
ggplot(subset(of.df, lead != "lead time 7"), aes(x=fcst, y=obs, colour=year)) + 
  geom_point() +
  geom_smooth(method='lm') + 
  facet_wrap( ~ lead, ncol=2)

barplot(rbind(apply(fcst.mon.mn, 1, function(x) lm(x ~ seq(x))$coef[2]),
              apply(obs.mon, 1, function(x) lm(x ~ seq(x))$coef[2])),
        beside=TRUE)

fcst.detrend <- apply(fcst.mon.mn, 1, function(x) x - lm(x ~ seq(x))$fit)
obs.detrend <- apply(obs.mon, 1, function(x) x - lm (x ~ seq(x))$fit)
f.corr <- diag(cor(t(obs.mon), t(fcst.mon.mn)))
f.corr2 <- diag(cor(t(obs.mon[,1:30]), t(fcst.mon.mn[,1:30])))
fdetrend.corr <- diag(cor(obs.detrend, fcst.detrend))
barplot(rbind(f.corr, f.corr2), beside=TRUE)



fcst.seas <- apply(fcst.mon, 2:3, filter, rep(1/3,3))
obs.seas <- apply(obs.mon, 2, filter, rep(1/3,3))
fcst.seas.mn <- rowMeans(fcst.seas[,,1:15], dims=2)
f.corr <- diag(cor(t(obs.seas), t(fcst.seas.mn)))
barplot(f.corr)






## get dataset for comparison with ECOMS-UDG
lon0 <- 13.9
lat0 <- 61

fpath <- "/store/msclim/sysclim/ecmwf/system4/daily/deg075/Tavg"
ffiles <- list.files(fpath, pattern="^....0501_55_eobs.nc", full=TRUE)

nc <- nc_open(ffiles[1])
lon <- nc$dim$lon$vals
lat <- nc$dim$lat$vals
loni <- which.min((lon - lon0)**2)
lati <- which.min((lat - lat0)**2)
nc_close(nc)

fcst2 <- array(NA, c(214, length(ffiles), 51))
for (i in seq(ffiles)){
  nc <- nc_open(ffiles[i])
  fcst2[,i,] <- t(ncvar_get(nc, "mean2t24", start=c(loni, lati, 1, 1), count=c(1,1,-1,-1)))
  nc_close(nc)
}
obs2 <- obs

save(obs2, fcst2, file='/store/msclim/bhendj/tmp/system4_forecasts.Rdata')


