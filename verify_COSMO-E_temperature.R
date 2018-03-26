library(myhelpers)
library(tidyverse)
library(parallel)
library(easyVerification)

rfile <- "/store/msclim/bhendj/verification/Rdata/t2m_COSMO-E.Rdata"

ofile <- "/store/msclim/bhendj/ml_ce_zhang/full/observations/meteoswiss_t2m_20151001-20171114.nc"
ffs <- list.files("/store/msclim/bhendj/verification/COSMO-E", pattern=".nc", recursive=TRUE, full.names=TRUE)
fcst <- read_ncdf(ffs, 'T', n.cores=20)
ftime <- with_tz(attr(fcst, 'time'), 'UTC')
ftime <- array(ftime - ftime[1], dim(fcst)[3:4]) + ftime[1]

otmp <- read_ncdf(ofile, 't2m')

nc <- nc_open(ofile)
ll <- list()
for (nn in setdiff(names(nc$var), 't2m')){
  ll[[nn]] <- ncvar_get(nc, nn)
}
meta <- as.data.frame(ll)
nc_close(nc)

# range(apply(abs(fcst - rep(obs[,match(attr(fcst, 'time'), attr(obs, 'time'))], each=21) - 273.15), 2:3, mean, na.rm=T))
# range(apply(abs(fcst - rep(obs[,match(attr(fcst, 'time'), attr(obs, 'time')) - 1], each=21) - 273.15), 2:3, mean, na.rm=T))
# range(apply(abs(fcst - rep(obs[,match(attr(fcst, 'time'), attr(obs, 'time')) + 1], each=21) - 273.15), 2:3, mean, na.rm=T))
# range(apply(abs(fcst - rep(obs[,match(attr(fcst, 'time'), attr(obs, 'time')) + 2], each=21) - 273.15), 2:3, mean, na.rm=T))

obs <- array(otmp[,match(attr(fcst, 'time'), attr(otmp, 'time')) + 1], dim(fcst)[-1])

## change to degree Celsius
fcst <- fcst - 273.15

## reorder arrays for ease of computing
if (length(dim(fcst)) == 4){
  fcst <- aperm(array(fcst, c(21, 144, 121, 2, 62)), c(2,4,3,5,1))
  obs <- aperm(array(obs, c(144, 121, 2, 62)), c(1,3,2,4))  
  ftime <- aperm(array(ftime - ftime[1], c(121, 2, 62)), c(2,1,3)) + ftime[1]
}

bias <- apply(array(fcst - c(obs), c(144, 2, 121, 31, 2, 21)), c(1,2,3,5), mean, na.rm=T)
fdeb <- fcst - c(bias[,,,rep(2:1, each=31)])


## recalibration function
## compute the recalibration (single model)
recal <- function(fcst, obs, fcst.out = fcst, 
                  type = c('ensNGR', 'lm', "debias", 'NGR'),
                  crpsfun=FairCrps, initpars = NULL, 
                  maxit = 1e5){
  type <- match.arg(type)
  fmn <- rowMeans(fcst)
  fsd <- apply(fcst, 1, sd)
  fanom <- (fcst - fmn) / fsd
  fmn.out <- rowMeans(fcst.out)
  fsd.out <- apply(fcst.out, 1, sd)
  fanom.out <- (fcst.out - fmn.out) / fsd.out
  
  if (type == 'debias'){
    bias <- mean(fmn) - mean(obs)
    fout <- fcst.out - bias
  } else if (type == 'lm') {
    flm <- lm(obs ~ fcst, data.frame(obs = obs, fcst = fmn)) 
    plm <- predict(flm, newdata = data.frame(fcst=fmn.out),
                   interval='prediction', level= pnorm(1) - pnorm(-1))
    fout <- plm[,1] + -apply(plm[,1:2, drop=F], 1, diff) * fanom.out
    attr(fout, 'par') <- coef(flm)
  } else {
    scalefun <- function(pars, fmn, fsd, fanom){
      pars[1] + pars[2]*fmn + sqrt(pars[3]**2 + pars[4]**2 * fsd**2) * fanom
    }  
    
    ## switch for function
    if (type == 'ensNGR'){
      optfun <- function(pars, fmn, fsd, fanom, obs){
        mean(crpsfun(scalefun(pars, fmn, fsd, fanom), obs))
      }
    } else if (type == 'NGR'){
      optfun <- function(pars, fmn, fsd, fanom, obs){
        mean(GaussCrps(pars[1] + pars[2]*fmn, sqrt(pars[3]**2 + pars[4]**2*fsd**2), 
                       obs))
      }
    }
    if (is.null(initpars)) {
      initpars <- rep(NA, 4)
    }
    if (any(is.na(initpars))){
      flm <- lm(obs ~ fmn)
      suppressWarnings(plm <- -apply(predict(flm, interval='prediction', level = pnorm(1) - pnorm(-1))[,1:2], 1, diff))
      flm2 <- lm(plm ~ fsd)
      initpars[is.na(initpars)] <- c(coef(flm), coef(flm2))[is.na(initpars)]
    }
    if (maxit == 0){
      opars <- initpars
    } else {
      opars <- optim(initpars, optfun,
                     fmn=fmn, fsd=fsd, fanom=fanom, obs=obs,
                     control = list(maxit=maxit), method = 'BFGS')$par      
    }
    fout <- scalefun(opars, fmn.out, fsd.out, fanom.out)
    attr(fout, 'par') <- opars
  } 
  return(fout)
}




frecal <- fcst * NA
crecal <- fcst[,,,,1]*NA
for (si in 1:144){
  print(si)
  for (ii in 1:2){
    for (li in 1:121){
      for (ci in 1:2){
        ind <- intersect(1:31 + (ci - 1)*31, which(!is.na(obs[si,ii,li,])))
        if (length(ind) > 20){   
          ind2 <- setdiff(1:62, ind)
          frecal[si,ii,li,ind2,] <- recal(fcst[si,ii,li,ind,], obs[si,ii,li,ind], fcst.out=fcst[si,ii,li,ind2,], maxit=0)
          crecal[si,ii,li,ind2] <- predict(lm(obs ~ fcst, data.frame(obs = obs[si,ii,li,ind], fcst=fcst[si,ii,li,ind,1])),
                                           newdata=data.frame(fcst=fcst[si,ii,li,ind2,1]))
        }
      }
    }
  }
}



## compute the CRPS
crps <- veriApply("FairCrps", fcst, obs, parallel=T, maxncpus=20)
crps.deb <- veriApply("FairCrps", fcst - c(bias), obs, parallel=TRUE, maxncpus=20)
crps.recal <- veriApply("FairCrps", frecal, obs, parallel=TRUE, maxncpus=20)
mae <- abs(rowMeans(fcst, dims=4) - obs)
mae.deb <- abs(rowMeans(fdeb, dims=4) - obs)
mae.recal <- abs(rowMeans(frecal, dims=4) - obs)
mae0 <- abs(fcst[,,,,1] - obs)
mae0.deb <- abs(fcst[,,,,1] - obs  - c(rowMeans(fcst[,,,,1] - obs, dims=3, na.rm=T)))
mae0.recal <- abs(crecal - obs)
fmed <- apply(fcst, 1:4, function(y) sort(y)[11])
maemed <- abs(fmed - obs)
maemed.deb <- abs(fmed - obs  - c(rowMeans(fmed - obs, dims=3, na.rm=T)))
fopti <- apply(apply(abs(fcst - c(obs)), c(1,2,4,5), mean, na.rm=T), 1:3, function(y) if (all(!is.na(y))) which.min(y) else NA)
fopti.deb <- apply(apply(abs(fdeb - c(obs)), c(1,2,4,5), mean, na.rm=T), 1:3, function(y) if (all(!is.na(y))) which.min(y) else NA)
fopt <- fopt.deb <- fmed*NA
for (i in 1:144){
  for (j in 1:2){
    for (k in 1:62){
      fopt[i,j,,k] <- fcst[i,j,,k,fopti[i,j,k]]
      fopt.deb[i,j,,k] <- fdeb[i,j,,k,fopti.deb[i,j,k]]
    }
  }
}
maeopt <- abs(fopt - obs)
maeopt.deb <- abs(fopt.deb - obs)

plot(apply(crps.deb, 3, mean, na.rm=T), lwd=2, type='l', ylim=c(1, 2.5))
lines(apply(mae.deb, 3, mean, na.rm=T), lwd=2, col=2)
lines(apply(mae0.deb, 3, mean, na.rm=T), lwd=2, col=3)
lines(apply(maemed.deb, 3, mean, na.rm=T), lwd=2, col=4)
lines(apply(maeopt.deb, 3, mean, na.rm=T), lwd=2, col=5)


scrps1 <- 1 - apply(mae0.recal, c(1,2,4), mean) / apply(mae0.deb, c(1,2,4), mean)
scrps2 <- 1 - apply(crps.recal, c(1,2,4), mean) / apply(crps.deb, c(1,2,4), mean)

mi <- which(sqrt(scrps1) + sqrt(scrps2) > max(sqrt(scrps1) + sqrt(scrps2) - .2, na.rm=T), arr.ind=TRUE)

casei <- 11
si <- mi[casei,1]
ii <- mi[casei,2]
fi <- mi[casei,3]

png("~/tmp/cosmo-e_veri1.png", width=7, height=5, units='in', res=200)
xdate <- ftime[ii,,fi]
par(las=1, cex.axis=0.83, cex.lab=1,
    mgp=c(2, 0.5, 0), tcl=0.3, mar=c(2,3,3,0.5), mfrow=c(1,1), oma=rep(0, 4))
plot(xdate, fdeb[si,ii,,fi,1], type='l', lwd=2, lty=2, 
     ylim=range(13, crecal[si,ii,,fi], fdeb[si,ii,,fi,1], obs[si,ii,,fi]), 
     col=grey(0.3), 
     xlab='', xaxt='n',
     ylab = 'Temperatur (deg. C)', 
     main = paste("COSMO-E Vorhersage,", meta$name[si]))
axis.POSIXct(1, at = xdate[seq(1,121,24)], format = "%d.%m.%y")
lines(xdate, crecal[si,ii,,fi], lwd=2, col=grey(0.3))
lines(xdate, obs[si,ii,,fi], col=2, lwd=2)
legend("topright", 
       paste0(c("COSMO-E control,\t\t", "COSMO-E control PP,\t"),  " MAE = ", 
              round(c(mean(mae0.deb[si,ii,,fi]), mean(mae0.recal[si,ii,,fi])), 2)), 
       lwd=2, lty=2:1, bty='n')
dev.off()

png("~/tmp/cosmo-e_veri2.png", width=7, height=5, units='in', res=200)
par(las=1, cex.axis=0.83, cex.lab=1,
    mgp=c(2, 0.5, 0), tcl=0.3, 
    mar=rep(0.2, 4), 
    oma=c(1.8,2.8,2.8,0.3),mfrow=c(2,1))
matplot(x=xdate, frecal[si,ii,,fi,21:1], 
        type='l', lwd=2, lty=1, col=grey(0.8), 
        ylim=range(13, frecal[si,ii,,fi,], fdeb[si,ii,,fi,], obs[si,ii,,fi]), 
        xlab='', xaxt='n',
        ylab = 'Temperatur (deg. C)')
lines(xdate, obs[si,ii,,fi], col=2, lwd=2)
matplot(x=xdate, fdeb[si,ii,,fi,], 
        type='l', lwd=2, lty=1, col=grey(0.8), 
        ylim=range(13, frecal[si,ii,,fi,], fdeb[si,ii,,fi,], obs[si,ii,,fi]), 
        xlab='', xaxt='n',
        ylab = 'Temperatur (deg. C)')
axis.POSIXct(1, at = xdate[seq(1,121,24)], format = "%d.%m.%y")
lines(xdate, obs[si,ii,,fi], col=2, lwd=2)
mtext(side=3, outer=TRUE, line=0.8, cex=1/0.8333, font=2,
      text = paste("COSMO-E Vorhersage,", meta$name[si]))
mtext(side=2, outer=TRUE, line=1.8, cex=1, font=1,
      text = 'Temperatur (deg. C)', las=3)
dev.off()

png("~/tmp/cosmo-e_veri3.png", width=7, height=5, units='in', res=200)
par(las=1, cex.axis=0.83, cex.lab=1,
    mgp=c(2, 0.5, 0), tcl=0.3, 
    mar=rep(0.2, 4), 
    oma=c(1.8,2.8,2.8,0.3),mfrow=c(2,1))
matplot(x=xdate, frecal[si,ii,,fi,21:1], 
        type='l', lwd=2, lty=1, col=grey(0.8), 
        ylim=range(13, frecal[si,ii,,fi,], fdeb[si,ii,,fi,], obs[si,ii,,fi]), 
        xlab='', xaxt='n',
        ylab = 'Temperatur (deg. C)')
lines(xdate, crecal[si,ii,,fi], lwd=2, col=grey(0.3))
lines(xdate, obs[si,ii,,fi], col=2, lwd=2)
legend("bottomright", 
       paste0(c("COSMO-E control PP,\t\tMAE", "COSMO-E PP,\t\t\tCRPS"),  " = ", 
              round(c(mean(mae0.recal[si,ii,,fi]),
                      mean(crps.recal[si,ii,,fi])), 2)), 
       lwd=2, lty=1, col=c(grey(0.3),grey(0.8)), bty='n')
matplot(x=xdate, fdeb[si,ii,,fi,], 
        type='l', lwd=2, lty=1, col=grey(0.8), 
        ylim=range(13, frecal[si,ii,,fi,], fdeb[si,ii,,fi,], obs[si,ii,,fi]), 
        xlab='', xaxt='n',
        ylab = 'Temperatur (deg. C)')
axis.POSIXct(1, at = xdate[seq(1,121,24)], format = "%d.%m.%y")
lines(xdate, fdeb[si,ii,,fi,1], lwd=2, col=grey(0.3))
lines(xdate, obs[si,ii,,fi], col=2, lwd=2)
legend("bottomright", 
       paste0(c("COSMO-E control,\t\tMAE", "COSMO-E,\t\t\tCRPS"),  " = ", 
              round(c(mean(mae0.deb[si,ii,,fi]),
                      mean(crps.deb[si,ii,,fi])), 2)), 
       lwd=2, lty=1, col=c(grey(0.3),grey(0.8)), bty='n')
mtext(side=3, outer=TRUE, line=0.8, cex=1/0.8333, font=2,
      text = paste("COSMO-E Vorhersage,", meta$name[si]))
mtext(side=2, outer=TRUE, line=1.8, cex=1, font=1,
      text = 'Temperatur (deg. C)', las=3)
dev.off()



png("~/tmp/cosmo-e_veri2.png", width=7, height=5, units='in', res=200)
par(las=1, cex.axis=0.83, cex.lab=1,
    mgp=c(2, 0.5, 0), tcl=0.3, mar=c(2,3,3,0.5))
matplot(x=xdate, fdeb[si,ii,,fi,21:1], 
        type='l', lwd=2, lty=1, col=c(rep(grey(0.8), 20), grey(0.3)), 
        ylim=range(fcst[si,ii,,fi,1], fdeb[si,ii,,fi,], obs[si,ii,,fi]), 
        xlab='', xaxt='n',
        ylab = 'Temperatur (deg. C)', 
        main = paste("COSMO-E Vorhersage,", meta$name[si]))
axis.POSIXct(1, at = xdate[seq(1,121,24)], format = "%d.%m.%y")
lines(xdate, fcst[si,ii,,fi,1], lwd=2, col=grey(0.3), lty=2)
lines(xdate, obs[si,ii,,fi], col=2, lwd=2)
legend("topright", 
       paste0(c("COSMO-E control,\t\tMAE", "COSMO-E control PP,\tMAE", "COSMO-E PP,\t\tCRPS"),  " = ", 
              round(c(mean(mae0[si,ii,,fi]), mean(mae0.deb[si,ii,,fi]), mean(crps.deb[si,ii,,fi])), 2)), 
       lwd=2, lty=c(2:1,1), col=c(rep(grey(0.3), 2),grey(0.8)), bty='n')
dev.off()


## find a case with exceptionally low temperature
mmi <- unique(rbind(which(obs[si,ii,,] > 9.6, arr=T), which(fdeb[si,ii,,,] > 9.6, arr=T)[,1:2]))
mmi <- mmi[mmi[,1] > 30 & mmi[,1] < 60,]
ffi <- lapply(1:nrow(mmi), function(i) list(fcst = fcst[si,ii,mmi[i,1], mmi[i,2],] + 7, 
                                            fdeb = fdeb[si,ii,mmi[i,1], mmi[i,2],], 
                                            obs=obs[si,ii,mmi[i,1], mmi[i,2]]))

png("~/tmp/forecasters_dilemma.png", width=7, height=4, units='in', res=200)
par(las=1, cex.axis=0.83, cex.lab=1,
    mgp=c(2, 0.5, 0), tcl=0.3, 
    mar=c(0.5, 3, 0.5, 0.5))
boxplot(at = seq_along(ffi) - 0.15, lapply(ffi, function(x) x$fcst), 
        ylim=range(ffi), xlim=c(0.5, length(ffi) + 0.5), xaxs='i',  
        col="#55550022", border="#555500", lty=1, range=0, 
        xaxt = 'n', outline=FALSE, boxwex = 0.3, ylab="Temperature (deg. C)")
boxplot(at = seq_along(ffi) + 0.15, lapply(ffi, function(x) x$fdeb), 
        add = TRUE, col="#00555522", border="#005555", lty=1, range=0, 
        axes=FALSE, xlab='', ylab='', outline=FALSE, boxwex=0.3)
points(sapply(ffi, function(x) x$obs), pch=21, cex=1.4, bg=2, lwd=2)
abline(v=seq(1.5, length(ffi) - 0.5), lty=3, col='grey')
abline(v=4.5, lty=2, lwd=2)
dev.off()


png("~/tmp/local.png", width=3, height=3, units="in", res=200)
ind <- seq(-3,5,0.01)
plot(ind, dnorm(ind, mean=-1, sd=0.5), type='l', lwd=1, axes=FALSE, 
     xlab='', ylab='', col='darkred', yaxs='i')
lines(ind, dnorm(ind, mean=2, sd=2), col='darkblue', lwd=1)
abline(v=-0.03, lwd=1)
dev.off()


