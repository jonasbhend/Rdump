library(ncdf4)
library(myhelpers)

date <- 20170720

fcfile <- paste0("/store/mch/msopr/ifs-ens-ext-exchange/", date, "/IFS_EXT_FC_", date,"_T_2M.nc")

zlon <- 8.55
zlat <- 47.4

nc <- nc_open(fcfile)
i <- which.min((nc$dim$x_1$vals - zlon)**2)
j <- which.min((nc$dim$y_1$vals - zlat)**2)
t2m <- ncvar_get(nc, 'T_2M', start=c(i,j,1,1), count=c(1,1,-1,-1)) - 273.15
ntime <- nc_time(nc)
nc_close(nc)

t2m <- t2m[,1:32]


t2m.mn <- t(apply(t2m, 1, filter, rep(1/7,7)))
t2m.mn[,1:7] <- NA
t2m.quant <- apply(t2m.mn, 2, quantile, c(0, 0.25, 0.5, 0.75, 1), na.rm=T)

xind <- as.Date(as.character(date), '%Y%m%d') + 1:32

png("~/tmp/mofc1.png", width=7, height=5, units='in', res=300)
par(cex.axis=0.83, cex.lab=1, mar=c(4,4,3,0.1), mgp=c(2, 0.7, 0))
matplot(x=xind, t(t2m), type='l', lwd=2, col="#55555555", lty=1,
        las = 1, 
        ylab="Temperatur in deg. C",
        xlab="Datum",
        main="Monatsvorhersage für Zürich",
        ylim = range(t2m[,1:32]), 
        xaxt='n')
axis.Date(1, at=xind[seq(1,32,7)], format='%d.%m.')
grid(nx = NA, ny=NULL)
dev.off()

png("~/tmp/mofc2.png", width=7, height=5, units='in', res=300)
par(cex.axis=0.83, cex.lab=1, mar=c(4,4,3,0.1), mgp=c(2, 0.7, 0))
matplot(x=xind, t(t2m), type='l', lwd=2, col="#cccccc55", lty=1,
        las = 1, 
        ylab="Temperatur in deg. C",
        xlab="Datum",
        main="Monatsvorhersage für Zürich",
        ylim = range(t2m[,1:32]), 
        xaxt='n')
axis.Date(1, at=xind[seq(1,32,7)], format='%d.%m.')
grid(nx = NA, ny=NULL)
matplot(xind[seq(8,32,7)],t(t2m.mn[,seq(8, 32, 7)]), 
        type='l', lwd=2, col='#33333355', lty=1, add=T, pch=16)
dev.off()

for (j in 1:4){
  png(paste0("~/tmp/mofc", j + 2, ".png"), width=7, height=5, units='in', res=300)
  par(cex.axis=0.83, cex.lab=1, mar=c(4,4,3,0.1), mgp=c(2, 0.7, 0))
  matplot(x=xind, t(t2m), type='l', lwd=2, col="#cccccc55", lty=1,
          las = 1, 
          ylab="Temperatur in deg. C",
          xlab="Datum",
          main="Monatsvorhersage für Zürich",
          ylim = range(t2m[,1:32]), 
          xaxt='n')
  axis.Date(1, at=xind[seq(1,32,7)], format='%d.%m.')
  grid(nx = NA, ny=NULL)
  for (i in 1:j){
    segments(x0=xind[i*7 + 1], y0=t2m.quant[1,i*7 + 1], y1=t2m.quant[5,i*7 + 1], 
             lwd=5, lend=3, col=hcl((i- 0.5)/10*360, l=90, c=80))  
    segments(x0=xind[i*7 + 1], y0=t2m.quant[2,i*7 + 1], y1=t2m.quant[4,i*7 + 1], 
             lwd=100, lend=3, col=hcl((i- 0.5)/10*360, l=70, c=80))  
    segments(x0=xind[i*7 + 1], y0=t2m.quant[3,i*7 + 1] - 0.1, y1=t2m.quant[3,i*7 + 1] + 0.1, 
             lwd=100, lend=3, col=hcl((i- 0.5)/10*360, l=50, c=80))  
  }
  dev.off()
}

