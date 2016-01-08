library(ncdf4)
library(geocors)
library(maps)

nc <- nc_open("/store/msclim/bhendj/EUPORIAS//grids/eobs0.44_lsm.nc")
rlon <- nc$dim$rlon$vals
rlat <- nc$dim$rlat$vals
lsm <- ncvar_get(nc, "FR_LAND")

lola <- geocors.trafo(rep(rlon, length(rlat)), 
                      rep(rlat, each=length(rlon)),
                      from.type='rotpol',
                      from.pars=list(plon=-162, plat=39.25),
                      to.type='lonlat')


uk <- map(region = 'UK', plot=FALSE)

library(sp)
ii <- uk$x > -10 & uk$y > 40
pip <- as.logical(point.in.polygon(lola$lon, lola$lat, uk$x[ii], uk$y[ii])) & lsm > 0.5

map(xlim=c(-10,10), ylim=c(40,60))
points(lat ~ lon, lola, pch=16, col=c(NA, 2)[pip*1 + 1])


setwd("/store/msclim/bhendj/EUPORIAS/grids")
mask.nc <- ncvar_def(name="mask", dim=nc$dim, units='1')
ncout <- nc_create(mask.nc, file="UK_mask_eobs0.44.nc")
ncvar_put(ncout, "mask", pip*1)
nc_close(ncout)
