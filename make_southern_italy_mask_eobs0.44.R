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

## southern Italy bounding box
lo <- c(13, 17, 20, 16, 13)
la <- c(40.5, 43, 39, 38, 40.5)

map("italy")
points(lat ~ lon, lola, pch=16, col=c(0,1)[(lsm > 0.5)*1 + 1])
lines(lo, la)

library(sp)
pip <- as.logical(point.in.polygon(lola$lon, lola$lat, lo, la)) & lsm > 0.5

points(lat ~ lon, lola, pch=16, col=c(NA, 2)[pip*1 + 1])

mask.nc <- ncvar_def(name="mask", dim=nc$dim, units='1')
ncout <- nc_create(mask.nc, file="southern_italy_mask_eobs0.44.nc")
ncvar_put(ncout, "mask", pip*1)
nc_close(ncout)
