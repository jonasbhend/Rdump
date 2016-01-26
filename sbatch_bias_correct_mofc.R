library(myhelpers)
library(biascorrection)

args <- commandArgs(TRUE)
filenum <- as.numeric(args[1])
method <- args[2]

print(filenum)
print(method)

## deparse method
crossval <- length(grep("crossval", method)) == 1
blocklength <- ifelse(crossval, as.numeric(gsub(".*crossval", "", method)), 1)
forward <- length(grep("forward", method)) == 1


## observations
obsfile <- '/store/msclim/bhendj/EUPORIAS/E-OBS/eobs0.44/daily/tas/tg_0.44deg_rot_v11.0.nc'
obsname <- "E-OBS"
ncobs <- nc_open(obsfile)
obs.time <- as.Date(nc_time(ncobs))
starti <- max(which(obs.time <= as.Date("1990-01-01")))
obs.time <- obs.time[starti:length(obs.time)]
otmp <- ncvar_get(ncobs, 'tg', start=c(1,1,starti)) + 273.15

## forecast files
mofcpath <- '/store/msclim/bhendj/MOFC/eobs0.44/daily/tas/none'
fcfiles <- list.files(mofcpath, full.names=T)
fcfile <- fcfiles[filenum]

## get data
nc <- nc_open(fcfile)
hdate <- as.Date(as.character(ncvar_get(nc, 'hdate')), format="%Y%m%d")
fcst <- ncvar_get(nc, 'tas')
fcst <- aperm(array(fcst, c(nrow(fcst), ncol(fcst), 20, 5, 32)), c(1,2,5,3,4))
fct <- as.Date(nc_time(nc))
fc.time <- outer(hdate[1:20], fct - fct[1], '+')
o.i <- obs.time %in% fc.time
obs <- array(otmp[,,o.i], dim(fcst)[-5])

## get year range
deb.years <- range(as.numeric(format(fc.time[,1], '%Y')))
## set-up method string
method.string <- paste(method, paste(deb.years, collapse='-'), obsname, sep="_")

## set up output directory
outdir <- gsub("none", method.string, mofcpath)
outfile <- gsub(".nc", paste0("_", method.string, ".nc"), basename(fcfile))

system(paste("mkdir -p", outdir))

## set up forecast and obs mask
mask <- apply(!is.na(obs), 1:2, all) & apply(!is.na(fcst), 1:2, all)

fcst.debias <- fcst*NA
for (i in 1:nrow(fcst)){
  print(i)
  for (j in 1:ncol(fcst)){
    if (mask[i,j]){
      fcst.debias[i,j,,,] <- debias(fcst=fcst[i,j,,,], obs=obs[i,j,,],
                                    method=gsub("-.*", "", method), 
                                    fc.time=fc.time, 
                                    crossval=crossval,
                                    blocklength=blocklength,
                                    forward=forward)
    }
  }
}

fcst.debias <- aperm(array(fcst.debias, c(dim(fcst)[1:3], prod(dim(fcst)[4:5]))), c(1,2,4,3))

## write output file
nc_write(fcfile, paste(outdir, outfile, sep='/'), "tas", fcst.debias)

print("succesfully reached end of script")

q(save = "no")
