library(myhelpers)
library(biascorrection)
library(devtools)

args <- commandArgs(TRUE)
if (length(args) == 0){
  args <- c("/store/msclim/bhendj/MOFC/eobs0.44/weekly/tas/none/20140508_tas_eobs0.44.nc", "smooth-crossval1")
}
fcfile <- args[1]
method <- args[2]

print(fcfile)
print(method)

## deparse method and file
crossval <- length(grep("crossval", method)) == 1
blocklength <- ifelse(crossval, as.numeric(gsub(".*crossval", "", method)), 1)
forward <- length(grep("forward", method)) == 1
strategy <- list(type=c("none", "crossval", "forward")[1 + 1*crossval + 2*forward],
                 blocklength=blocklength)

mofcpath <- dirname(fcfile)
fileparts <- strsplit(fcfile, '/')[[1]]
granularity <- fileparts[length(fileparts) - 3]

## observations
obsfile <- list.files(paste0('/store/msclim/bhendj/EUPORIAS/E-OBS/eobs0.44/', granularity, '/tas'),
                      pattern='.nc', full.names=TRUE)
obsname <- "E-OBS"
ncobs <- nc_open(obsfile)
obs.time <- as.Date(nc_time(ncobs))
starti <- max(which(obs.time <= as.Date("1990-01-01")))
obs.time <- obs.time[starti:length(obs.time)]
obsvar <- names(ncobs$var)[grep("tg|tas", names(ncobs$var))]
otmp <- ncvar_get(ncobs, obsvar, start=c(1,1,starti)) + 273.15


## get data
nc <- nc_open(fcfile)
hdate <- as.Date(as.character(ncvar_get(nc, 'hdate')), format="%Y%m%d")
fcst <- ncvar_get(nc, 'tas')
fcst <- aperm(array(fcst, c(nrow(fcst), ncol(fcst), 20, 5, dim(fcst)[4])), c(1,2,5,3,4))
fct <- as.Date(nc_time(nc))
fcorig <- as.Date(substr(basename(fcfile), 1, 8), format='%Y%m%d')
fc.time <- outer(hdate[1:20], fct - fcorig, '+')
o.i <- obs.time %in% fc.time
obs <- array(otmp[,,o.i], dim(fcst)[-5])

## get year range
deb.years <- range(as.numeric(format(fc.time[,1], '%Y')))
## set-up method string
method.string <- paste(method, obsname, sep="_")

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
                                    fc.time=t(fc.time), 
                                    strategy=strategy)
    }
  }
}

fcst.debias <- aperm(array(fcst.debias, c(dim(fcst)[1:3], prod(dim(fcst)[4:5]))), c(1,2,4,3))

## write output file
nc_write(fcfile, paste(outdir, outfile, sep='/'), "tas", fcst.debias)

## add documentation in attributes
tp <- session_info()
biasversion <- tp$packages[tp$packages$package == "biascorrection", ]
system(paste0("ncatted -h -a comment,global,o,c,'Bias corrected using method ", method, " against ", obsname, "\nfrom package ", biasversion$package, " version ", biasversion$version, " installed on ", biasversion$date, "\nusing hindcasts initialized in ", paste(deb.years, collapse='-'), "' ", paste(outdir, outfile, sep='/')))

print("successfully reached end of script")

q(save = "no")
