source("/users/bhendj/.Rprofile")
library(myhelpers)
library(easyVerification)
library(devtools)

args <- commandArgs(TRUE)
if (length(args) == 0){
  args <- c("/store/msclim/bhendj/MOFC/eobs0.44/daily/tas/none/20140508_tas_eobs0.44.nc")
}
fcfile <- args[1]

## deparse method and file
mofcpath <- dirname(fcfile)
method <- gsub(".*/", "", mofcpath)
fileparts <- strsplit(fcfile, '/')[[1]]
granularity <- fileparts[length(fileparts) - 3]

## set scores
scores <- rev(c("FairCrps", "FairRps", "EnsMe"))
scoretxt <- c(EnsMe="Difference with ensemble mean",
              FairRps="Fair ranked probability score for tercile category forecasts",
              FairCrps="Fair continuous ranked probability score")
scorefun <- function(score, fcst, obs, ...){
  if (score == "EnsMe"){
    out <- rowMeans(fcst, dims=length(dim(obs))) - obs 
  } else {
    out <- veriApply(score, fcst=fcst, obs=obs, prob=if(score == "FairRps") 1:2/3 else NULL)
  }
  return(out)
}

## observations
obsfile <- list.files(paste0('/store/msclim/bhendj/EUPORIAS/E-OBS/eobs0.44/', granularity, '/tas'),
                      pattern = ".nc", full.names=TRUE)
obsname <- "E-OBS"
ncobs <- nc_open(obsfile)
obs.time <- as.Date(nc_time(ncobs))
starti <- max(which(obs.time <= as.Date("1990-01-01")))
obs.time <- obs.time[starti:length(obs.time)]
otmp <- ncvar_get(ncobs, 'tg', start=c(1,1,starti)) + 273.15


## get data
nc <- nc_open(fcfile)
hdate <- as.Date(as.character(ncvar_get(nc, 'hdate')), format="%Y%m%d")
fcst <- ncvar_get(nc, 'tas')
fcst <- aperm(array(fcst, c(nrow(fcst), ncol(fcst), 20, 5, dim(fcst)[4])), c(1,2,5,3,4))
fct <- as.Date(nc_time(nc))
fc.time <- outer(hdate[1:20], fct - fct[1], '+')
o.i <- obs.time %in% fc.time
obs <- array(otmp[,,o.i], dim(fcst)[-5])

## set up output directory
outdir <- gsub("MOFC", "MOFC_scores", mofcpath)
system(paste("mkdir -p", outdir))
outfile <- paste("scores", basename(fcfile), sep='_')


## set up output file
hcstart.nc <- ncdim_def(name="hcstart", "days since 1980-01-01",
                      vals=as.numeric(fc.time[,1] - as.Date("1980-01-01")),
                      unlim=TRUE)
time.nc <- nc$dim$time
time.nc$unlim <- FALSE
dims.nc <- list(rlon=nc$dim$rlon, rlat=nc$dim$rlat, time=time.nc, hcstart=hcstart.nc)
vars.nc <- list(rotated_pole=nc$var[['rotated_pole']])
for (score in scores){
  vars.nc[[score]] <- ncvar_def(score, if (score == 'FairRps') "1" else "K", dim=dims.nc, 
                                missval=-1e20, longname=scoretxt[score])
}
ncout <- nc_create(paste(outdir, outfile, sep='/'), vars.nc)

## run scores
for (score in scores){
  print(score)
  ncvar_put(ncout, score, vals=scorefun(score, fcst=fcst, obs=obs))
}
nc_close(ncout)

## add documentation in attributes
tp <- session_info()
ev <- tp$packages[tp$packages$package == "easyVerification", ]
system(paste0("ncatted -h -a comment,global,a,c,'\nScores computed with ", ev$package, " version ", ev$version, " installed on ", ev$date, "' ", paste(outdir, outfile, sep='/')))

if (method == 'none'){
  outdir <- gsub("none", "refclim", outdir)
  system(paste0("mkdir -p ", outdir))
  outfile <- paste0("refclim_", outfile)
  
  nobs <- dim(obs)[4]
  obs.ref <- array(obs[,,,rep(1:nobs, each=nobs)], c(dim(obs), nobs))
  
  ncout <- nc_create(paste(outdir, outfile, sep='/'), vars.nc)
  
  ## run scores
  for (score in scores){
    print(score)
    ncvar_put(ncout, score, vals=scorefun(score, fcst=obs.ref, obs=obs))
  }
  nc_close(ncout)

  ## add documentation in attributes
  system(paste0("ncatted -h -a comment,global,a,c,'\nScores computed with ", ev$package, " version ", ev$version, " installed on ", ev$date, "' ", paste(outdir, outfile, sep='/')))
  
}

print("successfully reached end of script")

q(save = "no")
