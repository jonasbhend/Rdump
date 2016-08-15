## script to compute skill metrics from debiased or raw forecast metrics

## load packages
library(methods) ## for processing with Rscript
library(ncdf4) ## ncdf is deprecated but still supported
library(biascorrection)
## library(veri)
library(easyVerification)

nc_time <- function(nc){
  t.i <- grep('tim', names(nc$dim))
  t.units <- nc$dim[[t.i]]$units
  t.mul <- c(second=1, 
             minute=60, 
             hour=3600,
             day=1, 
             month=1, 
             year=1) ## special cases for days-years
  mul.i <- pmatch(substr(t.units,1,3), names(t.mul))
  if (length(mul.i) == 1){
    if (mul.i < 4){
      t.out <- as.POSIXct(gsub('.*since ', '', t.units), tz='UTC') + t.mul[mul.i]*nc$dim[[t.i]]$val
    } else if (mul.i == 4){
      ## t.out <- as.POSIXct(round(as.Date(gsub('.*since ', '', t.units), tz='UTC')) + round(nc$dim[[t.i]]$val), tz='UTC')
      t.ref <- as.POSIXct(gsub('.*since ', '', t.units), tz='UTC')
      t.out <- as.POSIXct(t.ref + 3600*24*nc$dim[[t.i]]$val, tz='UTC')
    }
  } else {
    stop('Time coordinate format not implemented yet')
  }
  return(t.out)
}

## initialise
scratchdir <- Sys.getenv('SCRATCH')

## read command line arguments:
args <- commandArgs(trailingOnly=TRUE)
## check if there are any command line arguments
if (length(args) == 0){
  args <- c('ecmwf-system4',
            'ERA-INT',
            'tas',
            'global2',
            'smooth_????-????_ERA-INT',
            '05',
            TRUE,
            TRUE,
            FALSE)
} else if (length(args) < 5){
  stop('Not enough command line arguments')
}

## disaggregate command line arguments
index <- args[1]
grid <- args[2]
obsname <- args[3]
method <- args[4]
initmon <- args[5]

if (length(args) >= 7){
  seasonals <- as.logical(args[6])
  ccrs <- as.logical(args[7])
} else {
  seasonals <- c(FALSE, TRUE)
  ccrs <- c(FALSE, TRUE)
}
if (length(args) == 8){
  detrends <- as.logical(args[8])
} else {
  detrends <- c(FALSE, TRUE)  
}

## stop for monthly forecasts and for detrended forecasts
if (all(!seasonals) | all(detrends)) q(save='no')

## replace placeholder for method string
dpath <- '/store/msclim/bhendj/EUPORIAS'
if (length(grep('????-????', method)) == 1 & method != 'none-forward'){
  mfiles <- system(paste0('ls ', dpath, '/', model, '/', grid, '/monthly/', index, '/', method, '/', index, '_????', initmon, '??_*.nc'), intern=TRUE)
  method <- sapply(strsplit(mfiles, '/'), function(x) x[length(x) - 1])[1]
}
stopifnot(is.character(method))

## specify the file paths
fpath <- paste(dpath, model, grid, 'monthly', index, method, sep='/')
opath <- paste(dpath, obsname, grid, 'monthly', index, sep='/')

## get forecast and observation files
fcfiles <- list.files(fpath, pattern=paste0('^', index, '_....', initmon, '.._.*_',grid,'_',method, '.nc'), full.name=TRUE)
if (method == "none-forward"){
  fpath <- paste(dpath, model, grid, 'monthly', index, 'none', sep='/')
  fcfiles <- list.files(fpath, pattern=paste0('^', index, '_....', initmon, '.._.*_',grid,'_none.nc'), full.name=TRUE)
}
obsfiles <- list.files(opath, pattern=paste0('^', index, '_'), full.name=TRUE)

## check whether there is at least one of each
if (length(fcfiles) < 1) stop('no forecast file found')
if (length(obsfiles) < 1) stop('no observation file found')

## only use most recent obsfile
details <- file.info(obsfiles)
details <- details[with(details, order(as.POSIXct(mtime), decreasing=TRUE)),]
obsfiles <- rownames(details)[1]
rm(details)

## enquire file properties (names and dims from forecast)
nc <- nc_open(fcfiles[1])
## get variable name in forecasts
## parfc <- names(nc$var)[names(nc$var) %in% varlist[[index]]]
parfc <- index
## get dimensions
fdims <- nc$var[[parfc]]$dim
fdimnames <- sapply(fdims, function(x) x$name)
names(fdims) <- fdimnames
if (length(fdimnames) == 3){
  nlon <- fdims[[grep('ncells', fdimnames)]]$len
  nlat <- 1
} else {
  nlon <- fdims[[grep('lon', fdimnames)]]$len
  nlat <- fdims[[grep('lat', fdimnames)]]$len  
}
nens <- fdims[[which(!fdimnames %in% grep('tim|lon|lat|ncells', fdimnames, value=TRUE))]]$len
## get more information on time (i.e. initialisation month)
time.i <- grep('tim', fdimnames)
ntime <- fdims[[time.i]]$len
fctime <- as.Date(nc_time(nc))
## number of months with all days present
initmn <- as.numeric(format(fctime[1], '%m'))
ncomplete <- sum(as.numeric(format(fctime, '%d')) >= 28)
nc_close(nc)

## name from observations
nc <- nc_open(obsfiles[1])
## parobs <- names(nc$var)[names(nc$var) %in% varlist[[index]]]
parobs <- index
obstime <- as.Date(nc_time(nc))
nc_close(nc)

## set up file connections
fc.con <- lapply(as.list(fcfiles), nc_open)
obs.con <- lapply(as.list(obsfiles), nc_open)

## get forecast times
fc.times <- lapply(fc.con, function(x) as.Date(nc_time(x)[1:ncomplete] ))
names(fc.times) <- sapply(fc.times, function(x) format(x[1], '%Y'))

## check whether forecasts are in obs
is.in.obs <- sapply(fc.times, function(x) all(x[1:ncomplete] %in% obstime))
if (mean(is.in.obs) < 0.5){
  obstime <- obstime + 1
  is.in.obs <- sapply(fc.times, function(x) all(x[1:ncomplete] %in% obstime))
}

## number of years
nyears <- length(fcfiles)
years <- as.numeric(names(fc.times))

## read in the observations (corresponding time steps)
if (length(obsfiles) == 1){
  obs <- array(NA, c(nlon, nlat, ncomplete, sum(is.in.obs)))
  fc.i <- sapply(fc.times[is.in.obs], function(x) format(x, '%Y%m%d') %in% unlist(sapply(obstime, format, '%Y%m%d')))
  obs.i <- sapply(obstime, function(x) format(x, '%Y%m%d') %in% unlist(sapply(fc.times[is.in.obs], format, '%Y%m%d')))
  if (length(obs.con[[1]]$dim) == 2){
    obs.tmp <- ncvar_get(obs.con[[1]], parobs)[,obs.i]
  } else {
    obs.tmp <- ncvar_get(obs.con[[1]], parobs)[,,obs.i]    
  }
  obs[fc.i[rep(1:nrow(fc.i), each=nlon*nlat), ]] <- as.vector(obs.tmp)
  obs <- aperm(obs, c(3,1,2,4))
  dimnames(obs)[[length(dim(obs))]] <- as.character(years[is.in.obs])
  oobs <- obs
  rm(obs)
} else {
  stop('Multiple observation files not implemented yet')
}

## fix units
if (index %in% c('tas', 'tasmax', 'tasmin') & max(oobs, na.rm=T) < 200){
  oobs <- oobs + 273.15
}

# read in all the data
ffcst <- array(NA, c(ncomplete, nlon, nlat, nens, nyears))
if (length(fc.con[[1]]$dim) == 3) {
  for (i in 1:nyears){
    ffcst[,,,,i] <- aperm(ncvar_get(fc.con[[i]], parfc, count=c(-1,-1,ncomplete)), c(3,1,2))
  }
} else {
  for (i in 1:nyears){
    ffcst[,,,,i] <- aperm(ncvar_get(fc.con[[i]], parfc, count=c(-1,-1,-1,ncomplete)), c(4,1,2,3))
  }  
}
rownames(ffcst) <- unique(sapply(fc.times, format, '%m'))[1:ncomplete]
dimnames(ffcst)[[5]] <- names(fc.times)

## get the grid data for display of forecasts
nc <- fc.con[[1]]
if (length(grep('eobs', grid)) == 1){
  lon <- nc$dim$rlon$vals
  lat <- nc$dim$rlat$vals
  plon <- ncatt_get(nc, 'rotated_pole', 'grid_north_pole_longitude')$val
  plat <- ncatt_get(nc, 'rotated_pole', 'grid_north_pole_latitude')$val
} else if (length(nc$dim) == 3) {
  lon <- ncvar_get(nc, 'lon')
  lat <- ncvar_get(nc, 'lat')
} else {
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals        
  plon <- plat <- NULL
}
rm(nc)

## close file connections (desperate debugging)
sapply(fc.con, nc_close)
sapply(obs.con, nc_close)
rm(fc.con, obs.con)

## get crossvalidation string
crossval <- length(grep('crossval', method)) == 1
suppressWarnings(nblock <- as.numeric(gsub('_.*', '', gsub('.*crossval', '', method))))
forward <- length(grep("forward", method)) == 1
stopifnot(sum(forward, crossval) <= 1)
strategy <- list(type=c("none", "crossval", "forward")[1 + 1*crossval + 2*forward], 
                 blocklength=nblock)

## loop over seasonal or monthly forecast skill
for (seasonal in seasonals){
  
  outdir <- paste(dpath, 'skill_scores', grid, if (seasonal) 'seasonal' else 'monthly', index, sep='/')
  
  if (seasonal){
    ## compute three-monthly aggregates
    fcst.seas <- array(NA, dim(ffcst) - c(2,0,0,0,0))
    obs.seas <- array(NA, dim(oobs) - c(2,0,0,0))
    for (i in 1:nrow(fcst.seas)){
      fcst.seas[i,,,,] <- colMeans(ffcst[i+0:2,,,,,drop=F], dims=1)
      obs.seas[i,,,] <- colMeans(oobs[i + 0:2,,,,drop=F], dims=1)
    }
    
    ## rearrange time for output
    fcst.seastimes <- sapply(fc.times, function(x) format(x[setdiff(1:ncomplete, 1:2)], '%Y-%m-%d'))
    fcst.seastimes <- array(0, dim(fcst.seastimes)) + as.Date(fcst.seastimes)  
  } else {
    fcst.seas <- ffcst[1:ncomplete,,,,,drop=F]
    obs.seas <- oobs[1:ncomplete,,,,drop=F]
    fcst.seastimes <- sapply(fc.times, function(x) format(x[1:ncomplete], '%Y-%m-%d'))
    fcst.seastimes <- array(0, dim(fcst.seastimes)) + as.Date(fcst.seastimes)  
  }
  
  ## mask forecasts and obs where more than 80% of the years are missing
  fcst.mask <- rowMeans(is.na(fcst.seas)*1,  dims=3) <= .2
  obs.mask <- rowMeans(is.na(obs.seas)*1,  dims=3) <= .2
  
  ## mask according to observation climatology
  if (index %in% c('HDD', 'CDD', 'HDDch')){
    ## new criterion at least 67% of forecasts in climatology are non-zero
    nclim <- min(30, dim(obs.seas)[4])
    obs.mask <- obs.mask & rowMeans(obs.seas[,,,1:nclim] > 0, dims=3, na.rm=T) > 0.67
    ##obs.mask <- obs.mask & rowMeans(obs.seas, dims=3, na.rm=T) > 0.1
    ## this amounts to ~3 HDD per month, or 9 HDD per season
  } else if (index %in% c('FD', 'SD', 'HD')){
    obs.mn <- rowMeans(obs.seas, dims=3, na.rm=T)
    obs.mask <- obs.mask & obs.mn > 0.05 & obs.mn < 0.95
  }
  
  ## like this, all years are masked out (and therefore the corresponding scores also)
  obs.seas[rep(!(obs.mask & fcst.mask), length=length(obs.seas))] <- NA
  fcst.seas[rep(!(obs.mask & fcst.mask), length=length(fcst.seas))] <- NA
  
  ## set up temporary directory
  tmpdir <- paste0(scratchdir, '/skill_', index, '_', grid, '_', method, '_initmon', initmon, '_', ceiling(runif(1, min=0, max=9999)))
  if (! file.exists(tmpdir)) dir.create(tmpdir)
  
  ## loop through raw / detrended forecasts and obs
  for (detrend in detrends){
    
    if (detrend){
      ## detrend obs
      obs.seas <- aperm(apply(obs.seas, 1:3, function(y){
        if (all(!is.na(y))){
          mean(y) + lm(y ~ seq(along=y))$res
        } else {
          NA*y
        }
      }), c(2,3,4,1))
      ## detrend forecasts
      fcst.seas <- aperm(apply(fcst.seas, 1:4, function(y){
        if(all(!is.na(y))){
          mean(y) + lm(y ~ seq(along=y))$res
        } else {
          NA * y
        }
      }), c(2,3,4,5,1))
    }
    
    ## loop through raw / CCR forecasts
    for (is.ccr in ccrs){
      
      if (is.ccr){
        fcst.ccr <- array(NA, dim(fcst.seas))
        for (loi in 1:nlon){
          for (lai in 1:nlat){
            if (any(!is.na(obs.seas[,loi,lai,]))){
              fcst.ccr[,loi,lai,,] <- aperm(debias(
                fcst=aperm(fcst.seas[,loi, lai,,which(years < 2011)], c(1,3,2)),
                obs=obs.seas[,loi,lai,which(years<2011)],
                fcst.out=aperm(fcst.seas[,loi,lai,,], c(1,3,2)),
                method='ccr', 
                strategy=strategy), c(1,3,2))
            } ## end of if on missing values in obs
          } ## end of loop on latitudes
        } ## end of loop on longitudes
        
        fcst.seas <- fcst.ccr
      }
      
      ## save tercile forecasts to file for later display
      fcst.prob <- apply(fcst.seas, 1:3, function(x) {
        if(all(!is.na(x))) {
          return(convert2prob(t(x), threshold=quantile(c(x[,years < 2011]), 1:2/3, type=8)) / sum(!is.na(x[,1]))*100) 
        } else {
          return(rep(NA, 3*ncol(x)))
        }})
      fcst.prob <- array(fcst.prob, c(nrow(fcst.prob)/3, 3, dim(fcst.prob)[-1]))
      ## set up output file path
      Rdatastem <- paste0(paste(index, grid, obsname, method, 'initmon', sep='_'), initmon)
      if (detrend) Rdatastem <- gsub(paste0(index, '_'), paste0(index, '_detrend_'), Rdatastem)
      if (is.ccr) Rdatastem <- gsub(paste0(index, '_'), paste0(index, '_CCR_'), Rdatastem)
      outRdatapath <- paste(dpath, 'forecasts', model, grid, if (seasonal) 'seasonal' else 'monthly', index, method, Rdatastem, sep='/')
      if (!file.exists(outRdatapath)) dir.create(outRdatapath, recur=TRUE)
      leads <- apply(outer(1:nrow(fcst.seas), 0:2, '+'), 1, paste, collapse='')
      for (yi in seq(years)){
        for (li in seq(leads)){
          fcst <- ffcst[li,,,,yi,drop=F]
          obs <- if (yi <= dim(obs.seas)[length(dim(obs.seas))]) obs.seas[li,,,yi,drop=F] else NULL
          ## tricky selection to make sure we select the right 
          prob <- fcst.prob[yi,,li,,,drop=F]
          fcst.time <- fcst.seastimes[li,yi]
          save(fcst, prob, fcst.time, obs, 
               file=paste0(outRdatapath, '/', 
                           Rdatastem, '_', years[yi], '_', leads[li], '.Rdata'))
        }
      }
      
    } ## end of loop on CCR
    
  } ## end of loop on detrend
  
  ## remove temporary directory
  system(paste0('rm -rf ', tmpdir))
  
} ## end of loop on seasonal/monthly skill

warnings()
print('Finished successfully')

q(save='no')
