## script to bias correct seasonal forecasts (currently ECMWF system4)
## 
## feature requests:
##   - arbitrary variables and grids
##   - various bias correction methods
##   - various observation datasets
##
## REVISIONS
## 23/09/2014: Changed processing to scratch for faster access time 
## 23/09/2014: Changed processnig to run for multiple latitudes in one chunk
## 


## load packages
library(methods) ## for processing with Rscript
library(ncdf4) ## ncdf is deprecated but still supported
library(biascorrection) ## for additional bias correction methods

## initialise
varlist <- list(tasmin=c('MN2T24', 'mn2t24', 'tmin', 'tasmin', 'tn'),
                tasmax=c('MX2T24', 'mx2t24', 'tmax', 'tasmax', 'tx'),
                tas=c('MEAN2T24', 'mean2t24', 'tas', 'tg', 't2m'),
                pr=c('TOT_PREC', 'pr', 'tp', 'rr'))
spandays <- 91 ## number of days used in loess smoothing (centered) 
## roughly one month for 'unbiased' monthly means
scratchdir <- Sys.getenv('SCRATCH')
comp.time <- list()

## additional functions (may be moved to package later on)
#fast std dev?
colSdColMeans <- function(x, na.rm=TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x)) 
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x*x, na.rm=na.rm) - (colMeans(x, na.rm=na.rm))^2
  return(sqrt(colVar * n/(n-1)))
}

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
      t.out <- as.POSIXct(as.Date(gsub('.*since ', '', t.units), tz='UTC') + round(nc$dim[[t.i]]$val), tz='UTC')
    }
  } else {
    stop('Time coordinate format not implemented yet')
  }
  return(t.out)
}

## function to quickly write data to a NetCDF file
## with the same dimensions as in 
nc_write <- function(nctemplate, file, varname, data, append=FALSE, ...){
  
  ## check if an accordingly named variable exists
  if (!any(names(nctemplate$var) == varname)) stop('Variable to write does not exist')
  
  ## open netcdf file (and close on exit of function)
  on.exit(nc_close(ncout))
  if (append & file.exists(file)){
    ncout <- nc_open(filename=file, write=TRUE)
  } else {
    ncout <- nc_create(filename=file, nctemplate$var)
  }
  
  ## write the data
  ncvar_put(ncout, varid=varname, vals=data, ...)
  
}

## read command line arguments:
## 1. variable name
## 2. debiasing method (see below)
## 3. start and end year for de-biasing
## 4. list of input files
## 5. list of observation files
## 6. output directory stem (see below for how this is expanded)
args <- commandArgs(trailingOnly=TRUE)

## check if there are any command line arguments
if (length(args) == 0 | mode(args) == 'function'){
  args <- c('tasmin',
            'smooth',
            1981, 2010,
            system('\\ls /store/msclim/bhendj/EUPORIAS/ecmwf-system4/eobs0.44/daily/Tmin/????1101_52_eobs0.44.nc', intern=TRUE),
            '/store/msclim/bhendj/EUPORIAS/E-OBS/daily/tn_0.44deg_rot_v10.0.nc',
            '/store/msclim/bhendj/EUPORIAS/ecmwf-system4/eobs0.44/daily')
} else if (length(args) < 7){
  stop('Not enough command line arguments')
}

## disaggregate command line arguments
varname <- args[1]
method <- args[2]
yrange <- as.numeric(args[3:4])
files <- args[-c(1:4, length(args))]
outdir <- args[length(args)]

## further disaggregate the input files (into obs and fcst files)
## based on similarity of directory tree
fbits <- strsplit(files, '/')
## common length of directory structure
nmin <- min(sapply(fbits, length)) - 1 ## do not include the file name
direqual <- sapply(fbits, function(x) sum(x[1:nmin] == fbits[[1]][1:nmin]))
fcfiles <- files[direqual == max(direqual)]
obsfiles <- files[direqual < max(direqual)]
if (length(grep('E-OBS', obsfiles)) > 0){
  obsname <- 'E-OBS'
} else {
  obsname <- 'ERA-INT'
}

## check whether there is at least one of each
if (length(fcfiles) < 1) stop('no forecast file found')
if (length(obsfiles) < 1) stop('no observation file found')

## set up temporary directory for output files
tmpdir <- paste0(scratchdir, '/bias_correction_', method, '_', paste(yrange, collapse='-'), '_', ceiling(runif(1)*1000))
if (!file.exists(tmpdir)) dir.create(tmpdir, recursive=TRUE)

## enquire file properties (names and dims from forecast)
nc <- nc_open(fcfiles[1])
## get variable name in forecasts
parfc <- names(nc$var)[names(nc$var) %in% varlist[[varname]]]
## get dimensions
fdims <- nc$var[[parfc]]$dim
fdimnames <- sapply(fdims, function(x) x$name)
names(fdims) <- fdimnames
nlon <- fdims[[grep('lon', fdimnames)]]$len
nlat <- fdims[[grep('lat', fdimnames)]]$len
nens <- fdims[[which(!fdimnames %in% grep('tim|lon|lat', fdimnames, value=TRUE))]]$len
## get more information on time (i.e. initialisation month)
time.i <- grep('tim', fdimnames)
ntime <- fdims[[time.i]]$len
## correct forecast time stamp
if (varname %in% c('tasmin', 'tasmax', 'tas')){
  time.offset <- -86400
} else {
  warning('Please check time offset in forecast data set')
}
fctime <- nc_time(nc) + time.offset
initmn <- as.numeric(format(fctime[1], '%m'))
tvec <- 1:ntime
nyears <- length(fcfiles)
nc_close(nc)

## name from observations
nc <- nc_open(obsfiles[1])
parobs <- names(nc$var)[names(nc$var) %in% varlist[[varname]]]
obstime <- nc_time(nc)
nc_close(nc)

## set up file connections
fc.con <- lapply(as.list(fcfiles), nc_open)
obs.con <- lapply(as.list(obsfiles), nc_open)

## get time axes
fc.times <- lapply(fc.con, function(x) nc_time(x) + time.offset)
names(fc.times) <- sapply(fc.times, function(x) format(x[1], '%Y'))
## check that all the years for debiasing are available (else change yrange)
if (any(! as.character(yrange[1]:yrange[2]) %in% names(fc.times))){
  warning('Not all the proposed years are available for debiasing')
  yrange <- range(seq(yrange[1], yrange[2])[as.character(yrange[1]:yrange[2]) %in% names(fc.times)])
  warning(paste('Reset reference period to', paste(yrange, collapse=' to ')))
}
## produce time array for debiasing function
## slightly involved as dates get transformed to numeric
## in most built-in operations
fc.timarr <- array(0, dim(sapply(fc.times, as.Date))) + 
  as.Date(sapply(fc.times, format, '%Y-%m-%d'))
colnames(fc.timarr) <- names(fc.times)


debias.years <- as.character(yrange[1]:yrange[2])
obs.times <- lapply(obs.con, nc_time)


## set up output directory
outpath <- paste(outdir, varname, paste(method, paste(yrange, collapse='-'), obsname, sep='_'), sep='/')
if (!file.exists(outpath)) dir.create(outpath, recursive=TRUE, mode='0775')

## change file names for output
outfiles <- gsub('\\.nc', paste0('_', method, '_', paste(yrange, collapse='-'), '_', obsname, '.nc'), sapply(strsplit(fcfiles, '/'), function(x) x[length(x)]))

## compute how many latitudes to read in at once
latchunksize <- floor(1e9 / (nlon*nens*ntime*nyears))

## loop on  latitude band
for (lati in seq(1, nlat, latchunksize)){
  
  ## get list of object names that are not deleted at end of iteration
  no.rm <- ls()
  
  print(paste0(round((lati -1)/nlat*100), '% done'))
  
  ## compute effective number of latitudes to read in
  nchunk <- length(lati:min(lati+latchunksize - 1, nlat))
  
  ## read in all the observations for the corresponding years
  obs <- array(NA, c(diff(yrange) + 1, nlon, nchunk, ntime))
  if (length(obs.con) == 1){
    obs.i <- range(which(obs.times[[1]] %in% unlist(fc.times[debias.years])))
    obs.tmp <- ncvar_get(obs.con[[1]], varid=parobs, start=c(1,lati,obs.i[1]), count=c(-1,nchunk,diff(obs.i) + 1))
    obs.tim <- obs.times[[1]][obs.i[1]:obs.i[2]] 
    if (nchunk == 1){
      for (fi in seq(along=debias.years)) obs[fi,,,] <- obs.tmp[,obs.tim %in% fc.times[[fi]]]  
    } else {
      for (fi in seq(along=debias.years)) obs[fi,,,] <- obs.tmp[,,obs.tim %in% fc.times[[fi]]]        
    }
    rm(obs.tmp)
    gc()
  } else {
    stop('Case with multiple obs not implemented yet')
  }
  
  print('Read forecasts')
  comp.time[['read']][[paste(lati)]] <- system.time({
    ## read in all the forecasts
    fcst <- array(NA, c(nyears, nlon, nchunk, nens, ntime))
    for (fi in seq(along=fc.con)){
      if (nchunk == 1){
        fcst[fi,,,,] <- ncvar_get(fc.con[[fi]], varid=parfc, start=c(1,lati,1,1), count=c(-1,nchunk,-1,-1))[,,1:ntime]
      } else {
        fcst[fi,,,,] <- ncvar_get(fc.con[[fi]], varid=parfc, start=c(1,lati,1,1), count=c(-1,nchunk,-1,-1))[,,,1:ntime]
      }
    }
    rownames(fcst) <- names(fc.times)
  }) ## system time for reading forecasts
  print(comp.time[['read']][[paste(lati)]])
  
  print('Compute bias')
  comp.time[['compute']][[paste(lati)]] <- system.time({
    ## check units of observations
    if (varname %in% c('tasmin', 'tasmax', 'tas')){
      if (max(obs, na.rm=T) < 100) obs <- obs + 273.15
    }
    mm <- gsub('-crossval.*$', '', method)
    strategy <- list(type=ifelse(length(grep('-crossval.*$', method)) == 1,
                                 "crossval", "none"))
    if (crossval){
      nblock <- as.numeric(gsub('.*-crossval', '', method))
      if (is.na(nblock)) nblock <- 1
      strategy$blocklength <- nblock
    }
    
    fcst.debias <- array(NA, dim(fcst))
    for (loi in 1:nlon){
      for (lai in 1:nchunk){
        if (all(!is.na(obs[,loi,lai,]))){
          fcst.debias[,loi,lai,,] <- aperm(debias(
            fcst=aperm(fcst[debias.years,loi, lai,,], c(3,1,2)),
            obs=t(obs[,loi,lai,]),
            method=mm,
            strategy=strategy,
            fcst.out=aperm(fcst[,loi,lai,,], c(3,1,2)),
            fc.time=fc.timarr[, debias.years],
            fcout.time=fc.timarr,
            span=min(1, spandays/ntime)), c(2,3,1))
        } ## end of if on missing values in obs
      } ## end of loop on latitudes
    } ## end of loop on longitudes
    
    }) ## system time for computation
  print(comp.time[['compute']][[paste(lati)]])
  
  print('Writing output files')
  comp.time[['write']][[paste(lati)]] <- system.time({
    
    ## write output to output file (directly to output file path due to size limits on scratch)
    for (fi in seq(along=outfiles)){
      nc_write(nctemplate=fc.con[[fi]], 
               file=paste(outpath, outfiles[fi], sep='/'), 
               varname=parfc, 
               data=fcst.debias[fi,,,,],
               append=TRUE,
               start=c(1,lati,1,1),
               count=c(-1,nchunk,-1,-1))
    } ## end of loop on years
  }) ## system time for writing
  print(comp.time[['write']][[paste(lati)]])
  
  ## delete all old objects that are not needed in next iteration (to free memory)
  rm(list=setdiff(ls(), no.rm))
  gc()
  
} ## end of loop on latitude chunks

## loop through files and convert to NetCDF4 (with maximal compression)
print('Compress output')
comp.time[['compress']][[1]] <- system.time({
  for (f in outfiles){
    system(paste0('mv ', outpath, '/', f, ' ', tmpdir, '/', f))
    system(paste0('nccopy -d9 ', tmpdir, '/', f, ' ', outpath, '/', f))
    file.remove(paste0(tmpdir, '/', f))
  }
  system(paste('rm -rf', tmpdir))
})


## write output with statistics
outstat <- sapply(comp.time, function(y) paste0(round(sum(sapply(y, function(x) x['elapsed']))/60, 1), ' min.'))
print(outstat)

q(save='no')

