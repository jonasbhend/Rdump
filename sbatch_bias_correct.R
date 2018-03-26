## script to bias correct seasonal forecasts (currently ECMWF system4)
## 
## REVISIONS
## 25/03/2015: Added support for station data
## 23/09/2014: Changed processnig to run for multiple latitudes in one chunk
## 12/02/2014: Changed processing to include all available forecasts 
##             And reset reference period to 1981-2010
##             Check for lead time and ensemble members in forecasts
## 


## load packages
library(methods) ## for processing with Rscript
library(myhelpers) ## ncdf4 and additional functionality
library(biascorrection) ## for additional bias correction methods

## initialise
varlist <- list(tasmin=c('MN2T24', 'mn2t24', 'tmin', 'tasmin', 'tn'),
                tasmax=c('MX2T24', 'mx2t24', 'tmax', 'tasmax', 'tx'),
                tas=c('MEAN2T24', 'mean2t24', 'tas', 'tg', 't2m'),
                pr=c('TOT_PREC', 'pr', 'tp', 'rr'), 
                dtr=c('dtr'))
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

# nc_time <- function(nc){
#   t.i <- grep('tim', names(nc$dim))
#   t.units <- nc$dim[[t.i]]$units
#   t.mul <- c(second=1, 
#              minute=60, 
#              hour=3600,
#              day=1, 
#              month=1, 
#              year=1) ## special cases for days-years
#   mul.i <- pmatch(substr(t.units,1,3), names(t.mul))
#   if (length(mul.i) == 1){
#     if (mul.i < 4){
#       t.out <- as.POSIXct(gsub('.*since ', '', t.units), tz='UTC') + t.mul[mul.i]*nc$dim[[t.i]]$val
#     } else if (mul.i == 4){
#       t.out <- as.POSIXct(as.Date(gsub('.*since ', '', t.units), tz='UTC') + round(nc$dim[[t.i]]$val), tz='UTC')
#     }
#   } else {
#     stop('Time coordinate format not implemented yet')
#   }
#   return(t.out)
# }


## read command line arguments:
args <- commandArgs(trailingOnly=TRUE)

## check if there are any command line arguments
if (length(args) == 0 | mode(args) == 'function'){
  args <- c('UCAN-PP-LR-15PC', 
            'E-OBS',
            'tasmin', 
            'UK', 
            'fastqqmap', 
            '11',
            '1981-2010')
} else if (length(args) < 5){
  stop('Not enough command line arguments')
}

## decide whether to continue from previous computation
pickup <- FALSE

## disaggregate command line arguments
model <- args[1]
obsname <- args[2]
varname <- args[3]
grid <- args[4]
method <- args[5]
initmon <- args[6]
if (length(args) > 6){
  if (args[7] == "????-????"){
    startyear <- 1981
    endyear <- 2010
  } else {
    syear <- strsplit(args[7], '-')
    startyear <- as.numeric(syear[[1]][1])
    endyear <- as.numeric(syear[[1]][2])
  }
} else {
  startyear <- 1981
  endyear <- 2010
}

## get file names
dpath <- '/store/msclim/bhendj/EUPORIAS'
fpath <- paste(dpath, model, grid, 'daily', varname, 'none', sep='/')
opath <- paste(dpath, obsname, grid, 'daily', varname, sep='/')
outdir <- paste(dpath, model, grid, 'daily', sep='/')

## get forecast and observation files
fcfiles <- list.files(fpath, pattern=paste0('^....', initmon, '.._.*_', grid), full.name=TRUE)
obsfiles <- list.files(opath, pattern=paste0('.nc$'), full.name=TRUE)
if (length(grep('timmean', obsfiles)) > 0) obsfiles <- obsfiles[-grep('timmean', obsfiles)]

## check whether there is at least one of each
if (length(fcfiles) < 1) stop('no forecast file found')
if (length(obsfiles) < 1) stop('no observation file found')


## set up temporary directory for output files
tmpdir <- paste0(scratchdir, '/bias_correction_', method, '_initmon', initmon, '_', ceiling(runif(1)*1000))
if (pickup) {
  tmpdir <- c(list.files(scratchdir, paste0('bias_correction_', method, '_initmon', initmon, '_'), full.names=TRUE), tmpdir)[1]
}

if (!file.exists(tmpdir)) dir.create(tmpdir, recursive=TRUE)

## first get the time dimension of all forecasts
## and match with observations
fc.con <- lapply(as.list(fcfiles), nc_open)
obs.con <- lapply(as.list(obsfiles), nc_open)

## correct forecast time stamp (offset by 6 hours)
if (model == 'ecmwf-system4' & varname %in% c('tasmin', 'tasmax', 'tas', 'dtr', 'pr')){
  time.offset <- -6*3600
} else {
  time.offset <- 0
}

## get time axes
fc.times <- lapply(fc.con, function(x) as.Date(nc_time(x) + time.offset))
names(fc.times) <- sapply(fc.times, function(x) format(x[1], '%Y'))

## get times in observation files
obs.times <-lapply(obs.con, function(x) as.Date(nc_time(x)))
otimes <- as.Date(unlist(sapply(obs.times, format, '%Y-%m-%d')))
stopifnot(sum(duplicated(otimes)) == 0)


## now exclude all forecast years that are not fully present in obs
is.in.obs <- sapply(fc.times, function(x) all(x %in% otimes)) 
## constrain years to 1981-, but use all years if less than 25 years are available
if (sum(is.in.obs) > 25) is.in.obs <- is.in.obs & names(fc.times) %in% startyear:endyear
## if (sum(is.in.obs) > 25) is.in.obs <- is.in.obs & names(fc.times) %in% 1981:2014


## exclude forecast years that are not in obs
## fcfiles <- fcfiles[is.in.obs]
## fc.times <- fc.times[is.in.obs]
## fc.con <- fc.con[is.in.obs]

nlead <- min(sapply(fc.times, length))


## produce time array for debiasing function
## slightly involved as dates get transformed to numeric
## in most built-in operations
fc.timarr <- array(0, dim(sapply(fc.times, function(x) as.Date(x[1:nlead])))) + 
  as.Date(sapply(fc.times, function(x) format(x[1:nlead], '%Y-%m-%d')))
colnames(fc.timarr) <- names(fc.times)

## derive range of years for calibration description
yrange <- range(as.numeric(names(fc.times)[is.in.obs]))
debias.years <- as.character(yrange[1]:yrange[2])

## get parameter names
parfcs <- sapply(fc.con, function(x) names(x$var)[names(x$var) %in% varlist[[varname]]])
## get dimensions
if (length(grep('lon', names(fc.con[[1]]$dim))) == 0){
  nlons <- sapply(fc.con, function(x) x$dim[[grep('ncells', names(x$dim))]]$len)
  nlats <- rep(1, length(nlons))
} else {
  nlons <- sapply(fc.con, function(x) x$dim[[grep('lon', names(x$dim))]]$len)
  nlats <- sapply(fc.con, function(x) x$dim[[grep('lat', names(x$dim))]]$len)  
}
ntimes <- sapply(fc.con, function(x) x$dim[[grep('tim', names(x$dim))]]$len)
stopifnot(nlons == nlons[1], nlats == nlats[1])
nenses <- sapply(fc.con, function(x) x$dim[[which(!names(x$dim) %in% grep('tim|lon|lat|ncells|nb2|bnds', names(x$dim), value=TRUE))]]$len)
## stopifnot(nenses == nenses[1])

## get extent of fcst array
nlon <- nlons[1]
nlat <- nlats[1]
maxnens <- max(nenses)
minnens <- min(nenses)
ntime <- nlead
nyears <- length(fcfiles)

## name from observations
parobs <- names(obs.con[[1]]$var)[names(obs.con[[1]]$var) %in% varlist[[varname]]]

## get units
ounits <- ncatt_get(obs.con[[1]], parobs, attname='units')


## set up output directory
outpath <- paste(outdir, varname, paste(method, paste(yrange, collapse='-'), obsname, sep='_'), sep='/')
if (!file.exists(outpath)) dir.create(outpath, recursive=TRUE, mode='0775')

## change file names for output
outfiles <- gsub('\\_none.nc', paste0('_', method, '_', paste(yrange, collapse='-'), '_', obsname, '.nc'), sapply(strsplit(fcfiles, '/'), function(x) x[length(x)]))

## compute how many latitudes to read in at once
latchunksize <- min(nlat,floor(6e8 / (nlon*maxnens*ntime*nyears)))

if (pickup){
  nctmp <- nc_open(paste(tmpdir, outfiles[length(fc.con)], sep='/'))
  ftmp <- ncvar_get(nctmp, parfcs[length(fc.con)],
                    start=c(1,1,1,1), count=c(1,-1,1,1))
  nc_close(nctmp)
  initlat <- min(sum(!is.na(ftmp)) + 1, nlat)
} else {
  initlat <- 1
}


## loop on  latitude band
for (lati in seq(initlat, nlat, latchunksize)){
  
  ## get list of object names that are not deleted at end of iteration
  no.rm <- ls()
  
  print(paste0(round((lati -1)/nlat*100), '% done'))
  
  ## compute effective number of latitudes to read in
  nchunk <- length(lati:min(lati+latchunksize - 1, nlat))
  
  ## read in all the observations for the corresponding years
  obs <- array(NA, c(diff(yrange) + 1, nlon, nchunk, ntime))
  if (length(obs.con) == 1){
    obs.i <- range(which(as.Date(obs.times[[1]]) %in% unlist(fc.times[debias.years])))
    obs.tim <- obs.times[[1]][obs.i[1]:obs.i[2]] 
    if (length(obs.con[[1]]$dim) == 2 & nchunk == 1){
      obs.tmp <- ncvar_get(obs.con[[1]], varid=parobs, start=c(1,obs.i[1]), count=c(nlon,diff(obs.i) + 1))
    } else {
      obs.tmp <- ncvar_get(obs.con[[1]], varid=parobs, start=c(1,lati,obs.i[1]), count=c(nlon,nchunk,diff(obs.i) + 1))      
    }
    if (length(dim(obs.tmp)) == 2){
      for (fi in seq(along=debias.years)) obs[fi,,,] <- obs.tmp[,obs.tim %in% fc.times[[debias.years[fi]]]] 
    } else {
      for (fi in seq(along=debias.years)) obs[fi,,,] <- obs.tmp[,,obs.tim %in% fc.times[[debias.years[fi]]]]        
    }
    rm(obs.tmp)
    gc()
  } else {
    stop('Case with multiple obs not implemented yet')
  }
  
  print('Read forecasts')
  comp.time[['read']][[paste(lati)]] <- system.time({
    ## read in all the forecasts
    fcst <- array(NA, c(length(fc.con), nlon, nchunk, maxnens, ntime))
    for (fi in seq(along=fc.con)){
      if (length(fc.con[[1]]$dim) == 3){
        fcst[fi,,,1:nenses[fi],] <- ncvar_get(fc.con[[fi]], varid=parfcs[fi],
                                 start=c(1,1,1), 
                                 count=c(nlon, nenses[fi], ntime))
      } else {
        fcst[fi,,,1:nenses[fi],] <- ncvar_get(fc.con[[fi]], varid=parfcs[fi], start=c(1,lati,1,1), count=c(nlon,nchunk,nenses[fi],ntime))
      }
    }
    rownames(fcst) <- names(fc.times)
  }) ## system time for reading forecasts
  print(comp.time[['read']][[paste(lati)]])
  
  ## clear memory
  gc()

  ## constrain values to non-negative for precipitation
  if (varname == 'pr'){
    warning(paste0(round(mean(fcst < 0, na.rm=T)*100, 1), '% negative values constrained to zero in input'))
    fcst <- pmax(fcst, 0)
  }
  
  ## clear memory
  gc()
  
  print('Compute bias')
  comp.time[['compute']][[paste(lati)]] <- system.time({
    ## check units of observations
    if (varname %in% c('tasmin', 'tasmax', 'tas')){
      if (max(obs, na.rm=T) < 100) obs <- obs + 273.15
    } else if (varname == "pr" & obsname == "E-OBS") {
      obs <- obs/1000
    }
    mm <- gsub("-forward$", "", gsub('-crossval.*$', '', method))
    crossval <- length(grep('-crossval.*$', method)) == 1
    forward <- length(grep("-forward$", method)) == 1
    block <- length(grep("-block.*$", method)) == 1
    stopifnot(sum(crossval, forward, block) <= 1)
    
    strategy <- list(type=c("none", "crossval", "block", "forward")[1 + crossval*1 + block*2 + forward*3])
    if (crossval | block){
      nblock <- as.numeric(gsub(paste0('.*-', ifelse(crossval, 'crossval', 'block')), '', method))
      if (is.na(nblock)) nblock <- 1
      strategy$blocklength <- nblock
    }
    
    fcst.debias <- array(NA, dim(fcst))
    for (loi in 1:nlon){
      for (lai in 1:nchunk){
        if (any(!is.na(obs[,loi,lai,]) & apply(!is.na(fcst[debias.years, loi, lai,1:minnens,]), c(1,3), all))){
          fcst.debias[,loi,lai,,] <- aperm(debias(
            fcst=aperm(fcst[debias.years,loi, lai,1:minnens,], c(3,1,2)),
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
    
    ## clear memory
    rm(fcst)
    gc()
    
    
    ## constrain values to non-negative for precipitation
    if (varname == 'pr'){
      warning(paste0(round(mean(fcst.debias < 0, na.rm=T)*100, 1), '% negative values constrained to zero in output'))
      fcst.debias <- pmax(fcst.debias, 0)
    }
    
    ## clear memory
    gc()
    
    }) ## system time for computation
  print(comp.time[['compute']][[paste(lati)]])
  
  print('Writing output files')
  comp.time[['write']][[paste(lati)]] <- system.time({
    
    ## write output to output file (directly to output file path due to size limits on scratch)
    if (length(fc.con[[1]]$dim) == 3 & nchunk == 1){
      starti <- c(1,1,1)
      counti <- c(-1,-1,-1)
    } else {
      starti <- c(1,lati,1,1)
      counti <- c(-1,nchunk,-1,dim(fcst.debias)[5])
    }
    for (fi in seq(along=outfiles)){
      nc_write(nctempfile=fcfiles[fi], 
               file=paste(tmpdir, outfiles[fi], sep='/'), 
               varname=parfcs[fi], 
               data=fcst.debias[fi,,,1:nenses[fi],],
               append=TRUE, 
               start=starti,
               count=counti)
      ## fix output units (will be in units of the observations after BC)
      if (ounits$hasatt){
        ncout <- nc_open(paste(tmpdir, outfiles[fi], sep='/'), write=TRUE)
        ncatt_put(ncout, parfcs[fi], attname='units', attval=ounits$value)
        nc_close(ncout)
      }
            
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
    if (grid != 'eobs0.44'){
      system(paste0('cp ', tmpdir, '/', f, ' ', outpath, '/', f))
    } else {
      system(paste0('nccopy -d9 ', tmpdir, '/', f, ' ', outpath, '/', f))
    }
    file.remove(paste0(tmpdir, '/', f))
  }
  system(paste('rm -rf', tmpdir))
})  

## write output with statistics
outstat <- sapply(comp.time, function(y) paste0(round(sum(sapply(y, function(x) x['elapsed']))/60, 1), ' min.'))
print(outstat)

q(save='no')

