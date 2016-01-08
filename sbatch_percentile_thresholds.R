## compute percentile thresholds for different input data
## 
## 

## load packages
library(myhelpers) ## for input output etc.

## initialise
varlist <- list(tasmin=c('MN2T24', 'mn2t24', 'tmin', 'tasmin', 'tn'),
                tasmax=c('MX2T24', 'mx2t24', 'tmax', 'tasmax', 'tx'),
                tas=c('MEAN2T24', 'mean2t24', 'tas', 'tg', 't2m'),
                pr=c('TOT_PREC', 'pr', 'tp', 'rr'), 
                dtr=c('dtr'))

scratchdir <- Sys.getenv('SCRATCH')



## read command line arguments:
## 1. variable name
## 2. grid
## 3. method
## 4. initialization
## 5. startyear
## 6. endyear
args <- commandArgs(trailingOnly=TRUE)

## check if there are any command line arguments
if (length(args) == 0 | mode(args) == 'function'){
  args <- c('tas', 
            'global2', 
            'none', 
            '11')
} else if (length(args) < 4){
  stop('Not enough command line arguments')
}

## disaggregate command line arguments
varname <- args[1]
grid <- args[2]
method <- args[3]
initmon <- args[4]
if (length(args) == 4){
  startyear <- 1981
  endyear <- 2010  
}

## get file names
dpath <- '/store/msclim/bhendj/EUPORIAS'
fpath <- paste(dpath, 'ecmwf-system4', grid, 'daily', varname, 'none', sep='/')
outdir <- paste(dpath, 'ecmwf-system4', grid, 'fx', sep='/')

## get forecast and observation files
fcfiles <- list.files(fpath, pattern=paste0('^....', initmon, '.._.*_', grid), full.name=TRUE)

## check whether there is at least 10 forecast files
if (length(fcfiles) < 10) stop('no forecast file found')


## set up temporary directory for output files
tmpdir <- paste0(scratchdir, '/percentile_', method, '_initmon', initmon, '_', ceiling(runif(1)*1000))
if (!file.exists(tmpdir)) dir.create(tmpdir, recursive=TRUE)

## first get the time dimension of all forecasts
fc.con <- lapply(as.list(fcfiles), nc_open)

## correct forecast time stamp (offset by 6 hours)
if (varname %in% c('tasmin', 'tasmax', 'tas', 'dtr', 'pr')){
  time.offset <- -6*3600
} else {
  warning('Please check time offset in forecast data set')
}

## get time axes
fc.times <- lapply(fc.con, function(x) as.Date(nc_time(x) + time.offset))
names(fc.times) <- sapply(fc.times, function(x) format(x[1], '%Y'))

## check whether years are in range
is.in.range <- names(fc.times) %in% paste(startyear:endyear)

## exclude forecast years that are not in obs
fcfiles <- fcfiles[is.in.range]
fc.times <- fc.times[is.in.range]
fc.con <- fc.con[is.in.range]

nlead <- min(sapply(fc.times, length))

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
nenses <- sapply(fc.con, function(x) x$dim[[which(!names(x$dim) %in% grep('tim|lon|lat|ncells', names(x$dim), value=TRUE))]]$len)
stopifnot(nenses == nenses[1])

## get extent of fcst array
nlon <- nlons[1]
nlat <- nlats[1]
nens <- nenses[1]
ntime <- nlead
nyears <- length(fcfiles)
yrange <- range(as.numeric(names(fc.times)))

## set up output directory
outpath <- paste(outdir, varname, paste('none', paste(yrange, collapse='-'), sep='_'), sep='/')
if (!file.exists(outpath)) dir.create(outpath, recursive=TRUE, mode='0775')

## compute how many latitudes to read in at once
latchunksize <- min(nlat,floor(1e9 / (nlon*nens*ntime*nyears)))

## loop on  latitude band
for (lati in seq(1, nlat, latchunksize)){
  
  ## get list of object names that are not deleted at end of iteration
  no.rm <- ls()
  
  print(paste0(round((lati -1)/nlat*100), '% done'))
  
  ## compute effective number of latitudes to read in
  nchunk <- length(lati:min(lati+latchunksize - 1, nlat))
  
  ## read in all the forecasts
  fcst <- array(NA, c(length(fc.con), nlon, nchunk, nens, ntime))
  for (fi in seq(along=fc.con)){
    if (length(fc.con[[1]]$dim) == 3){
      fcst[fi,,,,] <- ncvar_get(fc.con[[fi]], varid=parfcs[fi],
                                start=c(1,1,1), 
                                count=c(nlon, nens, ntime))
    } else {
      fcst[fi,,,,] <- ncvar_get(fc.con[[fi]], varid=parfcs[fi], start=c(1,lati,1,1), count=c(nlon,nchunk,nens,ntime))
    }
  }
  rownames(fcst) <- names(fc.times)
  
  ## set up array for processing of seasonal quantiles
  fc.monarr <- t(sapply(fc.times, format, '%m'))
  fc.seasarr <- array(as.numeric(fc.monarr) %/% 3 %% 4, dim(fc.monarr))
  seas.ti <- fc.seasarr[rep(1:nyears, nens), ]
  seas.ind <- as.character(seas.ti[1,])

  
  for (pctl in c(1,2,5,10,25,75,90,95,98,99)){
    
    ## write output to output file
    if (length(fc.con[[1]]$dim) == 3 & nchunk == 1){
      starti <- c(1,1,1)
      counti <- c(-1,-1,-1)
    } else {
      starti <- c(1,lati,1,1)
      counti <- c(-1,nchunk,-1,ntime)
    }
    fcst.quant <- apply(fcst, c(2,3,5), quantile, pctl/100, na.rm=T)   
    
    ## write to outfile
    outfile <- paste0("enspctl", formatC(pctl, width=2, flag='0'),
                      '_ecmwf-system4_', varname, '_', 
                      paste(yrange, collapse='-'), '.nc')
    nc_write(nctempfile=fcfiles[1], 
             file=paste(tmpdir, outfile, sep='/'), 
             varname=parfcs[1], 
             data=fcst.quant[,rep(1:nchunk, nens),],
             append=TRUE, 
             start=starti,
             count=counti)
    
    
    ## compute seasonal percentile boundaries
    fcst.seasquant <- aperm(apply(fcst, c(2,3), function(x){
      tapply(c(x), seas.ti, quantile, pctl/100, na.rm=T)
    })[seas.ind, , ], c(2,3,1))
    ## write to outfile
    outfile <- paste0("seasenspctl", formatC(pctl, width=2, flag='0'),
                      '_ecmwf-system4_', varname, '_', 
                      paste(yrange, collapse='-'), '.nc')
    nc_write(nctempfile=fcfiles[1], 
             file=paste(tmpdir, outfile, sep='/'), 
             varname=parfcs[1], 
             data=fcst.seasquant[,rep(1:nchunk, nens),],
             append=TRUE, 
             start=starti,
             count=counti)
    
    
    
  }
  
  ## delete all old objects that are not needed in next iteration (to free memory)
  rm(list=setdiff(ls(), no.rm))
  gc()
  
} ## end of loop on latitude chunks


q(save='no')

