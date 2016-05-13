## script to compute percentile based indices from forecasts

## load packages
library(myhelpers) ## ncdf4 and additional functionality
library(methods)

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


## read command line arguments:
args <- commandArgs(trailingOnly=TRUE)

## check if there are any command line arguments
if (length(args) == 0 | mode(args) == 'function'){
  args <- c('ecmwf-system4', 
            'tasmin', 
            'global2', 
            '11')
} else if (length(args) < 4){
  stop('Not enough command line arguments')
}

## disaggregate command line arguments
model <- args[1]
varname <- args[2]
grid <- args[3]
init <- args[4]
if (length(args) > 4) {
  percentiles <- as.numeric(args[5:length(args)])
} else {
  percentiles <- c(1, 2, 5, 10, 25, 75, 90, 95, 98, 99)
}

## get file names
dpath <- '/store/msclim/bhendj/EUPORIAS'
fcfiles <- get_fcfiles(model=model, grid=grid, index=varname, method='none', 
                       init=init, granularity='daily', source='euporias')[1:30]
outdir <- paste(dpath, model, grid, 'fx', varname, sep='/')
if (!file.exists(outdir)){
  dir.create(outdir, recursive=TRUE)
}
templatefile <- paste0(dpath, '/grids/pctl_template_', grid, '.nc')
if (!file.exists(templatefile)){
  stop("Template file does not exist as of yet")
}

## read in the land sea mask to get grid dimensions
lsm <- read_ncdf(paste0("/store/msclim/bhendj/EUPORIAS/grids/", grid, '_lsm.nc'))

## loop through latitudes
for (lai in 1:ncol(lsm)){

  ## read in the data
  comp.time[['read']][[paste(lai)]] <- system.time({
    fcst <- read_ncdf(fcfiles, index=varname, lai=lai)
  })
  
  ntim <- dim(fcst)[3]
  years <- array(as.numeric(format(attr(fcst, 'time'), '%Y')), dim(fcst)[3:4])[1,]
  yearrange <- paste(range(years), collapse='-')
  
  ## loop through the days of the year
  comp.time[['pctl']][[paste(lai)]] <- system.time({
    fcst.pctl <- apply(fcst, c(1,3), quantile, percentiles/100)
  })
  
  comp.time[['smooth']][[paste(lai)]] <- system.time({
    fcst.pctlsmooth <- apply(fcst.pctl, 1:2, function(x) loess(x ~ seq(x), span=91/length(x))$fit)
  })
  
  
  ## write the percentiles to the respective files
  comp.time[['write']][[paste(lai)]] <- system.time({
    for (pii in seq(along=percentiles)){
      outfile <- paste(paste0('pctl', formatC(percentiles[pii], width=2, flag=0)),
                       model, varname, yearrange, paste0('initmon', init, '.nc'), sep='_')
      
      nc_write(nctempfile=templatefile,
               file=paste0(outdir,'/', outfile),
               varname=NULL,
               data=t(fcst.pctlsmooth[,pii,]),
               append=(lai > 1), 
               start=c(1,lai,1), count=c(-1,1,ntim))
      
      ## fix time, names and attributes in last iteration
      if (lai == ncol(lsm)){
        nc <- nc_open(paste0(outdir, '/', outfile), write=TRUE)
        ncvar_put(nc, 'time', 1:nc$dim$time$len - 1)
        nvarname <- names(nc$var)[1]
        nc_close(nc)
        
        ## fix time units
        system(paste0("ncatted -h -a units,time,o,c,'days since ", 
                      attr(fcst,'time')[1], "' ", outdir, '/', outfile))
        
        ## fix variable name
        system(paste0("ncrename -h -v ", nvarname, ',', varname, ' ', outdir, '/', outfile))
        
      } # end of if on last latitude
      
    }  # end of latitude loop
  })
  
  
}



## say good bye
print("successfully reached the end of this script")
# q(save='no')
