## script to compute percentile based indices from observations

## load packages
library(myhelpers) ## ncdf4 and additional functionality

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
  args <- c('ERA-INT', 
            'tasmin', 
            'global2')
} else if (length(args) < 3){
  stop('Not enough command line arguments')
}

## disaggregate command line arguments
observation <- args[1]
varname <- args[2]
grid <- args[3]
if (length(args) > 3) {
  percentiles <- as.numeric(args[4:length(args)])
} else {
  percentiles <- c(1, 2, 5, 10, 25, 75, 90, 95, 98, 99)
}

## get file names
dpath <- '/store/msclim/bhendj/EUPORIAS'
fpath <- paste(dpath, observation, grid, 'daily', varname, sep='/')
obsfile <- list.files(fpath, full.names=TRUE)
outdir <- paste(dpath, observation, grid, 'fx', varname, sep='/')
if (!file.exists(outdir)){
  dir.create(outdir, recursive=TRUE)
}
templatefile <- paste0(dpath, '/grids/pctl_template_', grid, '.nc')
if (!file.exists(templatefile)){
  stop("Template file does not exist as of yet")
}

## read in the land sea mask to set up a pseudo-mask (all TRUE)
lsm <- read_ncdf(paste0("/store/msclim/bhendj/EUPORIAS/grids/", grid, '_lsm.nc'))

## read all the data in one go
obs <- read_ncdf(obsfile, index=varname, 
                 mask=!is.na(lsm),
                 tlim=as.Date(c("1981-01-01", "2010-12-31")))
ntim <- length(attr(obs, 'time'))

## set up array for percentiles
obs.pctl <- array(NA, c(length(percentiles), nrow(obs), 366))

## loop through the days of the year
system.time({
  ii <- which(format(attr(obs, 'time'), '%m%d') == '0101')
  for (i in 1:366){
    iinds <- sort(c(outer(ii - 1 + i, seq(-5,5), '+')))
    iinds[iinds < 1] <- iinds[iinds < 1] + ntim
    iinds[iinds > ntim] <- iinds[iinds > ntim] - ntim
    obs.pctl[,,i] <- apply(obs[,iinds], 1, quantile, percentiles/100)
  }  
})

## smooth the resulting estimate
system.time({
  obs.pctlsmooth <- apply(obs.pctl, 1:2, 
                          function(x) loess(rep(x,3) ~ seq(along=rep(x,3)), 
                                            span=91/length(x)/3)$fit[seq(x) + length(x)])  
})

## write the percentiles to the respective files
for (pii in seq(along=percentiles)){
  outfile <- paste(paste0('pctl', formatC(percentiles[pii], width=2, flag=0)),
                   observation, varname, '1981-2010.nc', sep='_')
  
  nc_write(nctempfile=templatefile,
           file=paste0(outdir,'/', outfile),
           varname=NULL,
           data=t(obs.pctlsmooth[,pii,]),
           append=FALSE)
  
  ## fix variable names
  nc <- nc_open(templatefile)
  nvarname <- names(nc$var)[1]
  nc_close(nc)
  system(paste0('ncrename -h -v ', nvarname, ',', varname, ' ', outdir, '/', outfile))
}

## say good bye
print("successfully reached the end of this script")
# q(save='no')

