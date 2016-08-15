## script to compute skill metrics from debiased or raw forecast metrics

## load packages
library(methods) ## for processing with Rscript
library(ncdf4) ## ncdf is deprecated but still supported
library(biascorrection)
## library(veri)
library(verification)
library(easyVerification)

CRPSdecomp <- function(ens, obs, ...){
  cdec <- try(crpsDecomposition(obs, ens)[1:3], silent=TRUE)
  if (class(cdec) == 'try-error'){
    cdec <- list(CRPS=NA, CRPSpot=NA, Reli=NA)
  }
  return(cdec)
}

SpearCorr <- function(ens, obs, ...){
  cor(rowMeans(ens), obs, method='spearman', ...)
}

EnsTrend <- function(ens=ens, obs=obs, ...){
  stopifnot(is.matrix(ens), is.vector(obs), length(obs) == nrow(ens))
  nfcst <- length(obs)
  ebeta <- var(ens, seq(1,nfcst)) / var(seq(1,nfcst))
  obeta <- var(obs, seq(1,nfcst)) / var(seq(1,nfcst))
  return(list(EnsTrend = mean(ebeta) - obeta, EnsTrend.sigma=sd(ebeta)))
}

EnsCond <- function(ens=ens, obs=obs, ...){
  stopifnot(is.matrix(ens), is.vector(obs), length(obs) == nrow(ens))
  ens.mn <- rowMeans(ens)
  ebeta <- var(ens, ens.mn) / var(ens.mn)
  obeta <- var(obs, ens.mn) / var(ens.mn)
  return(list(EnsCond=mean(ebeta) - obeta, EnsCond.sigma=sd(ebeta)))
}

EnsVarlog <- function(ens=ens, obs=obs, ...){
  stopifnot(is.matrix(ens), is.vector(obs), length(obs) == nrow(ens))
  obs.std <- sd(obs, na.rm=T)
  ens.std <- apply(ens, 2, sd, na.rm=T)
  return(list(EnsVarlog=log(mean(ens.std**2) / obs.std**2), EnsVarlog.sigma=sd(log(ens.std**2))))
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
      ## t.out <- as.POSIXct(round(as.Date(gsub('.*since ', '', t.units), tz='UTC')) + round(nc$dim[[t.i]]$val), tz='UTC')
      t.ref <- as.POSIXct(gsub('.*since ', '', t.units), tz='UTC')
      t.out <- as.POSIXct(t.ref + 3600*24*nc$dim[[t.i]]$val, tz='UTC')
    }
  } else {
    stop('Time coordinate format not implemented yet')
  }
  return(t.out)
}


## scores <- c("Ens2AFC", "FairRpss", "EnsCorr", "SpearCorr")

scorelist <- list(Ens2AFC='generalized discrimination score', 
                  EnsCorr='correlation', 
                  SpearCorr="Spearman's correlation",                
                  EnsCrpss='continuous ranked probability skill score', 
                  FairCrpss='fair continuous ranked probability skill score', 
                  FairCrps='mean fair continuous ranked probability score',
                  EnsMae='mean absolute error',
                  EnsMaess='mean absolute error skill score',
                  EnsMe='mean error (bias)',
                  EnsMse='mean squared error',
                  EnsMsess='mean squared error skill score',
                  EnsRmse='root mean squared error',
                  EnsRmsess='root mean sqaured error skill score',
                  EnsRps='ranked probability score',
                  EnsRpss='ranked probability skill score',
                  FairRpss='fair ranked probability skill score',
                  EnsSprErr='spread error ratio',
                  EnsTrend="mean error in linear trend",
                  EnsCond="conditional error",
                  EnsVarlog="log ratio of variances",
                  EnsStdfrac="fraction of standard deviation",
                  FairSprErr='fair spread error ratio',
                  EnsRocss.cat1='ROC area skill score (lower tercile)',
                  EnsRocss.cat2='ROC area skill score (middle tercile)',
                  EnsRocss.cat3='ROC area skill score (upper tercile)', 
                  CRPSdecomp.CRPS='mean continuous ranked probability score',
                  CRPSdecomp.CRPSpot='potential CRPS (Resolution - Uncertainty)', 
                  CRPSdecomp.Reli='reliability term of the CRPS')


## scorefunction for use with veriApply
scorefun <- function(f, o, score, prob=NULL, ...){
  if (length(grep('Rps', score)) == 1 | score == 'EnsRocss'){
    prob <-  c(1/3,2/3)    
    return(veriApply(score, f, o, ensdim=4, tdim=5, prob=prob))
  } else {
    return(veriApply(score, f, o, ensdim=4, tdim=5, prob=prob, ...))
  }
}


## initialise
scratchdir <- Sys.getenv('SCRATCH')

## read command line arguments:
args <- commandArgs(trailingOnly=TRUE)
## check if there are any command line arguments
if (length(args) == 0){
  args <- c('ecmwf-system4', 
            'ERA-INT',
            'ITV',
            'global2',
            'smooth-forward_1981-2014_ERA-INT',
            '11',
            'small',
            TRUE,
            TRUE,
            FALSE)
} else if (length(args) < 5){
  stop('Not enough command line arguments')
}

## disaggregate command line arguments
model <- args[1]
obsname <- args[2]
index <- args[3]
grid <- args[4]
method <- args[5]
initmon <- args[6]
if (length(args) > 6){
  whatscores <- as.character(args[7])
} else {
  whatscores <- 'standard'
}

if (length(args) >= 9){
  seasonals <- as.logical(args[8])
  indmethod <- args[9]
} else {
  seasonals <- c(FALSE, TRUE)
  indmethod <- 'none'
}
ind.cali <- (indmethod != 'none')

if (length(args) == 10){
  detrends <- as.logical(args[10])
} else {
  detrends <- c(FALSE, TRUE)  
}

## check scores
if (whatscores == 'full'){
  scores <- c('EnsTrend', 'EnsCond', 'EnsVarlog', 'FairCrps', "CRPSdecomp", 'Ens2AFC', 'FairRpss', 'FairCrpss', 'EnsCorr', 'EnsMe', 'EnsMse', 'EnsMae', 'EnsRmse', 'EnsRmsess', 'FairSprErr', 'EnsRocss')
} else if (whatscores == 'small') {
  scores <- c("FairCrpss", "EnsCorr", "FairSprErr", "CRPSdecomp")
} else if (whatscores == 'standard') {
  scores <- c('Ens2AFC', 'FairRpss', 'FairCrpss', 'EnsCorr', 'EnsMe', 'EnsMse', 'EnsMae', 'EnsMsess', 'FairSprErr', 'EnsRocss')
}


## replace placeholder for method string
dpath <- '/store/msclim/bhendj/EUPORIAS'
if (length(grep('????-????', method)) == 1 & method != 'none-forward'){
  mfiles <- system(paste0('ls ', dpath, '/', model,'/', grid, '/monthly/', index, '/', method, '/', index, '_????', initmon, '??_*.nc'), intern=TRUE)
  method <- sapply(strsplit(mfiles, '/'), function(x) x[length(x) - 1])[1]
}
stopifnot(is.character(method))

## specify the file paths
fpath <- paste(dpath, model, grid, 'monthly', index, method, sep='/')
opath <- paste(dpath, obsname, grid, 'monthly', index, sep='/')

## get forecast and observation files
fcfiles <- list.files(fpath, pattern=paste0('^', index, '_....', initmon, '.._.*_',grid,'_',method, '.nc'), full.name=TRUE)[1:34]
if (method == "none-forward"){
  fpath <- paste(dpath, model, grid, 'monthly', index, 'none', sep='/')
  fcfiles <- list.files(fpath, pattern=paste0('^', index, '_....', initmon, '.._.*_',grid,'_none.nc'), full.name=TRUE)
}
## dirty hack for paper
## fcfiles <- fcfiles[1:min(34, length(fcfiles))]
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
nens <- fdims[[which(!fdimnames %in% grep('tim|lon|lat|ncells|nb2', fdimnames, value=TRUE))]]$len
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
## hack to account for missing day in DWD-CCLM4-8-21
if (model == 'DWD-CCLM4-8-21'){
  is.in.obs <- sapply(fc.times, function(x) all(x[1:ncomplete] %in% obstime | x[1:ncomplete] %in% obstime -1))  
}
if (mean(is.in.obs) < 0.5){
  obstime <- obstime + 1
  is.in.obs <- sapply(fc.times, function(x) all(x[1:ncomplete] %in% obstime))
}

stopifnot(sum(!is.in.obs) < 3)

## reduce to overlapping set
fcfiles <- fcfiles[is.in.obs]
fc.con <- fc.con[is.in.obs]
fc.times <- fc.times[is.in.obs]

## number of years
nyears <- length(fcfiles)
years <- as.numeric(names(fc.times))
if(length(grep('none', method)) == 1) {
  myears <- if (length(years) > 30) years[years %in% 1981:2010] else years  
} else {
  mtmp <- strsplit(method, '_')[[1]]
  mtmp <- as.numeric(strsplit(mtmp[length(mtmp) - 1], '-')[[1]])
  myears <- mtmp[1]:mtmp[2]
} 


## get the forecast ensemble members
nens <- min(sapply(fc.con, function(x) x$dim[[which(! names(x$dim) %in% grep("tim|lon|lat|ncells|nb2", names(x$dim), value=TRUE))]]$len)) 

## read in the observations (corresponding time steps)
if (length(obsfiles) == 1){
  obs <- array(NA, c(nlon, nlat, ncomplete, nyears))
  fc.i <- sapply(fc.times, function(x) format(x, '%Y%m') %in% unlist(sapply(obstime, format, '%Y%m')) & seq(along=as.numeric(x)) %in% 1:ncomplete)
  obs.i <- sapply(obstime, function(x) format(x, '%Y%m') %in% unlist(sapply(fc.times, function(x) format(x[1:ncomplete], '%Y%m'))))    
  if (length(obs.con[[1]]$dim) == 2){
    obs.tmp <- ncvar_get(obs.con[[1]], parobs)[,obs.i]
  } else {
    obs.tmp <- ncvar_get(obs.con[[1]], parobs)[,,obs.i]    
  }
  obs[fc.i[rep(1:nrow(fc.i), each=nlon*nlat), ]] <- as.vector(obs.tmp)
  obs <- aperm(obs, c(3,1,2,4))
} else {
  stop('Multiple observation files not implemented yet')
}

## fix units
if (index %in% c('tas', 'tasmax', 'tasmin') & max(obs, na.rm=T) < 200){
  obs <- obs + 273.15
}

## read in all the data
ffcst <- array(NA, c(ncomplete, nlon, nlat, nens, nyears))
if (length(fc.con[[1]]$dim) == 3) {
  for (i in 1:nyears){
    ffcst[,,,,i] <- aperm(ncvar_get(fc.con[[i]], parfc, count=c(-1,nens,ncomplete)), c(3,1,2))
  }
} else {
  for (i in 1:nyears){
    ffcst[,,,,i] <- aperm(ncvar_get(fc.con[[i]], parfc, count=c(-1,-1,nens,ncomplete)), c(4,1,2,3))
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
forward <- length(grep("forward", method)) == 1
block <- length(grep("block", method)) == 1
stopifnot(sum(forward, crossval, block) <= 1)
strategy <- list(nfcst=length(fc.times),
                 type=c("none", 'crossval', 'block', 'forward')[1 + 1*crossval + 2*block + 3*forward])
if (crossval | block) {
  suppressWarnings(strategy$blocklength <- as.numeric(gsub('_.*', '', gsub('.*crossval', '', method))))
  if (is.na(strategy$blocklength)) strategy$blocklength <- 1
}
if (strategy$type != "none"){
  fcyears <- as.numeric(names(fc.times))
  syears <- strsplit(method, '_')
  syears <- as.numeric(strsplit(syears[[1]][length(syears[[1]]) - 1], '-')[[1]])
  strategy$indices <- which(fcyears %in% syears[1]:syears[2])
}
## loop over seasonal or monthly forecast skill
for (seasonal in seasonals){
  
  outdir <- paste(dpath, 'skill_scores', grid, if (seasonal) 'seasonal' else 'monthly', index, sep='/')
  
  if (seasonal){
    ## compute three-monthly aggregates
    ## compute three-monthly aggregates
    if (grid == 'EAF-22'){
      nred <- 0
    } else {
      nred <- 2
    }
    fcst.seas <- array(NA, dim(ffcst) - c(nred,0,0,0,0))
    obs.seas <- array(NA, dim(obs) - c(nred,0,0,0))
    for (i in 1:(nrow(ffcst) - 2)){
      fcst.seas[i,,,,] <- colMeans(ffcst[i+0:2,,,,,drop=F], dims=1)
      obs.seas[i,,,] <- colMeans(obs[i + 0:2,,,,drop=F], dims=1)
    }
    if (grid == 'EAF-22'){
      fcst.seas[nrow(fcst.seas),,,,] <- colMeans(ffcst[2:5,,,,,drop=F], dims=1)
      obs.seas[nrow(obs.seas),,,] <- colMeans(obs[2:5,,,,drop=F], dims=1)      
    }
    
    ## rearrange time for output
    fcst.seastimes <- sapply(fc.times, function(x) format(if (grid == 'EAF-22') x[setdiff(1:ncomplete, 1:2)][c(seq(1,ncomplete-2), NA, 3)] else x[setdiff(1:ncomplete, 1:2)], '%Y-%m-%d'))
    fcst.seastimes <- array(0, dim(fcst.seastimes)) + as.Date(fcst.seastimes)  
  } else {
    fcst.seas <- ffcst[1:ncomplete,,,,,drop=F]
    obs.seas <- obs[1:ncomplete,,,,drop=F]
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
  } else if (index == 'WDF') {
    nclim <- min(30, dim(obs.seas)[4])
    obs.mask <- obs.mask & rowMeans(obs.seas[,,,1:nclim] > 0, dims=3, na.rm=T) > 0.67 & rowMeans(obs.seas[,,,1:nclim] < 1, dims=3, na.rm=T) > 0.67
  }
  
  ## like this, all years are masked out (and therefore the corresponding scores also)
  obs.seas[rep(!(obs.mask & fcst.mask), length=length(obs.seas))] <- NA
  
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
    
    ## find out whether seasonal/monthly indices need calibration
    if (ind.cali){
      fcst.cali <- array(NA, dim(fcst.seas))
      for (loi in 1:nlon){
        for (lai in 1:nlat){
          if (any(!is.na(obs.seas[,loi,lai,]))){
            fcst.cali[,loi,lai,,] <- aperm(debias(
              fcst=aperm(fcst.seas[,loi, lai,,which(years %in% myears)], c(1,3,2)),
              obs=obs.seas[,loi,lai,which(years %in% myears)],
              fcst.out=aperm(fcst.seas[,loi,lai,,], c(1,3,2)),
              method=indmethod, 
              strategy=strategy, 
              type= if(forward | crossval) "prediction" else "calibration"), c(1,3,2))
          } ## end of if on missing values in obs
        } ## end of loop on latitudes
      } ## end of loop on longitudes
      
      fcst.seas <- fcst.cali
      gc()
    }
    
    ## get number of years that are analysed (reduced for forward calibration)
    nyears <- dim(obs.seas)[4]
    ## replace if not succesful
    ## yind <- seq(if (forward) 16 else 1, nyears)
    yind <- seq(1,nyears)
    
    #       ## only write out seasonally aggregated forecasts
    #       if (seasonal){
    #         ## save tercile forecasts to file for later display
    #         fcst.prob <- apply(fcst.seas[,,,,yind, drop=F], 1:3, function(x) if(all(!is.na(x))) convert2prob(t(x), prob=c(1/3,2/3)) else rep(NA, 3*ncol(x))) / nens * 100
    #         fcst.prob <- array(fcst.prob, c(nrow(fcst.prob)/3, 3, dim(fcst.prob)[-1]))
    #         ## set up output file path
    #         Rdatastem <- paste0(paste(index, grid, obsname, method, 'initmon', sep='_'), initmon)
    #         if (detrend) Rdatastem <- gsub(paste0(index, '_'), paste0(index, '_detrend_'), Rdatastem)
    #         if (is.ccr) Rdatastem <- gsub(paste0(index, '_'), paste0(index, '_CCR_'), Rdatastem)
    #         outRdatapath <- paste(dpath, 'forecasts', model, grid, if (seasonal) 'seasonal' else 'monthly', index, method, Rdatastem, sep='/')
    #         if (!file.exists(outRdatapath)) dir.create(outRdatapath, recur=TRUE)
    #         leads <- apply(outer(1:nrow(fcst.seas), 0:2, '+'), 1, paste, collapse='')
    #         for (yi in yind){
    #           for (li in seq(leads)){
    #             fcst <- ffcst[li,,,,yi,drop=F]
    #             ## tricky selection to make sure we select the right 
    #             ## forecast probabilities for the forward method
    #             prob <- fcst.prob[yi - min(yind) + 1,,li,,,drop=F]
    #             fcst.time <- fcst.seastimes[li,yi]
    #             save(fcst, prob, fcst.time,
    #                  file=paste0(outRdatapath, '/', 
    #                              Rdatastem, '_', years[yi], '_', leads[li], '.Rdata'))
    #           }
    #         }
    #         ## save(fcst.seas, fcst.prob, fcst.seastimes, lon, lat, plon, plat, file=outRdata)  
    #       }
    
    ## set up output file to write to
    nc <- nc_open(obsfiles[1])
    nc2 <- nc_open(fcfiles[1])
    if (length(nc$dim) == 2){
      ncells.nc <- nc$dim$ncells
    } else if (length(grep('eobs', grid)) == 1){
      lon.nc <- nc2$dim$rlon
      lat.nc <- nc2$dim$rlat    
    } else {
      lon.nc <- nc2$dim$lon
      lat.nc <- nc2$dim$lat
    }
    time.nc <- ncdim_def(name='time', 
                         units='days since 1980-01-01', 
                         vals=as.numeric(fcst.seastimes[, ncol(fcst.seastimes)] - as.Date('1980-01-01')),
                         unlim=TRUE, 
                         longname=if (seasonal) 'end date of 3-monthly seasons' else 'end date of month')
    if (length(nc$dim) == 2){
      dims.nc <- list(ncells.nc, time.nc)
    } else {
      dims.nc <- list(lon.nc, lat.nc, time.nc)      
    }
    ## create variables
    vars.nc <- list()
    ## check whether this is a rotated grid
    if (length(grep('eobs', grid)) == 1){
      vars.nc[['lon']] <- ncvar_def('lon', 'degrees_east', dims.nc[1:2], missval=-1e20, longname='longitude')
      vars.nc[['lat']] <- ncvar_def('lat', 'degrees_north', dims.nc[1:2], missval=-1e20, longname='latitude')
      vars.nc[['rotated_pole']] <- ncvar_def('rotated_pole', units='', dim=list(), missval=NULL, prec='char')    
    } else if (length(nc$dim) == 2){
      vars.nc[['lon']] <- nc$var$lon
      vars.nc[['lat']] <- nc$var$lat
    }
    for (score in scores){
      if (score == 'EnsRocss'){
        for (i in 1:3) {
          rscore <- paste0(score, '.cat', i)
          vars.nc[[rscore]] <- ncvar_def(rscore, '1', dims.nc, missval=-1e20, longname=scorelist[[rscore]])
          vars.nc[[paste0(rscore, '.sigma')]] <- ncvar_def(paste0(rscore, '.sigma'), '1',
                                                           dims.nc, missval=-1e20,
                                                           longname=scorelist[[rscore]])
        }
      } else if (score == 'CRPSdecomp'){
        for (cii in c('CRPS', 'CRPSpot', 'Reli')){
          rscore <- paste(score, cii, sep='.')
          vars.nc[[rscore]] <- ncvar_def(rscore, '1', dims.nc, missval=-1e20, 
                                         longname=scorelist[[rscore]])
        }
      } else {
        vars.nc[[score]] <- ncvar_def(score, '1', dims.nc, missval=-1e20, longname=scorelist[[score]])
        if (score %in% c(scores[grep('pss$', scores)], 'EnsTrend', 'EnsCond', 'EnsVarlog')){
          vars.nc[[paste(score, 'sigma', sep='.')]] <- ncvar_def(paste(score, 'sigma', sep='.'), '1', dims.nc, missval=-1e20, longname=scorelist[[score]])
        } 
      }
    }
    
    
    ## get additional information on content of files
    fc.name <- sapply(strsplit(fcfiles, '/')[1], function(x) x[grep('EUPORIAS', x) + 1])
    obs.name <- sapply(strsplit(obsfiles, '/')[1], function(x) x[grep('EUPORIAS', x) + 1])
    ofile <- paste0(
      paste(
        paste0(index, 
               ifelse(detrend, '_detrend', ''), 
               ifelse(ind.cali, paste0('_', toupper(indmethod)), '')), 
        method, fc.name, 'vs', obs.name, 
        paste(range(names(fc.times)), collapse='-'), 
        sep="_"), '_initmon', initmon, '.nc')
    outfile <- paste(tmpdir, ofile, sep='/')
    desc <- paste0('Skill in ', 
                   ifelse(detrend, 'detrended ', ''), 
                   ifelse(ind.cali, paste0('bias corrected (', indmethod, ')'), ''), 
                   ifelse(seasonal, 'seasonal ', 'monthly '), 
                   index, 
                   ' from ', 
                   ifelse(method == 'none', '', paste0('bias-corrected (', method, ') ')), 
                   gsub('-', ' ', toupper(fc.name)), 
                   ' forecasts verified against ', obs.name, 
                   ' for ', paste(range(names(fc.times)[yind]), collapse='-'))
    
    ## create netcdf file
    ncout <- nc_create(outfile, vars.nc)
    ## write variables
    if (length(nc$dim) == 2){
      ncvar_put(ncout, varid='lon', ncvar_get(nc, 'lon'))
      ncvar_put(ncout, varid='lat', ncvar_get(nc, 'lat'))
    } else if (length(grep('eobs', grid)) == 1){
      if ('Actual_longitude' %in% names(nc$var)){
        ncvar_put(ncout, varid='lon', ncvar_get(nc, 'Actual_longitude'))
        ncvar_put(ncout, varid='lat', ncvar_get(nc, 'Actual_latitude'))
      } else {
        library(geocors)
        plon <- ncatt_get(nc, 'rotated_pole', attname='grid_north_pole_longitude')$val
        plat <- ncatt_get(nc, 'rotated_pole', attname='grid_north_pole_longitude')$val    
        lola <- geocors.trafo(rep(nc$dim$rlon$vals, nc$dim$rlat$len), 
                              rep(nc$dim$rlat$vals, each=nc$dim$rlon$len),
                              from.type='rotpol', 
                              from.pars=list(plon=as.numeric(plon), plat=as.numeric(plat)),
                              to.type='lonlat')
        ncvar_put(ncout, varid='lon', lola$lon)
        ncvar_put(ncout, varid='lat', lola$lat)    
      }
      ## fill in the rotated grid coordinates
      ncatt_put(ncout, varid='rotated_pole', attname='grid_mapping_name', 
                attval='rotated_latitude_longitude', prec='text')
      ncatt_put(ncout, varid='rotated_pole', attname='grid_north_pole_longitude', 
                attval=-162, prec='double')
      ncatt_put(ncout, varid='rotated_pole', attname='grid_north_pole_latitude', 
                attval=39.25, prec='double')
    }
    ## add in a global attribute with description
    ncatt_put(ncout, varid=0, attname='description', attval=desc, prec='text')
    nc_close(ncout)
    
    ## compute scores and write to output file
    for (score in scores){
      ncout <- nc_open(outfile, write=TRUE)
      print(score)
      print(system.time( sfo <- scorefun(fcst.seas[,,,,yind, drop=F], obs.seas[,,,yind,drop=F], score, strategy=strategy)))
      if (is.list(sfo)){
        sfo <- lapply(sfo, function(x){
          x[x == -Inf] <- -9999
          x[abs(x) == Inf] <- NA
          return(x)
        })
        if (score == 'EnsRocss'){
          for (i in 1:3){
            ncvar_put(ncout, varid=paste0(score, '.cat', i), vals=aperm(sfo[[paste0('cat', i)]], c(2,3,1)))
            ncvar_put(ncout, varid=paste0(score, '.cat', i, '.sigma'), vals=aperm(sfo[[paste0('cat', i, '.sigma')]], c(2,3,1)))
          }
        } else if (score == 'CRPSdecomp') {
          for (nn in names(sfo)){
            ncvar_put(ncout, varid=paste(score, nn, sep='.'), 
                      vals=aperm(sfo[[nn]], c(2,3,1)))
          }
        } else {
          ncvar_put(ncout, varid=score, vals=aperm(sfo[[1]], c(2,3,1)))
          ncvar_put(ncout, varid=paste(score, 'sigma', sep='.'), vals=aperm(sfo[[2]], c(2,3,1)))    
        }
      } else {
        ## set negative infinity (division by zero) to large negative value
        sfo[sfo == -Inf] <- -9999
        sfo[abs(sfo) == Inf] <- NA
        ## write to output file
        if (length(sfo) == length(obs.seas[,,,yind])){
          sfo <- rowMeans(sfo, dims=3)
        }
        ncvar_put(ncout, varid=score, vals=aperm(sfo, c(2,3,1))) 
      }
      nc_close(ncout)      
      
      ## readjust grid specification for variables
      ## this is very important to happen AFTER the file connection to the NetCDF
      ## is closed from R - otherwise the output is shifted latitudinally
      if (length(nc$dim) == 2 | length(grep('eobs', grid)) == 1){
        if (length(nc$dim) == 2){
          atts <- c(gridtype='unstructured', coordinates="lon lat")
        } else {
          atts <- c(grid_mapping='rotated_pole')
        }
        for (attn in names(atts)){
          if (score == 'EnsRocss'){
            for (i in 1:3) system(paste0("ncatted -h -a ", attn, ",", score, ".cat", i, ",c,c,'",atts[attn],"' ", outfile))            
          } else if (is.list(sfo)){
            system(paste0("ncatted -h -a ", attn, ",", score, ".sigma,c,c,'",atts[attn],"' ", outfile))      
            system(paste0("ncatted -h -a ", attn, ",", score, ",c,c,'",atts[attn],"' ", outfile))      
          } else {
            system(paste0("ncatted -h -a ", attn, ",", score, ",c,c,'",atts[attn],"' ", outfile))      
          }                
        }
      }
      
      ## make sure we free as much memory as possible/reasonable
      rm(sfo, ncout)
      gc()
      
      
    } ## end of loop on scores
    
    ## copy temporary file to permanent location
    ## make sure the output directory exists
    if (!file.exists(outdir)) dir.create(outdir, recursive=TRUE)
    system(paste0('mv ', outfile, ' ', outdir, '/', ofile))
    
    
  } ## end of loop on detrend
  
  ## remove temporary directory
  system(paste0('rm -rf ', tmpdir))
  
} ## end of loop on seasonal/monthly skill

warnings()
print('Finished successfully')

q(save='no')
